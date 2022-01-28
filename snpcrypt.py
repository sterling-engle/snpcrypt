# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 17:02:18 2021

@author: Sterling Engle
@uid: 904341227

Developed and tested only on Windows 10 running Ubuntu 18.04.6 LTS
under Python 3.9.7.

usage: snpcrypt.py [-h] [-b] [-c COUNT] [-d DELIM] [-e] [-f FILE]
                   [-g GENKEYPATH] [-i] [-k KEYPATH] [-l LOG] [-p POS]
                   [-r REMOVE] [-s SNPS] [-t] [-v] [-w PASSWORD] [-x REGION]
                   [inputfile]

positional arguments:
  inputfile             BAM/SAM or SNPs VCF input file [1000-genomes-phase-
												3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a
												.20130502.genotypes.vcf.gz]

optional arguments:
  -h, --help            show this help message and exit
  -b, --bases           include DNA bases with SNPs instead of indices [false]
  -c COUNT, --count COUNT
                        count of samples to output [all]
  -d DELIM, --delim DELIM
                        delimiter between fields [tab]
  -e, --encrypt         encrypt selected SNPs to --file=path [false]
  -f FILE, --file FILE  encryption file [none]
  -g GENKEYPATH, --genkeypath GENKEYPATH
                        generate private and public key path [none]
  -i, --ids             include sample ids with SNPs [false]
  -k KEYPATH, --keypath KEYPATH
                        RSA private and/or public key file path [none]
  -l LOG, --log LOG     log output file path [none]
	-o OUTFILE, --outfile OUTFILE
												BAM/SAM output file [none]
  -p POS, --pos POS     base reference position comma-separated list [none]
  -r REMOVE, --remove REMOVE
                        remove --pos SNPs to this VCF file path [none]
  -s SNPS, --snps SNPS  number of SNPs to output [0, use -1 for all]
  -t, --verbose         enables verbose output [false]
  -v, --header          output BAM/SAM/VCF header [false]
  -w PASSWORD, --password PASSWORD
                        encrypted private key password [none]
  -x REGION, --region REGION
                        eXtract BAM region(s), e.g. '1:100-105,2:1234:5678' [none]

CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT
['alleles',
 'alts',
 'chrom',
 'contig',
 'copy',
 'filter',
 'format',
 'header',
 'id',
 'info',
 'pos',
 'qual',
 'ref',
 'rid',
 'rlen',
 'samples',
 'start',
 'stop',
 'translate']
"""

import os
import sys

import argparse  # command line parsing library
from pysam import AlignmentFile
from pysam import VariantFile
from cryptography.fernet import Fernet
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.asymmetric import utils
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives import hashes

# accepts a variable number of arguments
def printlog(*args, sep=" "):
	if sep == "tab":
		print(*args, sep='\t')
		if logFile is not None:
			print(*args, file=logFile, sep='\t')
	else:
		print(*args)
		if logFile is not None:
			print(*args, file=logFile)

# returns alternative alleles
def getAlts(alts):
	altsStr = ""
	if alts is not None:
		first = True
		for alt in alts:
			if first == True:
				altsStr += f"{alt}"
				first = False
			else:
				altsStr += f",{alt}"
	return altsStr

# returns filter
def getFilter(filters):
	filterStr = ""
	first = True
	for filter in filters:
		if first == True:
			filterStr += f"{filter}"
			first = False
		else:
			filterStr += f",{filter}"
	return filterStr

# returns info
def getInfo(info):
	infoStr = ""
	first = True
	infoKeys = info.keys()
	for key in infoKeys:
		value = info[key]
		if type(value) is tuple:
			value = value[0]
		if isinstance(value, float):  # remove trailing zeroes
			value = f"{value:.9f}".rstrip('0').rstrip('.')
		if first == True:
			infoStr += f"{key}={value}"
			first = False
		else:
			infoStr += f";{key}={value}"
	return infoStr

# returns format
def getFormat(formats):
	formatStr = ""
	first = True
	for format in formats:
		if first == True:
			formatStr += f"{format}"
			first = False
		else:
			formatStr += f",{format}"
	return formatStr

# returns the SNP for each sample separated by spaces
#
# ids=True outputs the sample ID
# bases=True outputs 4 bases first letter: Adenine, Thymine, Cytosine & Guanine;
#       otherwise the index in the VCF file is output
def getSamples(samples, count=sys.maxsize, ids=False, bases=False, sep="tab"):
	sampleStr = ""
	first = True
	if sep == "tab":
		sep = '\t'
	else:
		sep = ' '
	# sampleKeys = samples.keys()
	cnt = 0
	if count > 0:
		for ss, rec in samples.items():
			# value = list(samples.values()[0])
			if ids:
				value = ss + "="
			else:
				value = ""
			if bases:
				value = value + rec.alleles[0] + "|" + rec.alleles[1]
			elif len(rec.allele_indices) > 1:
				value = value + str(rec.allele_indices[0]) + "|" + str(rec.allele_indices[1])
			elif len(rec.allele_indices) == 1:
				value = value + str(rec.allele_indices[0])

			if first == True:
				sampleStr += f"{value}"
				first = False
			else:
				sampleStr += f"{sep}{value}"
			cnt += 1;
			if cnt >= count:
				break
	return sampleStr


#
# convert a variable list of arguments to bytes and return
#
def makeBytes(*args, sep="tab"):
	if sep == "tab":
		sep = '\t'
	else:
		sep = ' '
	first = True
	rv = bytes(0)
	for arg in args:
		if first:
			rv += bytes(arg, "ascii")
			first = False
		else:
			arg = sep + arg
			rv += bytes(arg, "ascii")
	return rv


#
# reads private and/or public RSA keys from supplied files and returns them
#
def readPrivatePublicKeys(privateFile, publicFile, password):
	private_key = None
	if password != None and os.path.isfile(privateFile):
		if verbose:
			printlog(f"  reading RSA private key from: {privateFile}")
		with open(privateFile, "rb") as key_file:
			private_key = serialization.load_pem_private_key(key_file.read(),
																											password=password)
		if verbose:
			printlog(f"        RSA private key object: {private_key}")

	public_key = None
	if os.path.isfile(publicFile):
		if verbose:
			printlog(f"   reading RSA public key from: {publicFile}")
		with open(publicFile, "rb") as key_file:
			public_key = serialization.load_pem_public_key(key_file.read())
		if verbose:
			printlog(f"         RSA public key object: {public_key}")
	return private_key, public_key


#
# generates and returns private and public RSA keys, after saving them in files
#
def genPrivatePublicKeys(privateFile, publicFile, password):
	if password == None:
		printlog(" Error: must provide private key password with --password option")
		return None, None
	if verbose:
		printlog(f" generating RSA private key to: {privateFile}")
	genPrivateFile = open(privateFile, "wb")  # write bytes
	private_key = rsa.generate_private_key(public_exponent=65537, key_size=4096)
	public_key = private_key.public_key()
	if password != None:  # Encrypted Private Key
		pem = private_key.private_bytes(encoding=serialization.Encoding.PEM,
																		format=serialization.PrivateFormat.PKCS8,
				encryption_algorithm=serialization.BestAvailableEncryption(password))
	else:
		pem = private_key.private_bytes(encoding=serialization.Encoding.PEM,
													format=serialization.PrivateFormat.TraditionalOpenSSL,
													encryption_algorithm=serialization.NoEncryption())
	pemLines = pem.splitlines()
	for line in pemLines:
		genPrivateFile.write(line)
		genPrivateFile.write(b"\n")
	genPrivateFile.close()

	if verbose:
		printlog(f"  generating RSA public key to: {publicFile}")
	genPublicFile = open(publicFile, "wb")  # write bytes
	pemPublic = public_key.public_bytes(encoding=serialization.Encoding.PEM,
												format=serialization.PublicFormat.SubjectPublicKeyInfo)
	pemPublicLines = pemPublic.splitlines()
	for line in pemPublicLines:
		genPublicFile.write(line)
		genPublicFile.write(b"\n")
	genPublicFile.close()
	# test loading of new private and public keys and return them
	return readPrivatePublicKeys(privateFile, publicFile, password)


# encrypt selected SNPs with RSA-protected symmetric key
def encryptSNPs(encryptFilePath, RSAkeyFile, public_key, encryptFile, posTuple,
								vcfInfile):
	SNPcount = 0
	if encryptFilePath != None:  # encrypt selected SNPs
		if verbose:
			printlog(f"   encrypting selected SNPs to: {encryptFilePath + '.vcf.crypt'}")
			printlog(f"     with RSA-protected key in: {encryptFilePath + '.vcf.key'}")
		with RSAkeyFile as ek:
			key = Fernet.generate_key()
			ciphertext = public_key.encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			ek.write(ciphertext)
			f = Fernet(key)
			with encryptFile as ef:
				for pos in posTuple:
					for rec in vcfInfile.fetch(region = f"21:{pos}-{pos}"):
						qual = rec.qual
						if isinstance(qual, float):  # remove trailing zeroes
							qual = f"{qual:.6f}".rstrip('0').rstrip('.')
						if rec.pos == int(pos):
							SNPbytes = makeBytes(rec.chrom, str(rec.pos), rec.id, rec.ref,
																			getAlts(rec.alts), qual,
																			getFilter(rec.filter), getInfo(rec.info),
																			getFormat(rec.format),
																			getSamples(rec.samples))
							ef.write(f.encrypt(SNPbytes))
							ef.write(b"\n")
							SNPcount += 1
				if verbose:
					printlog(f"{SNPcount} selected SNPs encrypted with unique symmetric "
									"key protected by RSA public key")
	return SNPcount


def decryptSNPs(encryptFilePath, RSAkeyFile, private_key, encryptFile, vcfOutfile):
	SNPcount = 0
	if encryptFilePath != None:  # decrypt selected SNPs with RSA-protected symmetric key
		if verbose:
			printlog(f" decrypting selected SNPs from: {encryptFilePath + '.vcf.crypt'}")
			printlog(f"     with RSA-protected key in: {encryptFilePath + '.vcf.key'}")
			printlog(f"                 decrypting to: {encryptFilePath + '.vcf'}")
		with RSAkeyFile as ek:
			ciphertext = ek.read(512)
			key = private_key.decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			if verbose:
				printlog(f"          Fernet symmetric key: {key}")
			f = Fernet(key)
			with encryptFile as ef:
				with vcfOutfile as vcf:
					while True:
						cryptLine = ef.readline()
						if cryptLine == b"":  # end of file reached
							break
						cryptLine = cryptLine[:-1]
						vcf.write(f.decrypt(cryptLine))
						vcf.write(b"\n")
						SNPcount += 1
		if verbose:
			printlog(f"{SNPcount} selected SNPs decrypted with unique symmetric "
							"key protected by RSA public key")
	return SNPcount


def encryptFernet(encryptFilePath, encryptFile, encryptKeys, posTuple, vcfInfile):
	SNPcount = 0
	if encryptFilePath != None:  # encrypt selected SNPs using Fernet
		if verbose:
			printlog(f"                 encrypting to: {encryptFilePath + '.vcf.SNPcrypt'}")
			printlog(f" individual SNP symmetric keys: {encryptFilePath + '.vcf.SNPkeys'}")
		with encryptFile as ef:
			with encryptKeys as ek:
				for pos in posTuple:
					for rec in vcfInfile.fetch(region = f"21:{pos}-{pos}"):
						qual = rec.qual
						if isinstance(qual, float):  # remove trailing zeroes
							qual = f"{qual:.6f}".rstrip('0').rstrip('.')
						if rec.pos == int(pos):
							key = Fernet.generate_key()
							ek.write(key)
							ek.write(b"\n")
							SNPbytes = makeBytes(rec.chrom, str(rec.pos), rec.id, rec.ref,
																getAlts(rec.alts), qual,
																getFilter(rec.filter), getInfo(rec.info),
																getFormat(rec.format),
																getSamples(rec.samples))
							f = Fernet(key)
							ef.write(f.encrypt(SNPbytes))
							ef.write(b"\n")
							SNPcount += 1
		if verbose:
			printlog(f"{SNPcount} selected SNPs encrypted using individual symmetric keys")
	return SNPcount


def decryptFernet(encryptFilePath, encryptFile, encryptKeys, vcfOutfile):
	SNPcount = 0
	if encryptFilePath != None:  # encrypt selected SNPs using Fernet
		if verbose:
			printlog(f"               decrypting from: {encryptFilePath + '.vcf.SNPcrypt'}")
			printlog(f" individual SNP symmetric keys: {encryptFilePath + '.vcf.SNPkeys'}")
			printlog(f"                 decrypting to: {encryptFilePath + '.vcf'}")
		with encryptFile as ef:
			with encryptKeys as ek:
				with vcfOutfile as vcf:
					while True:
						cryptLine = ef.readline()
						if cryptLine == b"":  # end of file reached
						  break
						cryptLine = cryptLine[:-1]
						key = (ek.readline())[:-1]
						if verbose:
							printlog(f"Fernet symmetric key: {key}")
						f = Fernet(key)
						vcf.write(f.decrypt(cryptLine))
						vcf.write(b"\n")
						SNPcount += 1
		if verbose:
			printlog(f"{SNPcount} SNPs decrypted using individual symmetric keys")
	return SNPcount


def removeSNPs(removeFilePath, posTuple, vcfInfile):
	SNPcopyCount = 0
	SNPcount = 0
	if removeFilePath != None:  # remove selected SNPs from VCF file
		removeFilePath = removeFilePath.strip("'")
		removeFilePathExt = f"{removeFilePath + '.noIDs.vcf'}"
		if verbose:
			printlog(f"   removing identity SNPs from: {removeFilePathExt}")
		"""
		[E::bcf_write] Unchecked error (2)
		Traceback (most recent call last):
			File "/mnt/g/My Drive/1-Sterling/1-UCLA/1-MSCS/CM226 - ML in Bioinformatics/project/data/chr21.py", line 508, in main
			  vcf_out.write(rec)
			File "pysam/libcbcf.pyx", line 4400, in pysam.libcbcf.VariantFile.write
			File "pysam/libcbcf.pyx", line 4437, in pysam.libcbcf.VariantFile.write
		OSError: [Errno 0] b'Success'
		vcf_out = VariantFile(removeFilePathExt, 'w', header=vcfInfile.header)
		Segmentation fault
		"""
		# store VCF header first
		vcf_out = VariantFile(removeFilePathExt, 'w', header=vcfInfile.header)
		vcf_out.close()
		vcfOutfile = open(removeFilePathExt, "ab")  # append VCF bytes
		for rec in vcfInfile.fetch():
			skip = False
			for pos in posTuple:
				if rec.pos == int(pos):
					skip = True
					SNPcount += 1
					if verbose:
						printlog(f"          removed identity SNP: {rec.id} "
										f"at pos {rec.pos} from output file")
					posList = list(posTuple)
					posList.remove(pos)  # remove matched identity SNP from posTuple
					posTuple = tuple(posList)
					break
			if skip == False:
				qual = rec.qual
				if isinstance(qual, float):  # remove trailing zeroes
					qual = f"{qual:.6f}".rstrip('0').rstrip('.')
				SNPbytes = makeBytes(rec.chrom, str(rec.pos), rec.id, rec.ref,
															getAlts(rec.alts), qual,
															getFilter(rec.filter), getInfo(rec.info),
															getFormat(rec.format),
															getSamples(rec.samples, count, ids, bases,
															sep=delim), sep=delim)
				vcfOutfile.write(SNPbytes)
				vcfOutfile.write(b"\n")
				SNPcopyCount += 1
		vcfOutfile.close()
	if verbose:
		printlog(f"{SNPcopyCount} non-identity SNPs output to VCF file and "
						f"{SNPcount} identity SNPs removed")
	return SNPcopyCount, SNPcount


def getArgs():
	ap = argparse.ArgumentParser()  # argument parser
	ap.add_argument("-b", "--bases", action="store_true",
									help="include DNA bases with SNPs instead of indices [false]")
	ap.add_argument("-c", "--count", type=int, default=-1,
									help="count of samples to output [all]")
	ap.add_argument("-d", "--delim", type=ascii, default="tab",
									help="delimiter between fields [tab]")
	ap.add_argument("-e", "--encrypt", action="store_true",
									help="encrypt selected SNPs to --file=path [false]")
	ap.add_argument("-f", "--file", type=ascii,
									help="encryption file [none]")
	ap.add_argument("-g", "--genkeypath", type=ascii,
									help="generate private and public key path [none]")
	ap.add_argument("-i", "--ids", action="store_true",
									help="include sample ids with SNPs [false]")
	ap.add_argument("-k", "--keypath", type=ascii,
									help="RSA private and/or public key file path [none]")
	ap.add_argument("-l", "--log", type=argparse.FileType('a'),
									help="log output file path [none]")
	ap.add_argument("-o", "--outfile", type=ascii,
									help="BAM/SAM output file [none]")
	ap.add_argument("-p", "--pos", type=ascii, default="",
									help="base reference position comma-separated list [none]")
	ap.add_argument("-r", "--remove", type=ascii,
									help="remove --pos SNPs to this VCF file path [none]")
	ap.add_argument("-s", "--snps", type=int, default=0,
									help="number of SNPs to output [0, use -1 for all]")
	ap.add_argument("-t", "--verbose", action="store_true",
									help="enables verbose output [false]")
	ap.add_argument("-v", "--header", action="store_true",
									help="output BAM/SAM/VCF header [false]")
	ap.add_argument("-w", "--password", type=ascii,
									help="encrypted private key password [none]")
	ap.add_argument("-x", "--region", type=ascii,
									help="eXtract BAM region(s), e.g. '1:100-105,2:1234:5678' [none]")
	ap.add_argument('inputfile', nargs='?', type=ascii,
									default = "1000-genomes-phase-3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
									help="BAM/SAM or SNPs VCF input file [1000-genomes-phase-3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz]")
	return vars(ap.parse_args())


def main():
	global logFile  # output file for printlog(s) [default: None]
	global verbose	# verbose output if True [default: False]

	args = getArgs()
	inputfile = args['inputfile']
	inputfile = args['inputfile'].strip("'")
	basename = os.path.basename(inputfile)
	split = os.path.splitext(basename)
	infiletype = None
	if len(split) > 1:
		infiletype = split[len(split) - 1]
	logFile = args['log']
	verbose = args['verbose']
	if verbose:
		printlog("                   "
						f"{os.path.basename(__file__)}: verbose mode on")
	bases = args['bases']
	count = args['count']
	if count == -1:
		count = sys.maxsize
	delim = args['delim'].strip("'")
	password=args['password']  # None if not set
	if password != None:
		password = bytes(password.strip("'"), "ascii")
	region = args['region']  # BAM/SAM file region to extract, if any, else None
	if region != None:
		regionTuple = tuple(region.strip("'").split(","))
	genKeyPath = args['genkeypath']  # None if not set
	keyPath = args['keypath']  # None if not set
	encrypt = args['encrypt']  # default: False (decrypt)
	encryptFilePath = args['file']  # None if not set
	if encryptFilePath != None:
		encryptFilePath = encryptFilePath.strip("'")
		if keyPath != None or genKeyPath != None:  # use RSA public keys
			if encrypt:  # symmetric encryption using public key to protect the key
				encryptFile = open(encryptFilePath + ".vcf.crypt", "wb") # write bytes
				RSAkeyFile = open(encryptFilePath + ".vcf.key", "wb")  # write bytes
			else:  # decrypt using private key to get symmetric decryption key
				encryptFile = open(encryptFilePath + ".vcf.crypt", "rb") # read bytes
				RSAkeyFile = open(encryptFilePath + ".vcf.key", "rb")  # read bytes
				vcfOutfile = open(encryptFilePath + ".vcf", "wb")  # write decrypted bytes
		else:  # perform symmetric per-SNP encryption using Fernet
			if encrypt:
				encryptFile = open(encryptFilePath + ".vcf.SNPcrypt", "wb") # write byte
				encryptKeys = open(encryptFilePath + ".vcf.SNPkeys", "wb")  # write byte
			else:
				encryptFile = open(encryptFilePath + ".vcf.SNPcrypt", "rb") # read bytes
				encryptKeys = open(encryptFilePath + ".vcf.SNPkeys", "rb")  # read bytes
				vcfOutfile = open(encryptFilePath + ".vcf", "wb")  # write decrypted bytes

	removeFilePath = args['remove']  # None if not set

	outfiletype = None
	outfile = args['outfile']  # BAM/SAM output file
	if outfile != None:
		outfile = args['outfile'].strip("'")
		basename = os.path.basename(outfile)
		split = os.path.splitext(basename)
		if len(split) > 1:
			outfiletype = split[len(split) - 1]

	if args['pos'] != "''":
		posTuple = tuple(args['pos'].strip("'").split(","))
	else:
		posTuple = ()
	ids = args['ids']
	snpLimit = args['snps']
	if snpLimit < 0:  # output all SNPs
		snpLimit = sys.maxsize

	private_key = None
	public_key = None

	if keyPath != None:  # read private and/or public RSA keys from files
		keyPath = keyPath.strip("'")
		keyPrivateFile = keyPath + ".key.private"
		keyPublicFile = keyPath + ".key.public"
		private_key, public_key = readPrivatePublicKeys(keyPrivateFile,
																										keyPublicFile,
																										password)

	# newly-generated keys are used when --keypath is also given
	if genKeyPath != None:  # generate private and public RSA key pair
		genKeyPath = genKeyPath.strip("'")
		genKeyPrivateFile = genKeyPath + ".key.private"
		genKeyPublicFile = genKeyPath + ".key.public"
		private_key, public_key = genPrivatePublicKeys(genKeyPrivateFile,
																									 genKeyPublicFile,
																									 password)
	# ".gz" is .vcf.gz compressed VCF file
	if infiletype == ".vcf" or infiletype == ".gz" or infiletype == ".bcf":
		try:
			vcfInfile = VariantFile(inputfile, mode='r', threads=4)
		except (FileNotFoundError, ValueError) as e:
			printlog(f"[VariantFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
			quit()
		if args['header']:
			if verbose:
				printlog(f"VCF file '{inputfile}' header:")
			printlog(vcfInfile.header)
		if encryptFilePath != None:  # encrypt or decrypt selected SNPs
			if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
				if encrypt:
					SNPcount = encryptSNPs(encryptFilePath, RSAkeyFile, public_key, encryptFile,
																posTuple, vcfInfile)
				else:  # decrypt is the default
					SNPcount = decryptSNPs(encryptFilePath, RSAkeyFile, private_key, encryptFile,
																vcfOutfile)
			elif encrypt:  # symmetric encryption using Fernet
				SNPcount = encryptFernet(encryptFilePath, encryptFile, encryptKeys, posTuple,
																vcfInfile)
			else:  # symmetric decryption using Fernet
				SNPcount = decryptFernet(encryptFilePath, encryptFile, encryptKeys, vcfOutfile)
		elif removeFilePath != None:
			SNPcopyCount, SNPcount = removeSNPs(removeFilePath, posTuple, vcfInfile)
		else:  # print plaintext SNPs
			SNPcount = 0
			for pos in posTuple:
				for rec in vcfInfile.fetch(region = f"21:{pos}-{pos}"):
					qual = rec.qual
					if isinstance(qual, float):  # remove trailing zeroes
						qual = f"{qual:.6f}".rstrip('0').rstrip('.')
					if rec.pos == int(pos):
						printlog (rec.chrom, rec.pos, rec.id, rec.ref, getAlts(rec.alts), qual,
											getFilter(rec.filter), getInfo(rec.info),
											getFormat(rec.format), getSamples(rec.samples, count, ids,
											bases, sep=delim), sep=delim)
						SNPcount += 1
			if verbose and SNPcount > 0:
				printlog(f"{SNPcount} selected SNPs output in plaintext")

		if snpLimit > 0:
			SNPcount = 0
			try:
				for rec in vcfInfile.fetch():
					try:
						qual = rec.qual
						if isinstance(qual, float):  # remove trailing zeroes
							qual = f"{qual:.6f}".rstrip('0').rstrip('.')
						printlog (rec.chrom, rec.pos, rec.id, rec.ref, getAlts(rec.alts), qual,
											getFilter(rec.filter), getInfo(rec.info), getFormat(rec.format),
											getSamples(rec.samples, count, ids, bases))
						SNPcount += 1
						if SNPcount >= snpLimit:
							if verbose:
								printlog(f"{SNPcount} SNPs output in plaintext")
							break
					except ValueError as e:
						printlog(f"[--snps={args['snps']}] ValueError ignored: {e}.")
			except OSError as o:
				printlog(f"[--snps={args['snps']}] OSError ignored: {o}.")

		vcfInfile.close()

	elif infiletype == ".bam":
		try:
			bamInfile = AlignmentFile(inputfile, mode='rb', threads=4)
			if args['header']:
				if verbose:
					printlog(f"BAM file '{inputfile}' header:")
				printlog(bamInfile.header)
			if region is None:
				bamIter = bamInfile.fetch()
				for x in bamIter:
					print(str(x))
					quit()
			else:  # support a list of ranges
				for reg in regionTuple:
					bamIter = bamInfile.fetch(region = reg)
					for x in bamIter:
						print(str(x))
				quit()
		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
			quit()
	elif infiletype == ".sam":
		try:
			samInfile = AlignmentFile(inputfile, mode='r', check_sq=False, threads=4)
			if args['header']:
				if verbose:
					printlog(f"SAM file '{inputfile}' header:")
				printlog(samInfile.header)
			samIter = samInfile.fetch()  # cannot extract regions from non-indexed SAM files
			for x in samIter:
				print(str(x))
			quit()
		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
			quit()
	elif infiletype == ".cram":
		try:
			vcfInfile = AlignmentFile(inputfile, mode='rc', threads=4)
		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
			quit()
	else:
		printlog(f"[file: '{inputfile}'] unrecognized file type.")
		printlog("program exiting.")

if __name__ == "__main__":
	main()
