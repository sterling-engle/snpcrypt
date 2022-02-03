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
  --delim DELIM					delimiter between fields [tab]
  -e, --encrypt         encrypt selected SNPs to --file=path [false]
  -f FILE, --file FILE  encryption file [none]
  -g GENKEYPATH, --genkeypath GENKEYPATH
                        generate private and public key path [none]
  -i, --ids             include sample ids with SNPs [false]
  -k KEYPATH, --keypath KEYPATH
                        RSA private and/or public key file path [none]
  -l LOG, --log LOG     log output file path [none]
	-m, --mask						mask BAM file short reads [false]
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
import pysam  # python lightweight wrapper of the htslib C-API
from pysam import AlignmentFile  # reads BAM/SAM files
from pysam import VariantFile  # reads VCF files
from cryptography.fernet import Fernet  # symmetric encryption
from cryptography.hazmat.primitives.asymmetric import rsa  # RSA asymmetric encryption
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.asymmetric import utils
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives import hashes

# accepts a variable number of arguments and prints to sys.stderr and logFile if open
def printlog(*args, sep=" "):
	if sep == "tab":
		print(*args, file=sys.stderr, sep='\t')
		if logFile is not None:
			print(*args, file=logFile, sep='\t')
	else:
		print(*args, file=sys.stderr)
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
			try:
				private_key = serialization.load_pem_private_key(key_file.read(),
																												password=password)
			except (ValueError) as e:
				printlog(f"RSA private key file {privateFile} password '{password}': {e}")
				printlog("program exiting.")
				quit()
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
def encryptSNPs(cryptFilePath, RSAkeyFile, public_key, cryptFile, posTuple,
								vcfInfile):
	SNPcount = 0
	if cryptFilePath != None:  # encrypt selected SNPs
		if verbose:
			printlog(f"   encrypting selected SNPs to: {cryptFilePath + '.vcf.crypt'}")
			printlog(f"     with RSA-protected key in: {cryptFilePath + '.vcf.key'}")
		with RSAkeyFile as ek:
			key = Fernet.generate_key()
			ciphertext = public_key.encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			ek.write(ciphertext)
			f = Fernet(key)
			with cryptFile as ef:
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

def decryptSNPs(cryptFilePath, RSAkeyFile, private_key, cryptFile, vcfOutfile):
	SNPcount = 0
	if cryptFilePath != None:  # decrypt selected SNPs with RSA-protected symmetric key
		if verbose:
			printlog(f" decrypting selected SNPs from: {cryptFilePath + '.vcf.crypt'}")
			printlog(f"     with RSA-protected key in: {cryptFilePath + '.vcf.key'}")
			printlog(f"                 decrypting to: {cryptFilePath + '.vcf'}")
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
			with cryptFile as ef:
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

def encryptFernet(cryptFilePath, cryptFile, encryptKeys, posTuple, vcfInfile):
	SNPcount = 0
	if cryptFilePath != None:  # encrypt selected SNPs using Fernet
		if verbose:
			printlog(f"                 encrypting to: {cryptFilePath + '.vcf.SNPcrypt'}")
			printlog(f" individual SNP symmetric keys: {cryptFilePath + '.vcf.SNPkeys'}")
		with cryptFile as ef:
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

def decryptFernet(cryptFilePath, cryptFile, encryptKeys, vcfOutfile):
	SNPcount = 0
	if cryptFilePath != None:  # encrypt selected SNPs using Fernet
		if verbose:
			printlog(f"               decrypting from: {cryptFilePath + '.vcf.SNPcrypt'}")
			printlog(f" individual SNP symmetric keys: {cryptFilePath + '.vcf.SNPkeys'}")
			printlog(f"                 decrypting to: {cryptFilePath + '.vcf'}")
		with cryptFile as ef:
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

#
# extract regions (masking if mask is True) in regionTuple from bamInfile to tempfile
# with writemode, sorting the output to outfile and indexing it if it is a BAM file
#
def extractRegions(mask, bamInfile, regionTuple, tempfile, writemode, outfile):
	with AlignmentFile(tempfile, mode=writemode, header=bamInfile.header) as outf:
		if mask:  # mask short reads?
			if verbose:
				printlog("masking short reads with N's")
			samline = []
			queryNames = []
			readsList = []
			refNames = []
			refPos = []
			numAlignedBases = []
			queryPos = []
			qualScores = []
			sequences = []
			for reg in regionTuple:
				bamIter = bamInfile.pileup(region = reg, truncate = True)
				for x in bamIter:  # x is of type pysam.PileupColumn
					queryNames.append(x.get_query_names())
					queryPos.append(x.get_query_positions())
					sequences.append(x.get_query_sequences())
					if verbose and trace:
						readsList.append(x.nsegments)
						refNames.append(x.reference_name)
						refPos.append(x.reference_pos + 1)  # note + 1
						numAlignedBases.append(x.get_num_aligned())
						qualScores.append(x.get_mapping_qualities())
					for p in x.pileups:
						samline.append(p.alignment.to_string())
			samlines = list(dict.fromkeys(samline))
			samitems = []
			for s in samlines:
				if verbose and trace:
					printlog(s)
				samitems.append(s.split('\t'))
			if verbose and trace:
				for s in samitems:
					printlog(s)
				printlog(f"query names: {queryNames}")
				printlog(f"number of reads: {readsList}")
				printlog(f"reference name: {refNames}")
				printlog(f"reference position: {refPos}")
				printlog(f"number of aligned bases: {numAlignedBases}")
				printlog(f"positions in read: {queryPos}")
				printlog(f"quality scores: {qualScores}")
				printlog(f"query sequences: {sequences}")
			for s in range(len(samitems)):  # replace sequence with N's to mask it
				samitems[s][9] = "N" * len(samitems[s][9])
			for q in range(len(queryNames)):
				for i in range(len(queryNames[q])):
					for s in range(len(samitems)):
						if samitems[s][0] == queryNames[q][i]:
							# first index is 0-based position
							samitems[s][9] = samitems[s][9][:queryPos[q][i]] + sequences[q][i].upper()  + samitems[s][9][queryPos[q][i] + 1:]
			if verbose and trace:
				for s in range(len(samitems)):
					printlog(samitems[s])
			for s in samitems:
				a = pysam.AlignedSegment()
				# could use fromstring(type cls, sam, AlignmentHeader header)
				a.query_name = s[0]
				a.flag = int(s[1])
				a.reference_id = int(s[2]) - 1  # subtract 1 for 0 origin
				a.reference_start = int(s[3]) - 1  # subtract 1 for 0 origin
				a.mapping_quality = int(s[4])
				a.cigarstring = s[5]
				a.next_reference_name = s[6]
				a.next_reference_start = int(s[7]) - 1  # subtract 1 for 0 origin
				a.template_length = int(s[8])
				a.query_sequence = s[9]
				a.query_qualities = pysam.qualitystring_to_array(s[10])
				tags = []
				for t in range(11, len(s)):
					tagpieces = s[t].split(":")
					tagtype = tagpieces[1]
					if tagtype == "i":
						tagpieces[1] = int(tagpieces[2])
					elif tagtype == "f":
						tagpieces[1] = float(tagpieces[2])
					else:
						tagpieces[1] = tagpieces[2]
					tagpieces[2] = tagtype
					tags.append(tagpieces)
				a.set_tags(tags)
				outf.write(a)
		else:
			for reg in regionTuple:
				bamIter = bamInfile.fetch(region = reg)
				# bamIter = bamInfile.pileup(region = reg)
				for x in bamIter:
					outf.write(x)
		outf.close()
	if verbose:
		printlog(f"sorting {outfile}")
	pysam.sort("-@", "4", "-o", outfile, tempfile)
	os.remove(tempfile)
	if writemode == "wb":  #  index BAM file
		if verbose:
			printlog(f"indexing {outfile}")
		pysam.index("-@", "4", outfile)
	return

# encrypt bamFile selected short read regions with RSA-protected symmetric key
def encryptRegions(bamFile, RSAkeyFile, public_key, cryptFile):
	if bamFile != None:
		if verbose:
			printlog(f"encrypting extracted regions to: {bamFile + '.crypt'}")
			printlog(f"      with RSA-protected key in: {bamFile + '.key'}")
		with RSAkeyFile as ek:
			key = Fernet.generate_key()
			ciphertext = public_key.encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			ek.write(ciphertext)
			f = Fernet(key)
			with open(bamFile, "rb") as bf:
				bamData = bf.read()
				with cryptFile as ef:
					ef.write(f.encrypt(bamData))
			os.remove(bamFile)
			if verbose:
				printlog(f"{bamFile} extracted short reads encrypted to {bamFile + '.crypt'}")
				printlog("  with unique symmetric key protected by RSA public key in "
								f"{bamFile + '.key'}")
	return

# decrypt bamFile selected short read regions with RSA-protected symmetric key
def decryptRegions(bamFile, RSAkeyFile, private_key, cryptFile):
	if private_key == None:
		printlog("error: --decrypt (-d) option requires RSA private key --password "
						"to be supplied.")
		quit()
	if bamFile != None:
		if verbose:
			printlog(f"decrypting extracted regions to: {bamFile} from {bamFile + '.crypt'}")
			printlog(f"     using RSA-protected key in: {bamFile + '.key'}")
		with RSAkeyFile as ek:
			ciphertext = ek.read(512)
			key = private_key.decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			if verbose:
				printlog(f"           Fernet symmetric key: {key}")
			f = Fernet(key)
			with cryptFile as ef:
				cryptBamData = ef.read()
				with open(bamFile, "wb") as bf:
					bf.write(f.decrypt(cryptBamData))
			if verbose:
				printlog(f"{bamFile} extracted short reads decrypted from {bamFile + '.crypt'}")
				printlog("  with unique symmetric key protected by RSA public key in "
								f"{bamFile + '.key'}")
	return

# encrypt bamFile selected short read regions with symmetric key
def encryptRegionsFernet(bamFile, FernetKeyFile, cryptFile):
	if bamFile != None:
		if verbose:
			printlog(f"encrypting extracted regions to: {bamFile + '.symcrypt'}")
			printlog(f"          with symmetric key in: {bamFile + '.symkey'}")
		with FernetKeyFile as fkf:
			key = Fernet.generate_key()
			fkf.write(key)
			fkf.write(b"\n")
			f = Fernet(key)
			with open(bamFile, "rb") as bf:
				bamData = bf.read()
				with cryptFile as ef:
					ef.write(f.encrypt(bamData))
			os.remove(bamFile)
			if verbose:
				printlog(f"{bamFile} extracted short reads encrypted to {bamFile + '.symcrypt'}")
				printlog(f"  with symmetric key in {bamFile + '.symkey'}")
	return

# decrypt bamFile selected short read regions with symmetric key
def decryptRegionsFernet(bamFile, FernetKeyFile, cryptFile):
	if bamFile != None:
		if verbose:
			printlog(f"    decrypting from: {bamFile + '.symcrypt'}")
			printlog(f"using symmetric key: {bamFile + '.symkey'}")
			printlog(f"      decrypting to: {bamFile}")
		with FernetKeyFile as fkf:
			key = (fkf.readline())[:-1]
			if verbose:
				printlog(f"Fernet symmetric key: {key}")
			f = Fernet(key)
			with cryptFile as ef:
				bamData = ef.read()
				with open(bamFile, "wb") as bf:
					bf.write(f.decrypt(bamData))
		if verbose:
			printlog(f"{bamFile + '.symcrypt'} extracted short reads decrypted to {bamFile}")
			printlog(f"  with symmetric key in {bamFile + '.symkey'}")
	return

# encrypt inputfile with RSA-protected symmetric key using public key
def encryptFile(inputfile, RSAkeyFile, public_key, cryptFile):
	if inputfile != None:
		if verbose:
			printlog(f"encrypting {inputfile} to: {inputfile + '.crypt'}")
			printlog(f"  with RSA-protected key in: {inputfile + '.key'}")
		with RSAkeyFile as ek:
			key = Fernet.generate_key()
			ciphertext = public_key.encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			ek.write(ciphertext)
			f = Fernet(key)
			with open(inputfile, "rb") as inf:
				fileData = inf.read()
				with cryptFile as ef:
					ef.write(f.encrypt(fileData))
			os.remove(inputfile)
			if verbose:
				printlog(f"{inputfile} encrypted to {inputfile + '.crypt'} then deleted")
				printlog("  with unique symmetric key protected by RSA public key in "
								f"{inputfile + '.key'}")
	return

# decrypt inputfile with RSA-protected symmetric key
def decryptFile(outputfile, RSAkeyFile, private_key, cryptFile):
	if private_key == None:
		printlog("error: --decrypt (-d) option requires RSA private key --password "
						"to be supplied.")
		quit()
	if outputfile != None:
		if verbose:
			printlog(f"             decrypting to: {outputfile} from {outputfile + '.crypt'}")
			printlog(f"using RSA-protected key in: {outputfile + '.key'}")
		with RSAkeyFile as ek:
			ciphertext = ek.read(512)
			key = private_key.decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			if verbose:
				printlog(f"           Fernet symmetric key: {key}")
			f = Fernet(key)
			with cryptFile as ef:
				cryptData = ef.read()
				with open(outputfile, "wb") as outf:
					outf.write(f.decrypt(cryptData))
			if verbose:
				printlog(f"{outputfile} decrypted from {outputfile + '.crypt'}")
				printlog("  with unique symmetric key protected by RSA public key in "
								f"{outputfile + '.key'}")
	return

# encrypt inputfile with symmetric key
def encryptFileFernet(inputfile, FernetKeyFile, cryptFile):
	if inputfile != None:
		if verbose:
			printlog(f"encrypting {inputfile} to: {inputfile + '.symcrypt'}")
			printlog(f"  with symmetric key in: {inputfile + '.symkey'}")
		with FernetKeyFile as fkf:
			key = Fernet.generate_key()
			fkf.write(key)
			fkf.write(b"\n")
			f = Fernet(key)
			with open(inputfile, "rb") as inf:
				fileData = inf.read()
				with cryptFile as ef:
					ef.write(f.encrypt(fileData))
			os.remove(inputfile)
			if verbose:
				printlog(f"{inputfile} encrypted to {inputfile + '.symcrypt'} then deleted")
				printlog(f"  with symmetric key in {inputfile + '.symkey'}")
	return

# decrypt bamFile selected short read regions with symmetric key
def decryptFileFernet(outputfile, FernetKeyFile, cryptFile):
	if outputfile != None:
		if verbose:
			printlog(f"    decrypting from: {outputfile + '.symcrypt'}")
			printlog(f"using symmetric key: {outputfile + '.symkey'}")
			printlog(f"      decrypting to: {outputfile}")
		with FernetKeyFile as fkf:
			key = (fkf.readline())[:-1]
			if verbose:
				printlog(f"Fernet symmetric key: {key}")
			f = Fernet(key)
			with cryptFile as ef:
				cryptedData = ef.read()
				with open(outputfile, "wb") as bf:
					bf.write(f.decrypt(cryptedData))
		if verbose:
			printlog(f"{outputfile + '.symcrypt'} decrypted to {outputfile}")
			printlog(f"  with symmetric key in {outputfile + '.symkey'}")
	return

#
# initialize and parse the command line arguments
#
def getArgs():
	ap = argparse.ArgumentParser()  # argument parser
	ap.add_argument("-b", "--bases", action="store_true",
									help="include DNA bases with SNPs instead of indices [false]")
	ap.add_argument("-c", "--count", type=int, default=-1,
									help="count of samples to output [all]")
	ap.add_argument("--delim", type=ascii, default="tab",
									help="delimiter between fields [tab]")
	ap.add_argument("-d", "--decrypt", action="store_true",
									help="decrypt selected SNPs to --file=path.vcf or "
									"extracted BAM regions to inputfile.bam from inputfile.bam.crypt [false]")
	ap.add_argument("-e", "--encrypt", action="store_true",
									help="encrypt selected SNPs to --file=path.vcf.crypt or "
									"extracted BAM region to --outfile=path.crypt [false]")
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
	ap.add_argument("-m", "--mask", action="store_true",
									help="mask BAM file short reads [false]")
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
	global trace		# detailed trace output if True

	trace = False
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
	decrypt = args['decrypt']  # default: False
	encrypt = args['encrypt']  # default: False
	if decrypt and encrypt:
		printlog("error: cannot use --decrypt and --encrypt options on the same command line")
		quit()
	cryptFilePath = args['file']  # None if not set
	if cryptFilePath != None:
		cryptFilePath = cryptFilePath.strip("'")
		if keyPath != None or genKeyPath != None:  # use RSA public keys
			if encrypt:  # symmetric encryption using public key to protect the key
				cryptFile = open(cryptFilePath + ".vcf.crypt", "wb") # write bytes
				RSAkeyFile = open(cryptFilePath + ".vcf.key", "wb")  # write bytes
			elif decrypt:  # decrypt using private key to get symmetric decryption key
				cryptFile = open(cryptFilePath + ".vcf.crypt", "rb") # read bytes
				RSAkeyFile = open(cryptFilePath + ".vcf.key", "rb")  # read bytes
				vcfOutfile = open(cryptFilePath + ".vcf", "wb")  # write decrypted bytes
			else:
				printlog("error: --file (-f) option requires either --decrypt (-d) "
								"or --encrypt (-e)")
				quit()
		else:  # perform symmetric per-SNP encryption using Fernet
			if encrypt:
				cryptFile = open(cryptFilePath + ".vcf.SNPcrypt", "wb") # write byte
				encryptKeys = open(cryptFilePath + ".vcf.SNPkeys", "wb")  # write byte
			elif decrypt:
				cryptFile = open(cryptFilePath + ".vcf.SNPcrypt", "rb") # read bytes
				encryptKeys = open(cryptFilePath + ".vcf.SNPkeys", "rb")  # read bytes
				vcfOutfile = open(cryptFilePath + ".vcf", "wb")  # write decrypted bytes
			else:
				printlog("error: --file (-f) option requires either --decrypt (-d) "
								"or --encrypt (-e)")
				quit()

	removeFilePath = args['remove']  # None if not set

	mask = args['mask']  # mask short reads when extracting from indexed BAM file
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
			print(vcfInfile.header)
		if cryptFilePath != None:  # encrypt or decrypt selected SNPs
			if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
				if encrypt:
					SNPcount = encryptSNPs(cryptFilePath, RSAkeyFile, public_key, cryptFile,
																posTuple, vcfInfile)
				elif decrypt:
					SNPcount = decryptSNPs(cryptFilePath, RSAkeyFile, private_key, cryptFile,
																vcfOutfile)
			elif encrypt:  # symmetric encryption using Fernet
				SNPcount = encryptFernet(cryptFilePath, cryptFile, encryptKeys, posTuple,
																vcfInfile)
			elif decrypt:  # symmetric decryption using Fernet
				SNPcount = decryptFernet(cryptFilePath, cryptFile, encryptKeys, vcfOutfile)
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

	elif infiletype == ".bam" and decrypt:
		if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
			RSAkeyFile = open(inputfile + ".key", "rb")  # read bytes
			cryptFile = open(inputfile + ".crypt", "rb")  # read bytes
			decryptRegions(inputfile, RSAkeyFile, private_key, cryptFile)
		else:  # symmetric decryption using Fernet without RSA-protected symmetric key
			FernetKeyFile = open(inputfile + ".symkey", "rb")  # read bytes
			cryptFile = open(inputfile + ".symcrypt", "rb")  # read bytes
			decryptRegionsFernet(inputfile, FernetKeyFile, cryptFile)

	elif infiletype == ".bam":
		writemode = None
		tempfile = None
		if outfiletype == ".bam":
			writemode = "wb"
			tempfile = ".temp.bam"
		elif outfiletype == ".sam":
			writemode = "w"
			tempfile = ".temp.sam"
		try:
			bamInfile = AlignmentFile(inputfile, mode='rb', threads=4)
			if args['header']:
				if verbose:
					printlog(f"BAM file '{inputfile}' header:")
				print(bamInfile.header)
			if region is None:
				bamIter = bamInfile.fetch()
				for x in bamIter:
					print(str(x))
			else:  # support a list of ranges
				if writemode != None:
					if verbose:
						printlog(f"extracting region(s): {regionTuple} from {inputfile} to {outfile}")
					extractRegions(mask, bamInfile, regionTuple, tempfile, writemode, outfile)
					if encrypt:
						if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
							RSAkeyFile = open(outfile + ".key", "wb")  # write bytes
							cryptFile = open(outfile + ".crypt", "wb")  # write bytes
							encryptRegions(outfile, RSAkeyFile, public_key, cryptFile)
						else:  # symmetric encryption using Fernet
							FernetKeyFile = open(outfile + ".symkey", "wb")  # write bytes
							cryptFile = open(outfile + ".symcrypt", "wb")  # write bytes
							encryptRegionsFernet(outfile, FernetKeyFile, cryptFile)
				else:
					if verbose:
						printlog(f"extracting unmasked region: {regionTuple} from {inputfile}",
										"to standard output")
					for reg in regionTuple:
						bamIter = bamInfile.fetch(region = reg)
						for x in bamIter:
							print(str(x))

		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
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

		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")
	elif infiletype == ".cram":
		try:
			cramInfile = AlignmentFile(inputfile, mode='rc', threads=4)
			if args['header']:
				if verbose:
					printlog(f"CRAM file '{inputfile}' header:")
				printlog(cramInfile.header)
			cramIter = cramInfile.fetch()
			for x in cramIter:
				print(str(x))
		except (FileNotFoundError, ValueError) as e:
			printlog(f"[AlignmentFile({inputfile})] error ignored: {e}.")
			printlog("program exiting.")

	elif encrypt:  # encrypt all other file types when requested
		if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
			RSAkeyFile = open(inputfile + ".key", "wb")  # write bytes
			cryptFile = open(inputfile + ".crypt", "wb")  # write bytes
			encryptFile(inputfile, RSAkeyFile, public_key, cryptFile)
		else:  # symmetric encryption using Fernet
			FernetKeyFile = open(inputfile + ".symkey", "wb")  # write bytes
			cryptFile = open(inputfile + ".symcrypt", "wb")  # write bytes
			encryptFileFernet(inputfile, FernetKeyFile, cryptFile)

	elif decrypt:  # decrypt all other file types when requested
		if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
			RSAkeyFile = open(inputfile + ".key", "rb")  # read bytes
			cryptFile = open(inputfile + ".crypt", "rb")  # read bytes
			decryptFile(inputfile, RSAkeyFile, private_key, cryptFile)
		else:  # symmetric decryption using Fernet without RSA-protected symmetric key
			FernetKeyFile = open(inputfile + ".symkey", "rb")  # read bytes
			cryptFile = open(inputfile + ".symcrypt", "rb")  # read bytes
			decryptFileFernet(inputfile, FernetKeyFile, cryptFile)

	else:
		printlog(f"[file: '{inputfile}'] unrecognized file type.")
		printlog("program exiting.")

if __name__ == "__main__":
	main()
