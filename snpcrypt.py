# -*- coding: utf-8 -*-
"""
@author: Sterling Engle
@uid: 904341227

Developed and tested only on Ubuntu 18.04.6 LTS under Python 3.9.7.

usage: snpcrypt.py [-h] [-b] [-c COUNT] [--delim DELIM] [-d] [-e] [-f FILE]
                   [-g GENKEYPATH] [-i] [-k KEYPATH] [-l LOG] [-m] [-n SNPS]
                   [-o OUTFILE] [-p POS] [-r REMOVE] [-s SAMPLES] [-t] [-v]
                   [-w PASSWORD] [-x REGION]
                   [inputfile]

positional arguments:
  inputfile             BAM/SAM, SNPs VCF/BCF, or other input file [1000-genomes-phase
                        -3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_inte
                        grated_v5a.20130502.genotypes.vcf.gz]

optional arguments:
  -h, --help            show this help message and exit
  -b, --bases           include DNA bases with SNPs instead of indices [false]
  -c COUNT, --count COUNT
                        count of samples to output [all]
  --delim DELIM         delimiter between fields [tab]
  -d, --decrypt         decrypt selected SNPs to --file=path.vcf or extracted
                        BAM regions to inputfile.bam from inputfile.bam.crypt
                        [false]
  -e, --encrypt         encrypt selected SNPs to --file=path.vcf.crypt or
                        extracted BAM region to --outfile=path.crypt [false]
  -f FILE, --file FILE  encryption file [none]
  -g GENKEYPATH, --genkeypath GENKEYPATH
                        generate private and public key path list [none]
  -i, --ids             include sample ids with SNPs [false]
  -k KEYPATH, --keypath KEYPATH
                        RSA private and/or public key file path list [none]
  -l LOG, --log LOG     log output file path [none]
  -m, --mask            mask BAM file short reads [false]
  -n SNPS, --snps SNPS  number of SNPs to output [0, use -1 for all]
  -o OUTFILE, --outfile OUTFILE
                        BAM/SAM output file [none]
  -p POS, --pos POS     base reference position(s) comma-separated list [none]
  -r REMOVE, --remove REMOVE
                        remove --pos SNPs to this VCF or BCF file path [none]
  -s SAMPLES, --samples SAMPLES
                        VCF/BCF sample(s) comma-separated list [all]
  -t, --verbose         enables verbose output [false]
  -v, --header          output BAM/SAM/VCF header [false]
  -w PASSWORD, --password PASSWORD
                        encrypted private key password list [none]
  -x REGION, --region REGION
                        eXtract BAM region(s), e.g. '1:100-105,2:1234:5678'
                        [none]

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
from pysam import AlignmentFile  # reads BAM and SAM files
from pysam import VariantFile  # reads VCF and BCF files
from cryptography.fernet import Fernet  # symmetric encryption
from cryptography.hazmat.primitives.asymmetric import rsa  # RSA asymmetric encryption
from cryptography.hazmat.primitives.asymmetric import padding
from cryptography.hazmat.primitives.asymmetric import utils
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives import hashes

# accepts a variable number of arguments and prints to sys.stdout and logFile if open
def printlogstdout(*args, sep=" ", end="\n"):
	if sep == "tab":
		print(*args, sep='\t', end=end)
		if logFile is not None:
			print(*args, file=logFile, sep='\t', end=end)
	else:
		print(*args, end=end)
		if logFile is not None:
			print(*args, file=logFile, end=end)

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
	if len(infoKeys) == 0:
		return "."  # indicates an empty field
	for key in infoKeys:
		try:
			keyvalue = info[key]
		except (ValueError) as e:
			printlog(f"[info['{key}']] ValueError ignored: {e}.")
			continue
		# printlog(f"{key}={keyvalue}; type is {type(keyvalue)}")
		if type(keyvalue) is tuple:
			value = ""
			firstval = True
			for val in keyvalue:
				if isinstance(val, float):  # remove trailing zeroes
					val = f"{val:.7f}".rstrip('0').rstrip('.')
				if firstval:
					firstval = False
					value += f"{val}"
				else:
					value += f",{val}"
		elif isinstance(keyvalue, float):  # remove trailing zeroes
			value = f"{keyvalue:.7f}".rstrip('0').rstrip('.')
		elif isinstance(keyvalue, bool):	# change True to "1" and False to "0"
			if keyvalue:
				value = "1"
			else:
				value = "0"
		else:
			value = keyvalue

		if first == True:
			if value == None:
				infoStr += f"{key}"
			else:
				infoStr += f"{key}={value}"
			first = False
		else:
			if value == None:
				infoStr += f";{key}"
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
			formatStr += f":{format}"
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
	cnt = 0
	if count > 0:
		for ss, rec in samples.items():
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
			skip = True
			for val in rec.itervalues():  # add any additional fields besides GT (genotype)
				if skip:
					skip = False
				else:
					if type(val) is str:
						value = value + ":" + val
					elif type(val) is tuple:
						firstval = True
						for v in val:
							if firstval:
								value = value + ":" + v
								firstval = False
							else:
								value = value + "," + v

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
				printlog(f"error: RSA private key file {privateFile} password '{password}': {e}")
				quit()
		if verbose and trace:
			printlog(f"        RSA private key object: {private_key}")

	public_key = None
	if os.path.isfile(publicFile):
		if verbose:
			printlog(f"   reading RSA public key from: {publicFile}")
		with open(publicFile, "rb") as key_file:
			public_key = serialization.load_pem_public_key(key_file.read())
		if verbose and trace:
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
def encryptSNPs(cryptFilePath, keyUsers, RSAkeyFiles, public_keys, cryptFile, posTuple,
								vcfInfile):
	SNPcount = 0
	if cryptFilePath != None:  # encrypt selected SNPs
		if verbose:
			printlog(f"   encrypting selected SNPs to: {cryptFilePath + '.vcf.crypt'}")
			for i in range(len(keyUsers)):
				printlog("     with RSA-protected key in: "
								f"{cryptFilePath + '.vcf.' + keyUsers[i] + '.key'}")
		key = Fernet.generate_key()
		for i in range(len(RSAkeyFiles)):
			with RSAkeyFiles[i] as ek:
				ciphertext = public_keys[i].encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
				ek.write(ciphertext)
		f = Fernet(key)
		with cryptFile as ef:
			ef.write(f.encrypt(bytes(str(vcfInfile.header), "ascii")))  # encrypted VCF header
			ef.write(b"\n")
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

def decryptSNPs(cryptFilePath, keyUsers, RSAkeyFiles, private_keys, cryptFile, vcfOutfile):
	if private_keys[0] == None:
		printlog("error: --decrypt (-d) option requires RSA private key --password "
						"to be supplied.")
		quit()
	SNPcount = -1  # do not count the header
	if cryptFilePath != None:  # decrypt selected SNPs with RSA-protected symmetric key
		if verbose:
			printlog(f" decrypting selected SNPs from: {cryptFilePath + '.vcf.crypt'}")
			printlog("     with RSA-protected key in: "
								f"{cryptFilePath + '.vcf.' + keyUsers[0] + '.key'}")
			printlog(f"                 decrypting to: {cryptFilePath + '.vcf'}")
		with RSAkeyFiles[0] as ek:
			ciphertext = ek.read(512)
			key = private_keys[0].decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			if verbose and trace:
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
						if SNPcount != -1:  # header already has a newline at the end
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
				key = Fernet.generate_key()
				ek.write(key)
				ek.write(b"\n")
				f = Fernet(key)
				ef.write(f.encrypt(bytes(str(vcfInfile.header), "ascii")))  # encrypted VCF header
				ef.write(b"\n")
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
	SNPcount = -1  # do not count the header
	if cryptFilePath != None:  # encrypt selected SNPs using Fernet without RSA
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
						if verbose and trace:
							printlog(f"Fernet symmetric key: {key}")
						f = Fernet(key)
						vcf.write(f.decrypt(cryptLine))
						if SNPcount != -1:  # header already has a newline at the end
							vcf.write(b"\n")
						SNPcount += 1
		if verbose:
			printlog(f"{SNPcount} SNPs decrypted using individual symmetric keys")
	return SNPcount

def removeSNPs(inputfile, removeFilePath, posTuple, vcfInfile):
	SNPcopyCount = 0
	SNPcount = 0
	if removeFilePath != None:  # remove selected SNPs from vcfInfile to this BCF or VCF file
		removeFilePath = removeFilePath.strip("'")
		if verbose:
			printlog(f"   removing identity SNPs from: {inputfile}")
			printlog(f"                            to: {removeFilePath}")
		"""
		This line fixes pysam error [E::vcf_format] in vcfOutfile.write() call
		[W::vcf_parse_format] FORMAT 'GT' at 21:9411239 is not defined in the header, assuming Type=String
		[E::vcf_format] Invalid BCF, the FORMAT tag id=28 at 21:9411239 not present in the header
		"""
		vcfInfile.header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
		bcf_or_vcf_out = VariantFile(removeFilePath, 'w', header=vcfInfile.header)
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
				bcf_or_vcf_out.write(rec)
				SNPcopyCount += 1
		bcf_or_vcf_out.close()
		vcfInfile.close()
	if verbose:
		printlog(f"{SNPcopyCount} non-identity SNPs output to {removeFilePath} and "
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
			samline = []  # list of BAM file input lines in SAM format: one per base in region
			queryNames = []
			readsList = []
			refNames = []
			refPos = []
			numAlignedBases = []
			queryPos = []
			qualScores = []
			sequences = []
			for reg in regionTuple:
				# perform pileup within the region
				bamIter = bamInfile.pileup(region = reg, truncate = True)
				for x in bamIter:  # x is of type pysam.PileupColumn
					queryNames.append(x.get_query_names())  # list of query names in region
					queryPos.append(x.get_query_positions())  # list of sequence positions in region
					sequences.append(x.get_query_sequences())  # list of sequences in the region
					if verbose and trace:
						readsList.append(x.nsegments)
						refNames.append(x.reference_name)
						refPos.append(x.reference_pos + 1)  # note + 1
						numAlignedBases.append(x.get_num_aligned())
						qualScores.append(x.get_mapping_qualities())
					for p in x.pileups:  # list of reads (pysam.PileupRead) aligned to this column
						# p.alignment is a pysam.AlignedSegment object of the aligned read
						# to_string() is a string representation of aligned segment in valid SAM format
						samline.append(p.alignment.to_string())
			samlines = list(dict.fromkeys(samline))  # removes duplicate SAM lines
			samitems = []
			for s in samlines:
				if verbose and trace:
					printlog(s)
				samitems.append(s.split('\t'))  # create a list of the SAM fields
			if verbose and trace:
				for s in samitems:
					printlog(s)
				printlog(f"query names: {queryNames}")
				printlog(f"number of reads: {readsList}")
				printlog(f"reference name: {refNames}")
				printlog(f"reference position: {refPos}")
				printlog(f"number of aligned bases: {numAlignedBases}")
				printlog(f"positions in read: {queryPos}")
				# printlog(f"quality scores: {qualScores}")
				printlog(f"query sequences: {sequences}")
			for s in range(len(samitems)):  # replace sequence with N's to mask all bases
				samitems[s][9] = "N" * len(samitems[s][9])
			for q in range(len(queryNames)):  # for each query name
				for i in range(len(queryNames[q])):  # for each matched sequence in that query
					for s in range(len(samitems)):  # look for matching query name
						if samitems[s][0] == queryNames[q][i]:  # does this line contain the query?
							# first index is 0-based position
							# before matched base + matched base + after matched base
							samitems[s][9] = samitems[s][9][:queryPos[q][i]] + \
																sequences[q][i].upper() + \
																samitems[s][9][queryPos[q][i] + 1:]  # insert matched base
			if verbose and trace:
				for s in range(len(samitems)):
					printlog(samitems[s])
			for s in samitems:
				a = pysam.AlignedSegment()  # get empty aligned segment for writing masked result
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
				tags = []  # for storing optional fields
				for t in range(11, len(s)):
					tagpieces = s[t].split(":")  # split subfields, e.g. AM:i:23
					tagtype = tagpieces[1]  # tag types are: A, i, f, Z, H, B
					if tagtype == "i":  # convert string to signed integer
						tagpieces[1] = int(tagpieces[2])
					elif tagtype == "f":  # convert string to single-precision floating point number
						tagpieces[1] = float(tagpieces[2])
					else:
						tagpieces[1] = tagpieces[2]
					tagpieces[2] = tagtype
					tags.append(tagpieces)  # expects list of tag, value, type
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
	pysam.sort("-@", "4", "-o", outfile, tempfile)  # sort alignments by leftmost coordinates
	os.remove(tempfile)  # remove temporary file
	if writemode == "wb":  #  index BAM file
		if verbose:
			printlog(f"indexing {outfile}")
		pysam.index("-@", "4", outfile)
	return

# encrypt bamFile selected short read regions with RSA-protected symmetric key
#
# bamFile is the name of the *.bam file to be encrypted, e.g. "test.bam"
# keyUsers contains the --keypath list, e.g. --keypath=user2,user3 ['user2', 'user3']
# RSAkeyFiles contains a list of those keyUsers opened bamFile.<user>.key files for
#   storing the RSA public key encrypted symmetric key that decrypts the encrypted bamFile.
#   For example: opened files ['test.bam.user2.key', 'test.bam.user3.key']
# public_keys contains a list of the keyUsers public keys
# cryptFile contains the opened file for writing the symmetric encrypted bamFile.
#   For example, opened file 'test.bam.crypt'.
def encryptRegions(bamFile, keyUsers, RSAkeyFiles, public_keys, cryptFile):
	if bamFile != None:
		if verbose:
			printlog(f"encrypting extracted regions to: {bamFile + '.crypt'}")
			for i in range(len(keyUsers)):
				printlog("     with unique symmetric key protected by RSA public key in: "
								f"{bamFile + '.' + keyUsers[i] + '.key'}")
		key = Fernet.generate_key()  # Fernet generates a random 32-byte symmetric key
		for i in range(len(RSAkeyFiles)):  # encrypt key with each public RSA key to its file
			with RSAkeyFiles[i] as ek:
				ciphertext = public_keys[i].encrypt(key,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))  # returns RSA public key encrypted symmetric key
				ek.write(ciphertext)  # save it in corresponding .key file, e.g. 'test.bam.user2.key'
		f = Fernet(key)  # get Fernet object for encryption using key
		with open(bamFile, "rb") as bf:  # open the .bam file for reading binary values
			bamData = bf.read()  # read the entire file to bamData
			with cryptFile as ef:  # cryptFile is encryption output file, e.g. 'test.bam.crypt'
				ef.write(f.encrypt(bamData))  # encrypt the entire file, saving it to cryptFile
		os.remove(bamFile)  # remove the .bam plaintext file
		if verbose:
			printlog(f"{bamFile} extracted short reads encrypted to {bamFile + '.crypt'}")
			for i in range(len(keyUsers)):
				printlog("     with unique symmetric key protected by RSA public key in "
								f"{bamFile + '.' + keyUsers[i] + '.key'}")
	return

# decrypt bamFile selected short read regions with RSA-protected symmetric key
#
# bamFile is the name of the *.bam file to be decrypted, e.g. "test.bam"
# keyUsers contains the --keypath list, e.g. --keypath=user2 ['user2']
# RSAkeyFiles contains a list of that keyUsers opened bamFile.<user>.key file which
#   stores the RSA public key encrypted symmetric key that decrypts the encrypted bamFile.
#   For example: opened file ['test.bam.user2.key']
# private_keys contains a list of the keyUsers private keys
# cryptFile contains the opened file for reading the symmetric encrypted bamFile.
#   For example, opened file 'test.bam.crypt'.
def decryptRegions(bamFile, keyUsers, RSAkeyFiles, private_keys, cryptFile):
	if private_keys[0] == None:
		printlog("error: --decrypt (-d) option requires RSA private key --password "
						"to be supplied.")
		quit()
	if bamFile != None:
		if verbose:
			printlog(f"decrypting extracted regions to: {bamFile} from {bamFile + '.crypt'}")
			printlog("     using unique symmetric key protected by RSA public key in: "
							f"{bamFile + '.' + keyUsers[0] + '.key'}")
		with RSAkeyFiles[0] as ek:  # opened file.bam.keyid.key or file.sam.keyid.key file
			ciphertext = ek.read(512)  # read RSA public key encrypted symmetric key
			key = private_keys[0].decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))  # returns RSA private key decrypted symmetric key
			if verbose and trace:
				printlog(f"           Fernet symmetric key: {key}")
		f = Fernet(key)  # get Fernet object for decryption using key
		with cryptFile as ef:
			cryptBamData = ef.read()
			with open(bamFile, "wb") as bf:
				bf.write(f.decrypt(cryptBamData))
		if verbose:
			printlog(f"{bamFile} extracted short reads decrypted from {bamFile + '.crypt'}")
			printlog("     using unique symmetric key protected by RSA public key in "
							f"{bamFile + '.' + keyUsers[0] + '.key'}")
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
			if verbose and trace:
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
def encryptFile(inputfile, keyUsers, RSAkeyFiles, public_keys, cryptFile):
	if inputfile != None:
		if verbose:
			printlog(f"encrypting {inputfile} to: {inputfile + '.crypt'}")
			for i in range(len(keyUsers)):
				printlog("     with unique symmetric key protected by RSA public key in: "
								f"{inputfile + '.' + keyUsers[i] + '.key'}")
		key = Fernet.generate_key()
		for i in range(len(RSAkeyFiles)):
			with RSAkeyFiles[i] as ek:
				ciphertext = public_keys[i].encrypt(key,
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
			for i in range(len(keyUsers)):
				printlog("     with unique symmetric key protected by RSA public key in "
								f"{inputfile + '.' + keyUsers[i] + '.key'}")
	return

# decrypt inputfile with RSA-protected symmetric key
def decryptFile(outputfile, keyUsers, RSAkeyFiles, private_keys, cryptFile):
	if private_keys[0] == None:
		printlog("error: --decrypt (-d) option requires RSA private key --password "
						"to be supplied.")
		quit()
	if outputfile != None:
		if verbose:
			printlog(f"             decrypting to: {outputfile} from {outputfile + '.crypt'}")
			printlog(f"using RSA-protected key in: {outputfile +  '.' + keyUsers[0] + '.key'}")
		with RSAkeyFiles[0] as ek:
			ciphertext = ek.read(512)
			key = private_keys[0].decrypt(
											ciphertext,
											padding.OAEP(mgf=padding.MGF1(algorithm=hashes.SHA256()),
											algorithm=hashes.SHA256(),
											label=None))
			if verbose and trace:
				printlog(f"           Fernet symmetric key: {key}")
			f = Fernet(key)
			with cryptFile as ef:
				cryptData = ef.read()
				with open(outputfile, "wb") as outf:
					outf.write(f.decrypt(cryptData))
			if verbose:
				printlog(f"{outputfile} decrypted from {outputfile + '.crypt'}")
				printlog("  using unique symmetric key protected by RSA public key in: "
								f"{outputfile +  '.' + keyUsers[0] + '.key'}")
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

# decrypt inputfile with symmetric key
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
									help="generate private and public key path list [none]")
	ap.add_argument("-i", "--ids", action="store_true",
									help="include sample ids with SNPs [false]")
	ap.add_argument("-k", "--keypath", type=ascii,
									help="RSA private and/or public key file path list [none]")
	ap.add_argument("-l", "--log", type=argparse.FileType('a'),
									help="log output file path [none]")
	ap.add_argument("-m", "--mask", action="store_true",
									help="mask BAM file short reads [false]")
	ap.add_argument("-n", "--snps", type=int, default=0,
									help="number of SNPs to output [0, use -1 for all]")
	ap.add_argument("-o", "--outfile", type=ascii,
									help="BAM/SAM output file [none]")
	ap.add_argument("-p", "--pos", type=ascii, default="",
									help="base reference position(s) comma-separated list [none]")
	ap.add_argument("-r", "--remove", type=ascii,
									help="remove --pos SNPs to this VCF or BCF file path [none]")
	ap.add_argument("-s", "--samples", type=ascii, default="",
									help="VCF/BCF sample(s) comma-separated list [all]")
	ap.add_argument("-t", "--verbose", action="store_true",
									help="enables verbose output [false]")
	ap.add_argument("-v", "--header", action="store_true",
									help="output BAM/SAM/VCF header [false]")
	ap.add_argument("-w", "--password", type=ascii,
									help="encrypted private key password list [none]")
	ap.add_argument("-x", "--region", type=ascii,
									help="eXtract BAM region(s), e.g. '1:100-105,2:1234:5678' [none]")
	ap.add_argument('inputfile', nargs='?', type=ascii,
									default = "1000-genomes-phase-3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
									help="BAM/SAM, SNPs VCF/BCF, or other input file [1000-genomes-phase-3_vcf-20150220_ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz]")
	return vars(ap.parse_args())

def main():
	global logFile  # output file for printlog(s) [default: None]
	global verbose	# verbose output if True [default: False]
	global trace		# detailed trace output if True

	trace = False
	args = getArgs()
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
	region = args['region']  # BAM/SAM file region to extract, if any, else None
	if region != None:
		regionTuple = tuple(region.strip("'").split(","))
	# support a list of RSA keypaths and passwords for generation and
	# loading public and/or private keys
	password=args['password']  # None if not set
	genKeyPath = args['genkeypath']  # None if not set
	if genKeyPath != None:
		if password == None:
			printlog("error: -g/--genkeypath argument also requires -w/--password argument")
			quit()
		genKeyTuple = tuple(genKeyPath.strip("'").split(","))
		keyUsers = list()  # used to generate .<keyUser>.key file name
		for i in range(len(genKeyTuple)):
			keyUsers.append(os.path.basename(genKeyTuple[i]))
	keyPath = args['keypath']  # None if not set
	if keyPath != None:
		if genKeyPath != None:
			printlog("error: only one of -g/--genkeypath or -k/--keypath arguments "
							"may be given at a time")
			quit()
		keyTuple = tuple(keyPath.strip("'").split(","))
		keyUsers = list()  # used to generate .<keyUser>.key file name
		for i in range(len(keyTuple)):
			keyUsers.append(os.path.basename(keyTuple[i]))
	if password != None:
		passwordList = list(password.strip("'").split(","))
		for i in range(len(passwordList)):
			passwordList[i] = bytes(passwordList[i], "ascii")
		passwordTuple = tuple(passwordList)
		if genKeyPath != None and len(genKeyTuple) != len(passwordTuple):
			printlog("error: -g/--genkeypath and -w/--password list size are not the same:")
			printlog(f"       genkeypath(s): {genKeyTuple} and password(s): {passwordTuple}")
			quit()
		elif keyPath != None and len(keyTuple) != len(passwordTuple):
			printlog("error: -k/--keypath and -w/--password list size are not the same:")
			printlog(f"       keypath(s): {keyTuple} and password(s): {passwordTuple}")
			quit()
		elif genKeyPath == None and keyPath == None:
			printlog("error: -w/--password argument also requires either -k/--keypath "
							"or -g/--genkeypath")
			quit()
	decrypt = args['decrypt']  # default: False
	encrypt = args['encrypt']  # default: False
	if decrypt and encrypt:
		printlog("error: cannot use --decrypt and --encrypt options on the same command line")
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
	if args['samples'] != "''":
		samplesTuple = tuple(args['samples'].strip("'").split(","))
	else:
		samplesTuple = None
	snpLimit = args['snps']
	if snpLimit < 0:  # output all SNPs
		snpLimit = sys.maxsize

	private_keys = list()
	public_keys = list()

	if keyPath != None:  # read private and/or public RSA keys from files
		for i in range(len(keyTuple)):
			keyPrivateFile = keyTuple[i] + ".key.private"
			keyPublicFile = keyTuple[i] + ".key.public"
			if "passwordTuple" in locals():
				priv_key, pub_key = readPrivatePublicKeys(keyPrivateFile,
																									keyPublicFile,
																									passwordTuple[i])
			else:
				priv_key, pub_key = readPrivatePublicKeys(keyPrivateFile,
																									keyPublicFile,
																									None)
			private_keys.append(priv_key)  # None if not retrieved
			public_keys.append(pub_key)  # None if not retrieved

	# newly-generated keys are used when -g/--genkeypath is given
	if genKeyPath != None:  # generate private and public RSA key pair
		for i in range(len(genKeyTuple)):
			genKeyPrivateFile = genKeyTuple[i] + ".key.private"
			genKeyPublicFile = genKeyTuple[i] + ".key.public"
			priv_key, pub_key = genPrivatePublicKeys(genKeyPrivateFile,
																							 genKeyPublicFile,
																							 passwordTuple[i])
			private_keys.append(priv_key)  # None if not generated
			public_keys.append(pub_key)  # None if not generated

	cryptFilePath = args['file']  # None if not set
	if cryptFilePath != None:
		cryptFilePath = cryptFilePath.strip("'")
		if keyPath != None or genKeyPath != None:  # use RSA public keys
			RSAkeyFiles = list()
			if encrypt:  # symmetric encryption using public key to protect the key
				cryptFile = open(cryptFilePath + ".vcf.crypt", "wb") # write bytes
				for i in range(len(keyUsers)):
					RSAkeyFiles.append(open(cryptFilePath + ".vcf." + keyUsers[i] + ".key", "wb"))
			elif decrypt:  # decrypt using private key to get symmetric decryption key
				cryptFile = open(cryptFilePath + ".vcf.crypt", "rb") # read bytes
				for i in range(len(keyUsers)):
					RSAkeyFiles.append(open(cryptFilePath + ".vcf." + keyUsers[i] + ".key", "rb"))
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

	# ".gz" is .vcf.gz compressed VCF file
	if infiletype == ".vcf" or infiletype == ".gz" or infiletype == ".bcf":
		try:
			vcfInfile = VariantFile(inputfile, mode='r', threads=4)
		except (FileNotFoundError, ValueError) as e:
			printlog(f"error: [VariantFile({inputfile})]: {e}.")
			quit()
		if samplesTuple != None:  # read only a subset of samples
			printlog(f"    reading only these samples: {samplesTuple}")
			try:
				vcfInfile.subset_samples(samplesTuple)
			except (ValueError) as e:
				printlog(f"error: --samples: {e}.")
				quit()
		if args['header']:
			if verbose:
				printlog(f"VCF file '{inputfile}' header:")
			printlogstdout(vcfInfile.header, end = "")
		if cryptFilePath != None:  # encrypt or decrypt selected SNPs
			if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
				if encrypt:
					SNPcount = encryptSNPs(cryptFilePath, keyUsers, RSAkeyFiles, public_keys,
																cryptFile, posTuple, vcfInfile)
				elif decrypt:
					SNPcount = decryptSNPs(cryptFilePath, keyUsers, RSAkeyFiles, private_keys,
																cryptFile, vcfOutfile)
			elif encrypt:  # symmetric encryption using Fernet
				SNPcount = encryptFernet(cryptFilePath, cryptFile, encryptKeys, posTuple,
																vcfInfile)
			elif decrypt:  # symmetric decryption using Fernet
				SNPcount = decryptFernet(cryptFilePath, cryptFile, encryptKeys, vcfOutfile)
		elif removeFilePath != None:
			SNPcopyCount, SNPcount = removeSNPs(inputfile, removeFilePath, posTuple, vcfInfile)
		else:  # print plaintext SNPs
			SNPcount = 0
			for pos in posTuple:
				for rec in vcfInfile.fetch(region = f"21:{pos}-{pos}"):
					qual = rec.qual
					if isinstance(qual, float):  # remove trailing zeroes
						qual = f"{qual:.6f}".rstrip('0').rstrip('.')
					if rec.pos == int(pos):
						printlogstdout (rec.chrom, rec.pos, rec.id, rec.ref, getAlts(rec.alts),
														qual, getFilter(rec.filter), getInfo(rec.info),
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
						printlogstdout (rec.chrom, rec.pos, rec.id, rec.ref, getAlts(rec.alts),
														qual, getFilter(rec.filter), getInfo(rec.info),
														getFormat(rec.format), getSamples(rec.samples, count, ids,
														bases, sep=delim), sep=delim)
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
			RSAkeyFiles = list()
			for i in range(len(keyUsers)):
				RSAkeyFiles.append(open(inputfile + "." + keyUsers[i] + ".key", "rb"))
			cryptFile = open(inputfile + ".crypt", "rb")  # read bytes
			decryptRegions(inputfile, keyUsers, RSAkeyFiles, private_keys, cryptFile)
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
				if outfiletype != None:
					with AlignmentFile(tempfile, mode=writemode, header=bamInfile.header) as outf:
						for x in bamIter:
							outf.write(x)
					if verbose:
						printlog(f"sorting {outfile}")
					pysam.sort("-@", "4", "-o", outfile, tempfile)  # sort alignments by leftmost coord.
					os.remove(tempfile)  # remove temporary file
					if writemode == "wb":  #  index BAM file
						if verbose:
							printlog(f"indexing {outfile}")
						pysam.index("-@", "4", outfile)
				else:
					for x in bamIter:
						print(str(x))
			else:  # support a list of ranges
				if writemode != None:
					if verbose:
						printlog(f"extracting region(s): {regionTuple} from {inputfile} to {outfile}")
					extractRegions(mask, bamInfile, regionTuple, tempfile, writemode, outfile)
					if encrypt:
						if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
							RSAkeyFiles = list()
							for i in range(len(keyUsers)):
								RSAkeyFiles.append(open(outfile + "." + keyUsers[i] + ".key", "wb"))
							cryptFile = open(outfile + ".crypt", "wb")  # write bytes
							encryptRegions(outfile, keyUsers, RSAkeyFiles, public_keys, cryptFile)
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
			printlog(f"error: [AlignmentFile({inputfile})]: {e}.")
	elif infiletype == ".sam":
		writemode = None
		tempfile = None
		if outfiletype == ".bam":
			writemode = "wb"
			tempfile = ".temp.bam"
		elif outfiletype == ".sam":
			writemode = "w"
			tempfile = ".temp.sam"
		try:
			samInfile = AlignmentFile(inputfile, mode='r', check_sq=False, threads=4)
			if args['header']:
				if verbose:
					printlog(f"SAM file '{inputfile}' header:")
				printlog(samInfile.header)
			samIter = samInfile.fetch()  # cannot extract regions from non-indexed SAM files
			if outfiletype != None:
				with AlignmentFile(tempfile, mode=writemode, header=samInfile.header) as outf:
					for x in samIter:
						outf.write(x)
				if verbose:
					printlog(f"sorting {outfile}")
				pysam.sort("-@", "4", "-o", outfile, tempfile)  # sort alignments by leftmost coord.
				os.remove(tempfile)  # remove temporary file
				if writemode == "wb":  #  index BAM file
					if verbose:
						printlog(f"indexing {outfile}")
					pysam.index("-@", "4", outfile)
			else:
				for x in samIter:
					print(str(x))

		except (FileNotFoundError, ValueError) as e:
			printlog(f"error: [AlignmentFile({inputfile})] error ignored: {e}.")
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
			printlog(f"error: [AlignmentFile({inputfile})]: {e}.")

	elif encrypt:  # encrypt all other file types when requested
		if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
			RSAkeyFiles = list()
			for i in range(len(keyUsers)):
				RSAkeyFiles.append(open(inputfile + "." + keyUsers[i] + ".key", "wb"))
			cryptFile = open(inputfile + ".crypt", "wb")  # write bytes
			encryptFile(inputfile, keyUsers, RSAkeyFiles, public_keys, cryptFile)
		else:  # symmetric encryption using Fernet
			FernetKeyFile = open(inputfile + ".symkey", "wb")  # write bytes
			cryptFile = open(inputfile + ".symcrypt", "wb")  # write bytes
			encryptFileFernet(inputfile, FernetKeyFile, cryptFile)

	elif decrypt:  # decrypt all other file types when requested
		if keyPath != None or genKeyPath != None:  # use RSA-protected symmetric key
			RSAkeyFiles = list()
			for i in range(len(keyUsers)):
				RSAkeyFiles.append(open(inputfile + "." + keyUsers[i] + ".key", "rb"))
			cryptFile = open(inputfile + ".crypt", "rb")  # read bytes
			decryptFile(inputfile, keyUsers, RSAkeyFiles, private_keys, cryptFile)
		else:  # symmetric decryption using Fernet without RSA-protected symmetric key
			FernetKeyFile = open(inputfile + ".symkey", "rb")  # read bytes
			cryptFile = open(inputfile + ".symcrypt", "rb")  # read bytes
			decryptFileFernet(inputfile, FernetKeyFile, cryptFile)

	else:
		printlog(f"error: [file: '{inputfile}'] unrecognized file type.")

if __name__ == "__main__":
	main()
