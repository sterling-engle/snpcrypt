# snpcrypt
                     UNIVERSITY OF CALIFORNIA LOS ANGELES
                        DEPARTMENT OF COMPUTER SCIENCE

Securing Masked Short Reads, Identification SNPs, Phenotypes, and Clinical Data
using Cryptography

Capstone Project for: MSCS - Faculty Advisor: Prof. Sriram Sankararaman
                         Date: February 19, 2022

Abstract---Many people are understandably concerned about protecting their genomic privacy.
For this reason, they are unwilling to donate their genomes for genomic research. 
Huang et al discovered 116 identification (ID) SNPs on the human genome that can
be used to uniquely identify every individual. My snpcrypt Python program extracts 
and encrypts these ID SNPs, and can return a subset of the sample genotypes. It 
encrypts masked extracted short reads of raw genomic data. It can encrypt any file 
containing phenotypes or clinical data. All files are encrypted and decrypted using 
a symmetric key that is protected by asymmetric RSA public key cryptography. Multiple 
RSA public keys are supported for all file types for situations where more than one 
individual needs to access the same file. Encryption enables this extremely sensitive 
data to be securely stored in a genomic databank and securely transmitted to researchers 
and clinicians. I validated ID SNP extraction, optional sample selection, and encryption 
using a 2,504-sample chromosome 21 variant call format (VCF) file from phase 3 of the 
Human Genome Project. I tested masked genomic sequence short read extraction and 
encryption with an indexed binary alignment map (BAM) file of one of the phase 1 human 
DNA samples. I verified error-free parsing of SAM and VCF versions 4.1, 4.2, and 4.3 
files using scripts to run snpcrypt against all the "passed" test files in the Samtools 
hts-specs repository. Synpcrypt.py source code is available in my public code repository:
https://github.com/sterling-engle/snpcrypt.
