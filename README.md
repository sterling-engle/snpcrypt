# snpcrypt
Securing Identification SNPs via Encryption with snpcrypt

Many people are understandably concerned about protecting their genomic privacy.
For this reason, they are unwilling to donate their genomes for genomic research. 
Previous efforts discovered 116 identification (ID) SNPs on the human genome 
that can be used to uniquely identify every individual.

I created the snpcrypt Python program to extract and encrypt these ID SNPs with
a symmetric key that is protected by RSA public key cryptography. This enables
this extremely sensitive data to be securely stored in a genomic databank and
securely transmitted to researchers. I validated snpcrypt using a 2,504-sample
chromosome 21 variant call format (VCF) file from phase 3 of the Human Genome
Project, measured its performance, and discuss ideas for future work in my paper.
