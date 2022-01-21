#!/bin/bash
#
# vcf_test.sh runs snpcrypt with every *.vcf VCF file in the current working directory
#
# @author: Sterling Engle
# @uid: 904341227
#
for f in *.vcf ; do
  echo "testing file:" $f ;
  python ~/capstone/snpcrypt.py --snps -1 $f ;
done
