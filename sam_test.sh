#!/bin/bash
#
# sam_test.sh runs snpcrypt with every *.sam SAM file in the current working directory
#
# @author: Sterling Engle
# @uid: 904341227
#
for f in *.sam ; do
  echo "testing file:" $f ;
  python ~/capstone/snpcrypt.py $f ;
done
