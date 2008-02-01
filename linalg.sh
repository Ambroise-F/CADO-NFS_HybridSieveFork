#!/bin/sh -
linalg=linalg

# parameters names from merge_linalg_sqrt.sh:
# $1 is $name.small
# $2 is $name.ker_raw
# $3 is $root (only used for bw)

mat=$1
ker=$2
root=$3

if false ; then
   echo "Calling Gauss"
   time $linalg/linalg $mat 1 > $ker
else
   # use Block-Wiedemann (recommended for large matrices)
   echo "Transposing the matrix"
   time $linalg/transpose $mat $mat.tr

   echo "Calling Block-Wiedemann"
   if [ ! -e $root.bw ] ; then
      mkdir $root.bw
   fi
   time $linalg/bw/doit.pl matrix=$mat.tr mn=64 vectoring=64 multisols=1 wdir=$root.bw

   echo "Converting dependencies to CADO format"
   # doit.pl puts the dependency file W in the directory where the matrix was
   d=`dirname $root`
   time $linalg/bw/mkbitstrings $d/W > $ker
fi
