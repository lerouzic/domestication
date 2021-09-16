#!/usr/bin/env bash

# Merge panels in a sigle PDF (horizontal) when
# called figXXA.pdf, figXXB.pdf, etc. 

# Requires pdfjam

for i in fig*A.pdf
do
	BN=`basename $i A.pdf`
	MYFILES=$BN?.pdf
	MYFILES_ARR=( $MYFILES )
	MYFILES_N=${#MYFILES_ARR[@]}
	pdfjam $MYFILES --nup ${MYFILES_N}x1 --landscape --outfile $BN.pdf
done
