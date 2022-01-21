#!/bin/sh

source ./env.sh

idir=.

do_het() {
	ml=$1
	stats=$2
	awk '{ print $2/($1+$3) }' $ml
}

for dir in $idir/{Cataponera_turdoides,Turdus_*}; do
	sm=$(basename $dir)
	ml=$dir/ANGSD/${sm}.saf.ml
	echo -en "${sm}\t"
	do_het $ml
done
