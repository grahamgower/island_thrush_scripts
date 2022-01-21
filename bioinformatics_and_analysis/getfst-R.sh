#!/bin/sh

for f in `find ANGSD_sfs -name \*.Reynolds.stats`; do
	spair=$(basename ${f%.stats})
	s1=$(echo $spair | cut -f1 -d\.)
	s2=$(echo $spair | cut -f2 -d\.)
	echo -en "$s1 $s2 "
	cat $f
done
