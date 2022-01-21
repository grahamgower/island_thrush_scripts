#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
bamlist1=$idir/island_thrush.bamlist.no-canecsens-no-outgroup
bamlist2=$idir/island_thrush.bamlist
bamlist=./bamlist.abbababa
export ref=$idir/P12064_102.fasta
export regions=$idir/autosome-angsd.txt
export angsd=$HOME/angsd/angsd

cat $bamlist1 > $bamlist
grep merula $bamlist2 >> $bamlist

export minInd=$(($(wc -l $bamlist | cut -f1 -d" ") / 2))

$angsd \
	-doAbbababa2 1 \
	-bam $bamlist \
	-rf $regions \
	-rmTrans 1 \
	-useLast 1 \
	-doCounts 1 \
	-enhance 1 \
	-minInd $minInd \
	-minQ 20 \
	-minMapQ 30 \
	\
	-ref $ref \
	-baq 2 \
	-C 50 \
	\
	-out island_thrush_abbababa \
	|| die "angsd bailed"

