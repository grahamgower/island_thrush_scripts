#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
bamlist=$idir/island_thrush.bamlist.dp10
export ref=$idir/P12064_102.fasta
export regions=$idir/autosome-angsd.txt
export angsd=$HOME/angsd/angsd

export minInd=$(($(wc -l $bamlist | cut -f1 -d" ") / 2))

$angsd \
	-bam $bamlist \
	-rf $regions \
	-ref $ref \
	-anc $ref \
	-nThreads 24 \
	\
	-minMapQ 30 \
	-minQ 20 \
	-baq 2 \
	-C 50 \
	-uniqueOnly 1 \
	-minInd $minInd \
	-noTrans 1 \
	\
	-GL 1 \
	-doMaf 1 \
	-doGlf 2 \
	-doMajorMinor 1 \
	-doSaf 1 \
	-SNP_pval 1e-6 \
	-out island_thrush_dp10 \
	|| die "angsd bailed"

