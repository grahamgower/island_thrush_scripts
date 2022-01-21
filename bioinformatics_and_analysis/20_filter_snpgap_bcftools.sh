#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
export bamlist=$idir/island_thrush.bamlist
export ref=$idir/P12064_102.fasta
export regions=$idir/autosome.txt

bcf_in=island_thrush_all.bcf
bcf_out=island_thrush_all_SnpGap.bcf

bcftools filter \
	--SnpGap 10 \
	-O b \
	-o $bcf_out \
	$bcf_in
