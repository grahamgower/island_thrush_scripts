#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
export bamlist=$idir/island_thrush.bamlist
export ref=$idir/P12064_102.fasta
export regions=$idir/autosome.txt


bcftools mpileup \
	-b $bamlist \
	-f $ref \
	-R $regions \
	-a AD \
	-Ou \
	--threads 8 \
	-d 100 \
| bcftools call \
	--threads 8 \
	-m \
	-v \
	-O b \
	-o island_thrush_all.bcf \
	|| die "bcftools died"

