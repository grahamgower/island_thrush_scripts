#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
export regions=$idir/autosome.txt

bcf=island_thrush_all_SnpGap_shorter_names.bcf

mkdir -p eig

for s in `cat $regions`; do
	echo ~/eig-utils/vcf2eig \
		-r $s \
		-F SnpGap \
		-j \
		-t \
		-l \
		-a 1000000 \
		-o eig/island_thrush_all_filter.${s} \
		$bcf
done | parallel -j 64

for s in `cat $regions`; do
	cat eig/island_thrush_all_filter.${s}.snp
done > island_thrush_all_filter.snp

for s in `cat $regions`; do
	cat eig/island_thrush_all_filter.${s}.geno
done > island_thrush_all_filter.geno

s=`head -n1 $regions`
ln eig/island_thrush_all_filter.${s}.ind island_thrush_all_filter.ind
