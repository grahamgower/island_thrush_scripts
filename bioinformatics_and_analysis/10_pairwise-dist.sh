#!/bin/sh

ngsDist=/home/srx907/ngsDist/ngsDist
pfx=island_thrush_subset

n_lines=$(zcat ${pfx}.beagle.gz | wc -l | cut -f1 -d" ")
n_sites=$((n_lines-1))

bamlist=/willerslev/datasets/Island_Thrush/island_thrush_subset.bamlist
cut -f5 -d\/ $bamlist > sample.labels
n_ind=$(wc -l $bamlist | cut -f1 -d" ")


$ngsDist \
	-verbose 1 \
	-geno ${pfx}.beagle.gz \
	-probs \
	-n_ind ${n_ind} \
	-n_sites ${n_sites} \
	-labels sample.labels \
	--n_boot_rep 100 \
	--boot_block_size 1500 \
	-o ${pfx}.dist \
	-n_threads 24
