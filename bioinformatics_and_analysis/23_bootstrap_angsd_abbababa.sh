#!/bin/sh

export NUMBA_NUM_THREADS=40
./bootstrap_angsd_abbababa.py \
	P12064_102.fasta.fai \
	autosome.txt \
	abbababa1.indlist \
	island_thrush_abbababa1.abbababa \
	> island_thrush_abbababa1.bootstrap2
