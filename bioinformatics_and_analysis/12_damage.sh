#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
bamlist=$idir/island_thrush.bamlist
export ref=$idir/P12064_102.fasta
export odir=condamage
export cddir=$HOME/condamage

mkdir -p $odir

do_damage() {
	bam=$1
	pfx=$(basename ${bam%%.*.bam})
	$cddir/condamage $bam $ref > $odir/${pfx}.txt
	$cddir/plot_condamage.py -o $odir/${pfx}.pdf $odir/${pfx}.txt
}
export -f do_damage

for bam in `cat $bamlist`; do
	echo do_damage $bam
done | parallel --nice 10
