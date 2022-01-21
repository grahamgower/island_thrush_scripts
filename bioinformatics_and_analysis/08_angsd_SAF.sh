#!/bin/sh

# http://www.popgen.dk/angsd/index.php/SFS_Estimation
# http://www.popgen.dk/angsd/index.php/Fst

source ./env.sh

export angsd=$HOME/angsd/angsd

idir=/willerslev/datasets/Island_Thrush/
bamlist=$idir/island_thrush.bamlist
export ref=$idir/P12064_102.fasta
export fai=${ref}.fai

export regions=$idir/autosome-angsd.txt
# Sites that are polymorphic across the species, excluding transitions. (maf2regions.sh)
export sites=/willerslev/scratch/grg/Island_Thrush_cleaned_Genomes/island_thrush_dp10.snppos.tab

do_saf() {
	sm=$1
	bam=$2
	dir=$3

	mkdir -p $dir
	cd $dir

	$angsd \
		-i $bam \
		-ref $ref \
		-anc $ref \
		-out $sm \
		-rf $regions \
		-sites $sites \
		\
		-minMapQ 30 \
		-minQ 20 \
		-baq 2 \
		-C 50 \
		-uniqueOnly 1 \
		-noTrans 1 \
		\
		-doSaf 1 \
		-GL 1 \
	|| die "$sm: angsd -doSaf 1"
}
export -f do_saf

for bam in `cat $bamlist`; do
	sm=$(basename ${bam%%.*bam})
	dir=$sm/ANGSD/
	echo do_saf $sm $bam $dir
done | parallel -j 72 --nice 10
