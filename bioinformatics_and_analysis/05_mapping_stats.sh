#!/bin/sh

source ./env.sh

#idir=/willerslev/datasets/Island_Thrush/
idir=.
fdir=/willerslev/scratch/grg/Island_Thrush_cleaned_Genomes
export ref=$fdir/P12064_102.fasta

do_stats() {
	bam=$1
	sm=$(basename ${bam%.rmdup.realigned.calmd.MQ30.bam})

	samtools flagstat $bam > ${bam}.flagstat.txt \
		|| die "$sm: samtools flagstat"

	samtools stats -r $ref $bam > ${bam}.stats.txt \
		|| die "$sm: samtools stats"

	samtools idxstats $bam > ${bam}.idxstats.txt \
		|| die "$sm: samtools idxstats"
}
export -f do_stats

#for dir in $idir/Turdus_*; do
#for dir in $idir/Cataponera_turdoides; do
for dir in $idir/Turdus_merula; do
	sm=$(basename $dir)
	in=$dir/MAPPED/${sm}.rmdup.realigned.calmd.MQ30.bam
	echo do_stats $in
done | parallel --nice 10
