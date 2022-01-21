#!/bin/sh

source ./env.sh

fdir=/willerslev/scratch/grg/Island_Thrush_cleaned_Genomes
export ref=$fdir/P12064_102.fasta

do_stats() {
	sm=$1
	bam=${sm}/MAPPED/${sm}.*.calmd.bam

	samtools flagstat $bam > ${bam}.flagstat.txt \
		|| die "$sm: samtools flagstat"

	samtools stats -r $ref $bam > ${bam}.stats.txt \
		|| die "$sm: samtools stats"

	samtools idxstats $bam > ${bam}.idxstats.txt \
		|| die "$sm: samtools idxstats"
}
export -f do_stats

#for sm in Turdus_*; do
#for sm in Cataponera_turdoides; do
for sm in Turdus_merula; do
	echo do_stats $sm
done \
 | parallel --nice 10 --progress
