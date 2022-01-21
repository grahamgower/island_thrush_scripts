#!/bin/sh

mt_size=16733 
genome_size=1067231262

do_collate() {
	size=$1
	shift
	printf "individual\tCoverage\tDepth\tFraglen\tGC\n"
	for f in $@; do
		sm=$(basename ${f%%.*.txt})
		printf "$sm\t"
		awk -v genome_size=$size \
		'
			/^SN	average length:/ {mean_fraglen=$4}
			/^COV/ {dp+=$3*$4; cov+=$4}
			$1~/^GC(F|L)/ && $3>gc_i {gc_i=$3; gc=$2}
			END {
				OFS="\t"
				mean_cov = cov/genome_size
				mean_dp = dp/genome_size
				print mean_cov, mean_dp, mean_fraglen, gc
			}
		' $f
	done
}

idir=/willerslev/datasets/Island_Thrush/
do_collate $genome_size \
	$idir/{Cataponera_turdoides,Turdus_*}/MAPPED/*calmd.MQ30.bam.stats.txt \
	> collated_mapping_stats.tsv
do_collate $mt_size \
	$idir/{Cataponera_turdoides,Turdus_*}/MAPPED_Mt/*calmd.MQ30.bam.stats.txt \
	> collated_mapping_stats_Mt.tsv
