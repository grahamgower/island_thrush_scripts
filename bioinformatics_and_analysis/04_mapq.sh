#!/bin/sh

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
idir=.

do_mq_filter() {
	in=$1
	out=${in%.bam}.MQ30.bam
	samtools view -q 30 -b -o $out $in 
	samtools index $out
}
export -f do_mq_filter

#for dir in $idir/Turdus_*; do
#for dir in $idir/Cataponera_turdoides; do
for dir in $idir/Turdus_merula; do
	sm=$(basename $dir)
	in=$dir/MAPPED/${sm}.rmdup.realigned.calmd.bam
	echo do_mq_filter $in
done | parallel --nice 10
