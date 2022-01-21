#!/bin/bash

source ./env.sh

fdir=/willerslev/scratch/fernando/Island_Thrush_cleaned_Genomes
ref=$fdir/P12064_102.fasta
bwathreads=32
sortthreads=4
export fdir
export ref
export bwathreads
export sortthreads

do_map1() {
	id=$1
	sm=$2
	lb=$3
	out=$4
	shift; shift; shift; shift
	in="$@"

	#echo "$id ## $sm ## $lb ## $out ## $in"
	#return

        bwa mem \
		-t $bwathreads \
		-M \
		-R "@RG\tID:${id}\tSM:${sm}\tLB:${lb}\tPL:ILLUMINA" \
		$ref \
		$in \
	 | samtools view \
		-Sbu \
		- \
	 | samtools sort \
		-O bam \
		-m 8G \
		-@ $sortthreads \
		- \
	 	> $out \
	|| die "$sm:$id: bwa mem|samtools view|samtools sort"

	samtools index $out || die "$sm:$id: samtools index"
}
export -f do_map1

do_map() {
    sm=$1
    echo $sm
    mkdir -p $sm/MAPPED
    libdir=$sm/CLEANED/6.complex_checked/*
    lb=$(echo $libdir | cut -d "/" -f 4)
    for dir in $libdir/*; do
        echo $dir
        id=$(echo $dir | cut -d "/" -f 5)
        # SINGLE READS
	do_map1 \
		$id \
		$sm \
		$lb \
        	$sm/MAPPED/${id}_U.bam \
        	$dir/${id}_U.fastq.gz || exit 1
        # PAIRED-END READS
        do_map1 \
		$id \
		$sm \
		$lb \
        	$sm/MAPPED/${id}_P.bam \
        	$dir/${id}_R1.fastq.gz \
        	$dir/${id}_R2.fastq.gz || exit 1
    done
    samtools merge \
	$sm/MAPPED/${sm}.bam \
	$sm/MAPPED/*_{P,U}.bam \
     || die "$sm: samtools merge"
    samtools index \
	$sm/MAPPED/${sm}.bam \
     || die "$sm: samtools index"
}
export -f do_map

#for species in Turdus*; do
#for species in Cataponera_turdoides; do
for species in Turdus_merula; do
	echo do_map $species
done \
 | parallel -j 5 --nice 10 --progress
