#!/bin/sh

source ./env.sh

idir=/willerslev/scratch/fernando/Island_Thrush_cleaned_Genomes/Turdus_poliocephalus_vanikorensis_AMNH214400
odir=AMNH214400
export odir

mkdir -p $odir

do_idx() {
	sam=$1
	pfx=$(basename ${sam%_U.sam.gz})
	bam=$odir/${pfx}.bam
	samtools view \
		-Sbu \
		$sam \
	 | samtools sort \
		-O bam \
		-m 8G \
		> $bam \
	|| die "$pfx: samtools view|samtools sort"
	samtools index $bam || "$pfx: samtools index"
	samtools idxstats $bam > $odir/${pfx}.idxstats \
		|| "$pfx: samtools idxstats"
}
export -f do_idx

for sam in $idir/MAPPED/*U.sam.gz; do
	echo do_idx $sam
done | parallel --nice 10
