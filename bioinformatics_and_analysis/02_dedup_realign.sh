#!/bin/sh

source ./env.sh

fdir=/willerslev/scratch/grg/Island_Thrush_cleaned_Genomes
export ref=$fdir/P12064_102.fasta
export realign_threads=4

if [ ! -f ${ref%.fasta}.fai -o ! -f ${ref}.fai ]; then
	samtools faidx $ref || die "samtools faidx"
fi
if [ ! -f ${ref%.fasta}.dict ]; then
	# Watch out for this hanging with older picard tools on RedHat!
	# Version 2.10.5 hangs, whereas version 2.18.26 works ok.
	$java -jar $picard CreateSequenceDictionary R=$ref \
		|| die "CreateSequenceDictionary"
fi

do_rmdup_realign() {
	sm=$1
	in=$sm/MAPPED/${sm}.bam
	out1=${in%.bam}.rmdup.bam
	out2=${out1%.bam}.realigned.bam
	out3=${out2%.bam}.calmd.bam
	metrics=${in%.bam}.MarkDuplicates.txt

	# MarkDuplicates will happily create a file with any name,
	# but then RealignerTargetCreator complains if it doesn't have a
	# specific filename extension.  WTF were GATK devs thinking?
	intervals=${in%.bam}.realign.intervals

	if bam_check -i ${out3}; then
		# nothing to do
		return 0
	fi


	if ! bam_check -i ${out1}; then
		# Avoid "Too Many Open Files" errors in MarkDuplicates
		max_file_handles=$(($(ulimit -n) - 100))

		$java -jar $picard MarkDuplicates \
			I=$in \
			O=$out1 \
			M=$metrics \
			REMOVE_DUPLICATES=true \
			CLEAR_DT=false \
			CREATE_INDEX=true \
			MAX_FILE_HANDLES=$max_file_handles \
		|| die "$sm: MarkDuplicates"
	fi

	if ! bam_check ${out2}; then
		# Watch out for this hanging with GATK < 3.8-1-0 on RedHat!
		$java -jar $gatk -T RealignerTargetCreator \
			-nt $realign_threads \
			-R $ref \
			-I $out1 \
			-o $intervals \
		|| die "$sm: RealignerTargetCreator"

		$java -jar $gatk -T IndelRealigner \
			-R $ref \
			--bam_compression 0 \
			--disable_bam_indexing \
			-targetIntervals $intervals \
			-I $out1 \
			-o $out2 \
		|| die "$sm: IndelRealigner"
	fi

	samtools calmd \
		-b \
		$out2 \
		$ref \
		> $out3 \
	|| die "$sm: samtools calmd"

	samtools index $out3 || die "$sm: samtools index"

	rm $out1 ${out1%.bam}.bai $out2 $intervals
}
export -f do_rmdup_realign


#for sm in Turdus*; do
#for sm in Cataponera_turdoides; do
for sm in Turdus_merula; do
	echo do_rmdup_realign $sm
done \
 | parallel -j 40 --nice 10 --progress
