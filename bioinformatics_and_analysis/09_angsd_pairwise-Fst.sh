#!/bin/sh

# http://www.popgen.dk/angsd/index.php/SFS_Estimation
# http://www.popgen.dk/angsd/index.php/Fst

source ./env.sh

export realSFS=$HOME/angsd/misc/realSFS

idir=/willerslev/datasets/Island_Thrush/
export odir=./ANGSD_sfs
mkdir -p $odir

bamlist=$idir/island_thrush.bamlist
bamlist_subset=$idir/island_thrush.bamlist.no-canecsens-no-outgroup

do_sfs() {
	sm1=$1
	sm2=$2

	saf1=$sm1/ANGSD/${sm1}.saf.idx
	saf2=$sm2/ANGSD/${sm2}.saf.idx
	pair=$odir/${sm1}.${sm2}
	reynolds=${pair}.Reynolds
	bhatia=${pair}.Bhatia
	prior=${pair}.ml

	# Calculate the 2D-SFS prior.
	$realSFS $saf1 $saf2 \
		-maxIter 100 \
		-fold 1 \
		> $prior \
		|| die "$sm1/$sm2: realSFS"

	# Calculate per-site Reynolds Fst.
	$realSFS fst index $saf1 $saf2 \
		-sfs $prior \
		-whichFst 0 \
		-fstout $reynolds \
		-fold 1 \
		|| die "$sm1/$sm2: realSFS fst index -whichFst 0"

	# Calculate per-site Bhatia Fst.
	$realSFS fst index $saf1 $saf2 \
		-sfs $prior \
		-whichFst 1 \
		-fstout $bhatia \
		-fold 1 \
		|| die "$sm1/$sm2: realSFS fst index -whichFst 1"

	# Calculate global Reynolds Fst.
	$realSFS fst stats ${reynolds}.fst.idx \
		-fold 1 \
		-whichFst 0 \
		> ${reynolds}.stats \
		|| die "$sm1/$sm2: realSFS fst stats -whichFst 0"

	# Calculate global Bhatia Fst.
	$realSFS fst stats ${bhatia}.fst.idx \
		-fold 1 \
		-whichFst 1 \
		> ${bhatia}.stats \
		|| die "$sm1/$sm2: realSFS fst stats -whichFst 1"
}
export -f do_sfs

# Print all combinations of the lines in the input file.
comb() {
	file=$1
	awk '
		{ a[NR] = $0 }
		END {
			for (i=1; i<=length(a); i++)
				for (j=i+1; j<=length(a); j++)
					print a[i], a[j]
		}
	' $file
}

comb $bamlist | while read f1 f2; do
	sm1=$(basename ${f1%%.*bam})
	sm2=$(basename ${f2%%.*bam})
	#echo $sm1
	#echo $sm2
	echo do_sfs $sm1 $sm2 $odir
done | parallel -j 24 --nice 10


# Collate pairwise distances.
python fst_matrix.py -r $bamlist > island_thrush_all.reynoldsfstdist
python fst_matrix.py -b $bamlist > island_thrush_all.bhatiafstdist
python fst_matrix.py -r $bamlist_subset > island_thrush_subset.reynoldsfstdist
python fst_matrix.py -b $bamlist_subset > island_thrush_subset.bhatiafstdist
