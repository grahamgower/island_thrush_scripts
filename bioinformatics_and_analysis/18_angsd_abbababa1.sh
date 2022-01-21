#!/bin/sh
# http://www.popgen.dk/angsd/index.php/Abbababa
# Output: Each lines represents a block with a chromsome name (Column 1),
# a start position (Column 2), an end postion (Column 3). The new columns
# are the counts of ABBA and BABA sites. For each combination of 3 individuals
# (H1,H2,H3) two columns are printed. These number served as input to the
# R script called jackKnife.R. This script will skip combinations of
# individuals if there is less than 3 blocks with data.
# Type "Rscript R/jackKnife.R" to see additional options.

source ./env.sh

idir=/willerslev/datasets/Island_Thrush/
bamlist1=$idir/island_thrush.bamlist.no-canecsens-no-outgroup
bamlist2=$idir/island_thrush.bamlist
bamlist=./bamlist.abbababa1
export ref=$idir/P12064_102.fasta
export regions=$idir/autosome-angsd.txt
export angsd=$HOME/angsd/angsd

cat $bamlist1 > $bamlist
grep merula $bamlist2 >> $bamlist

export minInd=$(($(wc -l $bamlist | cut -f1 -d" ") / 2))

$angsd \
	-doAbbababa 1 \
	-bam $bamlist \
	-rf $regions \
	-rmTrans 1 \
	-useLast 1 \
	-doCounts 1 \
	-enhance 1 \
	-minInd $minInd \
	-minQ 20 \
	-minMapQ 30 \
	\
	-ref $ref \
	-baq 2 \
	-C 50 \
	\
	-out island_thrush_abbababa1 \
	|| die "angsd bailed"

