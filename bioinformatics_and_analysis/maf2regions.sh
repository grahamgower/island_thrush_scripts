#!/bin/sh

maf=island_thrush_dp10.mafs.gz
snppos=island_thrush_dp10.snppos.tab

zcat $maf \
 | awk 'BEGIN {OFS="\t"} $3$4!="CT" && $3$4!="TC" && $3$4!="GA" && $3$4!="AG" && NR>1 {print $1, $2}' \
 > $snppos
