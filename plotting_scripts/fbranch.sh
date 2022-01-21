#!/bin/sh
# Make Fbranch plots (actually Dbranch), using Dsuite and a bunch of scripts.
# What a mess!!!

# Dsuite location.
Dsuite_dir=$HOME/src/Dsuite/
Dsuite=$Dsuite_dir/Build/Dsuite
dtools=./dtools-grg.py  # with my mods

# Newick tree with branch lengths removed.
# Same tree as: island_thrush_no-canescens.pi-NJ-SPR.tree.raxml.support.rooted
newick_tree=tree.for-fbranch.nwk
# sed 's/Turdus_poliocephalus/Tp/g' island_thrush_no-canescens.pi-NJ-SPR.tree.raxml.support.rooted > tree-dsuite.nwk
newick_tree2=tree-dsuite.nwk

# ABBA-BABA output from my eig-utils/abba-baba program, in ANGSD format.
angsd_abbababa=island_thrush_all_filter.abba-baba-eigutils.bootstrap
# ABBA-BABA output, formatted for Dsuite. This file won't be recognised by
# Dsuite unless the filename ends with "_tree.txt"
dsuite_abbababa=${angsd_abbababa}.Dsuite_tree.txt

# P-value threshold for significance using Holm-Bonferroni correction.
p_thres=$(python holm-bonferroni.py $angsd_abbababa | cut -f2 -d" ")

# Convert ABBA-BABA output to match Dsuite format.
#python angsd_abbababa_to_Dsuite.py $angsd_abbababa > $dsuite_abbababa

# Calcualte "Fbranch" stats from the tree and the D stats.
$Dsuite Fbranch \
    -p ${p_thres} \
    $newick_tree \
    $dsuite_abbababa \
    > Dbranch.txt


# Creates figure "Dbranch.pdf".
$dtools \
    -n Dbranch \
    --outgroup Turdus_merula \
    --metadata "sample_metadata-grg8.csv" \
    Dbranch.txt \
    $newick_tree2
    #--use_distances
