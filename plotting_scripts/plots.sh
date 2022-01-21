#!/bin/sh

metadata=sample_metadata-grg8.csv
individuals=ind.txt
covariance=island_thrush_subset.pcangsd-inbreed3-filt.cov
pi_tree=island_thrush_no-canescens.pi-NJ-SPR.tree.raxml.support
fst_tree=island_thrush_no-canescens.rFst-NJ-SPR.tree
#pi_tree_no_outgroup=island_thrush_subset.pi-NJ-SPR.tree.raxml.support
#fst_tree_no_outgroup=island_thrush_subset.rFst-NJ-SPR.tree

# Fig 1. pairwise distance phylogeny with map on right.
./phylo_with_map.py -b $metadata $pi_tree Fig1_pi_tree_and_map.pdf

# Fig S3. Fst version of Fig 1.
./phylo_with_map.py $metadata $fst_tree FigS3_fst_tree_and_map.pdf

# Fig 2. pairwise distance phylogeny with lines connecting to map.
# Coloured by heterozygosity.
./phylo_with_map2.py -H $metadata $pi_tree Fig2_pi_tree_with_lines_to_map_heterozygosity.pdf

# Fst version of Fig 2.
#./phylo_with_map2.py -H $metadata $fst_tree fst_tree_with_lines_to_map_heterozygosity.pdf

# Fig S5. PCA plots.
./plot_cov.py $covariance $individuals $metadata FigS5_pcangsd_PCA.pdf

# Fig S6. Genetic correlation (from PCAngsd covariance) heatmap, with dendrograms in the margins.
./plot_cov_heatmap.py $covariance $individuals $metadata FigS6_pcangsd_heatmap.pdf

# Fig S7. "ADMIXTURE" barplots.
./plot_admix.py $covariance $individuals $metadata FigS7_pcangsd_admix.pdf \
	npy/island_thrush_subset.*.admix.Q.npy

# Fig S8. "ADMIXTURE" barplots, for all individuals with depth >= 10 (includes outgroup).
./plot_admix.py island_thrush_dp10.pcangsd.cov ind_dp10.txt $metadata FigS8_pcangsd_admix_with_outgroups.pdf \
	npy/island_thrush_dp10.*.admix.Q.npy

# Genetic similarity (from PCAngsd covariance) versus geographic distance.
#./plot_cov_geodist.py $covariance $individuals $metadata pcangsd_covariance_vs_geodist.pdf

# Fig S9. Genetic similarity (pairwise) versus geographic distance.
./plot_geneticdist_geodist.py island_thrush_subset.dist $metadata FigS9_pidist_vs_geodist.pdf

sh fbranch.sh && mv Dbranch.pdf Fig4_Dbranch.pdf
