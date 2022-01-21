#!/bin/sh

source ./env.sh

fastme=$HOME/FastME/src/fastme
raxml_ng=/willerslev/software/raxml-ng

./submatrix.py \
	Turdus_poliocephalus_canescens_QMO19743 \
	< island_thrush_all.dist \
	> island_thrush_no-canescens.dist
./submatrix.py \
	Turdus_poliocephalus_canescens_QMO19743 \
	< island_thrush_all.reynoldsfstdist \
	> island_thrush_no-canescens.reynoldsfstdist
./submatrix.py \
	Turdus_poliocephalus_canescens_QMO19743 \
	Turdus_celaenops_NRM90188019 \
	Turdus_chrysolaus_orii_UWBM83174 \
	Turdus_feae_BMNH1905910913 \
	Turdus_obscurus_UWBM78352 \
	Turdus_pallidus_UWBM75229 \
	Turdus_niveiceps_NRM569330 \
	Turdus_merula \
	Cataponera_turdoides \
	< island_thrush_all.dist \
	> island_thrush_subset.dist
./submatrix.py \
	Turdus_poliocephalus_canescens_QMO19743 \
	Turdus_celaenops_NRM90188019 \
	Turdus_chrysolaus_orii_UWBM83174 \
	Turdus_feae_BMNH1905910913 \
	Turdus_obscurus_UWBM78352 \
	Turdus_pallidus_UWBM75229 \
	Turdus_niveiceps_NRM569330 \
	Turdus_merula \
	Cataponera_turdoides \
	< island_thrush_all.reynoldsfstdist \
	> island_thrush_subset.reynoldsfstdist

for pfx in island_thrush_no-canescens island_thrush_subset; do
	$fastme \
		-D 1 \
		-i ${pfx}.reynoldsfstdist \
		-o ${pfx}.rFst-NJ-SPR.tree \
		-m N \
		-s \
		|| die "fastme bailed"

	tree=${pfx}.pi-NJ-SPR.tree
	$fastme \
		-D 101 \
		-i ${pfx}.dist \
		-o ${tree} \
		-m N \
		-s \
		|| die "fastme bailed"

	head -n 1 ${tree} > ${tree}.main
	tail -n +2 ${tree} > ${tree}.boot
	$raxml_ng --support --tree ${tree}.main --bs-trees ${tree}.boot --prefix ${tree}
done
