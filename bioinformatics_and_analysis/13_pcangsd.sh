#!/bin/sh

if [ "$CONDA_DEFAULT_ENV" != "pcangsd" ]; then
	echo "Error: pcangsd conda env not loaded" >&2
	exit 1
fi

pcangsd="python $HOME/pcangsd/pcangsd.py"

ipfx=island_thrush_subset

for inbreed in 3; do
	opfx=${ipfx}.pcangsd-inbreed$inbreed

	nice -n 5 $pcangsd \
		-beagle ${ipfx}.beagle.gz \
		-o ${opfx} \
		-sites_save \
		-selection \
		-snp_weights \
		-inbreedSites \
		-inbreed $inbreed \
		-kinship \
		-threads 80

# Filter outlier sites using inbreeding coefficients from the last run.
	nice -n 5 $pcangsd \
		-beagle ${ipfx}.beagle.gz \
		-o ${opfx}-filt \
		-sites_save \
		-selection \
		-snp_weights \
		-hwe ${opfx}.lrt.sites.npy \
		-inbreedSites \
		-inbreed $inbreed \
		-kinship \
		-threads 80
done
