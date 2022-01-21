#!/bin/sh

if [ "$CONDA_DEFAULT_ENV" != "pcangsd" ]; then
	echo "Error: pcangsd conda env not loaded" >&2
	exit 1
fi

pcangsd="python $HOME/pcangsd/pcangsd.py"

ipfx=island_thrush_subset

for k in 2 3 4 5 6 7 8; do
	opfx=${ipfx}.pcangsd.admix-k${k}

	nice -n 5 $pcangsd \
		-beagle ${ipfx}.beagle.gz \
		-o ${opfx} \
		-admix \
		-admix_save \
		-admix_K $k \
		-hwe ${ipfx}.pcangsd-inbreed3.lrt.sites.npy \
		-threads 80
done
