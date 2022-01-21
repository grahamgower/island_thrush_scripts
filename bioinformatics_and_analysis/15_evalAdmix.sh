#!/bin/sh

if [ "$CONDA_DEFAULT_ENV" != "pcangsd" ]; then
	echo "Error: pcangsd conda env not loaded" >&2
	exit 1
fi

evalAdmix=$HOME/evalAdmix/evalAdmix
pfx=island_thrush_dp10

# convert npy file from PCAngsd to text file expected by evalAdmix
dump_npy()
{
	infile=$1
	outfile=$2
	transpose=$3
	if $transpose; then
		python -c "import numpy as np; np.savetxt('$outfile', np.load('$infile').T)"
	else
		python -c "import numpy as np; np.savetxt('$outfile', np.load('$infile'))"
	fi
}

for k in 2 3 4 5 6 7 8; do
	opfx=${pfx}.pcangsd.admix-k${k}

	#dump_npy ${opfx}.admix.F.npy ${opfx}.admix.F.txt true
	#dump_npy ${opfx}.admix.Q.npy ${opfx}.admix.Q.txt false

	nice $evalAdmix \
		-beagle ${pfx}.beagle.gz \
		-fname ${opfx}.admix.F.txt \
		-qname ${opfx}.admix.Q.txt \
		-o ${opfx}.evalAdmix \
		-P 80
done
