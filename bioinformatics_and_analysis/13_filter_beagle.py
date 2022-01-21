#!/usr/bin/env python3
# Beagle format is "marker allele1 allele2 Ind0 Ind0 Ind0 Ind1 Ind1 Ind1 ..."
# So we need the column numbers of individuals to remove them.

import sys

bamlist="/willerslev/datasets/Island_Thrush/island_thrush.bamlist"

# Top one is very low coverage, the rest are outgroups.
to_remove="""
Turdus_poliocephalus_canescens_QMO19743
Cataponera_turdoides
Turdus_celaenops_NRM90188019
Turdus_chrysolaus_orii_UWBM83174
Turdus_feae_BMNH1905910913
Turdus_merula
Turdus_niveiceps_NRM569330
Turdus_obscurus_UWBM78352
Turdus_pallidus_UWBM75229
"""

def parse_bamlist(filename):
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            fields = line.split("/")
            yield fields[4]

indlist = list(parse_bamlist(bamlist))
#print(*indlist, sep="\n", file=sys.stderr)

rm = set()
for ind in to_remove.split():
    j = indlist.index(ind)
    rm.add(j)

# indexes to keep
indexes = []
for i in range(3 + 3 * len(indlist)):
    if i >= 3:
        j = (i - 3) // 3
        if j in rm:
            continue
    indexes.append(i)

assert len(indexes) == 3 + 3*(len(indlist) - len(rm))

for line in sys.stdin:
    fields = line.split()
    print(*[fields[i] for i in indexes], sep="\t")
