#!/usr/bin/env python3

import sys
import itertools

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} indfile outgroup")
        exit(1)

    ind_file = sys.argv[1]
    outgroup_str = sys.argv[2]

    inds = set()
    outgroup = None

    with open(ind_file) as f:
        for line in f:
            ind, *_ = line.split()
            if outgroup_str in ind:
                assert outgroup is None
                outgroup = ind
            inds.add(ind)

    if outgroup is None:
        print("outgroup not gound")
        exit(1)

    for a, b, c in itertools.permutations(inds, 3):
        print(a, b, c, outgroup)
