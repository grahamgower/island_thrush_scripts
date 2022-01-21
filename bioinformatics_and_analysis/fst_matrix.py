#!/usr/bin/env python3
# Output a phylip/fastME compatible distance matrix.

import sys
import glob
import collections

def parse_bamlist(filename):
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            fields = line.split("/")
            yield fields[4]

if __name__ == "__main__":
    def usage():
        print(f"usage: {sys.argv[0]} [-r|-b] bamlist", file=sys.stderr)
        exit(1)

    if len(sys.argv) != 3:
        usage()

    if sys.argv[1] == "-r":
        which = "Reynolds"
    elif sys.argv[1] == "-b":
        which = "Bhatia"
    else:
        print(f"unknown first parameter '{sys.argv[1]}'", file=sys.stderr)
        usage()

    bamlist = sys.argv[2]
    taxa = list(parse_bamlist(bamlist))

    etoolong = False
    for i, taxon in enumerate(taxa):
        #if len(taxon) > 10:
        #    print("{}: name too long for stupid phylip".format(taxon),
        #        file=sys.stderr)
        if len(taxon) > 64:
            print("{}: name too long for FastME".format(taxon),
                file=sys.stderr)
            etoolong = True

    if etoolong:
        exit(1)

    print(len(taxa))

    mat = collections.defaultdict(dict)

    for i, t1 in enumerate(taxa):
        mat[t1][t1] = 0
        for t2 in taxa[i+1:]:
            #fn1 = "ANGSD_sfs/{}.{}.Reynolds.stats".format(t1, t2)
            #fn2 = "ANGSD_sfs/{}.{}.Reynolds.stats".format(t2, t1)
            fn1 = "ANGSD_sfs/{}.{}.{}.stats".format(t1, t2, which)
            fn2 = "ANGSD_sfs/{}.{}.{}.stats".format(t2, t1, which)
            try:
                f = open(fn1)
            except:
                f = open(fn2)
            weighted = float(f.read().split()[-1])
            f.close()

            if weighted < 0:
                print("negative Fst", t1, t2, weighted, file=sys.stderr)
                weighted = 0

            mat[t1][t2] = weighted
            mat[t2][t1] = weighted

    for t1 in taxa:
        print(t1, end="")
        for t2 in taxa:
            try:
                x =  mat[t1][t2]
            except IndexError:
                print(t1, file=sys.stderr)
                print(t2, file=sys.stderr)
                raise
            print("\t", x, sep="", end="")
        print()
