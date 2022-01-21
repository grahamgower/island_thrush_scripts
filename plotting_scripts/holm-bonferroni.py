#!/usr/bin/env python3
import sys

import numpy as np
import scipy.stats

def parse_angsd_abbababa(filename):
    # print("H1", "H2", "H3", "nABBA", "nBABA", "Dstat", "bootEst", "SE", "Z", sep="\t")
    return np.loadtxt(
        filename,
        usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8),
        skiprows=1,
        dtype=[
            ("h1", object),
            ("h2", object),
            ("h3", object),
            ("nABBA", float),
            ("nBABA", float),
            ("Dstat", float),
            ("bootEst", float),
            ("SE", float),
            ("Z", float),
        ],
        unpack=True,
    )

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} angsd.abbababa")
        exit(1)

    angsd_file = sys.argv[1]
    *_, Z = parse_angsd_abbababa(angsd_file)
    Z[np.where(Z < 0)] *= -1  # make negative values positve
    p = scipy.stats.norm.sf(Z)  # 1-tailed p-value
    p = np.sort(p)
    alpha = 0.05
    for i, pval in enumerate(p):
        thres = alpha / (len(p) - i)
        if pval >= thres:
            print(i, thres)
            break
