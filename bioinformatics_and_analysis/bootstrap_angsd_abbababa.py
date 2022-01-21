#!/usr/bin/env python3
import sys

import numpy as np
from numba import njit, prange


def parse_fai(filename, contig_filter):
    seqid, length = np.loadtxt(
        filename,
        usecols=(0, 1),
        dtype=[("seqid", object), ("length", int)],
        unpack=True,
    )
    contig_filter = set(contig_filter)
    idx = []
    for i, sid in enumerate(seqid):
        if sid in contig_filter:
            idx.append(i)
    return seqid[idx], length[idx]


def parse_list(filename):
    seqid = np.loadtxt(
        filename,
        dtype=object,
    )
    return seqid


def parse_angsd_abbababa(filename):
    with open(filename) as f:
        first_line = next(f)
    ncols = len(first_line.split())
    return np.loadtxt(
        filename,
        usecols=range(3, ncols),
        dtype=np.int32,
    )


def angsd_triplets(n):
    """
    Indices of population triplets used by ANGSD.
    This was translated from the ANGSD jackKnife.R code to avoid errors.
    """
    triplets = dict()
    i = 0
    for h3 in range(n):
        for h2 in range(n):
            if h2 == h3:
                continue
            for h1 in range(n):
                if h1 == h3:
                    continue
                if h1 >= h2:
                    continue
                triplets[i] = (h1, h2, h3)
                i += 1
    return triplets


@njit()
def boot1(abba, baba, weight):
    target_weight = weight.sum()
    wsum = 0
    nsum = 0
    dsum = 0
    while True:
        i = np.random.randint(len(weight))
        wsum += weight[i]
        if wsum > target_weight:
            break
        nsum += abba[i] - baba[i]
        dsum += abba[i] + baba[i]
    return nsum / dsum


@njit #(parallel=True)
def bootstrap(abba, baba, weight, seed, n=10_000):
    np.random.seed(seed)
    reps = np.empty(n)
    for i in range(n):
        reps[i] = boot1(abba, baba, weight)
    return np.mean(reps), np.std(reps, ddof=1)


@njit(parallel=True)
def do_bootstrap_all(data, seed=123):
    #print("H1", "H2", "H3", "nABBA", "nBABA", "Dstat", "bootEst", "SE", "Z", sep="\t")
    rng = np.random.seed(seed)
    seeds = np.random.randint(0, 2**31, size=n)
    results = np.zeros(shape=(n, 6))
    for i in prange(n):
        seed = seeds[i]
        abba = data[:,2*i]
        baba = data[:,2*i + 1]
        nABBA = np.sum(abba)
        nBABA = np.sum(baba)
        Dstat = (nABBA - nBABA) / (nABBA + nBABA)
        bootEst, SE = bootstrap(abba, baba, length, seed=seed)
        Z = bootEst / SE
        results[i, :] = np.array([nABBA, nBABA, Dstat, bootEst, SE, Z])
    return results

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"usage: {sys.argv[0]} fasta.fai contig.list ind.list angsd.abbababa")
        exit(1)

    fai_file = sys.argv[1]
    contig_list_file = sys.argv[2]
    ind_list_file = sys.argv[3]
    abbababa_file = sys.argv[4]

    contig_list = parse_list(contig_list_file)
    ind_list = parse_list(ind_list_file)
    seqid, length = parse_fai(fai_file, contig_filter=contig_list)
    data = parse_angsd_abbababa(abbababa_file)
    n = data.shape[1] // 2
    m = len(ind_list)
    ecols = m*(m-1)*(m-2) // 2
    if ecols != n:
        raise RuntimeError(
            f"abbababa file has 2*{n} cols, but {ind_list_file} has "
            f"{len(ind_list)} inds (expected 2*{ecols} cols)"
        )

    # neuter length for testing
    #length = length[:data.shape[0]]

    #print(data.shape, n, m)
    triplets = angsd_triplets(m)

    print("H1", "H2", "H3", "nABBA", "nBABA", "Dstat", "bootEst", "SE", "Z", sep="\t")
    for i, result in enumerate(do_bootstrap_all(data)):
        h1, h2, h3 = triplets[i]
        print(ind_list[h1], ind_list[h2], ind_list[h3], *result, sep="\t")
