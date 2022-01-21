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


def parse_eigutils_abbababa(filename):
    # chr     blockstart      P1      P2      P3      P4      AAAA    AAAB    AABA    ABAA    BAAA BBAA    ABBA    BABA    nsites  F4sum   Ddensum F4bc
    """
    np.loadtxt(
        filename,
        usecols=(0, 2, 3, 4, 12, 13, 15, 16),
        dtype=[
            ("chr", object),
            ("P1", object),
            ("P2", object),
            ("P3", object),
            ("ABBA", int),
            ("BABA", int),
            ("F4sum", float),
            ("Ddensum", float),
        ]
    )
    """
    data = dict()
    with open(filename) as f:
        next(f)  # skip header
        for line in f:
            fields = line.split()
            #chrom = fields[0]
            #blockstart = int(fields[1])
            p1 = fields[2]
            p2 = fields[3]
            p3 = fields[4]
            #p4 = fields[5]
            #AAAA = int(fields[6])
            #AAAB = int(fields[7])
            #AABA = int(fields[8])
            #ABAA = int(fields[9])
            #BAAA = int(fields[10])
            #BBAA = int(fields[11])
            ABBA = int(fields[12])
            BABA = int(fields[13])
            #nsites = int(fields[14])
            #F4sum = float(fields[14])
            #Ddensum = float(fields[15])

            key = (p1, p2, p3)
            if key not in data:
                data[key] = []
            #data[key].append((ABBA, BABA, F4sum, Ddensum))
            data[key].append((ABBA, BABA))

    for key, value in data.items():
        ABBA = np.array([v[0] for v in value], dtype=np.int32)
        BABA = np.array([v[1] for v in value], dtype=np.int32)
        #F4sum = np.array([v[2] for v in value])
        #Ddensum = np.array([v[3] for v in value])
        #data[key] = (ABBA, BABA, F4sum, Ddensum)
        data[key] = (ABBA, BABA)

    return data


@njit()
def boot1(abba, baba, weight):
    target_weight = weight.sum()
    wsum = 0
    nsum = 0
    dsum = 0
    while wsum < target_weight:
        i = np.random.randint(len(weight))
        wsum += weight[i]
        nsum += abba[i] - baba[i]
        dsum += abba[i] + baba[i]
    return nsum / dsum


@njit(parallel=True)
def bootstrap(abba, baba, weight, seed, n=10_000):
    np.random.seed(seed)
    reps = np.empty(n)
    for i in range(n):
        reps[i] = boot1(abba, baba, weight)
    return np.mean(reps), np.std(reps)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} fasta.fai contig.list eigutils.abbababa")
        exit(1)

    fai_file = sys.argv[1]
    contig_list_file = sys.argv[2]
    eigutils_file = sys.argv[3]

    contig_list = parse_list(contig_list_file)
    seqid, length = parse_fai(fai_file, contig_filter=contig_list)
    data = parse_eigutils_abbababa(eigutils_file)

    print("H1", "H2", "H3", "nABBA", "nBABA", "Dstat", "bootEst", "SE", "Z", sep="\t")
    rng = np.random.seed(123)
    seeds = np.random.randint(0, 2**31, size=len(data))
    for (key, (abba, baba)), seed in zip(data.items(), seeds):
        nABBA = np.sum(abba)
        nBABA = np.sum(baba)
        Dstat = (nABBA - nBABA) / (nABBA + nBABA)
        bootEst, SE = bootstrap(abba, baba, length, seed=seed)
        Z = bootEst / SE
        print(*key, nABBA, nBABA, Dstat, bootEst, SE, Z, sep="\t")
