#!/usr/bin/env python3

import sys

import numpy as np


def parse_dist_matrices(f):
    num_inds = -2
    state = 0
    indlist = []
    matrix = []
    for lineno, line in enumerate(f, 1):
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) == 1 and state == 0:
            num_inds = int(fields[0])
            state = 1
        elif len(fields) == num_inds+1:
            indlist.append(fields[0])
            matrix.append(fields[1:])
            if state == num_inds:
                yield indlist, np.array(matrix)
                state = 0
                indlist = []
                matrix = []
            else:
                state += 1
        else:
            raise RuntimeError(f"{lineno}: unexpected line")


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
            description="Remove individuals from distance matrices.")
    parser.add_argument("inds", nargs="+")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    args.inds = set(args.inds)

    with sys.stdin as f:
        for indlist, matrix in parse_dist_matrices(f):
            to_remove = [indlist.index(ind) for ind in args.inds]
            indlist2 = [ind for ind in indlist if ind not in args.inds]
            matrix = np.delete(matrix, to_remove, axis=0)
            matrix = np.delete(matrix, to_remove, axis=1)
            print()
            print(matrix.shape[0])
            for ind, row in zip(indlist2, matrix):
                print(ind, *row, sep="\t")
