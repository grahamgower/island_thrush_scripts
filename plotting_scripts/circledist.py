#!/usr/bin/env python

import sys
import csv
from numpy import sin, cos, arccos, pi, floor


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Print phylip/fastme distance matrix from lat/lon data."
    )
    parser.add_argument(
        "csv_fn", metavar="file.csv", help="input csv, with Lat/Long columns."
    )
    return parser.parse_args()


def parse_locations(fn):
    with open(fn) as f:
        reader = csv.DictReader(f)
        for row in reader:
            taxid = row["Species"]
            indid = row["Museum Number"]
            lat = float(row["Lat"])
            lon = float(row["Long"])
            yield taxid, indid, lat, lon


# Distance on the surface of the earth.
def circledist(x1, y1, x2, y2, radius=6371, degrees=True):
    if degrees:
        x1 *= pi / 180
        x2 *= pi / 180
        y1 *= pi / 180
        y2 *= pi / 180

    sigma = sin(y1) * sin(y2) + cos(y1) * cos(y2) * cos(x2 - x1)

    # Sigma should be in [-1,1], but floating point is not perfect.
    if sigma > 1:
        sigma = 1
    elif sigma < -1:
        sigma = -1

    return radius * arccos(sigma)


if __name__ == "__main__":
    args = parse_args()

    inds = []
    for tax, id, lat, lon in parse_locations(args.csv_fn):
        name = tax + "_" + id
        inds.append((name, lon, lat))

    print(len(inds))
    for name, lon1, lat1 in inds:
        print(name, end="")
        for _, lon2, lat2 in inds:
            d = circledist(lon1, lat1, lon2, lat2)
            # Round to nearest km.
            d = int(floor(d + 0.5))
            print("\t", d, sep="", end="")
        print()
