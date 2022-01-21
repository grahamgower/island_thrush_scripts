#!/usr/bin/env python3

import csv

import numpy as np
from sklearn import decomposition

import matplotlib
import matplotlib.cm

matplotlib.use("Agg")  # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

from circledist import circledist
import matplotlib_style


def parse_locations(fn):
    loc = {}
    with open(fn) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # taxid = row["Species"]
            indid = row["Museum Number"]
            try:
                lat = float(row["Lat"])
            except ValueError:
                lat = None
            try:
                lon = float(row["Long"])
            except ValueError:
                lon = None

            # Adjust the lat/lon to avoid negative numbers and give sensible
            # distances between individuals.
            if lat is not None:
                lat = lat + 180
            if lon is not None:
                lon = lon if lon > 0 else 360 + lon

            loc[indid] = (lat, lon)
    return loc


def parse_list(fn):
    l = []
    with open(fn) as f:
        for line in f:
            l.append(line.rstrip())
    return l


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot PCA from covariance file.")
    parser.add_argument(
        "cov_file", metavar="file.cov", help="Covariance matrix, as output by PCAngsd."
    )
    parser.add_argument(
        "ind_file",
        metavar="ind.txt",
        help="Order of individuals in the cov file, with IDs matching the csv.",
    )
    parser.add_argument(
        "csv_file", metavar="file.csv", help="Metadata, with IDs and Lat/Lon columns."
    )
    parser.add_argument("out_file", metavar="out.pdf", help="output plot")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    loc = parse_locations(args.csv_file)
    indlist = parse_list(args.ind_file)
    rmidx = [i for i, ind in enumerate(indlist) if ind not in loc]
    for i in rmidx:
        print(f"{indlist[i]} has no location")
    lats = [loc[ind][0] for ind in indlist if ind in loc]
    lons = [loc[ind][1] for ind in indlist if ind in loc]

    ref_lat, ref_lon = 30, 120
    dists = [
        circledist(ref_lon, ref_lat, loc[ind][1], loc[ind][0])
        for ind in indlist
        if ind in loc
    ]

    n_pcs = 6
    C = np.loadtxt(args.cov_file)
    pca = decomposition.PCA(n_components=n_pcs)
    pc = pca.fit_transform(C)

    pdf = PdfPages(args.out_file)
    fig_w, fig_h = plt.figaspect(9.0 / 16.0)

    cmap = matplotlib.cm.get_cmap("plasma")
    distnorm = matplotlib.colors.Normalize(vmin=np.min(dists), vmax=np.max(dists))

    for pc_i in range(n_pcs - 1):
        fig1 = plt.figure(figsize=(fig_w, fig_h))
        gs1 = gridspec.GridSpec(1, 1)
        ax1 = fig1.add_subplot(gs1[0])

        x = np.delete(pc[:, pc_i], rmidx)
        y = np.delete(pc[:, pc_i + 1], rmidx)

        ax1.scatter(
            x,
            y,
            s=50,
            marker="o",
            alpha=1,
            lw=1,
            # edgecolor=cmap(latnorm(lats)),
            # edgecolor=cmap(lonnorm(lons)),
            facecolor=cmap(distnorm(dists)),
            # facecolor="none",
        )
        for i in rmidx:
            ax1.scatter(pc[i, pc_i], pc[i, pc_i + 1], s=50, marker="x", c="black")
        ax1.set_xlabel(f"PC{pc_i+1}")
        ax1.set_ylabel(f"PC{pc_i+2}")

        cb = fig1.colorbar(matplotlib.cm.ScalarMappable(norm=distnorm, cmap=cmap))
        cb.ax.get_yaxis().labelpad = 15
        cb.ax.set_ylabel("Distance from 120$^\circ$E, 30$^\circ$N", rotation=270)

        fig1.tight_layout()
        pdf.savefig(figure=fig1)

    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])
    ax1.bar(list(range(1, n_pcs + 1)), pca.explained_variance_)
    ax1.set_xlabel("Principal component")
    ax1.set_ylabel("Percentage variance explained")
    ax1.set_title(
        "Scree plot (total variance explained: {:.2f}\%)".format(
            np.sum(pca.explained_variance_)
        )
    )
    fig1.tight_layout()
    pdf.savefig(figure=fig1)

    pdf.close()
