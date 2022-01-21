#!/usr/bin/env python3

import csv

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, leaders
from scipy.spatial.distance import squareform

import matplotlib
import matplotlib.cm

matplotlib.use("Agg")  # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits import axes_grid1

import matplotlib_style
from display_name import italicise_display_name


def parse_metadata(fn):
    namemap = {}
    with open(fn) as f:
        reader = csv.DictReader(f)
        for row in reader:
            #taxid = row["Species"]
            indid = row["Museum Number"]
            dispname = row["Display Name"]
            namemap[indid] = italicise_display_name(dispname)
    return namemap


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot PCA from covariance file.")
    parser.add_argument(
        "cov_file",
        metavar="file.cov",
        help="Covariance matrix, as output by PCAngsd.",
    )
    parser.add_argument(
        "ind_file",
        metavar="ind.txt",
        help="Order of individuals in the cov file, with IDs matching the csv.",
    )
    parser.add_argument(
        "csv_file",
        metavar="file.csv",
        help="Metadata, with IDs and Lat/Lon columns.",
    )
    parser.add_argument("out_file", metavar="out.pdf", help="output plot")
    return parser.parse_args()


def corrcov(C):
    """
    Return the correlation matrix for covariance matrix C.
    """
    n = C.shape[0]
    R = np.empty_like(C)
    sigma = np.sqrt(np.diag(C))
    for j in range(n):
        for k in range(n):
            R[j, k] = C[j, k] / (sigma[j] * sigma[k])
    return R


def cov2dist(C):
    """
    Convert covariance matrix into a distance matrix.
    """
    R = 1 - (1 + corrcov(C)) / 2
    np.fill_diagonal(R, 0)
    return R


if __name__ == "__main__":
    args = parse_args()
    ind = np.loadtxt(args.ind_file, dtype=str)
    namemap = parse_metadata(args.csv_file)

    n_pcs = 6
    C = np.loadtxt(args.cov_file)
    n = C.shape[0]
    R = corrcov(C)
    D = cov2dist(C)
    y = squareform(D)

    Z = linkage(y, "average", optimal_ordering=True)
    # Z = linkage(y, "weighted", optimal_ordering=True)
    groups = fcluster(Z, n_pcs, "maxclust")

    pdf = PdfPages(args.out_file)
    fig1, ax_heatmap = plt.subplots(figsize=3.5 * plt.figaspect(1.0))
    ax_heatmap.xaxis.tick_top()

    divider = axes_grid1.make_axes_locatable(ax_heatmap)
    ax_colorbar = divider.append_axes("right", 0.2, pad=0.05)

    for sp in ax_heatmap.spines.keys():
        ax_heatmap.spines[sp].set_visible(False)

    dendro = dendrogram(Z, distance_sort=True, no_plot=True)

    leaves = dendro["leaves"]
    rleaves = list(reversed(leaves))
    R = R[:, leaves]
    R = R[rleaves, :]

    xnames = [namemap[ind[i]] for i in leaves]
    ynames = [namemap[ind[i]] for i in rleaves]

    im = ax_heatmap.imshow(
        R, vmin=-1, vmax=1,
    )
    im.set_cmap("coolwarm")
    ax_heatmap.set_xticks(range(len(xnames)))
    ax_heatmap.set_xticklabels(xnames, rotation=90)
    ax_heatmap.set_yticks(range(len(ynames)))
    ax_heatmap.set_yticklabels(ynames)
    cbar = plt.colorbar(im, cax=ax_colorbar)
    cbar.ax.set_xlabel("Correlation", labelpad=10, size=14)
    cbar.set_ticks([-1, 0, 1])
    cbar.ax.tick_params(labelsize=14)

    fig1.tight_layout()
    pdf.savefig(figure=fig1)
    pdf.close()
