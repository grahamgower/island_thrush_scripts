#!/usr/bin/env python3

import operator
import csv
import random

import ete3
from plot_eteTree import plot_tree

import cartopy
import cartopy.crs as ccrs

import numpy as np
import matplotlib

matplotlib.use("Agg")  # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from circledist import circledist
import matplotlib_style
from display_name import italicise_display_name


def tree_reorder(tree, reverse=False):
    # Order child nodes by the number of decendents,
    for node in tree.traverse():
        if len(node.children) != 2:
            continue
        order = len(node.children[0].get_leaves()) > len(node.children[1].get_leaves())
        if reverse:
            if not order:
                node.swap_children()
        else:
            if order:
                node.swap_children()


def tree_order_by_geography(tree, location):
    # Order child nodes by the mean longitude of descendents
    """
    def mean_lat(nodes):
        lats = []
        for node in nodes:
            _, lat = location[node.name]
            lat + 180
            lats.append(lat)
        return np.mean(lats)
    """

    def mean_lon(nodes):
        lons = []
        for node in nodes:
            lon, _ = location[node.name]
            if lon < 0:
                lon += 360
            lons.append(lon)
        return np.mean(lons)

    for node in tree.traverse("postorder"):
        if len(node.children) != 2:
            continue
        order = mean_lon(node.children[0].get_leaves()) > mean_lon(
            node.children[1].get_leaves()
        )

        if order:
            node.swap_children()


def plot_map(ax, location, fgcolours, bgcolours, axis_ticklabels=False, jitter=1):
    ax.set_extent([90, 195, -36, 72], crs=ccrs.PlateCarree())

    ax.gridlines(
        xlocs=list(range(75, 250 + 1, 15)),
        draw_labels=False,
        color="black",
        alpha=0.2,
        linestyle="--",
    )

    if axis_ticklabels:
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        gl = ax.gridlines(
            xlocs=[90, 120, 150, 180],
            ylocs=[-20, 0, 20, 40, 60],
            draw_labels=True,
            color="black",
            alpha=0.2,
            linestyle="--",
        )
        gl.xlabels_top = False
        #gl.ylabels_left = False
        gl.xlines = False
        gl.ylines = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

    feature_scale = "110m"  # options are: 10m, 50m, 110m.
    linewidth = 0.5
    ax.add_feature(cartopy.feature.OCEAN.with_scale(feature_scale), linewidth=0)
    ax.add_feature(cartopy.feature.LAND.with_scale(feature_scale), linewidth=0)
    ax.coastlines(feature_scale, linewidth=linewidth)

    for taxid, (lon, lat) in location.items():
        if taxid not in bgcolours:
            continue
        bc = bgcolours[taxid]
        fc = fgcolours[taxid]
        if jitter:
            # jitter
            jit_dir = random.uniform(0, 2 * np.pi)
            jit_dist = random.uniform(0, jitter)
            lon += jit_dist * np.cos(jit_dir)
            lat += jit_dist * np.sin(jit_dir)
        ax.scatter(
            lon,
            lat,
            marker="o",
            transform=ccrs.PlateCarree(),
            facecolor=fc,
            edgecolor=bc,
            s=50,
            zorder=5,
            linewidth=linewidth,
        )


def parse_metadata(fn):
    location = dict()
    namemap = dict()
    with open(fn) as f:
        reader = csv.DictReader(f)
        for row in reader:
            taxid = row["Species"]
            indid = row["Museum Number"]
            dispname = row["Display Name"]
            lon = float(row["Long"])
            lat = float(row["Lat"])
            if indid:
                taxid = f"{taxid}_{indid}"
            location[taxid] = (lon, lat)
            namemap[taxid] = italicise_display_name(dispname)
    return location, namemap


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot phylogeny and a map showing the individuals."
    )
    parser.add_argument(
        "-j",
        "--jitter",
        default=0,
        type=float,
        help="randomly jitter map points by up to this many degrees "
        "in any direction [default=%(default)s]",
    )
    parser.add_argument(
        "-b",
        "--bootstrap",
        default=False,
        action="store_true",
        help="Show bootstrap support values.",
    )
    parser.add_argument(
        "csv_fn", metavar="file.csv", help="sample metadata, with Lat and Long columns."
    )
    parser.add_argument("tree_fn", metavar="tree.nwk", help="newick format input tree")
    parser.add_argument("pdf_fn", metavar="out.pdf", help="output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    location, namemap = parse_metadata(args.csv_fn)
    lon_centre = 140

    tree = ete3.Tree(args.tree_fn)

    if False:
        # midpoint rooting
        midpoint = tree.get_midpoint_outgroup()
        tree.set_outgroup(midpoint)
    if False:
        # root by outgroup clade
        out = [
            "Turdus_celaenops_NRM90188019",
            "Turdus_chrysolaus_orii_UWBM83174",
            "Turdus_feae_BMNH1905910913",
            "Turdus_obscurus_UWBM78352",
            "Turdus_pallidus_UWBM75229",
            "Turdus_niveiceps_NRM569330",
            "Turdus_merula",
            "Cataponera_turdoides",
        ]

        outgroup = tree.get_common_ancestor(out)
        tree.set_outgroup(outgroup)
    if True:
        tree.set_outgroup("Turdus_merula")

        # use less space for the branch leading to the ingroup
        outgroup, ingroup = tree.get_children()
        delta = 0.9 * outgroup.dist
        outgroup.dist += delta
        ingroup.dist -= delta
        # XXX: ete3 sets this to 1 by default
        ingroup.support = 100

    tree_reorder(tree)

    pdf = PdfPages(args.pdf_fn)
    # fig_w, fig_h = plt.figaspect((72+36)/(195-90))
    fig_w, fig_h = plt.figaspect(9.0 / 16.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))

    ax1 = fig1.add_subplot(121)
    ax1.set_position([0, 0, 0.34, 1])
    ax2 = fig1.add_subplot(
        122, projection=ccrs.PlateCarree(central_longitude=lon_centre), frame_on=True
    )
    ax2.set_position([0.5, 0.10, 0.5, 0.8])

    cmap = matplotlib.cm.get_cmap("plasma")
    ref_lat, ref_lon = 30, 120
    dists = {
        k: circledist(ref_lon, ref_lat, lon, lat) for k, (lat, lon) in location.items()
    }
    distnorm = matplotlib.colors.Normalize(
        vmin=np.min(list(dists.values())), vmax=np.max(list(dists.values()))
    )
    fgcolours = {}
    bgcolours = {}
    for n in tree.iter_leaves():
        style = n._get_style()
        taxid = n.name
        n.name = namemap[n.name]

        if taxid not in location:
            raise RuntimeError(f"{taxid} not in .csv")

        lon, lat = location[taxid]
        if taxid.startswith("Turdus_poliocephalus"):
            # colour by distance to lat/lon reference
            fgcolours[taxid] = cmap(distnorm(circledist(ref_lon, ref_lat, lon, lat)))
            bgcolours[taxid] = "black"
        else:
            fgcolours[taxid] = "darkgray"
            bgcolours[taxid] = "black"
            style["shape"] = None

        style["fgcolor"] = fgcolours[taxid]
        style["bgcolor"] = bgcolours[taxid]
        style["size"] = 3.5

    for n in tree.iter_descendants(strategy="postorder"):
        if not n.is_leaf() and n.support < 100:
            style = n._get_style()
            style["fgcolor"] = "red"
            style["bgcolor"] = "red"
            style["shape"] = None #"square"
            style["size"] = 0

    for n in tree.traverse():
        style = n._get_style()
        style["hz_line_width"] = 0.4
        style["vt_line_width"] = 0.4

    coords = plot_tree(tree, axe=ax1, font_size=4, bootstrap=args.bootstrap, align_names=True, align_colour="darkgray", scale_line_width=0.5)
    plot_map(ax2, location, fgcolours, bgcolours, jitter=args.jitter)

    # fig1.tight_layout()
    pdf.savefig(figure=fig1)
    pdf.close()
