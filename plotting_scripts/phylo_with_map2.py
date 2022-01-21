#!/usr/bin/env python3

# Do the same thing as in phylo-with-map.py, but without the outgroup clade.
# Also draw lines between phylogeny and map.

import operator
import csv

import ete3
from plot_eteTree import plot_tree2, to_coord

import cartopy
import cartopy.crs as ccrs

import numpy as np
import matplotlib

matplotlib.use("Agg")  # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection

from circledist import circledist
from phylo_with_map import parse_metadata, tree_order_by_geography

import matplotlib_style


def plot_map(ax, location, fgcolours, bgcolours, axis_ticklabels=False):
    ax.set_extent([90, 195, -36, 20.0001], crs=ccrs.PlateCarree())

    # ax.gridlines(xlocs=list(range(75,250+1,15)),
    #        draw_labels=False, color='black', alpha=0.2, linestyle='--')

    if axis_ticklabels:
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        gl = ax.gridlines(
            xlocs=[90, 120, 150, 180],
            ylocs=[-20, 0, 20],
            draw_labels=True,
            color="black",
            alpha=0.2,
            linestyle="--",
        )
        gl.xlabels_top = False
        gl.ylabels_left = False
        gl.xlines = False
        gl.ylines = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

    feature_scale = "110m"  # options are: 10m, 50m, 110m.
    linewidth = 0.5
    #ax.add_feature(cartopy.feature.OCEAN.with_scale(feature_scale), linewidth=0)
    ax.add_feature(cartopy.feature.LAND.with_scale(feature_scale), linewidth=0)
    ax.coastlines(feature_scale, linewidth=linewidth)

    for taxid, (lon, lat) in location.items():
        if taxid not in fgcolours:
            continue
        fc = fgcolours[taxid]
        bc = bgcolours[taxid]
        # fc = "red" if taxid.startswith("Turdus_poliocephalus") else "yellow"
        ax.scatter(
            lon,
            lat,
            transform=ccrs.PlateCarree(),
            facecolor=fc,
            edgecolor=bc,
            s=50,
            zorder=10,
            linewidth=linewidth,
        )

    ax.set_axis_off()
    # ax.outline_patch.set_visible(False)


def parse_het(fn):
    with open(fn) as f:
        reader = csv.DictReader(f)
        for row in reader:
            taxid = row["Species"]
            indid = row["Museum Number"]
            H = float(row["Heterozygosity"])
            if indid:
                taxid = f"{taxid}_{indid}"
            yield taxid, H


def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Plot individuals on a map.")
    parser.add_argument(
        "-H",
        "--heterozygosity",
        default=False,
        action="store_true",
        help="Colour by heterozygosity",
    )
    parser.add_argument(
        "-b",
        "--bootstrap",
        default=False,
        action="store_true",
        help="Show bootstrap support values.",
    )
    parser.add_argument(
        "csv_fn", metavar="file.csv", help="input csv, with Lat/Long columns."
    )
    parser.add_argument("tree_fn", metavar="tree.nwk", help="newick format input tree")
    parser.add_argument("pdf_fn", metavar="out.pdf", help="output file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    location, namemap = parse_metadata(args.csv_fn)
    location = {k: v for k, v in location.items() if k.startswith("Turdus_poliocephalus")}
    lon_centre = 140

    if args.heterozygosity:
        het = {tax: het for tax, het in parse_het(args.csv_fn)}

    tree = ete3.Tree(args.tree_fn)

    remove = set([
        "Turdus_celaenops_NRM90188019",
        "Turdus_chrysolaus_orii_UWBM83174",
        "Turdus_feae_BMNH1905910913",
        "Turdus_obscurus_UWBM78352",
        "Turdus_pallidus_UWBM75229",
        "Turdus_niveiceps_NRM569330",
        "Turdus_merula",
        "Cataponera_turdoides",
    ])
    tree.prune([n for n in tree.get_leaf_names() if n not in remove])

    if False:
        # midpoint rooting
        midpoint = tree.get_midpoint_outgroup()
        tree.set_outgroup(midpoint)
    if True:
        # pairwise-distance phylogeny outgroup
        tree.set_outgroup("Turdus_poliocephalus_mindorensis_AMNH575605")
        # use less space for the branch leading to the ingroup
        outgroup, ingroup = tree.get_children()
        delta = 0.9 * outgroup.dist
        outgroup.dist += delta
        ingroup.dist -= delta
        # XXX: ete3 sets this to 1 by default
        ingroup.support = 100
    if False:
        out = [
            "Turdus_poliocephalus_mindorensis_AMNH575605",
            "Turdus_poliocephalus_sspunknown_FMNH358378",
            "Turdus_poliocephalus_thomassoni_AMNH416872",
            "Turdus_poliocephalus_mayonensis_ZMUC02632",
        ]
        outgroup = tree.get_common_ancestor(out)
        tree.set_outgroup(outgroup)
    if False:
        # rFst phylogeny outgroup
        out = [
            "Turdus_poliocephalus_mindorensis_AMNH575605",
            "Turdus_poliocephalus_nigrorum_ZMUC138023",
        ]
        outgroup = tree.get_common_ancestor(out)
        tree.set_outgroup(outgroup)

    # tree_reorder(tree)
    tree_order_by_geography(tree, location)

    if args.heterozygosity:
        cmap = matplotlib.cm.get_cmap("coolwarm")
    else:
        cmap = matplotlib.cm.get_cmap("plasma")

    ref_lat, ref_lon = 30, 120
    dists = {
        taxid: circledist(ref_lon, ref_lat, lon, lat)
        for taxid, (lon, lat) in location.items()
    }
    distnorm = matplotlib.colors.Normalize(
        vmin=np.min(list(dists.values())), vmax=np.max(list(dists.values()))
    )
    if args.heterozygosity:
        list(het.values())
        hetlist = [het[node.name] for node in tree.get_leaves()]
        hetnorm = matplotlib.colors.Normalize(
            vmin=np.min(hetlist), vmax=np.max(hetlist)
        )
    fgcolours = {}
    bgcolours = {}
    for n in tree.iter_leaves():
        style = n._get_style()
        taxid = n.name
        n.name = namemap[n.name]
        n.taxid = taxid
        if not taxid.startswith("Turdus_poliocephalus"):
            raise RuntimeError(f"{taxid} isn't a T. poliocephalus subspecies")
        if args.heterozygosity:
            # colour by heterozygosity
            fgcolours[taxid] = cmap(hetnorm(het[taxid]))
            bgcolours[taxid] = "black"
        else:
            # colour by distance to reference point
            lon, lat = location[taxid]
            fgcolours[taxid] = cmap(distnorm(circledist(ref_lon, ref_lat, lon, lat)))
            bgcolours[taxid] = "black"

        style["fgcolor"] = fgcolours[taxid]
        style["bgcolor"] = bgcolours[taxid]
        style["size"] = 10

    if args.bootstrap:
        for n in tree.iter_descendants(strategy="postorder"):
            if not n.is_leaf() and n.support < 100:
                style = n._get_style()
                style["fgcolor"] = "red"
                style["bgcolor"] = "red"
                style["shape"] = "square"
                style["size"] = 3

    # Height of the tree, from root to most distant leaf node, in branch units.
    tree_height = max(n.get_distance(tree) for n in tree.iter_leaves())

    pdf = PdfPages(args.pdf_fn)
    fig1 = plt.figure(figsize=1.5*plt.figaspect(12.0 / 16.0))

    ax1 = fig1.add_subplot(211)

    # [left, bottom, width, height]
    plt1_xmin, plt1_ymin, plt1_width, plt1_height = -0.025, 0.42, 1.05, 0.30
    ax1.set_position([plt1_xmin, plt1_ymin, plt1_width, plt1_height])

    # coords = plot_tree2(tree, axe=ax1, name_offset=0.02, font_size=4)
    coords = plot_tree2(
        tree,
        axe=ax1,
        font_size=5,
        align_names=True,
        fig=fig1,
        bootstrap=args.bootstrap,
        # name_offset=-1.045,  # Fst phylogeny
        name_offset=-2 * tree_height,  # pairwise-distance phylogeny
        # name_offset=-(plt1_ymin + plt1_height) - 0.15,
        scale_line_width=0.5,
    )

    ax2 = fig1.add_subplot(
        212, projection=ccrs.PlateCarree(central_longitude=lon_centre), frame_on=True
    )
    plt2_xmin, plt2_ymin, plt2_width, plt2_height = 0, 0.005, 1, 0.4
    ax2.set_position([plt2_xmin, plt2_ymin, plt2_width, plt2_height])
    plot_map(ax2, location, fgcolours, bgcolours)

    # calling this updates some internal value(s) such that my transformations work!
    _ = ax2.get_position()

    ax3 = fig1.add_axes([0, 0, 1, 1], frame_on=False, zorder=20)
    ax3.set_xticks([])
    ax3.set_yticks([])

    lines = []
    linec = []

    ax1data_to_fig = ax1.transData + fig1.transFigure.inverted()
    ax2data_to_fig = (
        ccrs.Geodetic()._as_mpl_transform(ax2) + fig1.transFigure.inverted()
    )

    def mytrans(x, y, plt_xmin, plt_ymin, plt_width, plt_height):
        tx = plt_xmin + x * plt_width
        ty = plt_ymin + y * plt_height
        return tx, ty

    for node, (x, y) in coords.items():
        if not (node.name and node.name in tree):
            continue

        # get phylo coords
        phylox, phyloy = ax1data_to_fig.transform([x, y])

        taxid = node.taxid
        lon, lat = location[taxid]

        # get map coords
        mapx, mapy = ax2data_to_fig.transform([lon, lat])

        lines.append(((phylox, phyloy), (mapx, mapy)))
        linec.append(fgcolours[taxid])
        # ax3.scatter(phylox, phyloy, color='red', marker='x', s=50, transform=fig1.transFigure)
        # ax3.scatter(mapx, mapy, facecolor="none", edgecolor='red', marker='o', s=50, transform=fig1.transFigure)
        # print(taxid, lat, lon, phylox, phyloy, mapx, mapy)

    linecol = LineCollection(
        lines, colors=linec, lw=1, linestyle="-", alpha=1, transform=fig1.transFigure
    )
    ax3.add_collection(linecol)

    #    ax3.set_xlim([0,1])
    #    ax3.set_ylim([0,1])

    if args.heterozygosity:
        # add colourbar
        ax4 = fig1.add_axes([0, 0, 0.05, 0.4], frame_on=False, zorder=20)
        ax4.set_xticks([])
        ax4.set_yticks([])
        mappable = matplotlib.cm.ScalarMappable(norm=hetnorm, cmap=cmap)
        cbar = fig1.colorbar(
            mappable, ax=ax4, label="Heterozygosity", aspect=20, fraction=0.20
        )
        #cbar.ax.yaxis.set_label_position("left")
        het_label = cbar.ax.yaxis.get_label()
        het_label.set_rotation(270)
        cbar.ax.yaxis.labelpad = 20

    # fig1.tight_layout()
    pdf.savefig(figure=fig1)
    pdf.close()
