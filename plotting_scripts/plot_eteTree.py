## Modified from:
## https://gist.github.com/fransua/da703c3d2ba121903c0de5e976838b71
from itertools import chain

from matplotlib.collections import LineCollection
from matplotlib import markers
from matplotlib.path import Path

import numpy as np


def round_sig(x, sig=2):
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


def to_coord(x, y, xmin, xmax, ymin, ymax, plt_xmin, plt_ymin, plt_width, plt_height):
    x = (x - xmin) / (xmax - xmin) * plt_width + plt_xmin
    y = (y - ymin) / (ymax - ymin) * plt_height + plt_ymin
    return x, y


def plot_tree(
    tree,
    align_names=False,
    align_colour="red",
    name_offset=None,
    max_dist=None,
    font_size=9,
    axe=None,
    bootstrap=False,
    boot_colour="red",
    scale_line_width=1,
):
    """
    Plots a ete3.Tree object using matploltib.

    :param tree: ete Tree object
    :param False align_names: if True names will be aligned vertically
    :param None max_dist: if defined any branch longer than the given value will be
       reduced by this same value.
    :param None name_offset: offset relative to tips to write leaf_names. In bL scale
    :param 12 font_size: to write text
    :param None axe: a matploltib.Axes object on which the tree will be painted.

    :returns: a dictionary of node objects with their coordinates
    """

    if axe is None:
        axe = plt.subplot(111)

    def __draw_edge_nm(c, x):
        h = node_pos[c]
        hlinec.append(((x, h), (x + c.dist, h)))
        hlines.append(cstyle)
        return (x + c.dist, h)

    def __draw_edge_md(c, x):
        h = node_pos[c]
        if c in cut_edge:
            offset = max_x / 600.0
            hlinec.append(((x, h), (x + c.dist / 2 - offset, h)))
            hlines.append(cstyle)
            hlinec.append(((x + c.dist / 2 + offset, h), (x + c.dist, h)))
            hlines.append(cstyle)
            hlinec.append(
                ((x + c.dist / 2, h - 0.05), (x + c.dist / 2 - 2 * offset, h + 0.05))
            )
            hlines.append(cstyle)
            hlinec.append(
                ((x + c.dist / 2 + 2 * offset, h - 0.05), (x + c.dist / 2, h + 0.05))
            )
            hlines.append(cstyle)
            axe.text(
                x + c.dist / 2,
                h - 0.07,
                "+%g" % max_dist,
                va="top",
                ha="center",
                size=2.0 * font_size / 3,
            )
        else:
            hlinec.append(((x, h), (x + c.dist, h)))
            hlines.append(cstyle)
        return (x + c.dist, h)

    __draw_edge = __draw_edge_nm if max_dist is None else __draw_edge_md

    vlinec = []
    vlines = []
    hlinec = []
    hlines = []
    nodes = []
    nodex = []
    nodey = []
    ali_lines = []
    boot_lines = []

    # to align leaf names
    tree = tree.copy()
    max_x = max(n.get_distance(tree) for n in tree.iter_leaves())

    coords = {}
    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))
    node_list = tree.iter_descendants(strategy="postorder")
    node_list = chain(node_list, [tree])

    # reduce branch length
    cut_edge = set()
    if max_dist is not None:
        for n in tree.iter_descendants():
            if n.dist > max_dist:
                n.dist -= max_dist
                cut_edge.add(n)

    if name_offset is None:
        name_offset = max_x / 25.0
    # draw tree
    for n in node_list:
        style = n._get_style()
        x = np.sum(n2.dist for n2 in n.iter_ancestors()) + n.dist
        if n.is_leaf():
            y = node_pos[n]
            if align_names:
                axe.text(
                    max_x + name_offset,
                    y,
                    n.name,
                    va="center",
                    size=font_size,
                    c=style["bgcolor"],
                    fontweight="bold",
                    fontstretch="condensed",
                )
                text_x = max_x
                if style["shape"] is None:
                    text_x += name_offset
                ali_lines.append(((x, y), (text_x , y)))
            else:
                axe.text(
                    x + name_offset,
                    y,
                    n.name,
                    va="center",
                    size=font_size,
                    c=style["bgcolor"],
                    fontweight="bold",
                    fontstretch="condensed",
                )
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            # draw vertical line
            vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
            vlines.append(style)

            # draw horizontal lines
            for child in n.children:
                cstyle = child._get_style()
                coords[child] = __draw_edge(child, x)

            if bootstrap and n.support < 100:
                text_x = 0.082
                axe.text(
                    text_x,
                    y,
                    f"{int(n.support)}",
                    va="center",
                    size=font_size - 1,
                    c=style["bgcolor"],
                    fontweight="bold",
                    fontstretch="condensed",
                )
                boot_lines.append(((text_x + name_offset / 2, y), (x, y)))
            else:
                # don't draw shapes for internal nodes
                continue
        nodes.append(style)
        if align_names and n.is_leaf():
            nodex.append(max_x + name_offset/2)
        else:
            nodex.append(x)
        nodey.append(y)

    # draw root
    __draw_edge(tree, 0)

    lstyles = ["-", "--", ":"]
    hline_col = LineCollection(
        hlinec,
        colors=[l["hz_line_color"] for l in hlines],
        linestyle=[lstyles[l["hz_line_type"]] for l in hlines],
        linewidth=[(l["hz_line_width"] + 0.1) / 2 for l in hlines],
    )
    vline_col = LineCollection(
        vlinec,
        colors=[l["vt_line_color"] for l in vlines],
        linestyle=[lstyles[l["vt_line_type"]] for l in vlines],
        linewidth=[(l["vt_line_width"] + 0.1) / 2 for l in vlines],
    )
    ali_line_col = LineCollection(ali_lines, colors=align_colour, linestyle=":", linewidth=0.3)
    boot_line_col = LineCollection(boot_lines, colors=boot_colour, linestyle=":", linewidth=0.3)

    axe.add_collection(hline_col)
    axe.add_collection(vline_col)
    axe.add_collection(ali_line_col)
    axe.add_collection(boot_line_col)

    nshapes = dict((("circle", "o"), ("square", "s"), ("sphere", "o")))
    shapes = set(n["shape"] for n in nodes if n["shape"] is not None)
    for shape in shapes:
        indexes = [i for i, n in enumerate(nodes) if n["shape"] == shape]
        scat = axe.scatter(
            [nodex[i] for i in indexes],
            [nodey[i] for i in indexes],
            s=0,
            marker=nshapes.get(shape, shape),
        )
        scat.set_sizes([(nodes[i]["size"]) ** 2 / 2 for i in indexes])
        scat.set_facecolor([nodes[i]["fgcolor"] for i in indexes])
        # scat.set_edgecolor([nodes[i]['bgcolor'] for i in indexes])
        # scat.set_linewidth(.5)
        scat.set_zorder(10)

    # scale line
    xmin, xmax = axe.get_xlim()
    ymin, ymax = axe.get_ylim()
    diffy = ymax - ymin
    dist = round_sig((xmax - xmin) / 5, sig=1)
    ymin += diffy / 18.0
    xmin += dist / 3
    axe.plot([xmin, xmin + dist], [ymin, ymin], color="k", linewidth=scale_line_width)
    axe.plot([xmin, xmin], [ymin - diffy / 200.0, ymin + diffy / 200.0], color="k", linewidth=scale_line_width)
    axe.plot(
        [xmin + dist, xmin + dist],
        [ymin - diffy / 200.0, ymin + diffy / 200.0],
        color="k",
        linewidth=scale_line_width,
    )
    axe.text(
        (xmin + xmin + dist) / 2,
        ymin - diffy / 200.0,
        dist,
        va="top",
        ha="center",
        size=font_size,
    )
    axe.set_axis_off()
    return coords


def plot_tree2(
    tree,
    align_names=False,
    name_offset=None,
    max_dist=None,
    font_size=9,
    axe=None,
    rotate=True,
    vflip=True,
    hflip=True,
    text=True,
    fig=None,
    bootstrap=False,
    scale_line_width=1,
):
    """
    Plots a ete3.Tree object using matploltib.

    :param tree: ete Tree object
    :param False align_names: if True names will be aligned vertically
    :param None max_dist: if defined any branch longer than the given value will be
       reduced by this same value.
    :param None name_offset: offset relative to tips to write leaf_names. In bL scale
    :param 12 font_size: to write text
    :param None axe: a matploltib.Axes object on which the tree will be painted.
    :param kwargs: for tree edge drawing (matplotlib LineCollection)

    :returns: a dictionary of node objects with their coordinates
    """

    def trcoords(x, y):
        if rotate:
            x, y = y, x
        if vflip:
            y = -y
        if hflip:
            x = -x
        return x, y

    if rotate:
        textrot = -90
        textha = "center"
        if vflip:
            textva = "top"
        else:
            textva = "bottom"
    else:
        textrot = "horizontal"
        textva = "center"
        if hflip:
            textha = "right"
        else:
            textha = "left"

    if axe is None:
        axe = plt.subplot(111)

    def __draw_edge_nm(c, x):
        h = node_pos[c]
        hlinec.append(((x, h), (x + c.dist, h)))
        hlines.append(cstyle)
        return (x + c.dist, h)

    def __draw_edge_md(c, x):
        h = node_pos[c]
        if c in cut_edge:
            offset = max_x / 600.0
            hlinec.append(((x, h), (x + c.dist / 2 - offset, h)))
            hlines.append(cstyle)
            hlinec.append(((x + c.dist / 2 + offset, h), (x + c.dist, h)))
            hlines.append(cstyle)
            hlinec.append(
                ((x + c.dist / 2, h - 0.05), (x + c.dist / 2 - 2 * offset, h + 0.05))
            )
            hlines.append(cstyle)
            hlinec.append(
                ((x + c.dist / 2 + 2 * offset, h - 0.05), (x + c.dist / 2, h + 0.05))
            )
            hlines.append(cstyle)
            if text:
                axe.text(
                    *trcoords(x + c.dist / 2, h - 0.07),
                    "+%g" % max_dist,
                    va=textva,
                    ha="center",
                    rotation=textrot,
                    size=2.0 * font_size / 3,
                )
        else:
            hlinec.append(((x, h), (x + c.dist, h)))
            hlines.append(cstyle)
        return (x + c.dist, h)

    __draw_edge = __draw_edge_nm if max_dist is None else __draw_edge_md

    vlinec = []
    vlines = []
    hlinec = []
    hlines = []
    nodes = []
    nodex = []
    nodey = []
    ali_lines = []

    # to align leaf names
    tree = tree.copy()
    max_x = max(n.get_distance(tree) for n in tree.iter_leaves())

    coords = {}
    node_pos = dict((n2, i) for i, n2 in enumerate(tree.get_leaves()[::-1]))
    node_list = tree.iter_descendants(strategy="postorder")
    node_list = chain(node_list, [tree])

    inv = axe.transData.inverted()

    # reduce branch length
    cut_edge = set()
    if max_dist is not None:
        for n in tree.iter_descendants():
            if n.dist > max_dist:
                n.dist -= max_dist
                cut_edge.add(n)

    if name_offset is None:
        name_offset = max_x / 100.0
    # draw tree
    for n in node_list:
        style = n._get_style()
        x = np.sum(n2.dist for n2 in n.iter_ancestors()) + n.dist
        if n.is_leaf():
            y = node_pos[n]
            if text:
                if align_names:
                    t = axe.text(
                        *trcoords(max_x + name_offset, y),
                        n.name,
                        # va=textva, ha="right", rotation=textrot, size=font_size, fontweight="bold", fontstretch="condensed")
                        va=textva,
                        ha=textha,
                        rotation=textrot,
                        size=font_size,
                        fontweight="bold",
                        fontstretch="condensed",
                    )
                    if name_offset >= 0:
                        # bb = t.get_window_extent(fig.canvas.get_renderer())
                        # tx, ty = inv.transform([bb.x0, bb.height])
                        # ali_lines.append(((x, y), (max_x + name_offset +tx, y))) # How to get correct `tx'?
                        ali_lines.append(((x, y), (max_x + name_offset, y)))
                else:
                    # c=style["bgcolor"]
                    axe.text(
                        *trcoords(x + name_offset, y),
                        n.name,
                        va=textva,
                        ha=textha,
                        rotation=textrot,
                        size=font_size,
                        fontweight="bold",
                        fontstretch="condensed",
                    )
        else:
            y = np.mean([node_pos[n2] for n2 in n.children])
            node_pos[n] = y

            # draw vertical line
            vlinec.append(((x, node_pos[n.children[0]]), (x, node_pos[n.children[-1]])))
            vlines.append(style)

            # draw horizontal lines
            for child in n.children:
                cstyle = child._get_style()
                coords[child] = __draw_edge(child, x)

            if bootstrap and n.support < 100:
                axe.text(
                    *trcoords(0, y),
                    f"{int(n.support)}",
                    va=textva,
                    ha=textha,
                    rotation=textrot,
                    size=font_size - 1,
                    c=style["bgcolor"],
                    fontweight="bold",
                    fontstretch="condensed",
                )
                ali_lines.append(((0.05 * max_x, y), (x, y)))
            else:
                # don't draw shapes for internal nodes
                continue
        nodes.append(style)
        nodex.append(x)
        nodey.append(y)

    # draw root
    __draw_edge(tree, 0)

    hlinec = [(trcoords(*hc[0]), trcoords(*hc[1])) for hc in hlinec]
    vlinec = [(trcoords(*vc[0]), trcoords(*vc[1])) for vc in vlinec]
    ali_lines = [(trcoords(*al[0]), trcoords(*al[1])) for al in ali_lines]
    nodex, nodey = zip(*[trcoords(x, y) for x, y in zip(nodex, nodey)])

    lstyles = ["-", "--", ":"]
    hline_col = LineCollection(
        hlinec,
        colors=[l["hz_line_color"] for l in hlines],
        linestyle=[lstyles[l["hz_line_type"]] for l in hlines],
        linewidth=[(l["hz_line_width"] + 1.0) / 2 for l in hlines],
    )
    vline_col = LineCollection(
        vlinec,
        colors=[l["vt_line_color"] for l in vlines],
        linestyle=[lstyles[l["vt_line_type"]] for l in vlines],
        linewidth=[(l["vt_line_width"] + 1.0) / 2 for l in vlines],
    )
    ali_line_col = LineCollection(
        ali_lines, colors="darkgray", linestyle=":", linewidth=0.5
    )

    axe.add_collection(hline_col)
    axe.add_collection(vline_col)
    axe.add_collection(ali_line_col)

    nshapes = dict((("circle", "o"), ("square", "s"), ("sphere", "o")))
    shapes = set(n["shape"] for n in nodes)
    for shape in shapes:
        indexes = [i for i, n in enumerate(nodes) if n["shape"] == shape]
        scat = axe.scatter(
            [nodex[i] for i in indexes],
            [nodey[i] for i in indexes],
            s=0,
            marker=nshapes.get(shape, shape),
        )
        scat.set_sizes([(nodes[i]["size"]) ** 2 / 2 for i in indexes])
        scat.set_facecolor([nodes[i]["fgcolor"] for i in indexes])
        scat.set_edgecolor([nodes[i]["bgcolor"] for i in indexes])
        scat.set_linewidth(0.5)
        scat.set_zorder(10)

    # scale line
    if rotate:
        xmin, xmax = axe.get_xlim()
        ymin, ymax = axe.get_ylim()
        diffx = xmax - xmin
        dist = round_sig((ymax - ymin) / 5, sig=1)
        xmin += diffx / 40.
        ymax -= dist
        axe.plot([xmin, xmin], [ymax, ymax - dist], color="k", linewidth=scale_line_width)  # line
        axe.plot(
            [xmin - diffx / 200.0, xmin + diffx / 200.0], [ymax, ymax], color="k", linewidth=scale_line_width
        )  # whisker 1
        axe.plot(
            [xmin - diffx / 200.0, xmin + diffx / 200.0],
            [ymax - dist, ymax - dist],  # whisker 2
            color="k",
            linewidth=scale_line_width,
        )
        axe.text(
            xmin - diffx / 200,
            (ymax + ymax - dist) / 2,
            dist,
            va="center",
            ha="right",
            rotation=textrot,
            size=font_size,
        )
    else:
        xmin, xmax = axe.get_xlim()
        ymin, ymax = axe.get_ylim()
        diffy = ymax - ymin
        dist = round_sig((xmax - xmin) / 5, sig=1)
        xmin += (xmax - xmin) / 10.0
        ymin -= diffy / 100.0
        axe.plot([xmin, xmin + dist], [ymin, ymin], color="k", linewidth=scale_line_width)
        axe.plot([xmin, xmin], [ymin - diffy / 200.0, ymin + diffy / 200.0], color="k", linewidth=scale_line_width)
        axe.plot(
            [xmin + dist, xmin + dist],
            [ymin - diffy / 200.0, ymin + diffy / 200.0],
            color="k",
            linewidth=scale_line_width,
        )
        axe.text(
            (xmin + xmin + dist) / 2,
            ymin - diffy / 200.0,
            dist,
            va="top",
            ha="center",
            rotation=textrot,
            size=font_size,
        )
    axe.set_axis_off()

    coords2 = {k: trcoords(*v) for k, v in coords.items()}
    return coords2
