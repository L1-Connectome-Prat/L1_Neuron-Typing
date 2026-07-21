"""
Figure scaffolds for the per-lineage Diff and SU clustering panels.

Layouts are tuned by-hand to fit the manuscript figure book and assume
specific gridspec slots - changing the dimensions cascades into many
hard-coded position offsets.

"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from prats_helpers import place_panel_label

#------------------------------------------------------------------------------
__all__ = ["build_diff_layout", "build_su_layout", "build_sensory_layout"]


#------------------------------------------------------------------------------
def build_diff_layout(annot):
    """
    Per-lineage Diff figure layout - 5 anatomy columns, dendrogram + Volker
    bars, bubble map, and a 6x4 cluster anatomy grid.

    """
    fig = plt.figure(figsize = (22, 45))
    fig.subplots_adjust(hspace = 1, wspace = 1)
    fig.suptitle(f"{annot}", fontsize = 28, y = 0.91)

    gs = gridspec.GridSpec(45, 22)

    # Top anatomy strip - five columns of stacked anterior/dorsal pairs
    ax_a = [fig.add_subplot(gs[0:4, 0:4]),   fig.add_subplot(gs[4:8, 0:4])]
    ax_b = [fig.add_subplot(gs[0:4, 4:8]),   fig.add_subplot(gs[4:8, 4:8])]
    ax_c = [fig.add_subplot(gs[0:4, 8:12]),  fig.add_subplot(gs[4:8, 8:12])]
    ax_d = [fig.add_subplot(gs[0:4, 12:16]), fig.add_subplot(gs[4:8, 12:16])]
    ax_h = [fig.add_subplot(gs[0:4, 16:20]), fig.add_subplot(gs[4:8, 16:20])]

    # Dendrogram + Volker bars
    ax_e1 = fig.add_subplot(gs[9:17, 0:8])
    ax_e2 = fig.add_subplot(gs[9:17, 10:11])
    ax_e2.set_xticks([])
    ax_e2.set_yticks([])
    ax_e2.set_facecolor("none")
    for spine in ax_e2.spines.values():
        spine.set_visible(False)
    ax_e = [ax_e1, ax_e2]

    ax_f = fig.add_subplot(gs[9:17, 14:22])

    # 6x4 cluster anatomy grid - rows alternate anterior / dorsal
    sp = 3
    ax_g = [[fig.add_subplot(gs[j:j+sp, i:i+sp]) for i in range(0, 21, sp)]
            for j in range(18, 42, sp)]
    ax_g1 = ax_g[0::2]
    ax_g2 = ax_g[1::2]

    # Align top anatomy strip to the bubble map's right edge
    target_x1 = ax_f.get_position().x1
    last_x1   = ax_h[0].get_position().x1
    sh = (target_x1 - last_x1) / 4
    for i, gr in enumerate([ax_a, ax_b, ax_c, ax_d, ax_h]):
        for pan in gr:
            pos = pan.get_position()
            pan.set_position([pos.x0 + (i * sh), pos.y0, pos.width, pos.height])

    # Same alignment for the cluster anatomy grid
    last_x1 = ax_g1[0][-1].get_position().x1
    sh = (target_x1 - last_x1) / (len(ax_g1[0]) - 1)
    for orient in [ax_g1, ax_g2]:
        for row in orient:
            for i, col in enumerate(row):
                pos = col.get_position()
                col.set_position([pos.x0 + (i * sh), pos.y0, pos.width, pos.height])

    # Progressive vertical spacers between cluster-anatomy rows
    for i, (top, bottom) in enumerate(zip(ax_g1, ax_g2)):
        if i in [1, 2, 3, 4, 5, 6]:
            shift = 0.015 * i
            for a in top + bottom:
                pos = a.get_position()
                a.set_position([pos.x0, pos.y0 - shift, pos.width, pos.height])

    # Panel labels
    place_panel_label(fig, ax_a[0], "Ai")
    place_panel_label(fig, ax_a[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_b[0], "Bi")
    place_panel_label(fig, ax_b[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_c[0], "Ci")
    place_panel_label(fig, ax_c[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_d[0], "Di")
    place_panel_label(fig, ax_d[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_h[0], "Ei")
    place_panel_label(fig, ax_h[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_e1,   "F")
    place_panel_label(fig, ax_f,    "G")
    place_panel_label(fig, ax_g[0][0], "H", shy = 0.012)

    return fig, ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g1, ax_g2, ax_h


#------------------------------------------------------------------------------
def build_su_layout(annot):
    """
    Per-lineage SU figure layout - one anatomy column, dendrogram, and a
    smaller cluster anatomy grid.

    """
    fig = plt.figure(figsize = (15, 25))
    fig.subplots_adjust(hspace = 1, wspace = 1)
    fig.suptitle(f"{annot}", fontsize = 28, y = 0.93)

    gs = gridspec.GridSpec(25, 15)

    ax_a = [fig.add_subplot(gs[0:4, 0:4]), fig.add_subplot(gs[4:8, 0:4])]
    ax_b = fig.add_subplot(gs[0:8, 5:13])

    sp = 3
    ax_g = [[fig.add_subplot(gs[j:j+sp, i:i+sp]) for i in range(0, 14, sp)]
            for j in range(10, 25, sp)]
    ax_g1 = ax_g[0::2]
    ax_g2 = ax_g[1::2]

    for i, (top, bottom) in enumerate(zip(ax_g1, ax_g2)):
        if i in [1, 2, 3, 4, 5, 6]:
            shift = 0.015 * i
            for a in top + bottom:
                pos = a.get_position()
                a.set_position([pos.x0, pos.y0 - shift, pos.width, pos.height])

    return fig, ax_a, ax_b, ax_g1, ax_g2


#------------------------------------------------------------------------------
def build_sensory_layout(annot):
    """
    Per-annotation sensory figure layout - three anatomy columns (clone, PSDs,
    T-Bars), dendrogram with Volker bars, and a cluster anatomy grid. No
    tract overlay, no SU panel, and no bubble map.

    """
    fig = plt.figure(figsize = (15, 30))
    fig.subplots_adjust(hspace = 1, wspace = 1)
    fig.suptitle(f"{annot}", fontsize = 28, y = 0.91)

    gs = gridspec.GridSpec(30, 15)

    # Top anatomy strip - three columns of stacked anterior/dorsal pairs
    ax_a = [fig.add_subplot(gs[0:4, 0:4]),  fig.add_subplot(gs[4:8, 0:4])]
    ax_b = [fig.add_subplot(gs[0:4, 4:8]),  fig.add_subplot(gs[4:8, 4:8])]
    ax_c = [fig.add_subplot(gs[0:4, 8:12]), fig.add_subplot(gs[4:8, 8:12])]

    # Dendrogram + Volker bars
    ax_d1 = fig.add_subplot(gs[9:17, 0:8])
    ax_d2 = fig.add_subplot(gs[9:17, 9:10])
    ax_d2.set_xticks([])
    ax_d2.set_yticks([])
    ax_d2.set_facecolor("none")
    for spine in ax_d2.spines.values():
        spine.set_visible(False)
    ax_d = [ax_d1, ax_d2]

    # Cluster anatomy grid - rows alternate anterior / dorsal
    sp = 3
    ax_g = [[fig.add_subplot(gs[j:j+sp, i:i+sp]) for i in range(0, 15, sp)]
            for j in range(18, 30, sp)]
    ax_g1 = ax_g[0::2]
    ax_g2 = ax_g[1::2]

    # Progressive vertical spacers between cluster-anatomy rows
    for i, (top, bottom) in enumerate(zip(ax_g1, ax_g2)):
        if i in [1, 2, 3]:
            shift = 0.015 * i
            for a in top + bottom:
                pos = a.get_position()
                a.set_position([pos.x0, pos.y0 - shift, pos.width, pos.height])

    # Panel labels
    place_panel_label(fig, ax_a[0], "Ai")
    place_panel_label(fig, ax_a[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_b[0], "Bi")
    place_panel_label(fig, ax_b[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_c[0], "Ci")
    place_panel_label(fig, ax_c[1], "ii", shy = 0.004)
    place_panel_label(fig, ax_d1,   "D")
    place_panel_label(fig, ax_g[0][0], "E", shy = 0.012)

    return fig, ax_a, ax_b, ax_c, ax_d, ax_g1, ax_g2
