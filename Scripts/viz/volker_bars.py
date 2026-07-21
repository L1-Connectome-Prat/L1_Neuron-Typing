"""
Vertical Volker-class bars drawn alongside the dendrogram leaves of a
lineage clustering figure.

"""
import colorsys

import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.patches as patches

#------------------------------------------------------------------------------
__all__ = ["plot_volkers_bars"]


#------------------------------------------------------------------------------
def _label_acronym(label, get_id = False):
    """
    Pull the acronym (and optional skeleton id) out of a stitched dendrogram label.

    """
    pieces  = label.get_text().split("_")
    acronym = pieces[2]
    skel_id = pieces[-1]
    return (acronym, skel_id) if get_id else acronym


#------------------------------------------------------------------------------
def _darken_color(hex_color, saturation_factor = 2, darkness_factor = 0.9):
    """
    Boost saturation and reduce lightness so a propagated color reads as
    distinct from the original.

    """
    rgb = mcolors.to_rgb(hex_color)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    new_s = min(s * saturation_factor, 1)
    new_l = max(l * darkness_factor, 0)
    return mcolors.to_hex(colorsys.hls_to_rgb(h, new_l, new_s))


#------------------------------------------------------------------------------
def _match_skeleton_colors(labels, color_dict):
    """
    Map dendrogram skeleton ids to the bar colors assigned to their acronym.

    """
    skel_dict = {}
    for lab in labels:
        acronym, skel_id = _label_acronym(lab, get_id = True)
        skel_dict[skel_id] = color_dict.get(acronym, "grey")
    return skel_dict


#------------------------------------------------------------------------------
def plot_volkers_bars(ax1, ax2, border = "auto",
                      prop_color = False,
                      palette    = "hls"):
    """
    Draw vertical group bars along the dendrogram leaves on ax2 to match the
    acronym sequence on ax1.

    With prop_color=True the bars inherit the dendrogram leaf color (used for
    the Diff figure). Otherwise a distinct seaborn palette is drawn.

    """
    ax2.set_ylim(ax1.get_ylim())
    labels = ax1.get_yticklabels()

    label_seq = [_label_acronym(lab) for lab in labels]
    unq_groups   = np.unique(label_seq)
    group_colors = sns.color_palette(palette, len(unq_groups)).as_hex()
    group_dict   = dict(zip(unq_groups, group_colors))

    # Find consecutive runs of identical acronyms
    starts = [0]
    ends   = []
    for i, (lab1, lab2) in enumerate(zip(label_seq, label_seq[1:])):
        if lab1 != lab2:
            ends.append(i)
            starts.append(i + 1)
    ends.append(len(label_seq) - 1)
    fragments = list(zip(starts, ends))

    # Auto-set bar half-height to half the spacing between leaves
    if isinstance(border, str):
        if border.lower() == "auto":
            border = labels[0].get_position()[1] * 0.65
        else:
            raise ValueError("Enter a valid border specification.")

    if prop_color:
        prop_color_dict = {}

    for frag in fragments:
        first, last = labels[frag[0]], labels[frag[1]]
        ymin = first.get_position()[1] - border
        ymax = last.get_position()[1]  + border

        # Skip leftover fragments where there's only one leaf
        if ymax - ymin == 2 * border:
            continue

        text = _label_acronym(first)

        if prop_color:
            existing = prop_color_dict.get(text)
            used     = set(prop_color_dict.values())
            if existing is None:
                color = first.get_color()
                if color in used:
                    color = _darken_color(color, 2, 0.75)
                prop_color_dict[text] = color
            else:
                color = existing
        else:
            color = group_dict[text]

        rect = patches.Rectangle(
            (0, ymin), 0.15, ymax - ymin,
            facecolor = color,
            edgecolor = None,
            alpha     = 0.65)
        ax2.add_patch(rect)

        ax2.text(0.25, ymin + ((ymax - ymin) / 2),
                 text, va = "center",
                 fontdict = {"color" : color, "size" : 12})

    if prop_color:
        skeleton_dict = _match_skeleton_colors(labels, prop_color_dict)
        return prop_color_dict, skeleton_dict

    return group_dict
