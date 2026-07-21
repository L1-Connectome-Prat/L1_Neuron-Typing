"""
Standardized 2D anatomy and synapse-density panels used by the typing
figures. Consistent axes limits keep panels comparable across lineages.

"""
import numpy as np
import navis

import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
__all__ = ["plot_lineage_anatomy",
           "synapse_plot"]

# Global axis limits in PK_L1CNSsym coordinates
ANTERIOR_XLIM = (-1000, 105000)
ANTERIOR_YLIM = (93722, -2737)
DORSAL_XLIM   = (-1000, 105000)
DORSAL_YLIM   = (4000, 100459)


#------------------------------------------------------------------------------
def _format_axes(ax, ylim, xlim, ylabel, xlabel = None):
    """
    Common axis cleanup.

    """
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(visible = False)
    ax.set_ylabel(ylabel, fontsize = 10)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize = 10)
    else:
        ax.set_xlabel("")


#------------------------------------------------------------------------------
def plot_lineage_anatomy(neurons, volumes, ax, **decor_kwargs):
    """
    Anterior + dorsal projections of a neuron set in two stacked axes.

    """
    defaults = dict(color     = "green",
                    alpha     = 0.75,
                    method    = "2d",
                    rasterize = True)
    kwargs = {**defaults, **decor_kwargs}

    navis.plot2d([neurons, volumes], ax = ax[0], view = ("x", "-y"), **kwargs)
    navis.plot2d([neurons, volumes], ax = ax[1], view = ("x",  "z"), **kwargs)

    _format_axes(ax[0], ANTERIOR_YLIM, ANTERIOR_XLIM,
                 r"$\longleftarrow$ Ventral | Dorsal $\longrightarrow$")
    _format_axes(ax[1], DORSAL_YLIM, DORSAL_XLIM,
                 r"$\longleftarrow$ Anterior | Posterior $\longrightarrow$",
                 r"$\longleftarrow$ Right | Left $\longrightarrow$")


#------------------------------------------------------------------------------
def synapse_plot(neurons, neuropil, ax,
                 cmap  = {"tbars" : "red", "psds" : "green"},
                 alpha = 0.1,
                 s     = 2):
    """
    PSD (top row) and TBar (bottom row) scatter overlay on the neuropil
    outline. ax is a 2x2 grid: rows = PSDs/TBars, cols = anterior/dorsal.

    """
    # Connectors as arrays
    tbars = neurons.connectors[neurons.connectors.type == 0][["x", "y", "z"]].values
    psds  = neurons.connectors[neurons.connectors.type == 1][["x", "y", "z"]].values

    # Neuropil outline behind every panel
    for row in range(2):
        navis.plot2d([neuropil], method = "2d", view = ("x", "-y"),
                     ax = ax[row][0], volume_outlines = True,
                     alpha = 1, color = "black", rasterize = True)
        navis.plot2d([neuropil], method = "2d", view = ("x", "z"),
                     ax = ax[row][1], volume_outlines = True,
                     alpha = 1, color = "black", rasterize = True)

    # PSDs (top row)
    ax[0][0].scatter(psds[:, 0], psds[:, 1], s = s, alpha = alpha, color = cmap["psds"])
    ax[0][1].scatter(psds[:, 0], psds[:, 2], s = s, alpha = alpha, color = cmap["psds"])

    # TBars (bottom row)
    ax[1][0].scatter(tbars[:, 0], tbars[:, 1], s = s, alpha = alpha, color = cmap["tbars"])
    ax[1][1].scatter(tbars[:, 0], tbars[:, 2], s = s, alpha = alpha, color = cmap["tbars"])

    # Apply consistent limits
    for row in range(2):
        _format_axes(ax[row][0], ANTERIOR_YLIM, ANTERIOR_XLIM, "")
        _format_axes(ax[row][1], DORSAL_YLIM,   DORSAL_XLIM,   "")

    # Titles + labels on the outer axes only
    ax[0][0].set_title("Inputs (PSDs)",  fontsize = 15)
    ax[1][0].set_title("Outputs (T-Bars)", fontsize = 15)

    ax[0][1].set_xlabel(r"$\longleftarrow$ Right | Left $\longrightarrow$",
                        fontsize = 10, labelpad = 11)
    ax[1][1].set_xlabel(r"$\longleftarrow$ Right | Left $\longrightarrow$",
                        fontsize = 10, labelpad = 11)
