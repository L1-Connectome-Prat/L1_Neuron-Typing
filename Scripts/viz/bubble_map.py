"""
Compartment "bubble" connectivity map drawn alongside each lineage clustering
figure. Nodes are L1 compartments at fixed canonical positions; edges are
within-Volker-class connections, optionally split into dendrite-onwards and
axo-axonic.

Requires a precomputed node layout pickle and either the full or
axon/dendrite-split innervation table from the innervation_acronyms scripts.

"""
import itertools as it
import pickle as pkl

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
__all__ = ["bubble_map", "load_bubble_assets"]


#------------------------------------------------------------------------------
def load_bubble_assets(node_dict_path,
                       paired_full_path,
                       paired_split_path):
    """
    Load the node-layout pickle and both innervation tables used by bubble_map.

    """
    with open(node_dict_path, "rb") as f:
        node_dict = pkl.load(f)

    paired_full   = pd.read_parquet(paired_full_path)
    paired_pieces = pd.read_parquet(paired_split_path)

    # Drop zero-weight rows once - the figures only ever plot positive weights
    paired_full   = paired_full.loc[paired_full["Scaled_Weight"]     > 0].reset_index(drop = True)
    paired_pieces = paired_pieces.loc[paired_pieces["Scaled_Weight"] > 0].reset_index(drop = True)

    return node_dict, paired_full, paired_pieces


#------------------------------------------------------------------------------
def _bubble_layout(G, node_dict, size_scale = 12500):
    """
    Compartment positions and bubble sizes from the precomputed node dict.

    """
    pos = dict(node_dict["position"])
    vol = node_dict["volume"]

    for node in pos.keys():
        G.add_node(node)

    node_sizes = [vol.get(node, 0.001) * size_scale for node in G.nodes()]
    return pos, vol, node_sizes


#------------------------------------------------------------------------------
def _place_connections(G, conn_df, color_dict):
    """
    Add weighted, color-tagged edges from a per-acronym connection table.

    """
    for acronym in np.unique(conn_df["Acronym"]):
        sub_df = conn_df.loc[conn_df["Acronym"] == acronym]
        color  = color_dict.get(acronym, "grey")
        for _, row in sub_df.iterrows():
            G.add_edge(row.Source_Compartment.removeprefix("PK_"),
                       row.Target_Compartment.removeprefix("PK_"),
                       weight = row.Scaled_Weight,
                       color  = color)


#------------------------------------------------------------------------------
def _threshold(conn_df, method, vis_threshold):
    """
    Apply one of three connection-thresholding strategies before plotting.

    """
    if method == "percent":
        cutoff = np.percentile(conn_df["Scaled_Weight"].values, 100 - vis_threshold)
        return conn_df.loc[conn_df["Scaled_Weight"] > cutoff]

    if method == "top":
        return conn_df.sort_values(by = "Scaled_Weight", ascending = False).head(vis_threshold)

    if method == "raw":
        return conn_df.loc[conn_df["Scaled_Weight"] > vis_threshold].reset_index(drop = True)

    raise ValueError("thresh_method must be one of: percent, top, raw.")


#------------------------------------------------------------------------------
def _isolate_axo_axonic(conn_df):
    """
    Split off axo-axonic connections; keep dendrite-onwards as the main set.

    """
    axo_df = conn_df.loc[
        (conn_df["Source_Piece"] == "axon") & (conn_df["Target_Piece"] == "axon")
        ].reset_index(drop = True)
    main_df = conn_df.loc[conn_df["Source_Piece"] != "axon"].reset_index(drop = True)
    return main_df, axo_df


#------------------------------------------------------------------------------
def _scale_widths(G, method = "sqrt", extra_widths = None):
    """
    Map raw edge weights to drawn line widths. Z-scoring optionally pools
    against another set so two passes share a scale.

    """
    raw_widths = np.array([data.get("weight", 0.1)
                           for _, _, data in G.edges(data = True)])

    if isinstance(extra_widths, pd.DataFrame):
        extra_widths = np.array(extra_widths["Scaled_Weight"])

    if method == "log":
        return np.log1p(raw_widths)

    if method == "sqrt":
        return np.sqrt(raw_widths)

    if method == "z":
        if extra_widths is not None:
            merged  = np.concatenate((raw_widths, np.array(extra_widths)))
            average = np.mean(merged)
            st_dev  = np.std(merged)
        else:
            average = np.mean(raw_widths)
            st_dev  = np.std(raw_widths)
        return (raw_widths - average) / st_dev

    return raw_widths


#------------------------------------------------------------------------------
def _build_block(conn_df, color_dict, node_dict):
    """
    Assemble the per-pass MultiDiGraph used by bubble_map.

    """
    G = nx.empty_graph(n = 0, create_using = nx.MultiDiGraph())
    pos, vol, node_sizes = _bubble_layout(G, node_dict)
    _place_connections(G, conn_df, color_dict)
    edge_colors = [data.get("color", "grey") for _, _, data in G.edges(data = True)]
    return G, pos, vol, node_sizes, edge_colors


#------------------------------------------------------------------------------
def bubble_map(annot, color_dict, node_dict,
               paired_full, paired_pieces,
               ax                   = None,
               hemisphere           = "Right",
               width_scaling_method = "sqrt",
               width_plotting_scale = 8,
               thresh_method        = "percent",
               vis_threshold        = 11,
               axon_dendrite        = False):
    """
    Compartment connectivity bubble for a single lineage.

    color_dict maps acronym -> color (typically the dendrogram-derived
    Volker palette). axon_dendrite=True draws the dendrite/axon-split table
    plus dashed axo-axonic overlays.

    """
    if width_scaling_method.lower() not in ("log", "sqrt", "z"):
        raise ValueError("width_scaling_method must be one of: log, sqrt, z.")

    created = False
    if ax is None:
        fig, ax = plt.subplots(figsize = (8, 8))
        created = True

    paired_df = paired_pieces if axon_dendrite else paired_full
    axo_df    = None

    # Per-lineage / hemisphere isolation - typing figures only show one side
    conn_df = paired_df.loc[
        (paired_df["Lineage"]    == annot) &
        (paired_df["Hemisphere"] == hemisphere)
        ].reset_index(drop = True)
    conn_df = conn_df.fillna(0)

    conn_df = _threshold(conn_df, thresh_method, vis_threshold)

    if axon_dendrite:
        conn_df, axo_df = _isolate_axo_axonic(conn_df)

    G, pos, vol, node_sizes, edge_colors = _build_block(conn_df, color_dict, node_dict)
    edge_widths = _scale_widths(G, width_scaling_method, axo_df) * width_plotting_scale

    nx.draw_networkx_nodes(G, pos,
                           node_color = "#808080bf",
                           node_size  = node_sizes,
                           ax         = ax)

    connectionstyles = [f"arc3,rad={r}" for r in it.accumulate([0.175] * 5)]

    nx.draw_networkx_edges(G, pos,
                           edge_color      = edge_colors,
                           width           = edge_widths,
                           node_size       = node_sizes,
                           connectionstyle = connectionstyles,
                           ax              = ax)

    # Only label the left-hemisphere bubbles to avoid duplication
    labels = {node : node.split("_")[0] for node in G.nodes() if "Left" in node}
    nx.draw_networkx_labels(G, pos, labels, ax = ax, font_size = 5)

    # Axo-axonic overlay (dashed)
    if axon_dendrite and len(axo_df) > 0:
        G, pos, vol, node_sizes, edge_colors = _build_block(axo_df, color_dict, node_dict)
        edge_widths = _scale_widths(G, width_scaling_method, conn_df) * width_plotting_scale

        edges = nx.draw_networkx_edges(G, pos,
                                       edge_color      = edge_colors,
                                       width           = edge_widths,
                                       node_size       = node_sizes,
                                       connectionstyle = connectionstyles,
                                       style           = "dashed",
                                       ax              = ax)
        for edge in edges:
            edge.set_linestyle((0, (10, 10)))

    if created:
        return fig, ax
