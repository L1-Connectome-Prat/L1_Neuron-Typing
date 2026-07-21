"""
Per-lineage morphology-based typing of differentiated (Diff) L1 neurons via
NBLAST + Ward clustering. Produces a multi-panel figure per lineage and
appends each lineage's cluster assignments to a shared parquet table.

SU neurons are loaded only as visual context. Lineages with no SU neurons
just leave that panel blank.

"""
import gc
import os
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import plotly
import seaborn as sns
sns.set_style("whitegrid")

# Connectome
import navis
import flybrains
from navis.nbl.smat import Lookup2d

# Analysis package
import connectome_analysis_claude as ca

# Local helpers
from Scripts.viz._paths          import (LINEAGES_NPY, nblast_matrix_path,
                                         NODE_BUBBLE_PKL,
                                         PAIRED_FULL_PARQUET, PAIRED_PIECES_PARQUET,
                                         MB_RIGHT_STL, MB_LEFT_STL,
                                         skeleton_folder)
from Scripts.viz.loaders         import load_lineage_neurons, load_su_neurons, load_lineage_tract
from Scripts.viz.anatomy         import plot_lineage_anatomy, synapse_plot
from Scripts.viz.volker_bars     import plot_volkers_bars
from Scripts.viz.bubble_map      import bubble_map, load_bubble_assets
from Scripts.viz.figure_layouts  import build_diff_layout
from Scripts.viz.clustering_io   import save_linkage, load_linkage

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
smat_val    = 250
dot_tangent = 12
plot_3d     = False
re_blast    = False
folder      = "Analysis_Outputs/Neuron_Clusters_Diff"

params = {"out_folder"  : folder,
          "dot_tangent" : dot_tangent,
          "smat_val"    : smat_val,
          "use_alpha"   : False,
          "score"       : "min",
          "method"      : "ward",
          "threshold"   : 0.5,
          "plot_3d"     : plot_3d,
          "re_blast"    : re_blast}

#------------------------------------------------------------------------------
# NBLAST scoring matrix
smat   = pd.read_csv(nblast_matrix_path(smat_val), index_col = 0)
smat2d = Lookup2d.from_dataframe(smat)

# Annotations to iterate over
annotations = np.load(LINEAGES_NPY)

#------------------------------------------------------------------------------
# Volumes
cns      = flybrains.PK_L1CNSsym.mesh
neuropil = flybrains.PK_L1CBNeuropilsym.mesh
cns.color      = (0.85, 0.85, 0.85, 0.2)
neuropil.color = (0, 0, 0, 1)

mbr = ca.stl_to_navis(MB_RIGHT_STL)
mbl = ca.stl_to_navis(MB_LEFT_STL)
mbr.color = (0.85, 0.85, 0.85, 0.65)
mbl.color = (0.85, 0.85, 0.85, 0.65)
mbr = navis.xform_brain(mbr, source = "PK_L1CNS", target = "PK_L1CNSsym")
mbl = navis.xform_brain(mbl, source = "PK_L1CNS", target = "PK_L1CNSsym")

volumes = [cns, mbr, mbl]

#------------------------------------------------------------------------------
# Bubble map assets - pre-loaded once
node_dict, paired_full, paired_pieces = load_bubble_assets(
    NODE_BUBBLE_PKL, PAIRED_FULL_PARQUET, PAIRED_PIECES_PARQUET)

#------------------------------------------------------------------------------
# Initialize the run folder, error log, and per-run table path
errors, error_path, table_name = ca.nblast_logging(params, format = "parquet")

# Per-lineage figure folder
lineage_pdf_dir = f"{folder}/Lineages"
os.makedirs(lineage_pdf_dir, exist_ok = True)

#------------------------------------------------------------------------------
def generate_clusters(neuron_collection, annot,
                      su_neurons = None,
                      tract      = None):
    """
    Run (or reload) the NBLAST clustering for a lineage and draw the
    multi-panel figure.

    """
    annot_name   = annot.replace("/", "")
    linkage_path = f"{folder}/Linkage/{annot_name}.h5"

    if re_blast:
        nbl_dist, nbl_vec, Z = ca.prats_nblast_pipeline(
                                    neuron_collection,
                                    dot_tangent = dot_tangent,
                                    resample    = False,
                                    smat        = smat2d,
                                    use_alpha   = False,
                                    score       = "min",
                                    n_cores     = 32,
                                    method      = "ward")

        info_dict = ca.run_multihemispheric_clustering(
                        nbl_dist, Z, neuron_collection[2],
                        t_adj     = 0.001,
                        threshold = 0.5)
        info_dict["Z"]        = Z
        info_dict["nbl_dist"] = nbl_dist
        info_dict["threshold"] = 0.5
    else:
        info_dict = load_linkage(linkage_path)
        Z        = info_dict["Z"]
        nbl_dist = info_dict["nbl_dist"]

    #--------------------------------------------------------------------------
    # Cluster assignment table
    cluster_df = ca.assemble_annotations(info_dict, neuron_collection[2], annot)

    if re_blast:
        group_dict = dict(zip(cluster_df.bodyId.values.tolist(),
                              cluster_df.prat_cluster.values.tolist()))
    else:
        group_dict = dict(zip(neuron_collection[2].id, neuron_collection[2].cluster))

    # Cluster -> stitched acronym title (depends on the rename_lr_neurons format)
    cluster_df["titles"] = cluster_df.apply(lambda row:
                                f"{row['name'].split('_')[1][:2]}_{row['name'].split('_')[2]}",
                                axis = 1)

    titles_by_key = defaultdict(set)
    for key, title in zip(cluster_df.prat_cluster, cluster_df.titles):
        titles_by_key[key].add(title)

    acronym_dict = {key : " &\n".join(sorted(s)) for key, s in titles_by_key.items()}

    #--------------------------------------------------------------------------
    # Figure scaffold
    fig, ax_a, ax_b, ax_c, ax_d, ax_e, ax_f, ax_g1, ax_g2, ax_h = build_diff_layout(annot)
    right_neurons = neuron_collection[2][neuron_collection[2].hemisphere == "Right"]

    # Lineage clone (column A)
    plot_lineage_anatomy(right_neurons, volumes, ax_a)
    ax_a[0].set_title("Lineage Clone", fontsize = 15)

    # Synapse panels (columns B, C)
    synapse_plot(right_neurons, neuropil, ax = np.vstack((ax_b, ax_c)))

    # Dendrogram with Volker bars
    ax_e1, ax_e2 = ax_e
    color_dict = ca.build_dendrogram(
                        Z, nbl_dist,
                        info_dict["optimal_t"],
                        info_dict["cluster_labels"],
                        names   = True,
                        neurons = neuron_collection[2],
                        ax      = ax_e1)

    ax_e2.set_ylim(ax_e1.get_ylim())
    ax_e1.set_title("Morphological Clustering")
    ax_e1.set_ylabel("")
    ax_e1.set_xlim(ax_e1.get_xlim()[0], 0.75)

    volker_dict, volker_skel_dict = plot_volkers_bars(ax_e1, ax_e2, prop_color = True)

    ax_e1.invert_yaxis()
    ax_e2.invert_yaxis()

    # Reposition the bars right next to the dendrogram leaves
    pos_e1    = ax_e1.get_tightbbox()
    pos_e2    = ax_e2.get_position()
    shifted_x = pos_e1.x1 + 8
    ratio     = fig.get_figwidth() * fig.dpi
    ax_e2.set_position([shifted_x / ratio, pos_e2.y0, pos_e2.width, pos_e2.height])
    ax_e1.set_xlabel("")

    # Tract (column D)
    if len(tract):
        plot_lineage_anatomy(tract, volumes, ax_d, color = volker_skel_dict)
        plot_lineage_anatomy(tract, volumes, ax_h, color = "grey", alpha = 0.25)

    ax_d[0].set_title("Lineage Tract", fontsize = 15)
    ax_d[0].set_xlabel("")
    ax_d[0].set_ylabel("")
    ax_d[1].set_ylabel("")

    # SU overlay (column H) - skipped silently if no SU neurons exist
    if su_neurons and len(su_neurons):
        plot_lineage_anatomy(su_neurons, None, ax_h, color = "blue")
        ax_h[0].set_title("SU Neurons", fontsize = 15)
        ax_h[0].set_ylabel("")
        ax_h[0].set_xlabel("")
        ax_h[1].set_ylabel("")

    #--------------------------------------------------------------------------
    # Bubble map
    bubble_map(annot, volker_dict, node_dict, paired_full, paired_pieces,
               ax                   = ax_f,
               thresh_method        = "top",
               vis_threshold        = 50,
               axon_dendrite        = True,
               width_scaling_method = "sqrt",
               width_plotting_scale = 9)

    ax_f.set_xlim(-27.5, 27.5)
    ax_f.set_ylim(-23.9125, 62.1625)
    ax_f.set_title("Neuron Class Connectivity")
    ax_f.text(-13.75, -21.5, "Right Hemisphere", va = "center", ha = "center")
    ax_f.text( 13.75, -21.5, "Left Hemisphere",  va = "center", ha = "center")
    ax_f.axvline(0, color = "black", linestyle = "--", linewidth = 0.5, zorder = 1)
    ax_f.grid(False)

    #--------------------------------------------------------------------------
    # Cluster anatomy grid
    ca.make_cluster_subplots(info_dict, right_neurons, volumes,
                             color_dict = color_dict,
                             ax         = ax_g1,
                             view       = ("x", "-y"))

    ca.make_cluster_subplots(info_dict, right_neurons, volumes,
                             color_dict = color_dict,
                             ax         = ax_g2,
                             view       = ("x", "z"))

    for row in ax_g1:
        for col in row:
            col.set_xlim(-1000, 105000)
            col.set_ylim(93722, -2737)

    for row in ax_g2:
        for col in row:
            col.set_xlim(-1000, 105000)
            col.set_ylim(4000, 100459)
            col.set_title("")

    cluster_df["color"]  = cluster_df["bodyId"].map(color_dict)
    cluster_color_dict   = dict(zip(cluster_df.prat_cluster, cluster_df.color))
    ax_g1_flat = [a for row in ax_g1 for a in row]
    for key, val in acronym_dict.items():
        ax_g1_flat[key - 1].set_title(val, fontsize = 13,
                                      color = cluster_color_dict[key])

    #--------------------------------------------------------------------------
    if plot_3d:
        figure_3d = navis.plot3d(
            [right_neurons, volumes],
            color        = color_dict,
            legend_group = group_dict,
            width        = 1600,
            height       = 1600,
            inline       = False,
            backend      = "plotly")
        plotly.offline.plot(figure_3d,
                            filename  = f"{folder}/3D_View/{annot_name}.html",
                            auto_open = False)

    cluster_df.drop(labels = ["titles", "color"], axis = 1, inplace = True)

    if re_blast:
        save_linkage(info_dict, linkage_path)

    return fig, cluster_df


#------------------------------------------------------------------------------
def run(annot):
    """
    Per-lineage end-to-end: load, cluster, plot, append, save, free memory.

    """
    re_annot = annot.replace("/", "").replace(" ", "")

    neuron_collection = load_lineage_neurons(annot, category = "Diff")
    su_neurons        = load_su_neurons(annot, hemisphere = "Right")
    tract             = load_lineage_tract(annot, hemisphere = "Right")

    fig, cluster_df = generate_clusters(neuron_collection, annot,
                                        su_neurons = su_neurons,
                                        tract      = tract)

    ca.append_annotations_to_table(cluster_df, table_name, mode = "parquet")
    fig.savefig(f"{lineage_pdf_dir}/{re_annot}.pdf", bbox_inches = "tight", dpi = 400)

    plt.cla()
    plt.clf()
    plt.close(fig)
    del neuron_collection, fig, cluster_df
    gc.collect()


#------------------------------------------------------------------------------
def has_neurons(annot, category):
    """
    Skip the lineage when it has no skeletons of this category on disk.

    """
    return os.path.isdir(skeleton_folder(annot, category = category))


#------------------------------------------------------------------------------
# Iterate over all lineages. Lineages without Diff neurons are silently
# skipped; only genuine pipeline failures get logged.
for annot in annotations:
    if not has_neurons(annot, "Diff"):
        continue

    print(f"\n\tWorking on {annot}")

    try:
        run(annot)

    except Exception as e:
        errors.append(annot)
        with open(error_path, "a") as f:
            f.write(f"\n{annot}\n\t{e}")

    gc.collect()

#------------------------------------------------------------------------------
if not errors:
    print("\nNo Errors Encountered!")
    print("बढ़िया!\n")
else:
    print(f"\nError Compilation ({len(errors)})\n")

print("हो गया दोस्तों!")
