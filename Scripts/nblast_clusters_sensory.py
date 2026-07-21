"""
Per-annotation morphology-based typing of sensory input bundles via NBLAST +
Ward clustering. Same flow as the Diff script, with the added step of cutting
out the sensory bundle tract from each neuron before clustering so that
the bundle geometry doesn't dominate the NBLAST distances.

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
from Scripts.viz._paths          import (ANNOTATIONS_NPY, nblast_matrix_path,
                                         MB_RIGHT_STL, MB_LEFT_STL,
                                         skeleton_folder)
from Scripts.viz.loaders         import load_sensory_neurons
from Scripts.viz.anatomy         import plot_lineage_anatomy, synapse_plot
from Scripts.viz.volker_bars     import plot_volkers_bars
from Scripts.viz.figure_layouts  import build_sensory_layout
from Scripts.viz.clustering_io   import save_linkage, load_linkage

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
smat_val    = 250
dot_tangent = 12
plot_3d     = False
re_blast    = False
folder      = "Analysis_Outputs/Neuron_Clusters_Sensory"

params = {"out_folder"  : folder,
          "dot_tangent" : dot_tangent,
          "smat_val"    : smat_val,
          "use_alpha"   : False,
          "score"       : "min",
          "method"      : "ward",
          "threshold"   : 0.5,
          "plot_3d"     : plot_3d,
          "re_blast"    : re_blast,
          "category"    : "Sensory"}

#------------------------------------------------------------------------------
# NBLAST scoring matrix
smat   = pd.read_csv(nblast_matrix_path(smat_val), index_col = 0)
smat2d = Lookup2d.from_dataframe(smat)

# Sensory annotations - everything tagged with the "S-" prefix
annotations = np.load(ANNOTATIONS_NPY)
# Remove also the S-AN
annotations = [a for a in annotations if str(a).startswith("S-") and str(a) != "S-AN"]
annotations = sorted(annotations)

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
# Sensory bundle meshes - keep only the portions outside the CNS
SENSORY_BUNDLE_DIR = f"{os.environ['PROJECTS_HOME']}/L1_Lineages/Analysis_Outputs/Volumes/Sensory_Bundles"
bundle_files = ["S-GLN/S-GLN_Right_Mesh_smoothed.stl",
                "S-GLN/S-GLN_Left_Mesh_smoothed.stl",
                "S-AN/S-AN_Right_Mesh_smoothed.stl",
                "S-AN/S-AN_Left_Mesh_smoothed.stl"]

cut_bundles_cns = []
for fname in bundle_files:
    vx = ca.stl_to_navis(f"{SENSORY_BUNDLE_DIR}/{fname}")
    vx.fix_normals()
    vx.name = fname.split("/")[-1].split(".")[0]
    cut_bundles_cns.append(vx)

#------------------------------------------------------------------------------
# Initialize the run folder, error log, and per-run table path
errors, error_path, table_name = ca.nblast_logging(params, format = "parquet")

# Per-annotation figure folder
lineage_pdf_dir = f"{folder}/Lineages"
os.makedirs(lineage_pdf_dir, exist_ok = True)

#------------------------------------------------------------------------------
def generate_clusters(neuron_collection, annot):
    """
    Run (or reload) the NBLAST clustering for a sensory annotation and draw
    the multi-panel figure.

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
        info_dict["Z"]         = Z
        info_dict["nbl_dist"]  = nbl_dist
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

    cluster_df["titles"] = cluster_df.apply(lambda row:
                                f"{row['name'].split('_')[1][:2]}_{row['name'].split('_')[2]}",
                                axis = 1)

    titles_by_key = defaultdict(set)
    for key, title in zip(cluster_df.prat_cluster, cluster_df.titles):
        titles_by_key[key].add(title)

    acronym_dict = {key : " &\n".join(sorted(s)) for key, s in titles_by_key.items()}

    #--------------------------------------------------------------------------
    # Figure scaffold
    fig, ax_a, ax_b, ax_c, ax_d, ax_g1, ax_g2 = build_sensory_layout(annot)
    right_neurons = neuron_collection[2][neuron_collection[2].hemisphere == "Right"]

    # Sensory clone (column A)
    plot_lineage_anatomy(right_neurons, volumes, ax_a)
    ax_a[0].set_title("Sensory Clone", fontsize = 15)

    # Synapse panels (columns B = PSDs, C = T-Bars)
    synapse_plot(right_neurons, neuropil, ax = np.vstack((ax_b, ax_c)))

    # Dendrogram with Volker bars (column D)
    ax_d1, ax_d2 = ax_d
    color_dict = ca.build_dendrogram(
                        Z, nbl_dist,
                        info_dict["optimal_t"],
                        info_dict["cluster_labels"],
                        names   = True,
                        neurons = neuron_collection[2],
                        ax      = ax_d1)

    ax_d2.set_ylim(ax_d1.get_ylim())
    ax_d1.set_title("Morphological Clustering")
    ax_d1.set_ylabel("")
    ax_d1.set_xlim(ax_d1.get_xlim()[0], 0.75)

    volker_dict, volker_skel_dict = plot_volkers_bars(ax_d1, ax_d2, prop_color = True)

    ax_d1.invert_yaxis()
    ax_d2.invert_yaxis()

    pos_d1    = ax_d1.get_tightbbox()
    pos_d2    = ax_d2.get_position()
    shifted_x = pos_d1.x1 + 8
    ratio     = fig.get_figwidth() * fig.dpi
    ax_d2.set_position([shifted_x / ratio, pos_d2.y0, pos_d2.width, pos_d2.height])
    ax_d1.set_xlabel("")

    #--------------------------------------------------------------------------
    # Cluster anatomy grid (panel E)
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
    Per-annotation end-to-end: load, cut bundles, cluster, plot, append, save.

    """
    re_annot = annot.replace("/", "").replace(" ", "")

    neuron_collection = load_sensory_neurons(annot, bundles = cut_bundles_cns)

    fig, cluster_df = generate_clusters(neuron_collection, annot)

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
    Skip the annotation when it has no skeletons of this category on disk.

    """
    return os.path.isdir(skeleton_folder(annot, category = category))


#------------------------------------------------------------------------------
# Iterate over all sensory annotations. Annotations without Diff neurons are
# silently skipped; only genuine pipeline failures get logged.
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
