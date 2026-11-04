"""
Per-lineage morphology-based typing of Small Undifferentiated (SU) L1 neurons
via NBLAST + Ward clustering. Same flow as the Diff script but a slimmer
panel layout - no bubble map, Volker bars, or tract overlays since SU
neurons don't have meaningful synaptic connectivity.

Lineages without SU neurons are silently skipped.

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
import connectome_analysis as ca

# Local helpers
from Scripts.viz._paths        import (LINEAGES_NPY, nblast_matrix_path,
                                       MB_RIGHT_STL, MB_LEFT_STL,
                                       skeleton_folder)
from Scripts.viz.loaders        import load_lineage_neurons
from Scripts.viz.anatomy        import plot_lineage_anatomy
from Scripts.viz.figure_layouts import build_su_layout
from Scripts.viz.clustering_io  import save_linkage, load_linkage

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
smat_val    = 250
dot_tangent = 12
plot_3d     = False
re_blast    = False
folder      = "Analysis_Outputs/Neuron_Clusters_SU"

params = {"out_folder"  : folder,
          "dot_tangent" : dot_tangent,
          "smat_val"    : smat_val,
          "use_alpha"   : True,
          "score"       : "min",
          "method"      : "ward",
          "threshold"   : 0.33,
          "plot_3d"     : plot_3d,
          "re_blast"    : re_blast,
          "category"    : "SU"}

#------------------------------------------------------------------------------
smat   = pd.read_csv(nblast_matrix_path(smat_val), index_col = 0)
smat2d = Lookup2d.from_dataframe(smat)

annotations = np.load(LINEAGES_NPY)

#------------------------------------------------------------------------------
# Volumes
cns      = flybrains.PK_L1CNSsym.mesh
neuropil = flybrains.PK_L1CBNeuropilsym.mesh
cns.color      = (0.85, 0.85, 0.85, 0.2)
neuropil.color = (0, 0, 0, 1)

mbr = ca.stl_to_navis(MB_RIGHT_STL)
mbl = ca.stl_to_navis(MB_LEFT_STL)
mbr.color = (0.85, 0.85, 0.85, 0.6)
mbl.color = (0.85, 0.85, 0.85, 0.6)
mbr = navis.xform_brain(mbr, source = "PK_L1CNS", target = "PK_L1CNSsym")
mbl = navis.xform_brain(mbl, source = "PK_L1CNS", target = "PK_L1CNSsym")

volumes = [cns, mbr, mbl]

#------------------------------------------------------------------------------
errors, error_path, table_name = ca.nblast_logging(params, format = "parquet")

# Per-lineage figure folder
lineage_pdf_dir = f"{folder}/Lineages"
os.makedirs(lineage_pdf_dir, exist_ok = True)

#------------------------------------------------------------------------------
def has_neurons(annot, category):
    """
    Skip the lineage when it has no skeletons of this category on disk.

    """
    return os.path.isdir(skeleton_folder(annot, category = category))


#------------------------------------------------------------------------------
def generate_clusters(neuron_collection, annot):
    """
    Run (or reload) the NBLAST clustering for the SU set of a lineage.

    """
    annot_name   = annot.replace("/", "")
    linkage_path = f"{folder}/Linkage/{annot_name}.h5"

    if re_blast:
        nbl_dist, nbl_vec, Z = ca.prats_nblast_pipeline(
                                    neuron_collection,
                                    dot_tangent = dot_tangent,
                                    resample    = False,
                                    smat        = smat2d,
                                    use_alpha   = True,
                                    score       = "min",
                                    n_cores     = 32,
                                    method      = "ward")

        info_dict = ca.run_multihemispheric_clustering(
                        nbl_dist, Z, neuron_collection[2],
                        t_adj     = 0.001,
                        threshold = 0.33)
        info_dict["Z"]         = Z
        info_dict["nbl_dist"]  = nbl_dist
        info_dict["threshold"] = 0.33
    else:
        info_dict = load_linkage(linkage_path)
        Z        = info_dict["Z"]
        nbl_dist = info_dict["nbl_dist"]

    #--------------------------------------------------------------------------
    cluster_df = ca.assemble_annotations(info_dict, neuron_collection[2], annot)

    if re_blast:
        group_dict = dict(zip(cluster_df.bodyId.values.tolist(),
                              cluster_df.prat_cluster.values.tolist()))
    else:
        group_dict = dict(zip(neuron_collection[2].id, neuron_collection[2].cluster))

    cluster_df["titles"] = cluster_df.apply(lambda row:
                                f"{row['name'].split('_')[1][:3]}_{row['name'].split('_')[2]}",
                                axis = 1)

    titles_by_key = defaultdict(set)
    for key, title in zip(cluster_df.prat_cluster, cluster_df.titles):
        titles_by_key[key].add(title)
    acronym_dict = {key : " &\n".join(sorted(s)) for key, s in titles_by_key.items()}

    #--------------------------------------------------------------------------
    fig, ax_a, ax_b, ax_g1, ax_g2 = build_su_layout(annot)
    right_neurons = neuron_collection[2][neuron_collection[2].hemisphere == "Right"]

    plot_lineage_anatomy(right_neurons, volumes, ax_a)
    ax_a[0].set_title("SU Clone", fontsize = 15)

    color_dict = ca.build_dendrogram(
                        Z, nbl_dist,
                        info_dict["optimal_t"],
                        info_dict["cluster_labels"],
                        names   = True,
                        neurons = neuron_collection[2],
                        ax      = ax_b)

    ax_b.set_title("Morphological Clustering", fontsize = 15)
    ax_b.set_ylabel("")
    ax_b.set_xlabel("")
    ax_b.set_xlim(ax_b.get_xlim()[0], 0.75)
    ax_b.invert_yaxis()

    #--------------------------------------------------------------------------
    ca.make_cluster_subplots(info_dict, right_neurons, volumes,
                             color_dict = color_dict,
                             ax         = ax_g1,
                             view       = ("x", "-y"))

    ca.make_cluster_subplots(info_dict, right_neurons, volumes,
                             color_dict = color_dict,
                             ax         = ax_g2,
                             view       = ("x", "-z"))

    for row in ax_g1:
        for col in row:
            col.set_xlim(0, 105000)
            col.set_aspect("equal")

    for row in ax_g2:
        for col in row:
            col.set_ylim(110000, 0)
            col.set_aspect("equal")
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
    Per-lineage end-to-end for the SU set.

    """
    re_annot = annot.replace("/", "").replace(" ", "")

    neuron_collection = load_lineage_neurons(annot, category = "SU")

    fig, cluster_df = generate_clusters(neuron_collection, annot)

    ca.append_annotations_to_table(cluster_df, table_name, mode = "parquet")
    fig.savefig(f"{lineage_pdf_dir}/{re_annot}.pdf", bbox_inches = "tight", dpi = 400)

    plt.cla()
    plt.clf()
    plt.close(fig)
    del neuron_collection, fig, cluster_df
    gc.collect()


#------------------------------------------------------------------------------
print()

for annot in annotations:
    if not has_neurons(annot, "SU"):
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
