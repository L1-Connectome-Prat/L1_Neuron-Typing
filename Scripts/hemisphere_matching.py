"""
Best left-right neuron matches per lineage. Independent of any typing
scheme - this is a sanity check that mirrored NBLAST scores recover the
expected bilateral pairings.

Right vs mirrored-Left and Left vs mirrored-Right NBLAST scores are
averaged, then the Hungarian algorithm picks an optimal one-to-one match.

"""
import gc
import os
import shutil

import numpy as np
import pandas as pd

from scipy.optimize           import linear_sum_assignment
from scipy.cluster.hierarchy  import linkage, leaves_list

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks")

# Connectome
from navis.nbl.smat import Lookup2d

# Analysis package
import connectome_analysis_claude as ca
import flybrains 

# Local helpers
from Scripts.viz._paths  import (NEURON_INFO_CSV, nblast_matrix_path)
from Scripts.viz.loaders import load_lineage_neurons

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
smat_val    = 250
dot_tangent = 12
folder      = "Analysis_Outputs/Hemisphere_Matching"

params = {"out_folder"  : folder,
          "dot_tangent" : dot_tangent,
          "smat_val"    : smat_val,
          "use_alpha"   : True,
          "method"      : "average",
          "metric"      : "euclidean"}

#------------------------------------------------------------------------------
smat   = pd.read_csv(nblast_matrix_path(smat_val), index_col = 0)
smat2d = Lookup2d.from_dataframe(smat)

#------------------------------------------------------------------------------
# Lineage list comes from the L1_Skeletons annotation table - Diff lineages
# only. Sensory bundles (S-* prefix) and BolwigsAxon are excluded since they
# don't have meaningful bilateral pairs to recover.
df = pd.read_csv(NEURON_INFO_CSV).fillna("NA")
annotations = np.unique(df.loc[df["Lineage"] != "NA", "Lineage"]).astype(str)
annotations = np.char.rstrip(annotations)
annotations = [a for a in annotations
               if not a.startswith("S-") and a != "BolwigsAxon"]

#------------------------------------------------------------------------------
errors, error_path, table_name = ca.nblast_logging(params, format = "parquet")

# Per-lineage outputs
matrix_dir   = f"{folder}/Matrices"
lineage_dir  = f"{folder}/Lineages"
os.makedirs(matrix_dir,  exist_ok = True)
os.makedirs(lineage_dir, exist_ok = True)

#------------------------------------------------------------------------------
def average_lr_distances(result_r, result_l):
    """
    Reorder and average the right- and left-side NBLAST score matrices.

    """
    result_l = result_l.sort_index(ascending = True)
    result_l = result_l.reindex(sorted(result_l.columns), axis = 1)

    result_r = result_r.sort_index(ascending = True)
    result_r = result_r.reindex(sorted(result_r.columns), axis = 1)

    return (result_l + result_r) / 2


#------------------------------------------------------------------------------
def optimal_pairs(result):
    """
    Cluster-reorder the score matrix then use the Hungarian algorithm to
    select an optimal one-to-one matching of right -> left neurons.

    """
    row_linkage = linkage(result,   method = "average", metric = "euclidean")
    col_linkage = linkage(result.T, method = "average", metric = "euclidean")

    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)
    rx = result.iloc[row_order, col_order]

    cost_matrix      = -rx.values
    row_idx, col_idx = linear_sum_assignment(cost_matrix)

    matched_pairs = [(rx.index[r], rx.columns[c]) for r, c in zip(row_idx, col_idx)]
    match_df = pd.DataFrame(matched_pairs, columns = ["Right_ID", "Left_ID"])

    return match_df, rx, row_idx, col_idx


#------------------------------------------------------------------------------
def run_nblast(neuron_collection):
    """
    Score Right vs mirrored-Left and Left vs mirrored-Right.

    """
    nxr, nxl, _ = neuron_collection
    # Mirror copies keep their original .hemisphere tag, so split that way -
    # the previous len(nxr) // 2 split silently broke whenever the lineage had
    # unequal hemisphere counts.
    neurons_r  = nxr[nxr.hemisphere == "Right"]
    mirror_l2r = nxr[nxr.hemisphere == "Left"]
    neurons_l  = nxl[nxl.hemisphere == "Left"]
    mirror_r2l = nxl[nxl.hemisphere == "Right"]

    result_r = ca.run_source2target_nblast(neurons_r, mirror_l2r,
                                           dot_tangent = dot_tangent,
                                           smat        = smat2d,
                                           use_alpha   = True)

    result_l = ca.run_source2target_nblast(neurons_l, mirror_r2l,
                                           dot_tangent = dot_tangent,
                                           smat        = smat2d,
                                           use_alpha   = True)

    return result_r, result_l.T


#------------------------------------------------------------------------------
def plot_performance(rx, annot, row_idx, col_idx):
    """
    Heatmap of average NBLAST scores with stars on the Hungarian matches.

    """
    fig, ax = plt.subplots(figsize = (10, 10))

    sns.heatmap(rx, vmin = -0.08, vmax = 0.08, cmap = "vlag",
                square = True,
                xticklabels = True, yticklabels = True,
                cbar_kws    = {"shrink" : 0.5, "pad" : 0.015,
                               "anchor" : (1, 0.5),
                               "label"  : "NBLAST Score",
                               "aspect" : 25})

    for r, c in zip(row_idx, col_idx):
        ax.text(c + 0.5, r + 0.5, "❋",
                color = "black", ha = "center", va = "center", fontsize = 10)

    ax.set_xlabel("Left Neurons",  fontsize = 12)
    ax.set_ylabel("Right Neurons", fontsize = 12)
    plt.suptitle(annot, y = 0.94, fontsize = 20)

    return fig, ax


#------------------------------------------------------------------------------
def generate_matches(annot):
    """
    Load a lineage, score across hemispheres, and pick optimal pairs.

    """
    neuron_collection = load_lineage_neurons(annot, category = "Diff")
    nxr, nxl, _ = neuron_collection

    neurons_r = nxr[nxr.hemisphere == "Right"]
    neurons_l = nxl[nxl.hemisphere == "Left"]

    right_df = pd.DataFrame(zip(neurons_r.id, neurons_r.name),
                            columns = ["Right_ID", "Right_Name"])
    left_df  = pd.DataFrame(zip(neurons_l.id, neurons_l.name),
                            columns = ["Left_ID",  "Left_Name"])

    result_r, result_l = run_nblast(neuron_collection)
    result = average_lr_distances(result_r, result_l)

    match_df, rx, row_idx, col_idx = optimal_pairs(result)
    match_df = pd.merge(match_df, right_df)
    match_df = pd.merge(match_df, left_df)

    fig, _ = plot_performance(rx, annot, row_idx, col_idx)
    return match_df, fig, result


#------------------------------------------------------------------------------
def run(annot):
    """
    Per-lineage end-to-end matching, table append, matrix + figure save.

    """
    re_annot = annot.replace("/", "").replace(" ", "")

    match_df, fig, result = generate_matches(annot)
    ca.append_annotations_to_table(match_df, table_name, mode = "parquet")

    result.to_parquet(f"{matrix_dir}/{re_annot}.parquet", engine = "pyarrow")
    fig.savefig(f"{lineage_dir}/{re_annot}.pdf", bbox_inches = "tight")

    plt.clf()
    plt.cla()
    plt.close("all")
    gc.collect()


#------------------------------------------------------------------------------
for annot in annotations:
    print(f"\n\tWorking on {annot}")

    try:
        run(annot)

    except Exception as e:
        errors.append(annot)
        with open(error_path, "a") as f:
            f.write(f"\n{annot}\n\t{e}")

#------------------------------------------------------------------------------
# nblast_logging spins up Linkage / 3D_View folders for the clustering scripts -
# neither is produced here, so drop them.
for stale in ("Linkage", "3D_View"):
    shutil.rmtree(f"{folder}/{stale}", ignore_errors = True)

#------------------------------------------------------------------------------
if not errors:
    print("\nNo Errors Encountered!")
    print("बढ़िया!\n")
else:
    print(f"\nError Compilation ({len(errors)})\n")

print("हो गया दोस्तों!")
