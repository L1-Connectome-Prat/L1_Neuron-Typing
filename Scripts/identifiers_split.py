"""
L1 neuron identifiers built from axon-dendrite splits. For each neuron we
generate two compartment-based names (one for the dendrite, one for the
axon) and stitch them as `dendrite-|-axon`.

Each piece's name lists compartments that carry more than `frac_limit` of
that piece's connectors, sorted by fraction. Suffixes:

    .s   - synaptic input only      (PSDs > frac_limit)
    .b   - bouton output only       (TBars > frac_limit)
    .s.b - both input and output

Contralateral mentions are starred (``*``) so a name read on either side
makes the host hemisphere obvious.

Splits prefer Winding et al. (2023) tags (``mw axon split``) where present
and fall back to navis.split_axon_dendrite via synapse flow centrality.

Original annotations placed into CATMAID as PK_Identifier
These were generated using the archived version as of 2025-09-11. Any changes
to prod CATMAID will result in slight changes here. To get my (Prat's) assignments
just use the L1_NeuronInformation_PK.csv file in the L1_Skeletons repository.

This script, if run separately, requires that the neurons are at least annotated
by their hemisphere as are the compartments. This enables assignment of compartments
in the Identifiers as ipsi / contra(*).

NOTE : We use this same script to get the fragments for level 1 and level 2 compartments.
The Acronyms are defined using the details resolved at level 2.

Just a note to myself:
presynapses = TBars
postsynapses = PSDs

"""
import argparse
import os

import numpy as np
import pandas as pd

import navis
import pymaid
import connectome_analysis_claude as ca

from Scripts.viz._paths import NEURON_INFO_CSV, COMPARTMENTS_HOME

#------------------------------------------------------------------------------
# Run parameters
num_pieces = 100
frac_limit = 0.05
folder     = "Analysis_Outputs/Nomenclature"

#------------------------------------------------------------------------------
def load_compartments(args):
    """
    Load compartments from file.
    """
    compartment_file  = f"{COMPARTMENTS_HOME}/{args.level.capitalize()}.npy"
    compartment_names = np.load(compartment_file, allow_pickle = True)

    # Load each compartment
    compartments = []
    for comp in compartment_names:
        # Get the stl
        stl = ca.stl_to_navis(f"{COMPARTMENTS_HOME}/{comp}")
        # Fix the normals
        stl.fix_normals()
        # Get the name
        name = comp.split("/")[-1].replace(".stl", "")
        # Put name as property
        stl.name = name
        compartments.append(stl)

    return compartments

#------------------------------------------------------------------------------
def assemble_neurons():
    """
    Pull Diff skeletons fresh from CATMAID. Skipping the local SWC step
    avoids a sort-order bug where multiple connectors on the same node
    end up reordered.

    """
    neurons_info = pymaid.find_neurons(annotations = ["Kandimalla-etal.", "Diff"],
                                       intersect   = True)
    skeletons = pymaid.get_neurons(list(neurons_info.id))

    dendrites, axons = split_into_pieces(skeletons)

    # Total synapse counts recorded before splitting — used as denominators in piece_name
    syn_counter = pd.DataFrame(zip(skeletons.id,
                                   skeletons.n_presynapses,
                                   skeletons.n_postsynapses))
    syn_counter.columns   = ["bodyId", "TBars", "PSDs"]
    syn_counter["bodyId"] = syn_counter["bodyId"].astype(int)

    return skeletons, dendrites, axons, syn_counter


#------------------------------------------------------------------------------
def split_into_pieces(skeletons, metric = "synapse_flow_centrality"):
    """
    Split each skeleton into axon and dendrite. Prefer Winding cuts where
    tagged; fall back to navis.split_axon_dendrite.

    """
    dendrites = []
    axons     = []
    fallback_count = 0

    for skel in skeletons:
        try:
            cut_posn = skel.tags["mw axon split"][0]
            # navis.cut_skeleton returns (distal, proximal) — axon is distal
            axon, dendrite = navis.cut_skeleton(skel, cut_posn)
        except (KeyError, IndexError):
            fallback_count += 1
            try:
                pieces   = navis.split_axon_dendrite(skel, metric = metric)
                axon     = pieces[pieces.compartment == "axon"]
                dendrite = pieces[pieces.compartment == "dendrite"]
            except Exception:
                # Neuron is silently dropped from both lists if both strategies fail
                continue

        axons.append(axon)
        dendrites.append(dendrite)

    dendrites = navis.NeuronList(dendrites)
    axons     = navis.NeuronList(axons)

    print(f"Winding cuts not found for {fallback_count} neurons.")
    print(f"{1 - (fallback_count / len(skeletons)):.3f} consistent.")

    return dendrites, axons

#------------------------------------------------------------------------------
def syn_by_compartment(vol_dict, piece = "dendrite"):
    """
    Roll up per-piece synapses into a long-form table keyed by compartment.

    """
    # NeuronList attributes are vectorized — build one DataFrame per volume
    # instead of looping over individual neurons
    dfs = [pd.DataFrame({"Compartment" : vol,
                         "bodyId"      : neurons.id,
                         "Piece"       : piece,
                         "TBars"       : neurons.n_presynapses,
                         "PSDs"        : neurons.n_postsynapses})
           for vol, neurons in vol_dict.items()]

    return pd.concat(dfs, ignore_index = True)

#------------------------------------------------------------------------------
def piece_name(piece_df, org_info,
               num_pieces = num_pieces,
               frac_limit = frac_limit):
    """
    Build a single piece's name from its top compartments by fraction.

    """
    if piece_df is None:
        return ""

    # Get the info
    num_tbars  = org_info["TBars"].values[0]
    num_psds   = org_info["PSDs"].values[0]
    hemisphere = org_info["Hemisphere"].values[0]

    # Normalize to total counts
    piece_df["TBars"] = piece_df["TBars"] / num_tbars
    piece_df["PSDs"]  = piece_df["PSDs"]  / num_psds

    psd_df  = piece_df[piece_df["PSDs"]  > frac_limit].sort_values(by = "PSDs",  ascending = False)
    tbar_df = piece_df[piece_df["TBars"] > frac_limit].sort_values(by = "TBars", ascending = False)

    # Just get the compartments - I don't know why I did it like this :(
    psd_comps  = [c for c in psd_df["Compartment"]][:num_pieces]
    tbar_comps = [c for c in tbar_df["Compartment"]][:num_pieces]

    # Get the common and unique compartments
    common  = [c for c in psd_comps  if c in tbar_comps]
    exc_in  = [c for c in psd_comps  if c not in common]
    exc_out = [c for c in tbar_comps if c not in common]

    # Order: input-only (.s) → mixed (.s.b) → output-only (.b), each sorted by fraction
    name_list = [f"{c}.s"   for c in exc_in
                ] + [f"{c}.s.b" for c in common
                ] + [f"{c}.b"   for c in exc_out]

    name = "-".join(name_list)

    # Contralateral compartments get a star; ipsilateral suffix is dropped entirely
    if hemisphere == "Right":
        ipsi = "_Right"
        contra = "_Left"
    else:
        ipsi = "_Left"
        contra = "_Right"

    # Replace as needed
    name = name.replace(ipsi, "").replace(contra, "*")

    return name


#------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", choices = ["level1", "level2"], default = "level1",
                        help = "Compartment level to use for naming (default: level1).")
    args = parser.parse_args()

    # Compartment file name and tag move together; CATMAIDConnect loads via pickle
    tag          = args.level
    compartments = load_compartments(args)

    # Neuron assignment/information dictionary
    info_df = pd.read_csv(NEURON_INFO_CSV)
    info_df = info_df[["bodyId", "Lineage", "Hemisphere", "Acronym"]]

    # Instantiate CATMAIDConnect to get the skeletons
    larval_analysis = ca.CATMAIDConnect(
                            os.environ["L1_USERNAME"],
                            os.environ["L1_PASSWORD"],
                            os.environ["L1_TOKEN"],
                            os.environ["L1_SERVER"])

    #--------------------------------------------------------------------------
    # Get the skeletons from CATMAID and split them into dendrites and axons
    skeletons, dendrites, axons, syn_counter = assemble_neurons()
    # Merge information into fragment
    syn_counter = pd.merge(syn_counter, info_df, on = "bodyId", how = "left")

    # Dendrites
    dendrite_dict = navis.in_volume(dendrites,
                                    compartments,
                                    mode    = "IN",
                                    inplace = False)
    dendrites_df = syn_by_compartment(dendrite_dict, "dendrite")

    # Axons
    axon_dict = navis.in_volume(axons,
                                compartments,
                                mode    = "IN",
                                inplace = False)
    axons_df = syn_by_compartment(axon_dict, "axon")

    # Ints as bodyIds
    dendrites_df["bodyId"] = dendrites_df["bodyId"].astype(int)
    axons_df["bodyId"]     = axons_df["bodyId"].astype(int)

    #--------------------------------------------------------------------------
    # Pre-group once so the loop does O(1) lookups instead of full-frame scans
    den_groups = dict(tuple(dendrites_df.groupby("bodyId")))
    ax_groups  = dict(tuple(axons_df.groupby("bodyId")))
    neu_groups = dict(tuple(syn_counter.groupby("bodyId")))

    id_dict = {}
    for neuron in np.unique(syn_counter["bodyId"]):
        # Get the name fragments
        dendrite_name = piece_name(den_groups.get(neuron), 
                                   neu_groups.get(neuron), 
                                   num_pieces)
        axon_name     = piece_name(ax_groups.get(neuron), 
                                   neu_groups.get(neuron), 
                                   num_pieces)
        
        id_dict[neuron] = f"{dendrite_name}-|-{axon_name}"

    #--------------------------------------------------------------------------
    # Merge
    syn_counter["Prat_Identifier"] = syn_counter["bodyId"].map(id_dict)
    # Save the counter and information
    syn_counter.to_csv(f"{folder}/Split_Nomenclature_Identifiers_{tag}_{num_pieces}.csv",
                       index = False)
    syn_counter.to_parquet(f"{folder}/Split_Nomenclature_Identifiers_{tag}_{num_pieces}.parquet",
                           index  = False,
                           engine = "pyarrow")

    #--------------------------------------------------------------------------
    # Get the axon-dendrite individual compartment innervations
    innervation_df = pd.concat((dendrites_df, axons_df))
    # Drop compartments a fragment passed through but made no synapses in
    innervation_df = innervation_df.loc[(innervation_df["TBars"] != 0) | (innervation_df["PSDs"] != 0)]
    innervation_df = innervation_df.sort_values(by = ["bodyId", "Piece", "TBars", "PSDs"],
                                                ascending = [True, True, False, False]).reset_index(drop = True)
    innervation_df = innervation_df[["bodyId", "Piece", "Compartment", "TBars", "PSDs"]]

    # Save
    innervation_df.to_csv(f"{folder}/Split_Nomenclature_Fragments_{tag}.csv", 
                          index = False)
    innervation_df.to_parquet(f"{folder}/Split_Nomenclature_Fragments_{tag}.parquet",
                              index  = False,
                              engine = "pyarrow")

    #--------------------------------------------------------------------------
    print("हो गया दोस्तों!")