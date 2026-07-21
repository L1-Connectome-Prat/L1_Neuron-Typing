"""
Skeleton-loading helpers shared by the typing scripts. Reads from the
L1_Skeletons repository and applies the L1 mirror/symmetry transforms.

"""
import os

import navis

from ._paths import skeleton_folder, tract_folder

# Connectome package
import connectome_analysis_claude as ca

#------------------------------------------------------------------------------
__all__ = ["rename_lr_neurons",
           "load_lineage_neurons",
           "load_sensory_neurons",
           "load_su_neurons",
           "load_lineage_tract"]

# SWC node label -> connector type mapping baked into the L1_Skeletons writes
SWC_CONNECTOR_LABELS = {0 : 7, 1 : 8}


#------------------------------------------------------------------------------
def _as_neuronlist(neurons):
    """
    Force a navis result into a NeuronList - read_swc / xform_brain on a
    single skeleton can return a bare TreeNeuron, which then breaks `len()`.

    """
    if neurons is None:
        return navis.NeuronList([])
    if isinstance(neurons, navis.TreeNeuron):
        return navis.NeuronList([neurons])
    return neurons


#------------------------------------------------------------------------------
def rename_lr_neurons(neurons):
    """
    Stitch lineage / cluster / instance / acronym / hemisphere / id into a
    single neuron name used downstream for dendrogram labels.

    """
    for neuron in neurons:
        name = f"{neuron.lineage}_{neuron.cluster}{neuron.instance}_{neuron.acronym}_{neuron.hemisphere}_{neuron.id}"
        # Drop missing acronyms
        name = name.replace("NA", "")
        name = name.replace("__", "_NA_")
        neuron.name = name

    return neurons


#------------------------------------------------------------------------------
def load_lineage_neurons(annot, category = "Diff",
                         template = "PK_L1CNS",
                         symmetrize_to = "PK_L1CNSsym"):
    """
    Build the [nxr, nxl, all] neuron collection for a lineage annotation.

    Wraps connectome_analysis_claude.assemble_neuron_collection (forwarding
    connector_labels so the L1 SWC label-7/8 synapse markers populate
    .connectors), then names + symmetrizes for plotting.

    """
    folder = skeleton_folder(annot, category = category)

    if not os.path.isdir(folder):
        raise FileNotFoundError(f"No {category} skeletons for '{annot}' at {folder}")

    nxr, nxl, all_neurons = ca.assemble_neuron_collection(
                                folder,
                                template         = template,
                                connector_labels = SWC_CONNECTOR_LABELS)

    # Tag the lineage onto every neuron - mirrored copies inherit attributes
    for nl in (nxr, nxl, all_neurons):
        for n in nl:
            n.lineage = annot

    rename_lr_neurons(nxr)
    rename_lr_neurons(nxl)
    rename_lr_neurons(all_neurons)

    # Symmetrize the union (only used for plotting / synapse views)
    all_neurons = navis.xform_brain(all_neurons,
                                    source = template,
                                    target = symmetrize_to)

    return [nxr, nxl, all_neurons]


#------------------------------------------------------------------------------
def load_su_neurons(annot, hemisphere = "Right",
                    template = "PK_L1CNS",
                    symmetrize_to = "PK_L1CNSsym"):
    """
    Load the SU neurons for a lineage / hemisphere if any exist.

    Returns an empty NeuronList silently when the lineage has no SU neurons -
    not every lineage has them.

    """
    folder = f"{skeleton_folder(annot, category = 'SU')}/{hemisphere}"

    if not os.path.isdir(folder):
        return navis.NeuronList([])

    neurons = _as_neuronlist(navis.read_swc(
                folder, fmt = "{id}.swc",
                connector_labels = SWC_CONNECTOR_LABELS))

    if not len(neurons):
        return navis.NeuronList([])

    return _as_neuronlist(navis.xform_brain(
                neurons, source = template, target = symmetrize_to))

#------------------------------------------------------------------------------
def load_sensory_neurons(annot,
                         bundles,
                         template = "PK_L1CNS",
                         symmetrize_to = "PK_L1CNSsym"):
    """
    Build the [nxr, nxl, all] neuron collection for a sensory annotation.

    Mirrors the bundle-cut flow from the original notebook: read raw
    right/left in PK_L1CNS, cut the sensory bundle volumes BEFORE mirroring,
    drop stub fragments, then mirror the cut sets to build nxr / nxl. The
    symmetrized union is built from the RAW (uncut) skeletons so the
    figure panels show full neurons.

    """
    folder = skeleton_folder(annot, category = "Diff")

    if not os.path.isdir(folder):
        raise FileNotFoundError(f"No Sensory skeletons for '{annot}' at {folder}")

    right = _as_neuronlist(navis.read_swc(
                f"{folder}/Right/", fmt = "{id}.swc",
                connector_labels = SWC_CONNECTOR_LABELS))
    left  = _as_neuronlist(navis.read_swc(
                f"{folder}/Left/",  fmt = "{id}.swc",
                connector_labels = SWC_CONNECTOR_LABELS))

    if not len(right) or not len(left):
        raise FileNotFoundError(
            f"Empty or missing skeleton folder under '{folder}'. "
            f"Right: {len(right)} neurons, Left: {len(left)} neurons.")

    # Tag the lineage onto every neuron
    for n in right:
        n.lineage = annot
    for n in left:
        n.lineage = annot

    rename_lr_neurons(right)
    rename_lr_neurons(left)

    # NOTE : Since they are pre-clustered, we are removing the fragments that are unconnected
    # If you are trying to re-run it from scratch, you might want to remove this
    right = right[right.cluster != "NA"]
    left  = left[left.cluster   != "NA"]

    # Snapshot the uncut originals - the symmetrized union uses these
    raw_r = right.copy()
    raw_l = left.copy()

    # Cut sensory bundles BEFORE mirroring, in PK_L1CNS coordinates
    for bundle in bundles:
        navis.in_volume(right, mode = "OUT", volume = bundle,
                        prevent_fragments = True, inplace = True)
        navis.in_volume(left,  mode = "OUT", volume = bundle,
                        prevent_fragments = True, inplace = True)

    # Drop stub fragments and apply the same mask to the raw copies
    keep_r = right.n_nodes > 25
    keep_l = left.n_nodes  > 25
    right  = right[keep_r]
    left   = left[keep_l]
    raw_r  = raw_r[keep_r]
    raw_l  = raw_l[keep_l]

    # Mirror the cut neurons within PK_L1CNS to build the symmetric NBLAST sets
    mirror_l2r = navis.mirror_brain(left,  template = template, warp = True)
    mirror_r2l = navis.mirror_brain(right, template = template, warp = True)

    nxr = right + mirror_l2r
    nxl = left  + mirror_r2l

    # Symmetrized union built from raw (uncut) skeletons - plotting only
    all_neurons = navis.xform_brain(raw_r + raw_l,
                                    source = template,
                                    target = symmetrize_to)

    return [nxr, nxl, all_neurons]

#------------------------------------------------------------------------------
def load_lineage_tract(annot, hemisphere = "Right",
                       template = "PK_L1CNS",
                       symmetrize_to = "PK_L1CNSsym"):
    """
    Load the tract skeleton for a lineage / hemisphere if it exists.

    """
    folder = tract_folder(annot, hemisphere = hemisphere)

    if not os.path.isdir(folder):
        return navis.NeuronList([])

    tract = _as_neuronlist(navis.read_swc(folder, fmt = "{id}.swc"))

    if not len(tract):
        return navis.NeuronList([])

    return _as_neuronlist(navis.xform_brain(
                tract, source = template, target = symmetrize_to))


