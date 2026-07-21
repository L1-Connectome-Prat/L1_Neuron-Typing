"""
Upload typing/nomenclature annotations back into CATMAID.

The script reads a CSV with a neuron-id column and an info column, then
either uploads neuron-by-neuron (sequential) or once-per-info-value
(mass). Both attach the same meta-annotation so the upload is easy to
revoke later.

Edit the path / column / meta_annotation block at the bottom for each run.

We expect errors in:
1) DPLc2 : Does not have any differentiated neurons in the L1
2) DPMl3/4a : Has only one neuron per hemisphere (default match)

"""
import json
import os
from importlib.resources import files

import numpy as np
import pandas as pd

import pymaid
import connectome_analysis_claude as ca

#------------------------------------------------------------------------------
label_file = files("connectome_analysis_claude.connect.catmaid_helpers").joinpath("Prat_Labeling.json")
with open(str(label_file), "r") as f:
    prat_labels = json.load(f)

larval_analysis = ca.CATMAIDConnect(
                        os.environ["L1_USERNAME"],
                        os.environ["L1_PASSWORD"],
                        os.environ["L1_TOKEN"],
                        os.environ["L1_SERVER"],
                        custom_label = prat_labels)

#------------------------------------------------------------------------------
def _format_value(value, prefix = "g", width = 3):
    """
    Standard prefix + zero-padded value used for cluster / group ids.

    """
    return f"{prefix}{int(value):0{width}d}"


#------------------------------------------------------------------------------
def upload_info_seq(info_df, neu_id_col, info_col, meta_annotation,
                    prefix = "g", width = 3):
    """
    Upload the annotation neuron-by-neuron, then meta-annotate each value.

    """
    value_list = []
    for _, row in info_df.iterrows():
        value = _format_value(row[info_col], prefix = prefix, width = width)
        value_list.append(value)
        pymaid.add_annotations(row[neu_id_col], value)

    for value in value_list:
        pymaid.add_meta_annotations(value, meta_annotation)


#------------------------------------------------------------------------------
def upload_info_mass(info_df, neu_id_col, info_col, meta_annotation,
                     prefix = "g", width = 3):
    """
    Upload one bulk annotation per unique info value - faster than per-neuron
    when many neurons share the same value.

    """
    infos = np.unique(info_df[info_col].dropna()).astype(int)
    value_list = []

    for info in infos:
        ids   = list(info_df.loc[info_df[info_col] == info, neu_id_col])
        value = _format_value(info, prefix = prefix, width = width)
        value_list.append(value)
        pymaid.add_annotations(ids, value)

    for value in value_list:
        pymaid.add_meta_annotations(value, meta_annotation)


#------------------------------------------------------------------------------
def organize_info(info_path, neu_id_col, info_col, meta_annotation,
                  mode   = "seq",
                  prefix = "g",
                  width  = 3):
    """
    Read the source CSV and dispatch to the chosen upload mode.

    """
    info_df = pd.read_csv(info_path)

    if mode == "seq":
        upload_info_seq(info_df, neu_id_col, info_col, meta_annotation,
                        prefix = prefix, width = width)
    elif mode == "mass":
        upload_info_mass(info_df, neu_id_col, info_col, meta_annotation,
                         prefix = prefix, width = width)
    else:
        raise ValueError("mode must be one of: seq, mass.")


#------------------------------------------------------------------------------
# Edit per run
path            = "Analysis_Outputs/L1_NeuronInformation_PK.csv"
neu_id_col      = "bodyId"
info_col        = "ConnectivityGroup"
meta_annotation = "PK_ConnectivityGroup"
prefix          = "g"
width           = 3
mode            = "seq"

organize_info(path, neu_id_col, info_col, meta_annotation,
              mode   = mode,
              prefix = prefix,
              width  = width)

#------------------------------------------------------------------------------
print("हो गया दोस्तों!")