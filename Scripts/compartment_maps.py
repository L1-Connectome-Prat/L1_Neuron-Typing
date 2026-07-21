"""
Compartment-by-compartment connectivity maps for each lineage-hemisphere-acronym
unit. Reads directly from the Synapse-by-Acronym parquets and produces:

  - LineageAcronymCompartments_level1.parquet       : source × target compartment
                                                    cross-product with weights
  - LineageAcronymCompartments_Pieces_level1.parquet : same, split by axon/dendrite piece

The weight formula is identical to the original innervation_acronyms_full/split
pipeline:

    Ratio_PSDs  = Source_PSDs  / acronym_PSDs
    Ratio_TBars = Target_TBars / acronym_TBars
    Wraw        = Ratio_PSDs * Ratio_TBars
    Scale       = acronym_connectors / lineage_connectors
    Wscaled     = Wraw * Scale

"""
import pandas as pd

from Scripts.viz._paths import (PAIRED_FULL_PARQUET, PAIRED_PIECES_PARQUET)

#------------------------------------------------------------------------------
SYNAPSE_FOLDER = "Analysis_Outputs/Synapses"
LEVEL          = "level1"

#------------------------------------------------------------------------------

def build_scale(df):
    """
    Fraction of lineage-hemisphere connectors carried by each acronym.

    """
    lin_totals = (df.groupby(["Lineage", "Hemisphere"])["Count"]
                    .sum()
                    .reset_index()
                    .rename(columns={"Count": "Lineage_Connectors"}))

    acr_totals = (df.groupby(["Lineage", "Hemisphere", "Acronym"])["Count"]
                    .sum()
                    .reset_index()
                    .rename(columns={"Count": "Acronym_Connectors"}))

    scale_df = pd.merge(acr_totals, lin_totals,
                        on  = ["Lineage", "Hemisphere"],
                        how = "left")
    scale_df["Scale"] = scale_df["Acronym_Connectors"] / scale_df["Lineage_Connectors"]

    return scale_df[["Lineage", "Hemisphere", "Acronym", "Scale"]]


def build_acronym_pivot(df):
    """
    PSDs and TBars totals per lineage-hemisphere-acronym.

    """
    pivot = (df.groupby(["Lineage", "Hemisphere", "Acronym", "Connector"])["Count"]
               .sum()
               .unstack("Connector")
               .reset_index())
    pivot.columns.name = None
    return pivot


def attach_weights(paired, pivot, scale_df):
    """
    Merge acronym totals and scale, then compute all derived weight columns.

    """
    paired = pd.merge(paired, pivot,
                      on  = ["Lineage", "Hemisphere", "Acronym"],
                      how = "left")

    paired["Ratio_PSDs"]        = paired["Source_PSDs"]  / paired["PSDs"]
    paired["Ratio_TBars"]       = paired["Target_TBars"] / paired["TBars"]
    paired["Connection_Weight"] = paired["Ratio_PSDs"] * paired["Ratio_TBars"]

    paired = pd.merge(paired, scale_df,
                      on  = ["Lineage", "Hemisphere", "Acronym"],
                      how = "left")
    paired["Scaled_Weight"] = paired["Connection_Weight"] * paired["Scale"]

    return paired


#------------------------------------------------------------------------------
# Load the aggregated acronym parquet (no piece breakdown)
acr_df   = pd.read_parquet(f"{SYNAPSE_FOLDER}/Synapse-by-Acronym_{LEVEL}.parquet")
scale_df = build_scale(acr_df)
pivot    = build_acronym_pivot(acr_df)

#------------------------------------------------------------------------------
# Full compartment map (source compartment × target compartment)
inputs  = (acr_df[acr_df["Connector"] == "PSDs"]
             .groupby(["Lineage", "Hemisphere", "Acronym", "Compartment"])["Count"]
             .sum()
             .reset_index()
             .rename(columns={"Compartment": "Source_Compartment",
                               "Count":       "Source_PSDs"}))

outputs = (acr_df[acr_df["Connector"] == "TBars"]
             .groupby(["Lineage", "Hemisphere", "Acronym", "Compartment"])["Count"]
             .sum()
             .reset_index()
             .rename(columns={"Compartment": "Target_Compartment",
                               "Count":       "Target_TBars"}))

paired = pd.merge(inputs, outputs,
                  on  = ["Lineage", "Hemisphere", "Acronym"],
                  how = "outer")
paired = attach_weights(paired, pivot, scale_df)

col_order = [
    "Lineage", "Hemisphere", "Acronym",
    "Source_Compartment", "Source_PSDs",
    "Target_Compartment", "Target_TBars",
    "PSDs", "TBars",
    "Ratio_PSDs", "Ratio_TBars",
    "Connection_Weight", "Scale", "Scaled_Weight"
]
paired[col_order].to_parquet(PAIRED_FULL_PARQUET,
                             engine = "pyarrow",
                             index  = False)

#------------------------------------------------------------------------------
# Pieces compartment map (source piece+compartment × target piece+compartment)
frag_df = pd.read_parquet(f"{SYNAPSE_FOLDER}/Synapse-by-Acronym_{LEVEL}_Fragments.parquet")

inputs_p  = (frag_df[frag_df["Connector"] == "PSDs"]
               .groupby(["Lineage", "Hemisphere", "Acronym", "Piece", "Compartment"])["Count"]
               .sum()
               .reset_index()
               .rename(columns={"Piece":       "Source_Piece",
                                 "Compartment": "Source_Compartment",
                                 "Count":       "Source_PSDs"}))

outputs_p = (frag_df[frag_df["Connector"] == "TBars"]
               .groupby(["Lineage", "Hemisphere", "Acronym", "Piece", "Compartment"])["Count"]
               .sum()
               .reset_index()
               .rename(columns={"Piece":       "Target_Piece",
                                 "Compartment": "Target_Compartment",
                                 "Count":       "Target_TBars"}))

paired_p = pd.merge(inputs_p, outputs_p,
                    on  = ["Lineage", "Hemisphere", "Acronym"],
                    how = "outer")
paired_p = attach_weights(paired_p, pivot, scale_df)

col_order_p = [
    "Lineage", "Hemisphere", "Acronym",
    "Source_Piece", "Source_Compartment", "Source_PSDs",
    "Target_Piece", "Target_Compartment", "Target_TBars",
    "PSDs", "TBars",
    "Ratio_PSDs", "Ratio_TBars",
    "Connection_Weight", "Scale", "Scaled_Weight"
]
paired_p[col_order_p].to_parquet(PAIRED_PIECES_PARQUET,
                               engine = "pyarrow",
                               index  = False)

#------------------------------------------------------------------------------
print("हो गया दोस्तों!")
