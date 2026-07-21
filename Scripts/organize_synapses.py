"""
Organizing the synapses assignments by neuron and acronyms at different levels.

Generated from the Identifiers_Split script.

This script:
- Reads neuron metadata (bodyId, Lineage, Hemisphere, Acronym)
- Loads fragmented nomenclature Parquet files for level1/level2
- Melts TBars/PSDs into long format
- Merges metadata
- Produces hierarchical aggregates (Neuron / Acronym / Lineage)
- Saves both fragmented (by Piece) and full (Piece-aggregated) Parquet files

"""

import pandas as pd
from Scripts.viz._paths import NEURON_INFO_CSV

#------------------------------------------------------------------------------
# Output
FOLDER = "Analysis_Outputs/Synapses"
#------------------------------------------------------------------------------
def groupby_piece(df):
    """
    Get the columns and grouby without the piece.

    """
    cols = df.columns.drop(["Piece", "Count"]).tolist()
    return df.groupby(cols)["Count"].sum().reset_index()

def organize_by_level(df, level):
    """
    Organize the dataframe hierarchically.

    """
    # Copy
    copy_df = df.copy()
        
    # Grouping by Acronym
    acronym_df = copy_df.groupby([
        "Lineage",
        "Acronym",
        "Piece",
        "Hemisphere", 
        "Compartment", 
        "Connector"
        ])["Count"].sum().reset_index()

    # Grouping by lineage
    lineage_df = copy_df.groupby([
        "Lineage",
        "Piece",
        "Hemisphere", 
        "Compartment", 
        "Connector"
        ])["Count"].sum().reset_index()

    # Group by piece 
    neuron_full_df = groupby_piece(copy_df)
    acronym_full_df = groupby_piece(acronym_df)
    lineage_full_df = groupby_piece(lineage_df)

    # Save the Fragmented dfs
    copy_df.to_parquet(f"{FOLDER}/Synapse-by-Neuron_{level}_Fragments.parquet",
                        engine = "pyarrow", 
                        index = False)
    acronym_df.to_parquet(f"{FOLDER}/Synapse-by-Acronym_{level}_Fragments.parquet", 
                        engine = "pyarrow", 
                        index = False)
    lineage_df.to_parquet(f"{FOLDER}/Synapse-by-Lineage_{level}_Fragments.parquet", 
                        engine = "pyarrow", 
                        index = False)

    # Save the full dfs
    neuron_full_df.to_parquet(f"{FOLDER}/Synapse-by-Neuron_{level}.parquet",
                             engine = "pyarrow", 
                             index = False)
    acronym_full_df.to_parquet(f"{FOLDER}/Synapse-by-Acronym_{level}.parquet",
                             engine = "pyarrow", 
                             index = False)
    lineage_full_df.to_parquet(f"{FOLDER}/Synapse-by-Lineage_{level}.parquet",
                             engine = "pyarrow", 
                             index = False)

#------------------------------------------------------------------------------
# Read in the neuron information
info_df = pd.read_csv(NEURON_INFO_CSV)
info_df = info_df[["bodyId", "Lineage", "Hemisphere", "Acronym"]]


# Loop through both the levels
levels = ["level1", "level2"]

# Run for each level
for level in levels:
    # Status 
    print(f"Organizing synapses for {level}...\n")

    # Load the dataframe
    syn_df = pd.read_parquet(f"Analysis_Outputs/Nomenclature/Split_Nomenclature_Fragments_{level}.parquet")
    # Melt to organize it correctly
    syn_df = syn_df.melt(
        id_vars     = ["bodyId", "Piece", "Compartment"],
        value_vars  = ["TBars", "PSDs"],
        var_name    = "Connector",
        value_name  = "Count"
    )

    # Merge with the info dataframe
    syn_df = pd.merge(syn_df, info_df, on = "bodyId", how = "left")

    # Orgnaize as neuron, acronym, lineage 
    organize_by_level(syn_df, level)

#------------------------------------------------------------------------------
print("हो गया दोस्तों!")