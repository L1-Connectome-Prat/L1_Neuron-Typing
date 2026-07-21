"""
Cross-repo path resolution. All paths flow from $PROJECTS_HOME so this
repository never carries a copy of annotations, skeletons, or scoring matrices.

"""
import os

#------------------------------------------------------------------------------
PROJECTS_HOME = os.environ["PROJECTS_HOME"]

# L1_Skeletons - per-neuron annotations and pre-organized skeletons
SKELETONS_HOME    = f"{PROJECTS_HOME}/L1_Skeletons"
SKELETONS_FOLDER  = f"{SKELETONS_HOME}/Skeletons"
DATA_SUMMARY      = f"{SKELETONS_HOME}/Data_Summary"
NEURON_INFO_CSV   = f"{DATA_SUMMARY}/L1_NeuronInformation_PK.csv"
LINEAGES_NPY      = f"{DATA_SUMMARY}/L1_Annotations-Lineages_PK.npy"
ANNOTATIONS_NPY   = f"{DATA_SUMMARY}/L1_Annotations_PK.npy"

# L1_EM-Transforms - NBLAST scoring matrices
NBLAST_MATRIX_DIR = f"{PROJECTS_HOME}/L1_EM-Transforms/NBLAST/Matrices"

# L1_Compartments - mushroom body and central brain meshes
COMPARTMENTS_HOME = f"{PROJECTS_HOME}/L1_Compartment-Atlas"
MB_RIGHT_STL      = f"{COMPARTMENTS_HOME}/Central-Brain_Compartments/MB/MB_Right.stl"
MB_LEFT_STL       = f"{COMPARTMENTS_HOME}/Central-Brain_Compartments/MB/MB_Left.stl"

# Local layout asset - bubble map node positions / volumes
NODE_BUBBLE_PKL   = "Data_Innervation/NodeBubble.pkl"

# Compartment connectivity tables produced by compartment_maps.py
PAIRED_FULL_PARQUET   = "Data_Innervation/LineageAcronymCompartments_level1.parquet"
PAIRED_PIECES_PARQUET = "Data_Innervation/LineageAcronymCompartments_Pieces_level1.parquet"


#------------------------------------------------------------------------------
def nblast_matrix_path(smat_val = 250):
    """
    Return the trained NBLAST scoring matrix path for a given dot tangent.

    """
    return f"{NBLAST_MATRIX_DIR}/NBLAST_smat_{smat_val}.csv"


#------------------------------------------------------------------------------
def skeleton_folder(annot, category = "Diff"):
    """
    Path to the per-lineage skeleton folder. Caller checks existence.

    """
    annot_name = annot.replace("/", "")
    return f"{SKELETONS_FOLDER}/{category}_Skeletons/{annot_name}"


#------------------------------------------------------------------------------
def tract_folder(annot, hemisphere = "Right"):
    """
    Path to the per-lineage tract folder.

    """
    annot_name = annot.replace("/", "")
    return f"{SKELETONS_FOLDER}/Diff_Tracts/{annot_name}/{hemisphere}"
