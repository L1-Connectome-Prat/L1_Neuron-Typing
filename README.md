# L1 Connectome Neuron Typing

Per-lineage morphology-based typing of differentiated and Small Undifferentiated
(SU) neurons in the first-instar *Drosophila* larval (L1) connectome, plus the
neuron-level nomenclature/identifier pipeline.

## Overview

Two parallel typing schemes are produced:

- **Morphological typing** — per-lineage NBLAST + Ward clustering of every
  Diff or SU neuron in the lineage, with multi-panel summary figures showing
  the lineage clone, synapse densities, dendrogram with Volker-class bars,
  compartment connectivity bubble map, and per-cluster anatomy grid.
- **Innervation-based naming** — every neuron's axon and dendrite are split
  (preferring Winding et al. tags, falling back to navis flow centrality);
  each piece is named by the compartments it innervates above a fractional
  threshold, then stitched as `dendrite-|-axon`.

Hemisphere matching is run as a sanity check: averaged left/right NBLAST
scores plus Hungarian assignment recover bilateral pairs independent of
the typing assignments.

## Pipeline

Scripts read paths via `os.environ["PROJECTS_HOME"]`. Annotations come
from `L1_Skeletons/Data_Summary/`, scoring matrices from
`L1_EM-Transforms/NBLAST/Matrices/`, compartment meshes from
`L1_Compartment-Atlas/`. Run sequentially within a category.

| #  | Script                          | Purpose                                                                                                                                                          |
|----|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1  | `identifiers_split.py`          | Pull Diff skeletons from CATMAID, split into axon/dendrite, name each piece by the compartments it innervates (`.s` input, `.b` output, `*` contralateral), stitch as `dendrite-\|-axon`. Writes `Split_Nomenclature_*` tables under `Analysis_Outputs/Nomenclature/`. |
| 2  | `organize_synapses.py`          | Read the fragment nomenclature parquets, melt TBars/PSDs into long format, and produce per-level hierarchical aggregates by neuron, acronym, and lineage. Writes six parquets per level to `Analysis_Outputs/Synapses/`. |
| 3  | `compartment_maps.py`           | Build source × target compartment connectivity maps per lineage-hemisphere-acronym unit from the Synapse-by-Acronym parquets. Writes `LineageAcronymCompartments_level1.parquet` and the axon/dendrite-split variant to `Data_Innervation/`. |
| 4  | `nblast_clusters_diff.py`       | Per-lineage NBLAST + Ward clustering of Diff neurons. Multi-panel figure per lineage; per-lineage cluster table appended to a single parquet. Includes SU neurons as visual context; optional 3D HTML plots. |
| 5  | `nblast_clusters_su.py`         | Same as #4 for SU neurons. Produces minimal figures (no bubble map or Volker bars). Lineages without SU neurons are silently skipped.                           |
| 6  | `nblast_clusters_sensory.py`    | Per-annotation clustering of sensory input bundles. Trims sensory bundle tracts from skeletons before NBLAST so bundle geometry does not dominate distances. Full multi-panel figures with synapse plots and Volker bars. |
| 7  | `hemisphere_matching.py`        | Independent sanity check — averaged L/R NBLAST + Hungarian assignment; per-lineage heatmap with matched-pair stars.                                             |
| 8  | `upload_neuron_info.py`         | Push typing/nomenclature annotations back into CATMAID with a meta-annotation tag.                                                                               |

Shared utilities live under `Scripts/viz/`:

- `_paths.py` — cross-repo paths derived from `$PROJECTS_HOME`.
- `loaders.py` — wrappers over `connectome_analysis_claude.assemble_neuron_collection` that add lineage tagging, naming, and symmetrization. Handles Diff, SU, and Sensory categories.
- `anatomy.py` — `plot_lineage_anatomy`, `synapse_plot` (consistent axes across lineages).
- `volker_bars.py` — vertical Volker-class bars alongside the dendrogram leaves.
- `bubble_map.py` — compartment bubble connectivity diagram.
- `figure_layouts.py` — gridspec layouts for the Diff, SU, and Sensory figures.
- `clustering_io.py` — `save_linkage` / `load_linkage` for the per-lineage HDF5 caches under `Linkage/`.

The clustering scripts cache each lineage's NBLAST distance matrix +
linkage + cluster labels as one HDF5 file per lineage
(`Analysis_Outputs/Neuron_Clusters_*/Linkage/<lineage>.h5`). HDF5 was
chosen over pickle so the caches are inspectable from the shell
(`h5dump`, `h5ls`) and survive a numpy/scipy major version bump.

## Directory Structure

```
L1_Neuron-Typing/
  Scripts/
    identifiers_split.py
    organize_synapses.py
    compartment_maps.py
    nblast_clusters_diff.py
    nblast_clusters_su.py
    nblast_clusters_sensory.py
    hemisphere_matching.py
    upload_neuron_info.py
    viz/                                # Shared helpers (paths, loaders, plots)
  Notebooks/                            # Per-figure exploration / stats notebooks
  Data_Innervation/                     # Compartment connectivity tables + bubble map assets
    NodeBubble.pkl                      # Bubble-map node positions/volumes
    LineageAcronymCompartments_level1.parquet
    LineageAcronymCompartments_Pieces_level1.parquet
  Analysis_Outputs/
    Nomenclature/                       # Identifier + fragment tables
    Synapses/                           # Synapse aggregates by neuron/acronym/lineage
    Neuron_Clusters_Diff/               # Per-lineage Diff clustering
      Linkage/                          # Cached linkage matrices (gitignored)
      3D_View/                          # Interactive plotly HTMLs (gitignored)
    Neuron_Clusters_SU/                 # Per-lineage SU clustering
    Neuron_Clusters_Sensory/            # Per-annotation sensory clustering
    Hemisphere_Matching/                # L/R matching tables + heatmaps
```

## Reproducing

1. Run `identifiers_split.py` → produces fragment tables under `Analysis_Outputs/Nomenclature/`.
2. Run `organize_synapses.py` → aggregates synapse counts to `Analysis_Outputs/Synapses/`.
3. Run `compartment_maps.py` → writes compartment connectivity tables to `Data_Innervation/`.
4. Run `nblast_clusters_diff.py` (set `re_blast = True` on first pass to compute and cache linkages; subsequent runs reload from HDF5).
5. Run `nblast_clusters_su.py` for SU neurons in lineages that have them.
6. Run `nblast_clusters_sensory.py` for sensory input bundle typing.
7. Optional: `hemisphere_matching.py` for the bilateral sanity check.
8. Optional: `upload_neuron_info.py` to push annotations back to CATMAID.
