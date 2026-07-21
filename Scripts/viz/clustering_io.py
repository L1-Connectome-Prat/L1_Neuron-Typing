"""
Persistent storage for per-lineage NBLAST clustering results.

One HDF5 file per lineage under Analysis_Outputs/<run_folder>/Linkage/. Files
are inspectable from the shell with `h5dump`/`h5ls` and survive a numpy /
scipy version bump that would break a pickle.

Layout
------
    /Z                           - Ward linkage matrix (n-1, 4)
    /nbl_dist/values             - distance matrix as a 2D float array
    /nbl_dist/index              - matching row/col labels (string)
    /cluster_labels              - 1D int array of cluster ids per neuron
    /neuron_ids                  - 1D string array of neuron ids
    root attrs:
        optimal_t                - distance cut threshold
        threshold                - input threshold passed to clustering

Use the two top-level functions; the module is intentionally small.

"""
import numpy as np
import pandas as pd
import h5py

#------------------------------------------------------------------------------
__all__ = ["save_linkage", "load_linkage"]


#------------------------------------------------------------------------------
def save_linkage(info_dict, path):
    """
    Persist a clustering info_dict + linkage + distance matrix to HDF5.

    Expects info_dict with keys: Z, nbl_dist, cluster_labels, neuron_ids,
    optimal_t (and optionally threshold).

    """
    nbl_dist = info_dict["nbl_dist"]

    with h5py.File(path, "w") as f:
        f.create_dataset("Z", data = np.asarray(info_dict["Z"]))

        grp = f.create_group("nbl_dist")
        grp.create_dataset("values", data = nbl_dist.values)
        grp.create_dataset("index",  data = np.asarray(nbl_dist.index, dtype = "S"))

        f.create_dataset("cluster_labels",
                         data = np.asarray(info_dict["cluster_labels"]))
        f.create_dataset("neuron_ids",
                         data = np.asarray(info_dict["neuron_ids"], dtype = "S"))

        f.attrs["optimal_t"] = float(info_dict["optimal_t"])
        if "threshold" in info_dict:
            f.attrs["threshold"] = float(info_dict["threshold"])


#------------------------------------------------------------------------------
def load_linkage(path):
    """
    Reload an info_dict previously written by save_linkage.

    """
    with h5py.File(path, "r") as f:
        Z = f["Z"][...]

        index    = np.asarray(f["nbl_dist/index"][...]).astype(str)
        values   = f["nbl_dist/values"][...]
        nbl_dist = pd.DataFrame(values, index = index, columns = index)

        cluster_labels = f["cluster_labels"][...]
        neuron_ids     = np.asarray(f["neuron_ids"][...]).astype(str)
        optimal_t      = float(f.attrs["optimal_t"])

        info_dict = {"Z"              : Z,
                     "nbl_dist"       : nbl_dist,
                     "cluster_labels" : cluster_labels,
                     "neuron_ids"     : neuron_ids,
                     "optimal_t"      : optimal_t}

        if "threshold" in f.attrs:
            info_dict["threshold"] = float(f.attrs["threshold"])

    return info_dict
