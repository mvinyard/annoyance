
__module_name__ = "_kNN_map_source_to_target.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import pandas as pd
import numpy as np
import anndata


# import local dependencies: --------------------------------------------------
from .._core._SpotifyAnnoy import SpotifyAnnoy as kNN


# Supporting functions: -------------------------------------------------------
def _add_source_target_exclude_col(
    adata,
    source_idx,
    target_idx,
):

    tmp = np.full(len(adata), "exclude")
    tmp[source_idx] = "source"
    tmp[target_idx] = "target"

    adata.obs["kNN_mapping"] = pd.Categorical(tmp)

    adata.obs["source"] = adata.obs["kNN_mapping"] == "source"
    adata.obs["target"] = adata.obs["kNN_mapping"] == "target"

    adata.obs.drop("kNN_mapping", axis=1, inplace=True)


def _source_target_query(adata, use_key):

    """build and query kNN index"""

    G = kNN(adata)
    G.build(subset="source")
    X_query = adata[adata.obs["target"]].obsm[use_key].toarray()

    return G.query(X_query)


def _merge_ref_map_values(adata, obs_key, X_nn):

    ref_values = adata[adata.obs["source"]].obs[obs_key].values
    map_values = ref_values[X_nn].mean(1)

    tmp = np.full(len(adata), None)
    tmp[adata.obs["source"]] = ref_values
    tmp[adata.obs["target"]] = map_values
    adata.obs["{}_mapped".format(obs_key)] = tmp

    adata.obs.drop(["source", "target"], axis=1, inplace=True)


# Main module function: -------------------------------------------------------
def map_source_to_target(
    adata: anndata.AnnData,
    source_idx: pd.Index,
    target_idx: pd.Index,
    obs_key: str,
    use_key="X_pca",
):

    """
    Using a kNN graph, map a value in a source set of cells to a target set of cells.

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    source_idx [ required ]
        type: pandas.Index

    target_idx [ required ]
        type: pandas.Index

    obs_key [ required ]
        type: str

    use_key
        type: str
        default: "X_pca"

    Returns:
    --------
    None, updates adata.obs with "{obs_key}_mapped"

    Notes:
    ------
    (1) annotates source/target in the adata.obs table
    (2) build (based on the source) a kNN graph and query (using the target)
    (3) merge mapped / references values
    """
    _add_source_target_exclude_col(adata, source_idx, target_idx)
    X_nn = _source_target_query(adata, use_key)
    _merge_ref_map_values(adata, obs_key, X_nn)
    
