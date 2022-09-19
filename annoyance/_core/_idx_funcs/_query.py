
__module_name__ = "_query.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import annoy
import numpy as np


# supporting function: --------------------------------------------------------
def _vector_query(
    annoy_idx: annoy.AnnoyIndex, query_vector: np.ndarray, n_neighbors: int
):
    """
    Parameters:
    -----------
    annoy_idx
        type: Annoy.AnnoyIndex

    query_vector
        type: np.ndarray

    n_neighbors
        type: int

    Returns:
    --------
    nn_idx
        Nearest neighbors of the passed vector from the annoy_idx.
        type: np.ndarray
    """
    return annoy_idx.get_nns_by_vector(query_vector, n_neighbors)


# primary function: -----------------------------------------------------------
def query(
    annoy_idx: annoy.AnnoyIndex, X_query: np.ndarray, n_neighbors: int = 20
) -> np.ndarray:

    """
    Query annoy index by vector.

    Parameters:
    -----------
    annoy_idx
        Spotify.AnnoyIndex of nearest neighbors built on X_ref.
        type: Annoy.AnnoyIndex

    X_query
        Vector of items x n_dim (on which the graph was built). To query
        the neigbor graph.
        type: np.ndarray

    n_neighbors
        Number of neighbors by which to query the graph.
        type: int
        default: 20

    Returns:
    --------
    nn_query
        item indices of the nearest neighbors to X_query (from X_ref)
        type: np.ndarray
    """

    nn_idx_list = []

    query_size = X_query.shape[0]
    for i in range(query_size):
        nn_idx = _vector_query(annoy_idx, X_query[i], n_neighbors)
        nn_idx_list.append(nn_idx)

    return np.array(nn_idx_list)