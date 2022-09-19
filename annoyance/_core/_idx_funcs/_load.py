
__module_name__ = "_load.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import annoy


# primary function: -----------------------------------------------------------
def load(
    path: str, n_features: int, metric: str = "euclidean"
) -> annoy.AnnoyIndex:

    """
    Load a previously saved AnnoyIndex.

    Parameters:
    -----------
    path
        file path  to annoy idx (.ann)

    n_features
        type: int

    metric
        type: str
        default: "euclidean"

    Returns:
    --------
    annoy_idx
        kNN graph
        type: annoy.AnnoyIndex

    Notes:
    ------
    """

    annoy_idx = annoy.AnnoyIndex(n_features, metric)
    annoy_idx.load(path)

    return annoy_idx