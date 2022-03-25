
__module_name__ = "_load_annoy_index.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
from annoy import AnnoyIndex


def _load_annoy_index(path, n_features, metric="euclidean"):
    
    """
    Load a pre-existing AnnoyIndex.
    
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
        AnnoyIndex
        type: annoy.Annoy
        
    Notes:
    ------
    """
    
    annoy_idx = AnnoyIndex(n_features, metric)
    annoy_idx.load(path)
    
    return annoy_idx
