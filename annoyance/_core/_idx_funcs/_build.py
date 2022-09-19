
__module_name__ = "_build.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import annoy


def build(
    x, metric: str = "euclidean", n_trees: int = 10
) -> annoy.AnnoyIndex:
    
    """
    Build a Spotify.annoy kNN index
    
    Parameters:
    -----------
    x
        input matrix
        type: numpy.ndarray
    
    metric
        distance metric
        type: str
        default: "euclidean"
    
    n_trees
        number of trees.
        type: int
        default: 10
        
    Returns:
    --------
    annoy_idx
        type: annoy.AnnoyIndex
    
    Notes:
    ------
    (1) Loops through each item. Appends n_dim feature vector to AnnoyIndex.
    """    

    annoy_idx = annoy.AnnoyIndex(x.shape[1], metric)

    for i in range(x.shape[0]):
        annoy_idx.add_item(i, x[i])

    annoy_idx.build(n_trees)

    return annoy_idx
