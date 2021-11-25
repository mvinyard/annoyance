
# _build_annoy_index.py
__module_name__ = "_build_annoy_index.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)


# import packages #
# --------------- #
import annoy


def _build_annoy_index(x_array, metric, n_trees):

    """
    Build an AnnoyIndex using Spotify.annoy.
    
    Parameters:
    -----------
    x_array
        input array
        type: numpy.ndarray
    
    metric
        metric passed to annoy. "euclidean" is typical. 
        type: str
    
    n_trees
        number of trees in the AnnoyIndex.
        type: int
        
    Returns:
    --------
    annoy_idx
        type: annoy.AnnoyIndex
    
    Notes:
    ------
    (1) Loops through each cell, appending the n_dimensional feature vector to the AnnoyIndex. 
    """

    annoy_idx = annoy.AnnoyIndex(x_array.shape[1], metric)

    for cell in range(x_array.shape[0]):
        annoy_idx.add_item(cell, x_array[cell])

    annoy_idx.build(n_trees)

    return annoy_idx
