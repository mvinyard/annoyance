import anndata
import pandas as pd
import numpy as np
import time


def _count_values(col: pd.Series) -> dict:
    return col.value_counts().to_dict()

def _max_count(col: pd.Series) -> str:
    return col.value_counts().idxmax()

def count_obs_neighbors(
    adata: anndata.AnnData, query_result: np.ndarray, obs_key: str, max_only=False,
) -> (list([dict, dict, ...]), float):

    """
    Queried kNN returns index of nearest neighbors. Query the adata.obs
    table using this index and count returned values (of a column).
    
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    query_result
        result of SpotifyAnnoy.query()
        type: np.ndarray

    obs_key
        type: str
        
    Returns:
    --------
    value_counts
        Value counts for each set of neighbors for a given input.
    
    query_time
        Time taken to perform the query.
        type: float
    
    
    Notes:
    ------
    (1) This implementation is much faster than a loop. Querying AnnData
        is very slow when done many times;best to do it only once.
    """

    query_time_start = time.time()

    nn_adata = adata[query_result.flatten()]
    # nn_idx = np.repeat(range(len(query_result)), 20)
    query_df = pd.DataFrame(
        nn_adata.obs[obs_key].to_numpy().reshape(query_result.shape).T
    )
    
    if max_only:
        value_counts = [
            _max_count(query_df[i]) for i in query_df.columns
        ] # list of values
    else:
        value_counts = [
            _count_values(query_df[i]) for i in query_df.columns
        ]  # list of dicts

    query_time = time.time() - query_time_start

    del nn_adata

    return value_counts, query_time