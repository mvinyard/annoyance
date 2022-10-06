
__module_name__ = "_SpotifyAnnoy.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import os, inspect
import numpy as np
import annoy
import anndata


# import local dependencies: --------------------------------------------------
from . import idx
from annoyance._interactive_funcs import *


# Main module class: ----------------------------------------------------------
class SpotifyAnnoy:
    def __init__(
        self,
        adata: anndata.AnnData,
        use_key: str = "X_pca",
        obs_key: str = None,
        n_neighbors: int = 20,
        n_trees: int = 10,
        metric: str = "euclidean",
    ):
        """
        Initialize the interactive module for building kNN graphs with AnnData
        using Spotify's super-fast Annoy Library.
        
        adata
            AnnData object for single-cell data
            type: anndata.AnnData
        
        use_key
            Key that specififes which matrix of AnnData object should be used
            for input items.
            Must be "X" or contained in adata.obsm_keys()
            type: str

        obs_key
            type: str
        
        n_neighbors [ optional, if defined above ]
            Number of neighbors by which to query the graph
            type: int
            default: None 
        
        n_trees
            number of trees 
            type: int
            default: 10
            
        metric
            distance metric; options: "angular", "euclidean", "manhattan", "hamming", or "dot".
            type: str
            default: "euclidean"
        """

        self._register_params(locals())
        self.X_ref = use_X(self._adata, use_key=self._use_key)

    def _register_params(self, local_params):

        sig = inspect.signature(SpotifyAnnoy)
        for key, value in sig.parameters.items():
            self.__setattr__("_{}".format(key), local_params[value.name])

    def build(
        self,
        subset: str = None,
        X_ref: np.ndarray = None,
        metric: str = "euclidean",
        n_trees: int = 10,
    )->annoy.AnnoyIndex:
        
        """
        Build a Spotify.annoy kNN index

        Parameters:
        -----------
        subset
            Graph will be built on the subset of data wherein: adata[adata.obs[col] == True]
            type: str
            default: None
            
        X_ref [ optional if defined during __init__ ]
            input matrix
            type: numpy.ndarray

        metric [ optional if defined during __init__ ]
            distance metric; options: "angular", "euclidean", "manhattan", "hamming", or "dot".
            type: str
            default: "euclidean"

        n_trees [ optional if defined during __init__ ]
            number of trees 
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

        # update values if necessary --------------------------------------------------
        if X_ref:
            self.X_ref = X_ref
            
        elif subset:
            adata_subset = self._adata[self._adata.obs[subset]]
            self.X_ref = use_X(adata_subset, use_key=self._use_key)
            
        if metric != self._metric:
            self._metric = metric
            
        if n_trees != self._n_trees:
            self._n_trees = n_trees        

        self.idx = idx.build(self.X_ref, metric=self._metric, n_trees=self._n_trees)

    def query(
        self, X_query: np.ndarray, n_neighbors = None, return_query: bool = True
    ) -> np.ndarray:
        
        
        """
        Query annoy index by vector - typically of unseen data.

        Parameters:
        -----------
        X_query
            Vector of items x n_dim (on which the graph was built). To query
            the neigbor graph.
            type: np.ndarray

        n_neighbors [ optional, if defined above ]
            Number of neighbors by which to query the graph
            type: int
            default: None 
            
        return_query
            Return the np.ndarray outside of the class
            type: bool
            default: True

        Returns:
        --------
        nn_query
            item indices of the nearest neighbors to X_query (from X_ref)
            type: np.ndarray
        """
        
        if not n_neighbors:
            n_neighbors = self._n_neighbors            
        
        self._nn_query = idx.query(self.idx, X_query, n_neighbors)
        if return_query:
            return self._nn_query

    def count(self, obs_key: str, max_only: bool = False) -> (list([dict, dict, ...]), float):
        
        """
        Queried kNN returns index of nearest neighbors. Query the adata.obs
        table using this index and count returned values (of a column).

        Parameters:
        -----------
        obs_key
            type: str

        Pre-defined arguments:
        ----------------------
        adata
            AnnData objects. Calls self._adata
            type: anndata.AnnData

        query_result
            result of SpotifyAnnoy.query(). Calls self._nn_query
            type: np.ndarray
        
        Returns:
        --------
        value_counts
            Value counts for each set of neighbors for a given input.
            Updated as: self._query_value_counts

        query_time
            Time taken to perform the query.
            type: float
            Updated as self._count_time


        Notes:
        ------
        (1) This implementation is much faster than a loop. Querying AnnData
            is very slow when done many times;best to do it only once.
        """
        
        self._query_value_counts, self._count_time = count_obs_neighbors(
            self._adata, self._nn_query, obs_key, max_only,
        )
        
    def load(
        self, path: str, n_features: int = None, metric: str = "euclidean"
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
            Options: "angular", "euclidean", "manhattan", "hamming", or "dot".
            type: str
            default: "euclidean"

        Returns:
        --------
        annoy_idx
            Updates existing `self.idx` with saved kNN graph.
            type: annoy.AnnoyIndex

        Notes:
        ------
        (1) `self._n_items` which is also return as a method of the annoy_idx corresponds to the number 
        of cells on which classifier was trained.
        """

        if not n_features:
            n_features = self.X_ref.shape[1]

        self.idx = idx.load(path, n_features, metric)

    def save(self, save_dir: str = "./", silent: bool = False)->None:
        
        """
        Save AnnoyIndex to /path/to/save_dir/annoyance/AnnoyIndex.ann

        Parameters:
        -----------
        save_dir
            Parant save directory path
            type: str
            default: "./"

        silent
            Silences update messages.
            type: bool
            default: False

        Returns:
        --------
        None, saves .ann file to disk.
        """
        
        n_dim = self.X_ref.shape[1]
        idx.save(self.idx, n_dim, self._n_neighbors, save_dir, silent)
