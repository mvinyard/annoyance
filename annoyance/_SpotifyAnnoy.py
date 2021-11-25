
# _SpotifyAnnoy.py
__module_name__ = "_SpotifyAnnoy.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)

# import packages #
# --------------- #
import os, annoy


# import local modules #
# -------------------- #
from ._annoy_functions._evaluate_annoy import _evaluate_annoy
from ._annoy_functions._build_annoy_index import _build_annoy_index
from ._annoy_functions._predict_on_test_data import _predict_on_test_data

from ._utility_functions._flexible_multilevel_mkdir import _flexible_multilevel_mkdir
from ._utility_functions._DataSplit import _DataSplit

class _SpotifyAnnoy:
    def __init__(self):

        """
        Notes:
        ------
        The following notes were taken directly from https://github.com/spotify/annoy

        AnnoyIndex(f, metric) returns a new index that's read-write and stores vector
        of f dimensions. Metric can be "angular", "euclidean", "manhattan", "hamming",
        or "dot".

        a.add_item(i, v) adds item i (any nonnegative integer) with vector v. Note
        that it will allocate memory for max(i)+1 items.

        a.build(n_trees, n_jobs=-1) builds a forest of n_trees trees. More trees gives
        higher precision when querying. After calling build, no more items can be
        added. n_jobs specifies the number of threads used to build the trees.
        n_jobs=-1 uses all available CPU cores.

        a.save(fn, prefault=False) saves the index to disk and loads it (see next
        function). After saving, no more items can be added.

        a.load(fn, prefault=False) loads (mmaps) an index from disk. If prefault is set
        to True, it will pre-read the entire file into memory (using mmap with
        MAP_POPULATE). Default is False.
        """

        self.data = _DataSplit()

    def split_data(
        self, adata, annot_col, labels, on=False, other="other", silent=False
    ):
        
        """
        Split the data on labels into test and training sets. 
        
        Parameters:
        -----------
        adata
            type: AnnData
        
        annot_col
            column in `adata.obs`
            type: str
        
        labels
            Categories in `annot_col`
            type: list
            
        on
            String indicating one item within an adata.obsm_keys()
            default: False
            type: bool or str
        
        other
            default: "other"
            type: str
        
        silent
            default: False
            type: bool
        
        """
        
        self.data.mask(adata, annot_col, labels, other, silent)
        self.X = self.data.split(adata, on)
        self.y_mask = self.data.y_mask

    def prebuild(
        self,
        n_trees=10,
        n_neighbors=20,
        metric="euclidean",
        silent=False,
        return_f1=False,
    ):

        """
        Build and test the model on split train/testing data. 
        
        Parameters:
        -----------
        n_trees
            default: 10
            type: int
        
        n_neighbors
            default: 20
            type: int
        
        metric
            default: "euclidean"
            type: str
        
        silent
            default: False
            type: bool
            
        return_f1
            default: False
            type: bool
        
        """

        self.n_trees = n_trees
        self.n_neighbors = n_neighbors
        self.metric = metric

        self.pre_annoy_idx = _build_annoy_index(
            self.data.x_train, self.metric, self.n_trees
        )

        f1 = _evaluate_annoy(
            self.pre_annoy_idx,
            self.data.x_test,
            self.data.y_train,
            self.data.y_test,
            self.n_neighbors,
            silent,
            return_f1,
        )

        if return_f1:
            return f1

    def build(self):

        """
        Builds 
        
        Parameters:
        -----------
        All parameters were passed in the prebuild phase. 
        
        Returns:
        --------
        None
        """

        self.annoy_idx = _build_annoy_index(self.data.X, self.metric, self.n_trees)

    def save(self, outdir=""):

        """
        Save AnnoyIndex to files (.ann and .txt). 
        
        Parameters:
        -----------
        outdir
            directory path
            default: ""
            type: str
        
        Returns:
        --------
        None
        """

        self._file_basename = (
            "annoyance.annoy_idx.{}_dims.{}_neighbors.{}_trees".format(
                self.data.X.shape[1], self.n_neighbors, self.n_trees
            )
        )

        self._ann_filename = "{}.ann".format(self._file_basename)
        self._txt_filename = "{}.txt".format(self._file_basename)

        self._outdir = os.path.join(outdir, "annoyance_index")
        _flexible_multilevel_mkdir(self._outdir)
        self._ann_path = os.path.join(self._outdir, self._ann_filename)
        self._txt_path = os.path.join(self._outdir, self._txt_filename)

        self.annoy_idx.save(self._ann_path)
        with open(self._txt_path, "w") as f:
            for y_val in self.y_mask.values:
                f.write(y_val + "\n")

        print(
            "Annoyance annoy index saved to:\n\n\t{}\n\t{}".format(
                self._ann_path, self._txt_path
            )
        )