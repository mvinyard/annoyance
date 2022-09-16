
__module_name__ = "_binominal_test_train_assignment.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import anndata
import numpy as np


# primary function: -----------------------------------------------------------
def binominal_test_train_assignment(adata: anndata.AnnData, p: float = 0.5) -> None:

    """
    Binomially assign test/train split.

    Parameters:
    -----------
    adata
        AnnData object
        type: anndata.AnnData

    p
        probability of training-set selection.
        type: float
        default: 0.5

    Returns:
    --------
    None
        adata is modified in-place, updated with adata.obs['kNN_test'] and adata.obs['kNN_train']
    """

    obs_df = adata.obs.copy()

    train = np.random.binomial(1, p, len(obs_df)).astype(bool)
    test = np.invert(train)

    obs_df["kNN_train"] = train
    obs_df["kNN_test"]  = test

    adata.obs = obs_df