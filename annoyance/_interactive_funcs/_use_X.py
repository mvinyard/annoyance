
__module_name__ = "_use_X.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -------------------------------------------------------
import licorice_font
import numpy


# supporting functions: --------------------------------------------------
def _is_numpy_array(x):
    return x.__class__ is numpy.ndarray


def _as_numpy_array(x):
    if _is_numpy_array(x):
        return x
    else:
        return x.toarray()

def _print_help_message(adata, use_key):
    msg = " - [{}] | {} is not a suitable array key for this dataset. {}:\n".format(
        licorice_font.font_format("NOTE", ["BLUE"]),
        licorice_font.font_format(use_key, ["BOLD", "RED"]),
        licorice_font.font_format("Choose from", ["BOLD"]),
    )
    print(msg)
    available = [licorice_font.font_format(key, ["BOLD"]) for key in adata.obsm_keys()]
    licorice_font.underline("Counts matrix:", ["BOLD", "GREEN"], n_newline=0)
    print(" - {}\n".format(licorice_font.font_format("X", ["BOLD"])))
    licorice_font.underline("adata.obsm_keys:", ["BOLD", "GREEN"], n_newline=0)
    for obsm_key in available:
        print(" - {}".format(obsm_key))
        
# main function: fetch data (X_ref): -------------------------------------
def use_X(adata, use_key="X"):

    """
    Return data from AnnData as numpy array (if not already).

    Parameters:
    -----------
    adata
        AnnData object
        type: anndata._core.anndata.AnnData

    use_key
        "X" or key in adata.obsm_keys(), adata.obs_keys(),
        type: str
        default: "X"

    Returns:
    --------
    X
        Formatted data matrix.
        type: np.ndarray
    """

    if use_key == "X":
        return _as_numpy_array(adata.X)

    elif use_key in adata.obsm_keys():
        return _as_numpy_array(adata.obsm[use_key])

    else:
        _print_help_message(adata, use_key)
