
__module_name__ = "_create_label_mask.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -------------------------------------------------------
import licorice_font
import pandas as pd
from anndata import AnnData


# supporting functions: --------------------------------------------------
def _print_label_mask_report(y_mask):
    msg = licorice_font.font_format("Label mask created", ["BOLD"])
    print("{}:\n\n{}".format(msg, y_mask.value_counts()))


# main function: ---------------------------------------------------------
def create_label_mask(
    adata: AnnData, col: str, labels: list, other: str = "other", silent: bool = False
) -> pd.Series:

    """

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    col
        type: str

    labels
        type: list of str

    other
        default: "other"

    silent
        type: bool
        default: False

    Returns:
    --------
    y_mask
        type: pandas.core.series.Series
    """

    # isolate col
    y = adata.obs[col]

    # add the new "other" category
    y_ = y.cat.add_categories(other)

    # rename categories outside of labels as others
    y_ = y_.mask(~y_.isin(labels), other=other)

    # remove the now empty categories that were renamed as others
    y_mask = y_.cat.remove_unused_categories()

    if not silent:
        _print_label_mask_report(y_mask)

    return y_mask