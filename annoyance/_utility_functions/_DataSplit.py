
# _DataSplit.py
__module_name__ = "_DataSplit.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)

# import packages #
# --------------- #
from sklearn.model_selection import train_test_split


def _create_label_mask(adata, annot_col, labels, other="other", silent=False):

    """"""

    y = adata.obs[annot_col]

    y_mask = (
        y.cat.add_categories(other)
        .mask(~y.isin(labels), other=other)
        .cat.remove_unused_categories()
    )

    if not silent:
        print("label mask created:\n\n{}".format(y_mask.value_counts()))

    return y_mask

class _DataSplit:
    def __init__(self):

        """"""

        self.y_mask = None

    def mask(self, adata, annot_col, labels, other="other", silent=False):

        """"""

        self.y_mask = _create_label_mask(adata, annot_col, labels, other, silent)

    def split(self, adata, on=False, y_mask=None):

        if y_mask != None:
            self.y_mask = y_mask

        if on:
            assert on in adata.obsm_keys(), print("Must be in adata.obsm_keys()")
            self.X = _X = adata.obsm[on]
        else:
            self.X = _X = adata.X
        [x_train, x_test, y_train, y_test] = train_test_split(
            _X, self.y_mask, stratify=self.y_mask
        )

        self.x_train, self.x_test = x_train, x_test
        self.y_train, self.y_test = y_train, y_test