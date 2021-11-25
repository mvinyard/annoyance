
# _evaluate_annoy.py
__module_name__ = "_evaluate_annoy.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(
    [
        "vinyard@g.harvard.edu",
    ]
)

# import packages #
# --------------- #
from sklearn.metrics import classification_report
from sklearn.metrics import f1_score

# import local modules #
# -------------------- #
from ._predict_on_test_data import _predict_on_test_data

def _evaluate_annoy(annoy_idx, x_test, y_train, y_test, n_neighbors, silent, return_f1):

    y_predicted = _predict_on_test_data(annoy_idx, x_test, y_train, n_neighbors)
    idx_build_report = classification_report(y_test, y_predicted)
    f1_metric = f1_score(y_test, y_predicted, average="macro")

    if not silent:
        print(idx_build_report)
        print("f1 score: {}".format(f1_metric))
    if return_f1:
        return f1_metric
