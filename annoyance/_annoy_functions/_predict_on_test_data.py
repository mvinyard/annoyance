
# _predict_on_test_data.py
__module_name__ = "_predict_on_test_data.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages #
# --------------- #
from collections import Counter


def _predict_on_test_data(annoy_idx, x_test, y_train, n_neighbors):

    y_predicted = []
    for cell in range(x_test.shape[0]):
        nearest_neighbors = y_train.iloc[
            annoy_idx.get_nns_by_vector(x_test[cell], n_neighbors)
        ]
        nearest_neighbors = Counter(nearest_neighbors.values).most_common(2)
        y_predicted.append(nearest_neighbors[0][0])

    return y_predicted