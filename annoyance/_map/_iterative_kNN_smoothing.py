import anndata
import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import vinplots
import matplotlib.pyplot as plt


from .._core._SpotifyAnnoy import SpotifyAnnoy as kNN
from ._annotate_adata_bool_idx import add_bool_idx

def _plot_smoothing(adata, score_idx, scores, n_iters):

    X_umap = adata[score_idx].obsm["X_umap"]

    nplots = n_iters + 1
    fig, axes = vinplots.quick_UMAP(nplots=nplots, ncols=nplots)

    c_idx = np.argsort(scores[0])
    axes[0].scatter(
        X_umap[:, 0], X_umap[:, 1], c="lightgrey", s=2, rasterized=True, zorder=0
    )
    axes[0].scatter(
        X_umap[c_idx, 0],
        X_umap[c_idx, 1],
        c=scores[0][c_idx],
        rasterized=True,
        zorder=1,
    )
    axes[0].set_title("Mapped (unsmoothed)")

    for i in range(1, nplots):
        axes[i].scatter(
            X_umap[:, 0],
            X_umap[:, 1],
            c="lightgrey",
            rasterized=True,
            zorder=0,
        )
        c_idx = np.argsort(scores[i])
        axes[i].scatter(
            X_umap[c_idx, 0],
            X_umap[c_idx, 1],
            c=scores[i][c_idx],
            rasterized=True,
            zorder=1,
        )
        axes[i].set_title("iter: {}".format(i))
    plt.show()


def smooth(
    adata: anndata.AnnData,
    score: pd.Series,
    plot=False,
    n_iters: int = 5,
    n_neighbors: int = 20,
    func: bool = False,
) -> np.ndarray:

    idx, score = score.index.astype(int), score.values.astype(float)
    curr = np.zeros(score.shape)
    scores, prev = [score], score

    add_bool_idx(adata, idx, key_added="tmp_bool", silent=True)

    kNN_Graph = kNN(adata)
    kNN_Graph.build(subset="tmp_bool")

    X_query = adata[idx].obsm["X_pca"]

    for _ in tqdm(range(n_iters)):
        X_nn = kNN_Graph.query(X_query)
        curr = prev[X_nn]
        if func:
            curr = func(curr)
        else:
            curr = curr.mean(1)
        prev = curr
        scores.append(curr)

    if plot:
        _plot_smoothing(adata, idx, scores, n_iters)

    return curr