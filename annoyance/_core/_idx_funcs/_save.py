
__module_name__ = "_save.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: ------------------------------------------------------------
import licorice_font
import annoy
import pydk
import os


# supporting functions: -------------------------------------------------------
def _save_message(ann_fname):

    note = licorice_font.font_format("NOTE", ["BLUE"])
    msg = licorice_font.font_format("Saving annoy idx to", ["BOLD"])
    print(" - [{}] | {}:\n\n\t{}".format(note, msg, ann_fname))


def _mk_idx_fname(dim, n_neighbors, n_trees):
    msg = "annoyance.annoy_idx.{}_dim.{}_neighbors.{}_trees"
    return msg.format(dim, n_neighbors, n_trees)


# main function: --------------------------------------------------------------
def save(
    idx: annoy.AnnoyIndex,
    n_dim: int,
    n_neighbors: int,
    save_dir: str = "./",
    silent=False,
) -> None:
    
    """
    Save AnnoyIndex to /path/to/save_dir/annoyance/AnnoyIndex.ann

    Parameters:
    -----------
    annoy_idx
        Spotify.AnnoyIndex of nearest neighbors built on X_ref.
        type: Annoy.AnnoyIndex
    
    n_dim
        type: int
    
    n_neighbors
        Number of neighbors by which to query the graph.
        type: int
        default: 20

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

    annoy_save_dir = os.path.join(save_dir, "annoyance")
    pydk.mkdir_flex(annoy_save_dir)

    base_fname = _mk_idx_fname(n_dim, n_neighbors, idx.get_n_trees())
    annoy_idx_base_path = os.path.join(annoy_save_dir, base_fname)
    ann_fname = "{}.ann".format(annoy_idx_base_path)
    idx.save(ann_fname)
    if not silent:
        _save_message(ann_fname)