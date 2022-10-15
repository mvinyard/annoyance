import licorice_font
import numpy as np

def _print_key_added(key_added):

    note = licorice_font.font_format("NOTE", ["BLUE"])
    key_print = licorice_font.font_format(key_added, ["BOLD"])
    print(" - [{}] | obs_key added: {}".format(note, key_print))
    
    
    
def add_bool_idx(adata, idx, key_added="bool_idx", silent=False):
    
#     _annotate_adata_bool_idx

    tmp = np.zeros(len(adata))
    tmp[idx] = 1
    adata.obs[key_added] = tmp.astype(bool)

    if not silent:
        _print_key_added(key_added)
