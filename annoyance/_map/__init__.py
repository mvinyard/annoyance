
__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])



# specify version: -----------------------------------------------------------------------
__version__ = "0.0.18"

# Modules and functions: ------------------------------------------------------
from ._kNN_map_source_to_target import map_source_to_target
from ._iterative_kNN_smoothing import smooth
from ._annotate_adata_bool_idx import add_bool_idx