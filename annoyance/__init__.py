
__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# VERSION: --------------------------------------------------------------------
__version__ = "0.0.18"


# Core module and functions: --------------------------------------------------
from ._core._SpotifyAnnoy import SpotifyAnnoy as kNN
from ._core._idx_funcs import *


# Data formatting and other interactive modules: ------------------------------
from ._interactive_funcs import *


# Data mapping functions: -----------------------------------------------------
from ._map._kNN_map_source_to_target import map_source_to_target
from ._map._iterative_kNN_smoothing import smooth