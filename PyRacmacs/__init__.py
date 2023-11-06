from rpy2.robjects.packages import importr
Racmacs = importr('Racmacs')

from .structures import (RacMap, RacOptimizerOptions, RacMergeOptions,
                         RacViewerOptions)

from . import plot_lib
from . import model_evaluation_lib

from .optimization import (make_map_from_table, optimize_map, relax_map)

from .io import read_racmap, write_racmap

from .bootstrap import (triangulation_blobs, bootstrap_map, bootstrap_blobs,
                        analyse_blobs)

from .visual import view

from .comparison import (procrustes_maps, procrustes_data, realign,
                         apply_plotspec, rotate_map, translate_map,
                         reflect_map, piecewise_procrustes)

from .merger import merge_maps, merge_tables

from .validation import dimension_test, check_hemisphering
