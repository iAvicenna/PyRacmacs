#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 18:38:05 2023

@author: avicenna
"""

import numpy as np
import time
from . import (Racmacs, RacMap, RacOptimizerOptions)
from .model_evaluation_lib.local_utils import DimensionTestSummaryStatistics
from rpy2 import robjects as ro
from .io import capture_r_output
from .utils.conversion import odict_to_dict


def dimension_test(
        racmap: RacMap,
        dimensions: list=None,
        test_proportion: float=0.1,
        minimum_column_basis:  str="none",
        fixed_column_bases: list=None,
        number_of_optimizations: int=1000,
        validations_per_dimension: int=100,
        options: list=None,
        verbose = True
        ):


    if dimensions is None:
        dimensions = ro.FloatVector([1.0,2.0,3.0,4.0,5.0])
    else:
        dimensions = ro.FloatVector([float(x) for x in dimensions]) # Racmacs seems to require float for this

    if fixed_column_bases is None:
        fixed_column_bases = ro.FloatVector([np.nan for _ in range(racmap.num_sera)])
    else:
        fixed_column_bases = ro.FloatVector(fixed_column_bases)

    if options is not None:
        options = ro.StrVector(options)
    else:
        options = ro.ListVector([])


    stoud, stderr = capture_r_output(verbose, False)

    t0 = time.time()

    dimension_test_result_R = Racmacs.dimensionTestMap(
        racmap._acmap_R,
        dimensions = dimensions,
        test_proportion = test_proportion,
        minimum_column_basis = minimum_column_basis,
        fixed_column_bases = fixed_column_bases,
        number_of_optimizations = number_of_optimizations,
        replicates_per_dimension = validations_per_dimension,
        options = options
        )

    t1 = time.time()
    print(f'\n{t1-t0:.2f} seconds.\n')

    dimension_test_result_R = np.reshape(
        np.array(dimension_test_result_R[1:]),
        (5,len(dimensions))
        ).T


    dimension_test_result = DimensionTestSummaryStatistics(dimensions,
                                                           validations_per_dimension,
                                                           dimension_test_result_R)


    return dimension_test_result

def check_hemisphering(racmap, optimization_number=0, grid_spacing=0.25,
                       stress_lim=0.1, options=None):

  if options is None:
    options = {}

  if isinstance(options,dict):
      options = RacOptimizerOptions(**options)


  result_R = Racmacs.checkHemisphering(racmap._acmap_R, optimization_number+1,
                                     grid_spacing, stress_lim,
                                     options.options_R)


  with ro.conversion.localconverter(ro.default_converter + ro.pandas2ri.converter):
    result = ro.conversion.get_conversion().rpy2py(result_R)


  ag_diagnostics=\
    list(result["optimizations"].items())[optimization_number][1]["ag_diagnostics"]

  ag_names = racmap.ag_names
  sr_names = racmap.sr_names

  ag_diagnostics = {ag:odict_to_dict(val[1]) for ag,val in zip(ag_names, ag_diagnostics.items())}

  sr_diagnostics=\
    list(result["optimizations"].items())[optimization_number][1]["sr_diagnostics"]

  sr_diagnostics = {sr:odict_to_dict(val[1]) for sr,val in zip(sr_names, sr_diagnostics.items())}

  return sr_diagnostics, ag_diagnostics
