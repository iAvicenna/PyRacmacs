#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 18:38:05 2023

@author: avicenna
"""

import numpy as np
import time
import itertools as it

from . import (Racmacs, RacMap, RacOptimizerOptions, make_map_from_table,
               procrustes_data)
from .model_evaluation_lib.local_utils import DimensionTestSummaryStatistics
from rpy2 import robjects as ro
from .io import capture_r_output
from .utils.conversion import odict_to_dict


def leave_pair_out_test(
    table,
    dilution_stepsize,
    number_of_dimensions,
    number_of_optimizations,
    serum_list=None,
    antigen_list=None,
    remove_sera=True,
    remove_antigens = True,
    include_thresholded=False,
    optimization_args: dict=None,
    ):

  if optimization_args is None:
    optimization_args = {}




  fullmap=\
    make_map_from_table(table, dilution_stepsize=dilution_stepsize,
                        number_of_dimensions=number_of_dimensions,
                        number_of_optimizations=number_of_optimizations)

  dif = fullmap.predicted_log_titers() - fullmap.log_titer_table.values

  if include_thresholded:
    I = fullmap.measured_titer_locations
  else:
    I = fullmap.detectable_titer_locations

  serum_to_diagnostics = {None:{"dif_std":np.nanstd(dif[I]), "red_chi2":fullmap.get_red_chi2(1),
                                "total_rmsd":0
                                }}
  antigen_to_diagnostics = {None:{"dif_std":np.nanstd(dif[I]), "red_chi2":fullmap.get_red_chi2(1),
                                  "total_rmsd":0}}
  if remove_sera:

    if serum_list is None:
      serum_list = table.columns

    serum_pairs = it.product(serum_list, serum_list)


    for serum0,serum1 in serum_pairs:

      if serum0<=serum1:
        continue

      cols = [s for s in table.columns if s not in [serum0, serum1]]
      subtable = table.loc[:, cols]

      submap=\
        make_map_from_table(subtable, dilution_stepsize=dilution_stepsize,
                            number_of_dimensions=number_of_dimensions,
                            number_of_optimizations=number_of_optimizations)

      dif = submap.predicted_log_titers() - submap.log_titer_table.values

      if include_thresholded:
        I = submap.measured_titer_locations
      else:
        I = submap.detectable_titer_locations


      proc_data = procrustes_data(submap, fullmap)

      serum_to_diagnostics[(serum0, serum1)] = {}

      serum_to_diagnostics[(serum0, serum1)]["dif_std"] = np.nanstd(dif[I])
      serum_to_diagnostics[(serum0, serum1)]["red_chi2"] = submap.get_red_chi2(1)
      serum_to_diagnostics[(serum0, serum1)]["total_rmsd"] = proc_data["total_rmsd"]



  if remove_antigens:

    if antigen_list is None:
      antigen_list = table.index

    antigen_pairs = it.product(antigen_list, antigen_list)

    for antigen0,antigen1 in antigen_pairs:

      index = [a for a in table.index if a not in [antigen0, antigen1]]
      subtable = table.loc[index, :]

      submap=\
        make_map_from_table(subtable, dilution_stepsize=dilution_stepsize,
                            number_of_dimensions=number_of_dimensions,
                            number_of_optimizations=number_of_optimizations)


      dif = submap.predicted_log_titers() - submap.log_titer_table.values
      if include_thresholded:
        I = submap.measured_titer_locations
      else:
        I = submap.detectable_titer_locations

      antigen_to_diagnostics[(antigen0, antigen1)] = {}

      proc_data = procrustes_data(submap, fullmap)

      antigen_to_diagnostics[(antigen0, antigen1)]["dif_std"] = np.nanstd(dif[I])
      antigen_to_diagnostics[(antigen0, antigen1)]["red_chi2"] = submap.get_red_chi2(1)
      antigen_to_diagnostics[(antigen0, antigen1)]["total_rmsd"] = proc_data["total_rmsd"]

  diagnostics = {"serum":serum_to_diagnostics,
                 "antigen":antigen_to_diagnostics}

  return diagnostics




def leave_one_out_test(
    table,
    dilution_stepsize,
    number_of_dimensions,
    number_of_optimizations,
    serum_list=None,
    antigen_list=None,
    remove_sera=True,
    remove_antigens = True,
    include_thresholded=False,
    optimization_args: dict=None,
    ):

  if optimization_args is None:
    optimization_args = {}


  fullmap=\
    make_map_from_table(table, dilution_stepsize=dilution_stepsize,
                        number_of_dimensions=number_of_dimensions,
                        number_of_optimizations=number_of_optimizations)

  dif = fullmap.predicted_log_titers() - fullmap.log_titer_table.values

  if include_thresholded:
    I = fullmap.measured_titer_locations
  else:
    I = fullmap.detectable_titer_locations

  serum_to_diagnostics = {None:{"dif_std":np.nanstd(dif[I]), "red_chi2":fullmap.get_red_chi2(1),
                                "total_rmsd":0
                                }}
  antigen_to_diagnostics = {None:{"dif_std":np.nanstd(dif[I]), "red_chi2":fullmap.get_red_chi2(1),
                                  "total_rmsd":0}}

  if remove_sera:

    if serum_list is None:
      serum_list = table.columns

    for serum in serum_list:

      cols = [s for s in table.columns if s != serum]
      subtable = table.loc[:, cols]

      submap=\
        make_map_from_table(subtable, dilution_stepsize=dilution_stepsize,
                            number_of_dimensions=number_of_dimensions,
                            number_of_optimizations=number_of_optimizations)

      dif = submap.predicted_log_titers() - submap.log_titer_table.values

      if include_thresholded:
        I = submap.measured_titer_locations
      else:
        I = submap.detectable_titer_locations

      serum_to_diagnostics[serum] = {}

      proc_data = procrustes_data(submap, fullmap)


      serum_to_diagnostics[serum]["dif_std"] = np.nanstd(dif[I])
      serum_to_diagnostics[serum]["red_chi2"] = submap.get_red_chi2(1)
      serum_to_diagnostics[serum]["total_rmsd"] = proc_data["total_rmsd"]



  if remove_antigens:

    if antigen_list is None:
      antigen_list = table.index

    for antigen in antigen_list:

      index = [a for a in table.index if a != antigen]
      subtable = table.loc[index, :]

      submap=\
        make_map_from_table(subtable, dilution_stepsize=dilution_stepsize,
                            number_of_dimensions=number_of_dimensions,
                            number_of_optimizations=number_of_optimizations)


      dif = submap.predicted_log_titers() - submap.log_titer_table.values
      if include_thresholded:
        I = submap.measured_titer_locations
      else:
        I = submap.detectable_titer_locations

      antigen_to_diagnostics[antigen] = {}

      proc_data = procrustes_data(submap, fullmap)

      antigen_to_diagnostics[antigen]["dif_std"] = np.nanstd(dif[I])
      antigen_to_diagnostics[antigen]["red_chi2"] = submap.get_red_chi2(1)
      antigen_to_diagnostics[antigen]["total_rmsd"] = proc_data["total_rmsd"]


  diagnostics = {"serum":serum_to_diagnostics,
                 "antigen":antigen_to_diagnostics}

  return diagnostics



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
