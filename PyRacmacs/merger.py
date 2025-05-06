#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:43:06 2022

@author: avicenna
"""

import time
from . import RacOptimizerOptions, RacMergeOptions, Racmacs, RacMap
from .io import capture_r_output


def merge_tables(tables:list, number_of_dimensions:int, dilution_step_size:int,
                 sd_limit:float=None, merge_args:dict=None):


    if len(tables)==1:
        return tables[0]

    if merge_args is None:
        merge_args = {}

    maps = [RacMap(titer_table=table) for table in tables]

    return merge_maps(maps, number_of_dimensions, dilution_step_size, sd_limit,
                      **merge_args).titer_table


def merge_maps(maps:list, number_of_dimensions:int, dilution_step_size:int,
               sd_limit:float=None, number_of_optimizations:int=1000, method="table",
               minimum_column_basis:str="none", optimizer_options=None,
               merge_options=None, verbose=True, merge_report=False):

    '''
    method should be table, reoptimized-merge, frozen-overlay, relaxed-overlay,
    frozen-merge
    '''

    if merge_options is None:
        merge_options = {}
    if optimizer_options is None:
        optimizer_options = {}

    if isinstance(optimizer_options,dict):
        optimizer_options = RacOptimizerOptions(**optimizer_options)
        optimizer_options.display_progress = verbose

    if isinstance(merge_options,dict):
        merge_options = RacMergeOptions(**merge_options)
        merge_options.dilution_step_size = dilution_step_size
        if sd_limit is not None:
          merge_options.sd_limit = sd_limit


    stdout, stderr = capture_r_output(verbose, True)

    t0 = time.time()

    merged_map_R =\
    Racmacs.mergeMaps(*[acmap._acmap_R for acmap in maps],
                      method=method, number_of_dimensions=number_of_dimensions,
                      number_of_optimizations=number_of_optimizations,
                      minimum_column_basis=minimum_column_basis,
                      optimizer_options = optimizer_options.options_R,
                      merge_options = merge_options.options_R,
                      verbose=verbose
                      )

    t1 = time.time()

    if verbose:
      print(f"Took {t1-t0:.2f} secs.\n")

    if merge_report:
        report = Racmacs.mergeReport(merged_map_R)

        return RacMap(merged_map_R), report

    else:

        return RacMap(merged_map_R)
