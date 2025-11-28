#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 21:40:33 2022

@author: Sina Tureli
"""

import pandas as pd

from rpy2 import robjects as ro

from . import Racmacs, RacMap, RacOptimizerOptions
from .io import capture_r_output
from .utils import regulations, conversion
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
base = importr('base')

rNULL = ro.rinterface.NULL
check = regulations.check


class OptimizationParameters(dict):

    def __init__(self):

        self['ag_names'] = None
        self['sr_names'] = None
        self['number_of_dimensions'] = 2
        self['minimum_column_basis'] = "none"
        self['fixed_column_bases'] = None
        self['number_of_optimizations'] = 1000
        self['titer_weights'] = None
        self['sort_optimizations'] = True
        self['check_convergence'] = True
        self['verbose'] = True
        self['options'] = None

    def update(self, new_dict):

        check(0, self, new_dict)

        for key in new_dict:
            self[key] = new_dict[key]


def make_map_from_table(
    titer_table: pd.DataFrame, dilution_stepsize: int, ag_names: list=None,
    sr_names: list=None, number_of_dimensions: int=2,
    minimum_column_basis: str="none", number_of_optimizations: int=1000,
    fixed_column_bases: list=None, titer_weights: list=None,
    sort_optimizations: bool=True, check_convergence: bool=True,
    verbose: bool=True, options:RacOptimizerOptions=None
    ):

    stdout, stderr = capture_r_output(verbose, True)

    _check_table(titer_table, number_of_dimensions)

    if ag_names is None:
        ag_names = ro.StrVector(list(titer_table.index))
    if sr_names is None:
        sr_names = ro.StrVector(list(titer_table.columns))

    if options is None:
        options = {}


    if fixed_column_bases is None:
        fixed_column_bases = rNULL
    else:
        fixed_column_bases = ro.FloatVector(fixed_column_bases)

    if titer_weights is None:
        titer_weights = rNULL
    else:
        assert titer_weights.shape==titer_table.shape,\
          "titer_table and titer_weights must have the same shape"

        ro.numpy2ri.activate()
        titer_weights = ro.r.matrix(titer_weights, nrow=titer_weights.shape[0],
                                    ncol=titer_weights.shape[1])


    if isinstance(options,dict):
        options = RacOptimizerOptions(**options)

    with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
        table_R = ro.conversion.get_conversion().py2rpy(titer_table)


    acmap_R = Racmacs.acmap(ag_names, sr_names, table_R)
    acmap_R = conversion.set_dilution_stepsize(acmap_R, dilution_stepsize)

    optimized_map_R = Racmacs.optimizeMap(acmap_R, number_of_dimensions=number_of_dimensions,
                                          number_of_optimizations=number_of_optimizations,
                                          minimum_column_basis=minimum_column_basis,
                                          fixed_column_bases=fixed_column_bases,
                                          titer_weights=titer_weights,
                                          sort_optimizations=sort_optimizations,
                                          check_convergence=check_convergence,
                                          verbose=verbose,
                                          options=options.options_R
                                          )

    optimized_racmap = RacMap(optimized_map_R)

    return optimized_racmap


def optimize_map(racmap, number_of_dimensions: int=2,
                 minimum_column_basis: str="none", number_of_optimizations: int=1000,
                 fixed_column_bases: list=None, titer_weights: list=None,
                 sort_optimizations: bool=True, check_convergence: bool=True,
                 verbose: bool=True, options:RacOptimizerOptions=None):

    _check_table(racmap.titer_table, number_of_dimensions)

    if options is None:
        options = {}

    if fixed_column_bases is None:
        fixed_column_bases = rNULL
    else:
        fixed_column_bases = ro.FloatVector(fixed_column_bases)

    if titer_weights is None:
        titer_weights = rNULL
    else:
        assert titer_weights.shape==racmap.shape,\
          "titer_table and titer_weights must have the same shape"

        ro.numpy2ri.activate()
        titer_weights = ro.r.matrix(titer_weights, nrow=titer_weights.shape[0],
                                    ncol=titer_weights.shape[1])

    if isinstance(options,dict):
        options = RacOptimizerOptions(**options)


    stdout, stderr = capture_r_output(verbose, True)

    optimized_map_R = Racmacs.optimizeMap(racmap._acmap_R, number_of_dimensions=number_of_dimensions,
                                          number_of_optimizations=number_of_optimizations,
                                          minimum_column_basis=minimum_column_basis,
                                          fixed_column_bases=fixed_column_bases,
                                          titer_weights=titer_weights,
                                          sort_optimizations=sort_optimizations,
                                          check_convergence=check_convergence,
                                          verbose=verbose,
                                          options=options.options_R
                                          )
    ro.numpy2ri.deactivate()
    optimized_racmap = RacMap(optimized_map_R)

    return optimized_racmap


def relax_map(acmap, optimization_number: int=0, fixed_antigens: bool=False,
              fixed_sera: bool=False, titer_weights: list=None,
              verbose=True,  options:RacOptimizerOptions=None):

    _check_table(acmap.titer_table, acmap.number_of_dimensions)

    if titer_weights is None:
        titer_weights = rNULL
    else:
        ro.numpy2ri.activate()
        titer_weights = ro.r.matrix(titer_weights, nrow=titer_weights.shape[0],
                                    ncol=titer_weights.shape[1])

    if options is None:
        options = {}

    if isinstance(options,dict):
        options = RacOptimizerOptions(**options)


    stdout, stderr = capture_r_output(verbose, True)

    relaxed_map_R = Racmacs.relaxMap(acmap._acmap_R,
                                     optimization_number=optimization_number+1,
                                     fixed_antigens=fixed_antigens,
                                     fixed_sera=fixed_sera,
                                     titer_weights=titer_weights,
                                     options=options.options_R)
    ro.numpy2ri.deactivate()

    relaxed_map = RacMap(relaxed_map_R)

    return relaxed_map


def _check_table(titer_table, dimensions):

    s1,s2 = titer_table.shape

    removed_antigens = ''
    poorly_coordinated_antigens = ''

    for i in range(s1):

        titers = titer_table.values[i,:].flatten().tolist()
        I = [x for x in titers if x!='*' and '<' not in str(x)]

        if len(I) == dimensions:
            poorly_coordinated_antigens += titer_table.index[i]+'\n'
        elif len(I)<dimensions:
            removed_antigens += titer_table.index[i]+'\n'

    if len(poorly_coordinated_antigens)>0:
        print(f'Following antigens have only {dimensions} detectable titrations. '
              f'Their position will be uncertain.')
        print(poorly_coordinated_antigens)
    if len(removed_antigens)>0:
        print(f'Following antigens have less than {dimensions} detectable titrations. '
              f'Their position will be set to NaN.')
        print(removed_antigens)


    removed_sera = ''
    poorly_coordinated_sera = ''
    for i in range(s2):

        titers = titer_table.values[:,i].flatten().tolist()

        I = [x for x in titers if x!='*' and '<' not in str(x)]

        if len(I) == dimensions:
            poorly_coordinated_sera += titer_table.columns[i]+'\n'
        elif len(I)<dimensions:
            removed_sera += titer_table.columns[i]+'\n'

    if len(poorly_coordinated_sera)>0:
        print(f'Following sera have only {dimensions} measured titrations. '
              f'Their position will be uncertain.')
        print(poorly_coordinated_sera)
    if len(removed_sera)>0:
        print(f'Following sera have less than {dimensions} measured titrations. '
              f'Their position will be set to NaN.')
        print(removed_sera)
