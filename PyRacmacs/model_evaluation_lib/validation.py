#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:01:54 2022

@author: Sina Tureli
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random
import sys

from rpy2 import robjects as ro

from .local_utils import DimensionTestSummaryStatistics, CrossValidationResultsTable
from . import local_utils
from .. import (Racmacs, optimization, RacMap)
from .. import io
from ..utils import regulations, plotting, math_tools, conversion

check = regulations.check
OptimizationParameters = optimization.OptimizationParameters

rNULL = ro.rinterface.NULL


def dimensional_cross_validation(
        racmap: RacMap,
        dimensions: list,
        test_proportion: float=0.1,
        validations_per_dimension: int=100,
        seed: int=None,
        optimization_parameters = None,
        ):

    if seed is None:
        seed = int(str(time.time()).replace('.',''))

    if optimization_parameters is None:
        optimization_parameters = OptimizationParameters()
    else:
        default_optimization_parameters = OptimizationParameters()
        default_optimization_parameters.update(optimization_parameters)
        optimization_parameters = default_optimization_parameters

    optimization_parameters['display_progress'] = False

    tests_per_validation = int(racmap.number_of_well_coordinated_titers*test_proportion)

    summary_statistics = DimensionTestSummaryStatistics(dimensions, validations_per_dimension)
    dimension_to_results_table = CrossValidationResultsTable(dimensions, validations_per_dimension,
                                                tests_per_validation)


    progress_bar = io.PyProgressBar(total_iterations=len(dimensions)*validations_per_dimension)

    print((f'Performing dimensional cross-validation for dimensions {dimensions} '
           f'and {validations_per_dimension} validations_per_dimension')
          )

    random_seed_generator = np.random.default_rng(seed)

    for dimension in dimensions:

        optimization_parameters['number_of_dimensions'] = dimension

        dimension_seed = random_seed_generator.integers(0,sys.maxsize)

        validation_to_titers, validation_to_test_indices = cross_validate(racmap,
                                                                          validations_per_dimension,
                                                                          test_proportion, dimension_seed,
                                                                          optimization_parameters,
                                                                          progress_bar)

        dimension_to_results_table.update_dimension_result(dimension, validation_to_titers,
                                              validation_to_test_indices)

        summary_statistics = local_utils.summary_statistics_for_result_tables(dimension_to_results_table,
                                                                  validations_per_dimension)


    return dimension_to_results_table, summary_statistics



def cross_validate(racmap:RacMap, number_of_validations, test_proportion,
                   seed, optimization_parameters, progress_bar=None):

    random_generator = random.Random()
    random_generator.seed(int(seed))

    validation_to_titers={}
    validation_to_test_indices ={}

    for validation_number in range(number_of_validations):

        validation_to_titers[validation_number] = {}

        masked_titer_table, masked_indices = local_utils.generate_masked_titer_table(racmap,
                                                                         test_proportion,
                                                                         random_generator)

        masked_map = optimization.make_map_from_dataframe(masked_titer_table,
                                                          racmap.dilution_stepsize,
                                                          **optimization_parameters)

        predicted_titers = masked_map.predicted_titers()[masked_indices]
        expected_titers = racmap.titer_table.values[masked_indices]

        validation_to_titers[validation_number]['predicted'] = predicted_titers
        validation_to_titers[validation_number]['expected'] = expected_titers

        validation_to_test_indices[validation_number] = racmap.ravel_multi_index(masked_indices)

        progress_bar.grow()

    return validation_to_titers, validation_to_test_indices





def run_dimension_test_map_R(
        racmap: RacMap,
        dimensions: list=None,
        test_proportion: float=0.1,
        minimum_column_basis:  str="none",
        fixed_column_bases: list=None,
        number_of_optimizations: int=1000,
        validations_per_dimension: int=100,
        options: list=None,
        display_progress = True
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


    stoud, stderr = io.capture_r_output(display_progress)

    t0 = time.time()

    dimension_test_results_R = Racmacs.runDimensionTestMap(
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

    dimension_test_results = conversion.recurse_r_tree(dimension_test_results_R)

    dimension_to_results_table = convert_dimension_test_results_to_tables(
        racmap,
        dimension_test_results
        )

    summary_statistics = local_utils.summary_statistics_for_result_tables(dimension_to_results_table,
                                                                          validations_per_dimension)

    return dimension_to_results_table, summary_statistics


def convert_dimension_test_results_to_tables(racmap:RacMap, dimension_test_results):

    validations_per_dimension = len(dimension_test_results['results'])
    tests_per_validation = dimension_test_results['results'][0]['test_indices'].size

    dimensions = list(dimension_test_results['results'][0]['coords'].keys())

    dimension_to_result_tables = CrossValidationResultsTable(dimensions,
                                                             validations_per_dimension,
                                                             tests_per_validation)

    titers = dimension_test_results['titers'] # this the titer table flattened in col major from R

    log_titers = [np.log2(conversion.titer_str_to_num(x)/10) for x in titers]

    for ind_run,run_result in enumerate(dimension_test_results['results'].values()):

        row_range = range(ind_run*tests_per_validation,
                          (ind_run+1)*tests_per_validation)

        test_indices = [int(x)-1 for x in run_result['test_indices'].flatten()]  #  R indices start from 1
        test_titers = [titers[x] for x in test_indices ]
        test_log_titers = [log_titers[x] for x in test_indices ]

        test_multi_indices = racmap.unravel_index(test_indices,order='col')   # R indices are col major so
        test_indices = racmap.ravel_multi_index(test_multi_indices,order='row') # convert them to row major

        for dim_ind,dimension in enumerate(dimensions):

            predicted_log_titers = run_result['predictions'][dim_ind]
            residuals = local_utils.calculate_residuals(test_titers, predicted_log_titers)
            detectability = [not('<' in x or '>' in x) for x in test_titers]

            dimension_to_result_tables[dimension].iloc[row_range,0] = [ind_run for _ in row_range]
            dimension_to_result_tables[dimension].iloc[row_range,1] = test_indices
            dimension_to_result_tables[dimension].iloc[row_range,2] = predicted_log_titers
            dimension_to_result_tables[dimension].iloc[row_range,3] = test_titers
            dimension_to_result_tables[dimension].iloc[row_range,4] = test_log_titers
            dimension_to_result_tables[dimension].iloc[row_range,5] = residuals
            dimension_to_result_tables[dimension].iloc[row_range,6] = detectability

    return dimension_to_result_tables


def plot_dimension_test_result(dimension_test_result: pd.DataFrame,
                               include_nondetectable=True,
                               include_variation=True,
                               ax=None):

    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=(5, 5))
    else:
        fig = ax.get_figure()

    mean_rmse_detectable = dimension_test_result.loc[:,'mean_rmse_detectable'].values
    var_rmse_detectable = dimension_test_result.loc[:,'var_rmse_detectable'].values

    mean_rmse_nondetectable = dimension_test_result.loc[:,'mean_rmse_nondetectable'].values
    var_rmse_nondetectable = dimension_test_result.loc[:,'var_rmse_nondetectable'].values

    legends = []
    xvalues = range(1,mean_rmse_detectable.size+1)

    if (not any(not np.isnan(x) and not np.isinf(x) for x in mean_rmse_detectable) and
        not any(not np.isnan(x) and not np.isinf(x) for x in mean_rmse_nondetectable)):

        raise ValueError('All values in dimension test result are either nan or inf.')

    if any(not np.isnan(x) and not np.isinf(x) for x in mean_rmse_detectable):

        ax.plot(xvalues, mean_rmse_detectable, color='black')

        if include_variation:
            ax.fill_between(xvalues, (mean_rmse_detectable-var_rmse_detectable),
                            (mean_rmse_detectable+var_rmse_detectable), color='black',
                            alpha=.1, label='_nolegend_')

        legends.append('mean_rmse_detectable')

    if (any(not np.isnan(x) and not np.isinf(x) for x in mean_rmse_nondetectable)
        and include_nondetectable):

        ax.plot(xvalues, mean_rmse_nondetectable, color='red')

        if include_variation:
            ax.fill_between(xvalues, (mean_rmse_nondetectable-var_rmse_nondetectable),
                            (mean_rmse_nondetectable+var_rmse_nondetectable), color='red',
                            alpha=.1, label='_nolegend_')

        legends.append('mean_rmse_nondetectable')

    ax.set_xlabel('Dimension')
    ax.set_xticks(xvalues)

    ax.set_ylabel('mean rmse')
    ax.legend(legends)

    plotting.add_grid(ax)

    return fig,ax


def plot_cross_validation_histogram(cross_validation_test_result,
                                    dimensions,
                                    density=False,
                                    fig_options=None,
                                    ax=None):

    if ax is None:
        if fig_options is None:
            fig_options = {}

        fig,ax = plt.subplots(1,1)
    else:
        fig = ax.get_figure()

    if density:
        ylabel = 'Density'
    else:
        ylabel = 'Count'

    fig,ax = plt.subplots(1,len(dimensions), figsize=(5*len(dimensions),5))

    if len(dimensions)==1:
        ax = [ax]

    axes_uniformizer = plotting.AxesUniformizer(equal_axes=False)

    for dim_ind,dimension in enumerate(dimensions):

        residuals = cross_validation_test_result[dimension]['residuals']
        s1= residuals.shape[0]

        detectable_titer_indices = np.argwhere(cross_validation_test_result[dimension]['detectables'].values)
        detectable_titer_indices = list(detectable_titer_indices.flatten())
        undetectable_titer_indices = set(range(s1)).difference(detectable_titer_indices)
        undetectable_titer_indices = list(sorted(undetectable_titer_indices))

        detectable_residuals = [residuals[x] for x in detectable_titer_indices]
        undetectable_residuals = [residuals[x] for x in undetectable_titer_indices]

        legends = []
        ax[dim_ind].set_title(f'Number of dimensions {dimension}')

        if (len(detectable_titer_indices)>0 and
            any(not np.isnan(x) for x in detectable_residuals)):

            ax[dim_ind].hist(detectable_residuals, color='black', edgecolor='black',
                             alpha=0.7, density=density, **fig_options)

            std, mean = math_tools.get_centered_data_std_and_mean(detectable_residuals)

            legends.append('Detectable Titers\n'
                           rf'$\mu$:{mean}, $\sigma$:{std}'
                           )

        if (len(undetectable_residuals)>0 and
            any(not np.isnan(x) for x in undetectable_residuals)):

            ax[dim_ind].hist(undetectable_residuals, color='red', edgecolor='black',
                             alpha=0.7, density=density, **fig_options)

            td, mean = math_tools.get_centered_data_std_and_mean(detectable_residuals)

            legends.append('Detectable Titers\n'
                           rf'$\mu$:{mean}, $\sigma$:{std}'
                           )

        ax[dim_ind].set_ylabel(ylabel)
        ax[dim_ind].set_xlabel('Measured Titer - Predicted Titer')
        ax[dim_ind].legend(legends)
        axes_uniformizer.add_axis(ax[dim_ind])

    axes_uniformizer.uniformize()

    return fig,ax
