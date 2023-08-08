#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 21:16:11 2022

@author: Sina Tureli
"""

import numpy as np
import pandas as pd
from ..utils.conversion import titer_str_to_num


class DimensionTestSummaryStatistics(pd.DataFrame):

    def __init__(self, dimensions, number_of_validations, data=None):

        index = [int(x) for x in dimensions]

        columns = ['mean_rmse_detectable',
                  'var_rmse_detectable',
                  'mean_rmse_nondetectable',
                  'var_rmse_nondetectable',
                  'number_of_validations'
                  ]

        shape = (len(index),len(columns))

        if data is None:
            data = np.zeros(shape)*np.nan

        super().__init__(data=data, index=index, columns=columns)

        self.index.name = 'dimension'

        self.iloc[:,4] = [int(number_of_validations) for _ in range(len(dimensions))]

    def update_detectable_result(self, detectable_titers_RMSEs, dimension):

        if not all(np.isnan(x) for x in  detectable_titers_RMSEs):
            self.loc[dimension,'mean_rmse_detectable'] = np.nanmean(detectable_titers_RMSEs)
            self.loc[dimension,'var_rmse_detectable'] = np.nanvar(detectable_titers_RMSEs)

    def update_nondetectable_result(self, nondetectable_titers_RMSEs, dimension):

        if not all(np.isnan(x) for x in  nondetectable_titers_RMSEs):
            self.loc[dimension,'mean_rmse_nondetectable'] = np.nanmean(nondetectable_titers_RMSEs)
            self.loc[dimension,'var_rmse_nondetectable'] = np.nanvar(nondetectable_titers_RMSEs)




class CrossValidationResultsTable(dict):

    def __init__(self, dimensions, validations_per_dimension, tests_per_validation):

        super().__init__()

        columns = ['run','test_indices','predicted_log_titers','expected_titers',
                   'expected_log_titers','residuals','detectables']

        N = validations_per_dimension*tests_per_validation

        for dimension in dimensions:

            self[dimension] = pd.DataFrame(np.empty((N,len(columns))),
                                            columns=columns)

    def update_dimension_result(self, dimension, validation_to_titers,
                            validation_to_test_indices):

        for validation_no in validation_to_titers:

            expected_titers = validation_to_titers[validation_no]['expected']
            expected_log_titers = [np.log2(titer_str_to_num(x)/10) for x in expected_titers]
            predicted_titers = validation_to_titers[validation_no]['predicted']
            predicted_log_titers = [np.log2(titer_str_to_num(x)/10) for x in predicted_titers]
            residuals = calculate_residuals(expected_titers, predicted_log_titers)

            detectability = ['<' not in x and '>' not in x for x in expected_titers]
            test_indices = [int(x) for x in validation_to_test_indices[validation_no]]

            N = len(test_indices)
            row_range = range(validation_no*N,(validation_no+1)*N)
            self[dimension].iloc[row_range,0] = validation_no
            self[dimension].iloc[row_range,1] = test_indices
            self[dimension].iloc[row_range,2] = predicted_log_titers
            self[dimension].iloc[row_range,3] = expected_titers
            self[dimension].iloc[row_range,4] = expected_log_titers
            self[dimension].iloc[row_range,5] = residuals
            self[dimension].iloc[row_range,6] = detectability


def generate_masked_titer_table(racmap, test_proportion, random_generator):


    titer_rows,titer_cols = racmap.well_coordinated_titer_locations
    n = int(racmap.number_of_well_coordinated_titers*test_proportion)

    indices = list(range(len(titer_rows)))

    random_generator.shuffle(indices)

    masked_indices = ([titer_rows[x] for x in indices[0:n]],
                      [titer_cols[x] for x in indices[0:n]])



    titer_table = racmap.titer_table.copy()

    titers = titer_table.values

    titers[masked_indices] = '*'

    masked_titer_table = pd.DataFrame(np.reshape(titers, titer_table.shape),
                                        index = titer_table.index,
                                        columns = titer_table.columns)

    return masked_titer_table, masked_indices


def summary_statistics_for_result_tables(dimension_to_result_tables, replicates_per_dimension):

    dimensions = list(dimension_to_result_tables.keys())
    summary_statistics = DimensionTestSummaryStatistics(dimensions,
                                                        replicates_per_dimension)


    for dimension in dimension_to_result_tables:

        result_table = dimension_to_result_tables[dimension]

        runs = sorted(list(set(result_table.iloc[:,0])))

        run_detectable_RMSEs = []
        run_nondetectable_RMSEs = []

        for run in runs:
            detectable_titer_indices = [ind for ind,x in
                                        enumerate(result_table.iloc[:,0].values.flatten())
                                        if x==run and result_table.iloc[ind,6]]

            nondetectable_titer_indices = [ind for ind,x in
                                           enumerate(result_table.iloc[:,0].values.flatten())
                                           if x==run and not result_table.iloc[ind,6]]

            if (len(detectable_titer_indices)>0 and not
                all(np.isnan(x) for x in result_table.loc[detectable_titer_indices,
                                                          'residuals'].values)):

                detectable_RMSE = np.sqrt(np.nanmean((result_table.loc[detectable_titer_indices,
                                                                       'residuals'].values)**2))

            else:
                detectable_RMSE = np.nan

            if (len(nondetectable_titer_indices)>0 and not
                all(np.isnan(x) for x in result_table.loc[nondetectable_titer_indices,
                                                          'residuals'].values)):

                nondetectable_RMSE = np.sqrt(np.nanmean((result_table.loc[nondetectable_titer_indices,
                                                                          'residuals'].values)**2))
            else:
                nondetectable_RMSE = np.nan

            run_detectable_RMSEs.append(detectable_RMSE)
            run_nondetectable_RMSEs.append(nondetectable_RMSE)

        summary_statistics.update_detectable_result(run_detectable_RMSEs, dimension)
        summary_statistics.update_nondetectable_result(run_nondetectable_RMSEs, dimension)

    return summary_statistics



def calculate_RMSEs_for_each_validation(validation_to_titers):

    titers_RMSEs = {'detectable':[], 'nondetectable':[]}

    for validation in validation_to_titers:

        predicted_titers = validation_to_titers[validation]['predicted']
        expected_titers = validation_to_titers[validation]['expected']

        detectable_titers_RMSE = calculate_detected_rmse(expected_titers, predicted_titers)
        nondetectable_titers_RMSE = calculate_undetected_rmse(expected_titers, predicted_titers)

        titers_RMSEs['detectable'].append(detectable_titers_RMSE)
        titers_RMSEs['nondetectable'].append(nondetectable_titers_RMSE)

    return titers_RMSEs


def calculate_detected_rmse(expected_titers, predicted_titers):

    detected_indices = [ind for ind,x in enumerate(expected_titers) if '<' not in x and '>' not in x]


    detected_expected_log_titers = np.array([np.log2(titer_str_to_num(x)/10) for x in
                           expected_titers[detected_indices]])

    detected_predicted_log_titers = np.array([np.log2(titer_str_to_num(x)/10) for x in
                            predicted_titers[detected_indices]])



    detected_se = (detected_expected_log_titers - detected_predicted_log_titers) ** 2

    if not all(np.isnan(x) for x in detected_se):
        return np.sqrt(np.nanmean(detected_se))
    else:
        return np.nan

def calculate_residuals(expected_titers, predicted_log_titers):

    residuals = []

    for titer_str1,log_titer2 in zip(expected_titers, predicted_log_titers):

        titer_num1 = titer_str_to_num(titer_str1)
        log_titer1 = np.log2(titer_num1/10)

        if '<' in titer_str1 and log_titer2<=log_titer1:
            residuals.append(0)
        elif '>' in titer_str1 and log_titer2>=log_titer1:
            residuals.append(0)
        else:
            residuals.append(log_titer1-log_titer2)

    return residuals

def calculate_undetected_rmse(expected_titers, predicted_titers):

    '''
    lod: limit of detection
    llod: lower lod
    ulod: upper lod
    '''

    llod_indices = [ind for ind,x in enumerate(expected_titers) if '<' in x]
    ulod_indices = [ind for ind,x in enumerate(expected_titers) if '>' in x]
    llod_se = np.array([])
    ulod_se = np.array([])


    if len(llod_indices)>0:
        llod_expected_log_titers = [np.log2(titer_str_to_num(x)/10) for x in
                               expected_titers[llod_indices]]

        llod_predicted_log_titers = [np.log2(titer_str_to_num(x)/10) for x in
                                predicted_titers[llod_indices]]

        llod_se = np.array(
            [0 if x<y else y-x for x,y in
             zip(llod_expected_log_titers,llod_predicted_log_titers)]
            )**2


    if len(ulod_indices)>0:

        ulod_expected_log_titers = np.array([np.log2(titer_str_to_num(x)/10) for x in
                               expected_titers[ulod_indices]])

        ulod_predicted_log_titers = np.array([np.log2(titer_str_to_num(x)/10) for x in
                                predicted_titers[ulod_indices]])

        ulod_se = np.array(
            [0 if x>y else y-x for x,y in
             zip(ulod_expected_log_titers,ulod_predicted_log_titers)]
           )**2


    undetected_se = np.append(llod_se,ulod_se)

    if undetected_se.size == 0:
        return np.nan
    else:
        if not all(np.isnan(x) for x in undetected_se):
            return np.sqrt(np.nanmean(undetected_se))
        else:
            return np.nan

def get_detectable_and_undetectable_log_titers(racmap):

    expected_log_titers = racmap.log_titer_table.values
    predicted_log_titers = racmap.predicted_log_titers

    detectable_titer_locations = racmap.detectable_titer_locations
    undetectable_titer_locations = racmap.undetectable_titer_locations

    detectable_expected_log_titers = expected_log_titers[detectable_titer_locations].flatten()
    detectable_predicted_log_titers = predicted_log_titers()[detectable_titer_locations].flatten()
    undetectable_expected_log_titers = expected_log_titers[undetectable_titer_locations].flatten()
    undetectable_predicted_log_titers = predicted_log_titers()[undetectable_titer_locations].flatten()


    return (detectable_expected_log_titers,
            detectable_predicted_log_titers,
            undetectable_expected_log_titers,
            undetectable_predicted_log_titers
            )


def get_detectable_and_undetectable_titers(racmap):

    expected_titers = racmap.titer_table.values
    predicted_titers = racmap.predicted_titers

    detectable_titer_locations = racmap.detectable_titer_locations
    undetectable_titer_locations = racmap.undetectable_titer_locations

    detectable_expected_titers = expected_titers[detectable_titer_locations].flatten()
    detectable_predicted_titers = predicted_titers()[detectable_titer_locations].flatten()
    undetectable_expected_titers = expected_titers[undetectable_titer_locations].flatten()
    undetectable_predicted_titers = predicted_titers()[undetectable_titer_locations].flatten()


    return (detectable_expected_titers,
            detectable_predicted_titers,
            undetectable_expected_titers,
            undetectable_predicted_titers
            )
