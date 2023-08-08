#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:06:35 2022

@author: avicenna
"""

import unittest
import PyRacmacs as pr
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

'''
These tests simply check that one call call the relevant
Racmacs functionalities without problems. These functionalities are
already tested for consistency within the Racmacs package so not
testing them here further
'''

def create_titer_table(dim, nAg, nSr):

    antigen_coords = np.random.rand(nAg,dim)*10
    serum_coords = np.random.rand(nSr, dim)*10

    distances = cdist(antigen_coords, serum_coords)
    log_titers = np.max(distances) - distances
    titers = 10*2**log_titers

    antigens = [f'Ag{i}' for i in range(nAg)]
    sera = [f'Sr{i}' for i in range(nSr)]

    return pd.DataFrame(titers, index=antigens, columns=sera)


def repeat_titration(titer_table, mean, std):

    repeat_titers = titer_table.values + np.random.normal(0,0.3,
                                                          titer_table.shape)

    return pd.DataFrame(repeat_titers, index=titer_table.index,
                        columns=titer_table.columns)


class MakeMapTests(unittest.TestCase):


    def test_make_map(self):


        for dim in [2,3,4]:

            nAg = np.random.randint(100,200)
            nSr = np.random.randint(50,80)

            titer_table = create_titer_table(dim, nAg, nSr)

            pr.make_map_from_table(titer_table, dilution_stepsize=1,
                                   verbose=False,
                                   number_of_dimensions=dim,
                                   number_of_optimizations=10)


    def test_make_map_with_options(self):


        for dim in [2,3,4]:

            nAg = np.random.randint(100,200)
            nSr = np.random.randint(50,80)

            titer_table = create_titer_table(dim, nAg, nSr)

            options1 = pr.RacOptimizerOptions(dim_annealing=True)
            options2 = pr.RacOptimizerOptions(maxit=1000)
            options3 = {'dim_annealing':True}

            pr.make_map_from_table(titer_table, dilution_stepsize=1,
                                   verbose=False,
                                   number_of_dimensions=dim,
                                   number_of_optimizations=10,
                                   options=options1)

            pr.make_map_from_table(titer_table, dilution_stepsize=1,
                                   verbose=False,
                                   number_of_dimensions=dim,
                                   number_of_optimizations=10,
                                   options=options2)

            pr.make_map_from_table(titer_table, dilution_stepsize=1,
                                   verbose=False,
                                   number_of_dimensions=dim,
                                   number_of_optimizations=10,
                                   options=options3)


    def test_optimize_map(self):


        for dim in [2,3,4]:

            nAg = np.random.randint(100,200)
            nSr = np.random.randint(50,80)

            titer_table = create_titer_table(dim, nAg, nSr)

            racmap = pr.RacMap(titer_table=titer_table)

            options1 = pr.RacOptimizerOptions(dim_annealing=True)
            options2 = pr.RacOptimizerOptions(maxit=1000)
            options3 = {'dim_annealing':True}

            pr.optimize_map(racmap, 2, number_of_optimizations=10,
                            verbose=False)

            pr.optimize_map(racmap, 2, number_of_optimizations=10,
                            options=options1, verbose=False)

            pr.optimize_map(racmap, 2, number_of_optimizations=10,
                            options=options2, verbose=False)

            pr.optimize_map(racmap, 2, number_of_optimizations=10,
                            options=options3, verbose=False)



class MergeMapsTests(unittest.TestCase):


    def test_merge_maps(self):

        for method in ["table", "reoptimized-merge", "incremental-merge",
                       "frozen-overlay", "relaxed-overlay", "frozen-merge"]:

            dim = 2
            nAg = np.random.randint(100,200)
            nSr = np.random.randint(50,80)

            titer_table1 = create_titer_table(dim, nAg, nSr)
            titer_table2 = repeat_titration(titer_table1, 0, 0.2)

            map1 = pr.RacMap(titer_table=titer_table1)
            map2 = pr.RacMap(titer_table=titer_table2)

            map1 = pr.optimize_map(map1, number_of_optimizations=10,
                                   verbose=False)
            map2 = pr.optimize_map(map2, number_of_optimizations=10,
                                   verbose=False)

            pr.merge_maps([map1, map2], 2, 1, 10, method=method,
                          verbose=False)

            merge_options = pr.RacMergeOptions(sd_limit=2)
            optimizer_options = pr.RacOptimizerOptions(dim_annealing=True,
                                                          maxit=1000)

            pr.merge_maps([map1, map2], 2, 1, 10, method=method,
                          merge_options=merge_options,
                          optimizer_options=optimizer_options,
                          verbose=False)


class TestValidation(unittest.TestCase):

  def test_dimension_test(self):

    nAg = 10
    nSr = 10

    titer_table = create_titer_table(2, nAg, nSr)

    racmap=\
    pr.make_map_from_table(titer_table, dilution_stepsize=1,
                           verbose=False,
                           number_of_dimensions=2,
                           number_of_optimizations=10)

    summary = pr.dimension_test(racmap, [1,2,3,4],
                                number_of_optimizations=10,
                                validations_per_dimension=100,
                                verbose=False
                                )



  def test_check_hemisphering(self):

    nAg = 10
    nSr = 10

    titer_table = create_titer_table(2, nAg, nSr)

    racmap=\
    pr.make_map_from_table(titer_table, dilution_stepsize=1,
                           verbose=False,
                           number_of_dimensions=2,
                           number_of_optimizations=10)

    data = pr.check_hemisphering(racmap)


class TestMapComparison(unittest.TestCase):

  def test_procrustes(self):

    nAg=10
    nSr=10

    titer_table = create_titer_table(2, nAg, nSr)

    racmap=\
    pr.make_map_from_table(titer_table, dilution_stepsize=1,
                           verbose=False,
                           number_of_dimensions=2,
                           number_of_optimizations=10)

    #comparing maps
    pr.procrustes_maps(racmap, racmap, 0, 5)
    pdata = pr.procrustes_data(racmap, racmap, 0, 5)

  def test_transformations(self):

    nAg=10
    nSr=10

    titer_table = create_titer_table(2, nAg, nSr)

    racmap=\
    pr.make_map_from_table(titer_table, dilution_stepsize=1,
                           verbose=False,
                           number_of_dimensions=2,
                           number_of_optimizations=10)

    #transformations of maps
    rotated_racmap = pr.rotate_map(racmap, axis="z", optimization_number=0, degrees=10)
    translated_racmap = pr.translate_map(racmap, translation=[1,1], optimization_number=0)
    reflected_racmap = pr.reflect_map(racmap)

    assert np.max((rotated_racmap.coordinates-racmap.coordinates).flatten())>1e-4
    assert np.max((translated_racmap.coordinates-racmap.coordinates).flatten())>1e-4
    assert np.max((reflected_racmap.coordinates-racmap.coordinates).flatten())>1e-4


def MakeMapTestSuite():

    suite = unittest.TestSuite()
    suite.addTest(MakeMapTests('test_make_map_with_options'))
    suite.addTest(MakeMapTests('test_make_map'))
    suite.addTest(MakeMapTests('test_optimize_map'))

    return suite


def MergeMapsTestSuite():

    suite = unittest.TestSuite()
    suite.addTest(MergeMapsTests('test_merge_maps'))

    return suite

def ValidationTestSuite():

    suite = unittest.TestSuite()
    suite.addTest(TestValidation('test_dimension_test'))
    suite.addTest(TestValidation('test_check_hemisphering'))

    return suite

def ComparisonTestSuite():

    suite = unittest.TestSuite()
    suite.addTest(TestMapComparison('test_procrustes'))
    suite.addTest(TestMapComparison('test_transformations'))

    return suite

if __name__ == '__main__':

    runner = unittest.TextTestRunner()
    runner.run(MakeMapTestSuite())
    runner.run(MergeMapsTestSuite())
    runner.run(ValidationTestSuite())
    runner.run(ComparisonTestSuite())
