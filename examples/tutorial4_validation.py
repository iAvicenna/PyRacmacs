#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 14:23:02 2023

@author: avicenna
"""

import PyRacmacs as pr

ndims = 2
racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")

#subset the map to make the process faster for demonstration purposes
racmap = racmap.subset_map(antigen_indices=range(30), serum_indices=range(10))
racmap = pr.optimize_map(racmap, ndims)

# method1: check hemisphering
ag_diagnostics, sr_diagnostics = pr.check_hemisphering(racmap)


# method2: dimensionality testing
test_summary = pr.dimension_test(racmap, dimensions=[1,2,3], number_of_optimizations=100)
