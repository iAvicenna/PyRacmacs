#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 22:28:39 2023

@author: avicenna
"""
import PyRacmacs as pr
import pandas as pd

# there are two ways to make maps

# 1- make_map_from_table: makes a map from a dataframe. dataframe is the python
#    equivalent of a table. row names become antigen names and column names
#    becomes serum names unless these are seperately specified as inputs.

titer_table = pd.read_csv("./tutorial_data/h3_table.csv", index_col=0, header=0)
racmap=\
pr.make_map_from_table(titer_table, dilution_stepsize=1, number_of_dimensions=2,
                       number_of_optimizations=100)

# 2- if you have already loaded a map and want to reoptimize it from scratch:
racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")
racmap = pr.optimize_map(racmap, number_of_optimizations=100)

# 3- if you have already loaded a map and want to relax it around its minimum
racmap = pr.relax_map(racmap)


#notes: RacOptimizer.options() can be supplied as the options argument
#which can be supplied as a dictionary using the same key words as in
#RacOptimizer.options() for instance
options = pr.RacOptimizerOptions(dim_annealing=True, num_cores=12)
racmap = pr.optimize_map(racmap, number_of_optimizations=100,
                          options=options)



#maps can be viewed with pr.view. similiar to optimization options,
#view options can be supplied as such:
options = pr.RacViewerOptions(opacity=1, show_errorlines=True,
                              translation=[0.2,0,0])
pr.view(racmap, options=options)

#you can also save maps and turn of automatic display
pr.view(racmap, export_path="./tutorial_data/h3map2004.html",
        display=False)

#pr.plot_lib provides other plotting utilities for instance interactive
#titer plots from the map data etc.
