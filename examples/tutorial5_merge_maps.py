#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:01:49 2023

@author: avicenna
"""
import PyRacmacs as pr
import pandas as pd
#some preparation

titer_table1 = pd.read_csv("./tutorial_data/h3map2004_table1.csv", header=0,
                           index_col=0)
titer_table2 = pd.read_csv("./tutorial_data/h3map2004_table2.csv", header=0,
                           index_col=0)

map1 = pr.read_racmap("./tutorial_data/h3map2004_map1.ace")
map2 = pr.read_racmap("./tutorial_data/h3map2004_map2.ace")


#there are several methods to merge maps all identical to how it is implemented
#in Racmacs apart from using Python like structures.

#method1 merge as table optimize from scratch
merged_titer_table = pr.merge_tables([titer_table1, titer_table2])
merged_map1 = pr.make_map_from_table(merged_titer_table, dilution_stepsize=1)

#equivalent to above
merged_map2 = pr.merge_maps([map1, map2], 2, 1, method="reoptimized-merge")

#other methods which can be used are frozen-overlay, relaxed-overlay,
#frozen-merge, see https://acorg.github.io/Racmacs/reference/mergeMaps.html
