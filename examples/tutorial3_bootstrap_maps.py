#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 23:37:19 2023

@author: avicenna
"""

import PyRacmacs as pr

ndims = 2
racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")

#subset the map to make the process faster for demonstration purposes
racmap = racmap.subset_map(antigen_indices=range(30), serum_indices=range(10))
racmap = pr.optimize_map(racmap, ndims)

#step 1: bootstrap
bootstrapped_map =\
  pr.bootstrap_map(racmap, method="bayesian", bootstrap_repeats=100,
                   optimizations_per_repeat=10)

#step 2: calculate blobs from bootstrap data
map_with_blobs = pr.bootstrap_blobs(bootstrapped_map, smoothing=6, grid_spacing=0.25) #coarse for speed

#step 3 :view the blobs
pr.view(map_with_blobs)

#you can directly extract bootstrap and blob data via the added functionalities
bootstrap = bootstrapped_map.get_bootstrap_data(0)
ag_blobs, sr_blobs = map_with_blobs.get_blobs()

#note that racmacs write function saves the bootstrap data but not the blob data
#so you can save and load bootstrapped_map but if you save blobs, any blob data
#will be lost and cant be reloaded.
pr.write_racmap(bootstrapped_map, "./tutorial_data/h3map2004_wbstrap.ace")
test_map = pr.read_racmap("./tutorial_data/h3map2004_wbstrap.ace")

try:
  bootstrap = test_map.get_bootstrap_data(0)
  print("bootstrap loaded succesfully\n")
except:
  print("bootstrap load failed\n")


pr.write_racmap(map_with_blobs, "./tutorial_data/h3map2004_wblobs.ace")
test_map = pr.read_racmap("./tutorial_data/h3map2004_wblobs.ace")


try:
  blobs = test_map.get_blobs()
  print("blobs loaded succesfully\n")
except:
  print("blobs load failed as expected\n")


#there are some functionalities that does not exist in Racmacs as well.
#after extracting the individual blob data (which is also available in Racmacs)
#one can calculate the blob areas/volumes and maximum sizes (also available in 3D
#and actually more useful in 3D since blobs are very messy to visualize in 3D)

ag_blobs, sr_blobs = map_with_blobs.get_blobs()

ag_blob_area_data, ag_blob_maxsize_data = pr.analyse_blobs(ag_blobs)
sr_blob_area_data, sr_blob_maxsize_data = pr.analyse_blobs(sr_blobs)

fig,ax = pr.plot_lib.bootstrap_bar(ag_blob_area_data, ndims) #bar plot of ag radii
fig,ax = pr.plot_lib.bootstrap_bar(sr_blob_area_data, ndims) #bar plot of sr radii

#you can also do a barplot of maxsizes but make sure to change the transformation
#as the default one calculates radius from volumes.
fig,ax = pr.plot_lib.bootstrap_bar(ag_blob_maxsize_data, ndims,
                                   transformation= lambda x:x/2)

colored_racmap = pr.plot_lib.bootstrap_agcolored_map(ag_blob_area_data, racmap) #map colored by antigen areas
pr.view(colored_racmap)

#one can also do the triangulation blobs and get the blob information
#exactly as in bootstrap_blobs and used them in the same plotting functionalities

racmap = pr.relax_map(racmap)
map_with_tblobs = pr.triangulation_blobs(racmap)
pr.view(map_with_tblobs)
ag_tblobs, sr_tblobs = map_with_tblobs.get_triangulation_blobs()

#some other related functionalities
