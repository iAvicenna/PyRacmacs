#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:14:45 2023

@author: avicenna
"""

import matplotlib.pyplot as plt
import PyRacmacs as pr
import numpy as np
from scipy.spatial.distance import cdist

racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")
racmap_3D = pr.read_racmap("./tutorial_data/h3map2004_3D.ace")

#transformations of maps
rotated_racmap = pr.rotate_map(racmap, axis="z", optimization_number=0, degrees=10)
translated_racmap = pr.translate_map(racmap, translation=[1,1], optimization_number=0)
reflected_racmap = pr.reflect_map(racmap)

#comparing maps
pr.procrustes_maps(racmap, racmap, 0, 5)
pdata = pr.procrustes_data(racmap, racmap, 0, 5)


#piecewise procrustes
result = pr.piecewise_procrustes(racmap, racmap_3D, 3, num_cpus=6,
                                 metric="sd")

fig,ax = plt.subplots(1, 1, figsize=(5,5))
ax.plot(np.array(result[-1])/(racmap.shape[0]+racmap.shape[1]))
ax.set_xlabel("Iteration")
ax.set_ylabel("Mean Squared Distance")
ax.grid("on", alpha=0.2)
