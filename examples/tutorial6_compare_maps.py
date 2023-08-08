#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:14:45 2023

@author: avicenna
"""

import PyRacmacs as pr


racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")

#transformations of maps
rotated_racmap = pr.rotate_map(racmap, axis="z", optimization_number=0, degrees=10)
translated_racmap = pr.translate_map(racmap, translation=[1,1], optimization_number=0)
reflected_racmap = pr.reflect_map(racmap)

#comparing maps
pr.procrustes_maps(racmap, racmap, 0, 5)
pdata = pr.procrustes_data(racmap, racmap, 0, 5)
