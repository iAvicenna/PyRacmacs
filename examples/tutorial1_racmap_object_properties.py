#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:07:14 2023

@author: avicenna
"""

import PyRacmacs as pr
import time
# the main object of the module is RacMap which a python equivalent of
# a map in Racmacs

racmap = pr.read_racmap("./tutorial_data/h3map2004.ace")
#or via class property
racmap = pr.RacMap.from_ace("./tutorial_data/h3map2004.ace")
#relax the map to make sure stress info is stored for later usage in the tutorial
racmap = pr.relax_map(racmap,0)


# racmap._acmap_R is the R object from which everything is extracted when
# requested. since each time a property is accessed it is converted from
# its R equivalent, it is a bit more expensive than directly accessing it
# so when using properties repeatly, you may want to make a copy. as
# an example:

t0 = time.time()

for i in range(500):

  racmap.ag_fills[0]

t1 = time.time()

ag_fills = racmap.ag_fills.copy()

for i in range(500):

  ag_fills[0]

t2 = time.time()

print(f"conversion access: {t1-t0:.4f}, direct access: {t2-t1:.4f}")

#this makes interfacing with Racmacs functions easier however care must
#be taken for repeated access to elements. one of the to dos is to write
#another object PacMap that will directly work with converted properties
#the following is still fine though:

t0 = time.time()
for i in range(10):
  for ag_fill in racmap.ag_fills:
    continue

t1 = time.time()

ag_fills = racmap.ag_fills
for i in range(10):
  for ag_fill in ag_fills:
    continue
t2 = time.time()

print(f"conversion access: {t1-t0:.4f}, direct access: {t2-t1:.4f}")


#many properties that exists in Racmacs as camelCase exists here as
#camel_case, for instance agFills->ag_fills, srOutlines->sr_outlines.
#like in Racmacs you can get and set using these:

print("before change " + racmap.ag_fills[0])
racmap.ag_fills = ["black"] + racmap.ag_fills[1:]
print("after change " + racmap.ag_fills[0])

#note that incorrect usage will not give a warning
racmap.ag_fills[0] = "white"
print("after wrong change still " + racmap.ag_fills[0])

#in Racmacs optimizations are indexed starting from 1, in PyRacmacs from 0

first_optimization = racmap.get_optimization(0)


#there are some functionalities that does not exists in Racmacs for instance
#reduced chi2 statistics:
print(racmap.get_red_chi2(1,0))

#as well as functionalities to get bootstrap data and blobs if they exist
#(they dont in this case, so see tutorial 3)
