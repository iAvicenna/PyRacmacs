#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 02:28:47 2023

@author: avicenna
"""

import numpy as np
import itertools as it
import trimesh
import time
from scipy.spatial.distance import cdist
from . import Racmacs
from scipy.sparse.csgraph import connected_components
from .optimization import RacOptimizerOptions
from .io import capture_r_output
from . import RacMap




def triangulation_blobs(racmap, optimization_number=0, stress_lim=1, grid_spacing=0.25,
                        antigens=True, sera=True, check_relaxation=True, options=None,
                        verbose=True):

  if options is None:
    options = {}

  if isinstance(options,dict):
      options = RacOptimizerOptions(**options)

  stdout, stderr = capture_r_output(verbose, False)

  if verbose:
    print("running triangulation blobs...", end="")

  t0 = time.time()
  map_with_tblobs_R = Racmacs.triangulationBlobs(racmap._acmap_R, optimization_number+1,
                                                 stress_lim, grid_spacing,
                                                 antigens, sera, check_relaxation,
                                                 options.options_R
                                                 )
  t1 = time.time()
  if verbose:
    print(f"triangulation blobs complete. Took {t1-t0:.2f} secs.\n")


  map_with_tblobs = RacMap(map_with_tblobs_R)

  return map_with_tblobs



def bootstrap_map(racmap, method, bootstrap_repeats=1000, bootstrap_ags=True,
                  bootstrap_sr=True, reoptimize=True, optimizations_per_repeat=100,
                  ag_noise_sd=0.7, titer_noise_sd=0.7,
                  options:RacOptimizerOptions=None,
                  verbose=True):


  if options is None:
      options = {}

  if isinstance(options,dict):
      options = RacOptimizerOptions(**options)

  stdout, stderr = capture_r_output(verbose, False)

  t0 = time.time()
  bootstrapped_map_R = Racmacs.bootstrapMap(racmap._acmap_R,
                                            method,
                                            bootstrap_repeats = bootstrap_repeats,
                                            bootstrap_ags = bootstrap_ags,
                                            bootstrap_sr = bootstrap_sr,
                                            reoptimize = reoptimize,
                                            optimizations_per_repeat = optimizations_per_repeat,
                                            ag_noise_sd = ag_noise_sd,
                                            titer_noise_sd = titer_noise_sd,
                                            options = options.options_R
                                            )
  t1 = time.time()
  if verbose:
    print(f" Bootstrap complete.\nTook {t1-t0:.2f} secs.\n")

  bootstrapped_map = RacMap(bootstrapped_map_R)

  return bootstrapped_map


def bootstrap_blobs(racmap, conf=0.68, smoothing=6, grid_spacing=0.25, antigens=True,
                    sera=True, method="ks", verbose=True):

    '''
    supplied map must be bootstrapped first via the bootstrap_map function
    '''

    stdout, stderr = capture_r_output(verbose, False)

    t0 = time.time()
    racmap_wblobs_R = Racmacs.bootstrapBlobs(racmap._acmap_R, conf, smoothing,
                                             grid_spacing, antigens, sera, method)

    racmap_wblobs = RacMap(racmap_wblobs_R)
    t1 = time.time()

    if verbose:
      print(f" blobs complete.\nTook {t1-t0:.2f} secs.\n")


    return racmap_wblobs




def analyse_blobs(blobs):

  if all([isinstance(blob, dict) and all(x in blob for x in ['vertices', 'faces', 'normals'])
          for blob in blobs.values()]):
    return _analyse_blobs_3D(blobs)
  elif all(isinstance(blob,list) for blob in blobs.values()):
    return _analyse_blobs_2D(blobs)
  else:
    raise ValueError("blobs should be a dict mapping antigen names to dictionaries "
                     "vertices,faces,normals in 3D and to a list of vertices in 2D.")



def _analyse_blobs_3D(blobs):

  blob_vol_data = {}
  blob_maxsize_data = {}

  for name in blobs:
    vertices = blobs[name]['vertices']
    faces = blobs[name]['faces']
    normals = blobs[name]['normals']

    I = list(sorted(set(faces.flatten())))
    blob_vol_data[name] = []
    blob_maxsize_data[name] = []

    vertices = vertices[I,:]
    normals = normals[I,:]
    faces=_reindex_faces(faces)

    if len(I)<4:
      blob_vol_data[name]=[0]
      blob_maxsize_data[name] = np.max(cdist(vertices,vertices))

      continue


    adj = _faces_to_adj(faces)

    blob_comps = connected_components(adj)
    component_levels = set(blob_comps[1])

    for level in component_levels:
      I = [ind for ind,cl in enumerate(blob_comps[1]) if cl==level]

      new_index = {i:ind for ind,i in enumerate(I)}

      subvertices = vertices[I,:]
      subnormals = normals[I,:]
      subfaces = [[new_index[x] for x in face] for face in faces if face[0] in I]

      max_dist = np.max(cdist(subvertices, subvertices).flatten())
      blob_maxsize_data[name].append(np.round(max_dist,3))

      if len(I)<4:
        blob_vol_data[name].append(0)
        continue

      mesh = trimesh.Trimesh(vertices=subvertices, faces=subfaces,
                             facenormal=subnormals)
      mesh.fix_normals()
      mesh.remove_degenerate_faces()
      mesh.remove_duplicate_faces()
      mesh.remove_infinite_values()
      mesh.remove_unreferenced_vertices()

      vol = mesh.volume

      if not(mesh.is_watertight and mesh.is_winding_consistent and\
                  np.isfinite(mesh.center_mass).all() and vol>0):
        vol = mesh.convex_hull.volume

      blob_vol_data[name].append(np.round(vol,3))


    I = _clean_blobs(blob_vol_data[name])
    blob_vol_data[name] = [blob_vol_data[name][i] for i in I]
    blob_maxsize_data[name] = [blob_maxsize_data[name][i] for i in I]

  return blob_vol_data, blob_maxsize_data



def _analyse_blobs_2D(blobs):

  blob_area_data = {key:[] for key in blobs}
  blob_maxsize_data = {key:[] for key in blobs}

  for name in blobs:
    blob_list = blobs[name]

    for vertices in blob_list:

      max_dist = np.max(cdist(vertices, vertices).flatten())
      blob_maxsize_data[name].append(np.round(max_dist,3))

      blob_area_data[name].append(_polygon_area(vertices[:,0],
                                               vertices[:,1]))


  return blob_area_data, blob_maxsize_data

def _polygon_area(x,y):
  '''
  shoelace formula
  https://en.wikipedia.org/wiki/Shoelace_formula
  '''


  correction = x[-1] * y[0] - y[-1]* x[0]
  main_area = np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:])
  return 0.5*np.abs(main_area + correction)

def _clean_blobs(volumes):

  if len(volumes)>1 and max(volumes)>1:
    I = [ind for ind,vol in enumerate(volumes) if vol>0.1]
  else:
    I = [ind for ind,vol in enumerate(volumes) if vol!=0]

  return I


def _reindex_faces(faces):

  I = sorted(set(faces.flatten()))
  index_map = {i:ind for ind,i in enumerate(I)}

  faces = [[index_map[x] for x in face] for face in faces]

  return np.array(faces)


def _faces_to_adj(faces):

  face_indices = sorted(list(set(faces.flatten())))
  nfaces = len(face_indices)
  adj = np.zeros((nfaces,nfaces))

  for face in faces:
    for i0,i1 in it.combinations(face, 2):
      adj[i0,i1]=1

  return adj


# def extract_triangulation_blobs(racmap):


#   Takasbank ® Altın Transfer Sistemi
