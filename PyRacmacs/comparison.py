#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 12:02:55 2022

@author: Sina Tureli
"""


import PyRacmacs as pr
import tempfile
import numpy as np
import tqdm
import warnings

from scipy.spatial.distance import cdist
from multiprocessing import Pool
from functools import partial
from . import Racmacs
from rpy2 import robjects as ro
rNULL = ro.rinterface.NULL

def apply_plotspec(
        target_map: pr.RacMap,
        source_map: pr.RacMap,
        ):

    target_map_R =\
        Racmacs.applyPlotspec(target_map._acmap_R, source_map._acmap_R)

    return pr.RacMap(target_map_R)


def realign(
            source_map: pr.RacMap,
            target_map: pr.RacMap,
            translation=True,
            scaling=False,
            sera=True,
            reflection=True
            ):

    '''
    If gives unexpected results, try setting map transformations to
    identity
    '''

    # for some reason Racmacs realignMap does not always work correctly
    # so realign is done within this package it might have to do with
    # how transformations are handled?
    # realigned_map_R = Racmacs.realignMap(map = source_map._acmap_R,
    #                                      target_map = target_map._acmap_R,
    #                                      translation = translation,
    #                                      scaling = scaling)
    # return pr.RacMap(realigned_map_R)


    if reflection:
      reflection="best"

    S = [] # source coordinates
    T = [] # target coordinates
    if len(source_map.ag_ids)!=0 and not all(x=='' for x in source_map.ag_ids)\
      and len(target_map.ag_ids)!=0 and not all(x=='' for x in target_map.ag_ids):
      source_ag = source_map.ag_ids
      target_ag = target_map.ag_ids
    else:
      source_ag = source_map.ag_names
      target_ag = target_map.ag_names

    common_ag = set(source_ag).intersection(target_ag)
    source_ag_coordinates = source_map.ag_coordinates
    target_ag_coordinates = target_map.ag_coordinates

    for ag in common_ag:
      i0 = source_ag.index(ag)
      i1 = target_ag.index(ag)
      S.append(source_ag_coordinates[i0,:])
      T.append(target_ag_coordinates[i1,:])

    if sera:
      if len(source_map.sr_ids)!=0 and len(target_map.sr_ids)!=0:
        source_sr = source_map.sr_ids
        target_sr = target_map.sr_ids
      else:
        source_sr = source_map.sr_names
        target_sr = target_map.sr_names

      common_sr = set(source_sr).intersection(target_sr)
      source_sr_coordinates = source_map.sr_coordinates
      target_sr_coordinates = target_map.sr_coordinates

      for sr in common_sr:
        i0 = source_sr.index(sr)
        i1 = target_sr.index(sr)
        S.append(source_sr_coordinates[i0,:])
        T.append(target_sr_coordinates[i1,:])


    I = ~np.isnan(np.array(S).sum(axis=1)) & ~np.isnan(np.array(T).sum(axis=1))
    Z,_,_,M = _procrustes(np.array(S)[I,:], np.array(T)[I,:], scaling=scaling,
                          reflection=reflection)
    source_map.coordinates = (M["scale"]*M["rotation"]@source_map.coordinates.T).T

    if translation:
      translation = np.nanmean(target_map.coordinates, axis=0) -\
        np.nanmean(source_map.coordinates, axis=0)

      source_map.coordinates = source_map.coordinates + translation

    return source_map


def procrustes_maps(
        source_map: pr.RacMap,
        target_map: pr.RacMap,
        source_optimization_number: int=0,
        target_optimization_number: int=0,
        compare_antigens = True,
        compare_sera = True,
        allow_translation = True,
        allow_scaling = False,
        export_path = None,
        display = True,
        options = None
        ):

    procrustes_R = Racmacs.procrustesMap(
        map = source_map._acmap_R,
        comparison_map = target_map._acmap_R,
        optimization_number = source_optimization_number+1,
        comparison_optimization_number = target_optimization_number+1,
        antigens = compare_antigens,
        sera = compare_sera,
        translation = allow_translation,
        scaling = allow_scaling
        )

    procrustes = pr.RacMap(procrustes_R)

    if export_path is None:
        with tempfile.NamedTemporaryFile(suffix='.html') as temp_fp:
            export_path = temp_fp.name

    pr.view(procrustes, export_path=export_path, display=display, options=options)


def procrustes_data(
        source_map: pr.RacMap,
        target_map: pr.RacMap,
        source_optimization_number: int=0,
        target_optimization_number: int=0,
        compare_antigens = True,
        compare_sera = True,
        allow_translation = True,
        allow_scaling = False,
        export_path = None
        ):


    procrustes_data_R = Racmacs.procrustesData(
        map = source_map._acmap_R,
        comparison_map = target_map._acmap_R,
        optimization_number = source_optimization_number+1,
        comparison_optimization_number = target_optimization_number+1,
        antigens = compare_antigens,
        sera = compare_sera,
        translation = allow_translation,
        scaling = allow_scaling
        )

    procrustes_data = {
        'ag_dists':list(procrustes_data_R[0]),
        'sr_dists':list(procrustes_data_R[1]),
        'ag_rmsd':list(procrustes_data_R[2])[0],
        'sr_rmsd':list(procrustes_data_R[3])[0],
        'total_rmsd':list(procrustes_data_R[4])[0]
        }

    return procrustes_data


def piecewise_procrustes_data(source_map: pr.RacMap, target_map: pr.RacMap,
                              compare_sera=True, compare_antigens=True,
                              labels=None):

    if compare_antigens and not compare_sera:
      X = source_map.ag_coordinates
      Y = target_map.ag_coordinates
    elif compare_sera and not compare_antigens:
      X = source_map.sr_coordinates
      Y = target_map.sr_coordinates
    elif compare_sera and compare_antigens:
      X = source_map.coordinates
      Y = target_map.coordinates
    else:
      raise ValueError("One of compare_antigens or compare_sera must be True.")

    s = X.shape[0]

    if labels is None:
      labels = [0 for _ in range(s)]

    label_levels = set(labels)

    max_ndims = max(X.shape[1], Y.shape[1])

    pc_data = np.zeros((X.shape[0], max_ndims))

    for label_level in label_levels:
        I = [ind for ind,label in enumerate(labels) if label==label_level]
        X_transform,_,_,_ = pr.comparison._procrustes(Y[I,:], X[I,:])
        Y_transform,_,_,_ = pr.comparison._procrustes(X[I,:], Y[I,:])

        pc_data[I,:] = X_transform - Y[I,:]



    return pc_data


def piecewise_procrustes(source_map: pr.RacMap, target_map: pr.RacMap,
                         npieces: int=2,
                         source_optimization_number: int=0,
                         target_optimization_number: int=0,
                         compare_antigens = True,
                         compare_sera = True,
                         scaling=False, reflection='best',
                         niterations=100000, threshold=2000, num_cpus=1,
                         ninitial_conditions=10, seeds=None,
                         hide_progress_bar=False, metric='sd'):

      '''
      metric can be sd, ad, rad for square distance, absolute distance, root
      absolute distance.

      works by assigning labels (npieces=number of labels) to antigens and
      then changing the partition iteratively so as to find the best partition
      of the antigens which would minimize the procrustes with the given number
      of pieces

      returns:
        - initial labels assigned to each cluster
        - final labels assigned to each cluster
        - initial stress
        - a history of how stress changes with each iteration.

      '''

      if seeds is None:
          seeds = np.random.SeedSequence().spawn(ninitial_conditions)
      else:
          ninitial_conditions = len(seeds)

      if compare_antigens and compare_sera:
        X = source_map.coordinates
        Y = target_map.coordinates
      elif compare_antigens:
        X = source_map.ag_coordinates
        Y = target_map.ag_coordinates
      elif compare_sera:
        X = source_map.sr_coordinates
        Y = target_map.sr_coordinates
      else:
        raise ValueError("one of compare_antigens or compare_sera must be "
                         "True.")

      partial_parfun = partial(_piecewise_procrustes, X, Y, npieces, scaling,
                        reflection, niterations, threshold, metric)

      if num_cpus>1:
        with Pool(num_cpus) as p:
            results = list(tqdm.tqdm(p.imap(partial_parfun, seeds),
                                     total=len(seeds), disable=hide_progress_bar))
      else:
          results = list([partial_parfun(seed) for seed in tqdm.tqdm(seeds)])


      I = np.argmin([result[-1][-1] for result in results])

      return results[I]


def _piecewise_procrustes(X, Y, npieces=2, scaling=False, reflection='best',
                         niterations=10000, threshold=2000, metric='ad', seed=None):

    if seed is None:
        seed = np.random.SeedSequence().spawn(1)[0]

    if X.shape[0] != Y.shape[0]:
        raise ValueError(f'X number of rows is {X.shape[0]}, Y number of rows is {Y.shape[0]} but '
                       'they should be equal.')

    if X.shape[0]<npieces*3:
        raise ValueError('Number of elements in the vector should be more than 3*npieces.')

    if X.shape[0]<npieces:
        raise ValueError('There are more clusters than points.')

    if X.shape[1] != Y.shape[1]:
        ncolsX = X.shape[1]
        ncolsY = Y.shape[1]

        if ncolsX>ncolsY:
            Y = np.concatenate([Y, np.zeros((Y.shape[0], ncolsX - ncolsY))],
                               axis=1)
        else:
            X = np.concatenate([X, np.zeros((X.shape[0], ncolsY - ncolsX))],
                               axis=1)

    npoints = X.shape[0]
    initial_dissimilarity = _total_dissimilarity(X, Y,
                                                 [0 for _ in range(npoints)], 1,
                                                 metric)
    rng = np.random.default_rng(seed)
    current_labels = []
    set_counter = 0
    dissimilarity_history = [initial_dissimilarity]

    while len(set(current_labels))<npieces:
        cluster_centers = X[rng.integers(0, npoints, npieces),:]

        distances = cdist(cluster_centers, X)
        initial_labels = np.argmin(distances,axis=0)

        current_labels = initial_labels.copy()

        if len(set(current_labels))==npieces:
            dissimilarity_history += [_total_dissimilarity(X, Y, current_labels, npieces,
                                                           metric)]

        set_counter += 1

        if set_counter>10*threshold:
          raise ValueError(f'Can not divide into {npieces} nonempty clusters')


    if npieces==1:
      return (initial_labels, current_labels, initial_dissimilarity,
              dissimilarity_history)

    for i in range(niterations):

        if (i+1)%threshold==0 and dissimilarity_history[-1] == dissimilarity_history[-threshold]:
            break

        perm_labels = rng.permutation(range(npieces))

        I = [ind for ind,label in enumerate(current_labels)
             if label==perm_labels[0]]

        if len(I)<3:
          dissimilarity_history.append(dissimilarity_history[-1])
          continue

        pt_ind = I[rng.integers(0,len(I))]

        new_labels = current_labels.copy()

        new_labels[pt_ind] = perm_labels[1]

        new_dissimilarity = _total_dissimilarity(X, Y, new_labels, npieces, metric)
        if new_dissimilarity < dissimilarity_history[-1]:

            dissimilarity_history.append(new_dissimilarity)
            current_labels = new_labels.copy()

        else:
            dissimilarity_history.append(dissimilarity_history[-1])


    return (initial_labels, current_labels, initial_dissimilarity,
            dissimilarity_history)


def _total_dissimilarity(X, Y, labels, n_groups, metric):

    total_dissimilarity = 0

    for i in range(n_groups):

        I = [ind for ind,label in enumerate(labels) if label==i]

        Xsub = X[I,:]
        Ysub = Y[I,:]

        try:
          Ytrans,_,_,_=\
            _procrustes(Xsub, Ysub, scaling=False)
        except Exception as e:
          Ytrans = Ysub
          print(e)

        if metric == 'sd':
            total_dissimilarity +=\
              np.sum(np.linalg.norm((Ytrans-Xsub),axis=1)**2)
        if metric == 'ad':
            total_dissimilarity +=\
              np.sum(np.abs(Ytrans-Xsub))
        if metric == 'rad':
            total_dissimilarity +=\
              np.sum(np.abs(Ytrans-Xsub)**(1/2))

    return total_dissimilarity


def _procrustes(X, Y, scaling=True, reflection='best'):
    """
    A port of MATLAB's `procrustes` function to Numpy.
    @author: https://stackoverflow.com/questions/18925181/procrustes-analysis-with-numpy

    X: target
    Y: transformed

    Procrustes analysis determines a linear transformation (translation,
    reflection, orthogonal rotation and scaling) of the points in Y to best
    conform them to the points in matrix X, using the sum of squared errors
    as the goodness of fit criterion.

        d, Z, [tform] = procrustes(X, Y)

    Inputs:
    ------------
    X, Y
        matrices of target and input coordinates. they must have equal
        numbers of  points (rows)

    scaling
        if False, the scaling component of the transformation is forced
        to 1

    reflection
        if 'best' (default), the transformation solution may or may not
        include a reflection component, depending on which fits the data
        best. setting reflection to True or False forces a solution with
        reflection or no reflection respectively.

    Outputs
    ------------
    d
        the residual sum of squared errors, normalized according to a
        measure of the scale of X, ((X - X.mean(0))**2).sum()

    Z
        the matrix of transformed Y-values

    tform
        a dict specifying the rotation, translation and scaling that
        maps X --> Y

    """

    sumX = np.sum(X,axis=1)
    sumY = np.sum(Y,axis=1)
    I = [ind for ind in range(X.shape[0]) if not np.isnan(sumX[ind])
         and not np.isnan(sumY[ind])]

    if len(I)==0:
        warnings.warn("Atleast one of the coordinates have all nan values. "
                      "Aborting procrustes")

        return Y, np.inf, np.sqrt(np.nanmean(np.linalg.norm((X - Y)**2, axis=1))),\
          {"scale":1, "rotation":np.eye(Y.shape[1],),
           "translation":np.zeros(Y.shape[1],)}

    if len(I)/X.shape[0]<0.5:
        print('Warning: More than half of the rows have a nan value in it')

    X_sub = X[I,:]
    Y_sub = Y[I,:]

    nx,mx = X_sub.shape
    ny,my = Y_sub.shape

    if my < mx:
        Y_sub = np.concatenate((Y_sub, np.zeros((ny, mx-my))), 1)
        my = mx
    elif my > mx:
        X_sub = np.concatenate((X_sub, np.zeros((nx, my-mx))), 1)
        mx = my


    muX = X_sub.mean(0)
    muY = Y_sub.mean(0)

    X0 = X_sub - muX
    Y0 = Y_sub - muY

    ssX = (X0**2.).sum()
    ssY = (Y0**2.).sum()

    # centred Frobenius norm
    normX = np.sqrt(ssX)
    normY = np.sqrt(ssY)

    if normX==0 or normY==0:
        warnings.warn("One of the vectors has only one element. "
                      "Aborting procrustes")
        #return Y, np.inf, np.inf, None
        return Y_sub, np.inf, np.sqrt(np.nanmean(np.linalg.norm((X - Y)**2, axis=1))),\
          {"scale":1, "rotation":np.eye(my,), "translation":np.zeros(my,)}


    # scale to equal (unit) norm
    X0 /= normX
    Y0 /= normY


    # optimum rotation matrix of Y
    A = np.dot(X0.T, Y0)

    U,s,Vt = np.linalg.svd(A,full_matrices=False)
    V = Vt.T
    T = np.dot(V, U.T)

    if reflection != 'best':

        # does the current solution use a reflection?
        have_reflection = np.linalg.det(T) < 0

        # if that's not what was specified, force another reflection
        if reflection != have_reflection:
            V[:,-1] *= -1

            s[-1] *= -1
            T = np.dot(V, U.T)

    traceTA = s.sum()

    if scaling:

        # optimum scaling of Y
        b = traceTA * normX / normY

        # standarised distance between X and b*Y*T + c
        standardized_distance = 1 - traceTA**2

        # transformed coords
        Z = normX*traceTA*np.dot(Y0, T) + muX

    else:
        b = 1
        standardized_distance = 1 + ssY/ssX - 2 * traceTA * normY / normX
        Z = normY*np.dot(Y0, T) + muX

    # transformation matrix
    if my < mx:
        T = T[:my,:]
    c = muX - b*np.dot(muY, T)

    #transformation values
    tform = {'rotation':T, 'scale':b, 'translation':c}
    rmsd = np.sqrt(np.nanmean(np.linalg.norm((Z - X_sub)**2, axis=1)))

    Y_transformed = np.zeros((Y.shape[0],Y_sub.shape[1]))

    for i in set(range(ny)).difference(I):
        Y_transformed[i,:] = Y_sub[i,:]
    for indi,i in enumerate(I):
        Y_transformed[i,:] = Z[indi,:]

    return Y_transformed, standardized_distance, rmsd, tform


'''
Transformations of maps affect the results in the following way
(correspondance with Sam Wilks):

They affect the coordinates you get when you call agCoords, etc.
so yes are accounted for when plotting maps.  For procrustes, this aligns
the maps by adjusting the transformations and the base coordinates aren't
changed.  When you do relaxmap the transformation is ignored and only
the base coordinates are optimised.
'''

def rotate_map(racmap:pr.RacMap, degrees:float, optimization_number:int=None,
               axis:str=None):

  if optimization_number is None:
    optimization_number = rNULL
  else:
    optimization_number += 1

  if axis is None:
    axis = rNULL

  rotated_map_R = Racmacs.rotateMap(racmap._acmap_R, degrees, axis,
                                    optimization_number)


  rotated_map = pr.RacMap(rotated_map_R)

  return rotated_map


def translate_map(racmap:pr.RacMap, translation, optimization_number:int=None):

  if optimization_number is None:
    optimization_number = rNULL
  else:
    optimization_number += 1


  translation_R = ro.FloatVector(translation)

  translated_map_R = pr.Racmacs.translateMap(racmap._acmap_R, translation_R,
                                             optimization_number)

  translated_map = pr.RacMap(translated_map_R)

  return translated_map


def reflect_map(racmap:pr.RacMap, axis=None, optimization_number:int=None):

  if optimization_number is None:
    optimization_number = rNULL
  else:
    optimization_number += 1

  if axis is None: axis="x"

  reflected_map_R = pr.Racmacs.reflectMap(racmap._acmap_R, axis,
                                          optimization_number)

  reflected_map = pr.RacMap(reflected_map_R)

  return reflected_map
