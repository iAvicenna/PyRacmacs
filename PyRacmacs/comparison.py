#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 12:02:55 2022

@author: Sina Tureli
"""


import PyRacmacs as pr
import tempfile

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
            map: pr.RacMap,
            target_map: pr.RacMap,
            translation=True,
            scaling=False
            ):

    realigned_map_R = Racmacs.realignMap(map = map._acmap_R,
                                         target_map = target_map._acmap_R,
                                         translation = translation,
                                         scaling = scaling
                                         )

    return pr.RacMap(realigned_map_R)


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
        'sr_dsts':list(procrustes_data_R[1]),
        'ag_rmsd':list(procrustes_data_R[2]),
        'sr_rmsd':list(procrustes_data_R[3]),
        'total_rmsd':list(procrustes_data_R[4])
        }

    return procrustes_data


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
