#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 21:40:57 2022

@author: Sina Tureli
"""

import numpy as np
import json
import rpy2
import rpy2.robjects as ro

from rpy2.robjects import r,  conversion, default_converter
from collections import OrderedDict
from rpy2.robjects.vectors import (DataFrame, BoolVector, FloatVector, FloatMatrix,
                                   IntVector, StrVector, ListVector, Matrix,
                                   StrMatrix)
from rpy2.robjects import pandas2ri, numpy2ri
rNULL = ro.rinterface.NULL

from .. import Racmacs


__ALL__ = ['convert_acmap_R_to_py','set_dilution_stepsize']

def titer_types(titers):

    titer_types = np.zeros(titers.shape)

    I1 = titers == '*'
    I2 = np.vectorize(lambda x: '<' in x)(titers)
    I3 = np.vectorize(lambda x: '>' in x)(titers)

    titer_types[I1] = 1
    titer_types[I2] = 2
    titer_types[I3] = 3

    return titer_types


def titer_str_to_num(titer_str):

    if titer_str in ['*','','nan',None]:
        return np.nan
    elif titer_str[0] in ['<','>']:
        return float(titer_str[1:])
    elif '/' in titer_str:
        titers = titer_str.split('/')
        if len(titers)==2:
            titer1,titer2 = titers
            return np.sqrt(float(titer1)*float(titer2))
        else:
            raise ValueError(f'Unknown format for titer {titer_str}')
    elif '?' in titer_str:
        titer_str = titer_str.replace('?','')
        return titer_str_to_num(titer_str)
    elif all(x.isnumeric() or x=='.' for x in titer_str):
        return float(titer_str)
    else:
        raise ValueError(f'Unknown format for titer {titer_str}')


def titers_num_to_str(num_titers, titer_types):

    s1,s2 = num_titers.shape
    str_titers = np.zeros((s1,s2),dtype=object)

    for i in range(s1):
        for j in range(s2):

            if titer_types[i,j] == 0:
                str_titers[i,j] = str(int(num_titers[i,j]))
            elif titer_types[i,j] == 1:
                str_titers[i,j] = '*'
            elif titer_types[i,j] == 2:
                str_titers[i,j] = '<' + str(int(num_titers[i,j]))
            elif titer_types[i,j] == 3:
                str_titers[i,j] = '>' + str(int(num_titers[i,j]))
            else:
                raise ValueError(f'titer_types can only contain 0 to 3 but '
                                 f'entry ({i},{j}) is {titer_types[i,j]}'
                                 )
    return str_titers

def odict_to_dict(odict):

  '''
  should be used with odicts which does not contain None keys or
  other odicts with None keys
  '''

  new_dict = {}

  for key in odict:

    if isinstance(odict[key], rpy2.rlike.container.OrdDict):
      new_dict[key] = odict_to_dict(odict[key])
    elif isinstance(odict[key], list):
      new_dict[key] = [odict_to_dict(val) for val in odict[key]]
    else:
      new_dict[key] = odict[key]

  return new_dict


def convert_options(option_dict):

    option_dict_R = {}

    with ro.conversion.localconverter(ro.default_converter + NONE_to_NA_converter):
        for key in option_dict:
            option_dict_R[key] = ro.conversion.py2rpy(option_dict[key])


    return option_dict_R


def recurse_r_tree(data):
    """
    step through an R object recursively and convert the types to python types as appropriate.
    Leaves will be converted to e.g. numpy arrays or lists as appropriate and the whole tree to a dictionary.

    author: user3324315
    link: https://stackoverflow.com/questions/24152160/converting-an-rpy2-listvector-to-a-python-dictionary

    """
    r_dict_types = [DataFrame, ListVector]
    r_array_types = [FloatVector, IntVector, Matrix, FloatMatrix, StrMatrix]
    r_list_types = [StrVector, BoolVector]

    if type(data) in r_dict_types:
        try:
            return OrderedDict(zip(data.names, [recurse_r_tree(elt) for elt in data]))
        except TypeError:
            return OrderedDict(zip(range(len(data)), [recurse_r_tree(elt) for elt in data]))

    elif type(data) in r_list_types:
        return [recurse_r_tree(elt) for elt in data]
    elif type(data) in r_array_types:
        return np.array(data)
    else:
        if hasattr(data, "rclass"):  # An unsupported r class
            raise KeyError('Could not proceed, type {} is not defined'
                           'to add support for this type, just add it to the imports '
                           'and to the appropriate type list above'.format(type(data)))
        else:
            return data  # We reached the end of recursion


def acmap_R_to_dict(acmap_R):

    acmap_R_json = Racmacs.as_json(acmap_R)
    acmap_dict = json.loads(recurse_r_tree(acmap_R_json)[0])

    return acmap_dict


def set_dilution_stepsize(acmap_R, dilution_stepsize):

    set_method = r("`dilutionStepsize<-`")
    return set_method(acmap_R, dilution_stepsize)


def return_list(val):
    if isinstance(val, rpy2.rinterface_lib.sexp.NULLType):
        return None
    else:
        with conversion.localconverter(default_converter + NA_to_NONE_converter
                                       + pandas2ri.converter + numpy2ri.converter
                                       + str_vec_to_str_converter) as cv:
            return [cv.rpy2py(x) for x in val]


NA_to_NONE_converter = conversion.Converter("NA to NONE")

def NA_to_NONE(na):
    return None

NONE_to_NA_converter = conversion.Converter("NONE to NA")

def NONE_to_NA(none):
    return rNULL


str_vec_to_str_converter = conversion.Converter('str_vector_converter')

def str_vector_to_str(str_vector):
    return ''.join([x for x in str_vector])


NA_to_NONE_converter.rpy2py.register(type(ro.NA_Integer),  NA_to_NONE)
NA_to_NONE_converter.rpy2py.register(type(ro.NA_Real),  NA_to_NONE)
NA_to_NONE_converter.rpy2py.register(type(ro.NA_Character),  NA_to_NONE)
NA_to_NONE_converter.rpy2py.register(type(ro.NULL),  NA_to_NONE)
NA_to_NONE_converter.rpy2py.register(type(ro.NA_Logical),  NA_to_NONE)
str_vec_to_str_converter.rpy2py.register(type(rpy2.robjects.vectors.StrVector('a')),
                                         str_vector_to_str)

NONE_to_NA_converter.rpy2py.register(type(None),  NA_to_NONE)
