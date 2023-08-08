#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:59:45 2022

@author: Sina Tureli
"""

def check(case, *args):
    
      if case==0:
            parameter_update_check(*args)
      if case==1:
            serum_group_existence_check(*args)
        

def error(case, *args):
    
     if case==0:
            parameter_update_error(*args)
     if case==1:
            serum_group_existence_error(*args)
        
def parameter_update_check(parameters, new_parameters):
    
    if not all(key in parameters for key in new_parameters):
        error(0, parameters, new_parameters)
    
def parameter_update_error(parameters, new_parameters):
    
    error_msg = (
        f'new parameters {set(new_parameters).difference(parameters)} are invalid'
        )
    
    raise ValueError(error_msg)
    

def serum_group_existence_error(racmap, sr_group_index):
    
    error_msg = (
        f'No such serum group with index {sr_group_index}. '
        f'Existing serum groups are {set(racmap.sr_groups)}'
        )
    
    raise ValueError(error_msg)
    

def serum_group_existence_check(racmap, sr_group_index):
    
    if sr_group_index not in set(racmap.sr_groups):
    
        error(1, racmap, sr_group_index)