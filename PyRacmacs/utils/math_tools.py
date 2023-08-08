#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:03:05 2022

@author: Sina Tureli
"""

import numpy as np

def get_centered_data_std_and_mean(data):
    
    mean = np.round(np.nanmean(data),2)
    std = np.round(np.nanstd(data-mean),2)
    
    return std,mean

    
def get_data_RMSE(ax, datax, datay):
    
    
    non_nan_datax, non_nan_datay = get_non_nan_data(datax, datay)

    rmse =  np.round(np.sqrt(((non_nan_datax- non_nan_datay)**2).mean()),2)

    return rmse
    
def get_data_correlation(ax, datax, datay):
    
    
    non_nan_datax, non_nan_datay = get_non_nan_data(datax, datay)
    
    correlation = np.round(np.corrcoef(non_nan_datax, non_nan_datay)[0,1],2)
    
    return correlation
    
    
def get_non_nan_data(datax,datay):
    
    non_nan_indices = [ind for ind,(x,y) in enumerate(zip(datax.flatten(),datay.flatten())) 
                       if not np.isnan(x) and not np.isnan(y)]
    
    non_nan_datax = datax[non_nan_indices]
    non_nan_datay = datay[non_nan_indices]
    
    return non_nan_datax, non_nan_datay
