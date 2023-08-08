#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 21:19:17 2022

@author: Sina Tureli
"""
import matplotlib.pyplot as plt
import numpy as np

from . import local_utils
from .. import RacMap
from ..utils import plotting

def plot_goodness_of_fit_scatter(
        racmap: RacMap,
        fig_options = None,
        lims = None,
        sr_group = None,
        ax = None
        ):
    
    if ax is None:
        if fig_options is None:
            fig_options = {}
        
        if 'figsize' not in fig_options:
            fig_options['figsize'] = (5,5)
    
        fig,ax = plt.subplots(1,1,**fig_options)
    else:
        fig = ax.get_figure()
    
    if sr_group is not None:
        selected_sr_indices = racmap.sr_indices_of_sr_group(sr_group)
        detectable_titer_color = plotting.convert_color_to_hex(racmap.srOutlines[selected_sr_indices[0]])
        undetectable_titer_color = 'white'
        serum_group_name = racmap.sr_group_levels[sr_group-1]
        title = f'Serum Group {serum_group_name}'
        alpha = 1
    else:
        selected_sr_indices = range(racmap.num_sera)
        detectable_titer_color = 'black'
        undetectable_titer_color = 'red'
        title = 'All Sera'
        alpha = 0.7

    subset_racmap = racmap.subset_map(serum_indices=selected_sr_indices)
    
    all_log_titers = local_utils.get_detectable_and_undetectable_log_titers(subset_racmap)
        
    detectable_expected_log_titers = all_log_titers[0]
    detectable_predicted_log_titers = all_log_titers[1]
    undetectable_expected_log_titers = all_log_titers[2]
    undetectable_predicted_log_titers = all_log_titers[3]
    
    legends = []
    
    if len(detectable_expected_log_titers)>0:
        ax.scatter(detectable_predicted_log_titers, 
                   detectable_expected_log_titers, 
                   c=detectable_titer_color, alpha=alpha, 
                   edgecolor='black',s=30)
       
        
        legends.append('Detectable Titers')

    if len(undetectable_expected_log_titers)>0:
        ax.scatter(undetectable_predicted_log_titers, 
                   undetectable_expected_log_titers, 
                   c=undetectable_titer_color, alpha=alpha,
                   edgecolor='red',s=30)
        
        
        
        legends.append('Detectable Titers')
    
    ax.set_title(title)
    ax.set_ylabel('Measured Titers')
    ax.set_xlabel('Fitted Titers')
    
    if lims is not None:
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    
    plotting.make_plot_square(ax)
    plotting.add_grid(ax)
    plotting.add_log_titer_ticklabels(ax)
    plotting.draw_diagonal(ax)

    return fig,ax


def plot_goodness_of_fit_histogram(
        racmap: RacMap,
        density=False,
        fig_options = None,
        lims = None,
        sr_group = None,
        ax = None
        ):
    
    if ax is None:
        if fig_options is None:
            fig_options = {}
    
        fig,ax = plt.subplots(1,1)
    else:
        fig = ax.get_figure()
        
    if density:
        ylabel = 'Density'
    else:
        ylabel = 'Count'
    
    if sr_group is not None:
        selected_sr_indices = racmap.sr_indices_of_sr_group(sr_group)
        detectable_titer_color = plotting.convert_color_to_hex(racmap.srOutlines[selected_sr_indices[0]])
        undetectable_titer_color = 'white'
        serum_group_name = racmap.sr_group_levels[sr_group-1]
        title = f'Serum Group {serum_group_name}'
        alpha = 1
    else:
        selected_sr_indices = range(racmap.num_sera)
        detectable_titer_color = 'black'
        undetectable_titer_color = 'red'
        title = 'All Sera'
        alpha = 0.7
    
    subset_racmap = racmap.subset_map(serum_indices=selected_sr_indices)
    
    all_titers = local_utils.get_detectable_and_undetectable_titers(subset_racmap)
    all_log_titers = local_utils.get_detectable_and_undetectable_log_titers(subset_racmap)
    
    detectable_expected_titers = all_titers[0]
    detectable_predicted_log_titers = all_log_titers[1]
    undetectable_expected_titers = all_titers[2]
    undetectable_predicted_log_titers = all_log_titers[3]    
    
    
    detectable_residuals = local_utils.calculate_residuals(detectable_expected_titers,
                                                           detectable_predicted_log_titers)
    
    undetectable_residuals = local_utils.calculate_residuals(undetectable_expected_titers,
                                                             undetectable_predicted_log_titers)
    
    legends = []

    if len(detectable_residuals)>0:
        ax.hist(detectable_residuals, color=detectable_titer_color, alpha=alpha, 
                   edgecolor='black', density=density, **fig_options)
       
        mean = np.nanmean(detectable_residuals)
        std = np.nanstd(detectable_residuals-mean)
        
        mean = np.round(mean,2)
        std = np.round(std,2)
        
        legend = ('Detectable Titers\n '
                  fr'($\mu$={mean:.2f}, '
                  fr'$\sigma$={std:.2f})'
                  )
        
        legends.append(legend)

    if len(undetectable_residuals)>0:
        ax.hist(undetectable_residuals, color=undetectable_titer_color, alpha=alpha,
                   edgecolor='red', density=density, **fig_options)
        
        mean = np.nanmean(undetectable_residuals)
        std = np.nanstd(undetectable_residuals-mean)
        
        mean = np.round(mean,2)
        std = np.round(std,2)
        
        legend = ('Undetectable Titers\n '
                  fr'($\mu$={mean:.2f}, '
                  fr'$\sigma$={std:.2f})'
                  )
        
        
        
        legends.append(legend)
        
        
    ax.legend(legends)
    ax.set_xlabel('Measured Titer - Fitted Titer')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
        
    return fig,ax