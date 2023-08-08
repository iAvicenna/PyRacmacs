#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 13:20:20 2022

@author: Sina Tureli
"""

import numpy as np
from . import colors

class AxesUniformizer():
    
    def __init__(self, equal_axes=True):
        
        self.axes = []
        self.equal_axes = equal_axes
        self.xlims = [np.inf, -np.inf]
        self.ylims = [np.inf, -np.inf]
        
    def add_axis(self, ax):
        
        self.axes.append(ax)
        
        ax_xlims,ax_ylims = ax_lims(ax)
        
        self.xlims = self.update_limits(self.xlims, ax_xlims)
        self.ylims = self.update_limits(self.ylims, ax_ylims)
        
        if self.equal_axes:
            self.xlims = self.update_limits(self.xlims, self.ylims)
            self.ylims = self.xlims
            
        
    def update_limits(self, limits, new_limits):
        
        limits = [min(limits[0],new_limits[0]),max(limits[1],new_limits[1])]
        
        return limits
        

    def uniformize(self):
        
        for ax_ind in range(len(self.axes)):
            
            ax = self.axes[ax_ind]
            
            ax.set_xlim(self.xlims)
            ax.set_ylim(self.ylims)
            
            if self.equal_axes:
                make_plot_square(ax)

    def add_log_titer_ticklabels(self):
        
        for ax_ind in range(len(self.axes)):
            add_log_titer_ticklabels(self.axes[ax_ind])
            
    def add_grid(self):
        
        for ax_ind in range(len(self.axes)):
            add_grid(self.axes[ax_ind])
            
    def draw_diagonal(self):
        
        for ax_ind in range(len(self.axes)):
            draw_diagonal(self.axes[ax_ind])
            
    def do_all(self):
        
        self.uniformize()
        self.add_log_titer_ticklabels()
        self.add_grid()
        self.draw_diagonal()
                

def ax_lims(ax):
    
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    
    return (x0,x1),(y0,y1)

def make_plot_square(ax):
    
    xlims,ylims = ax_lims(ax)
    
    max_val = max(xlims[1], ylims[1])
    min_val = min(xlims[0], ylims[0])
    val_lims = [min_val, max_val]
    
    ax.set_xlim(val_lims)
    ax.set_ylim(val_lims)
    ax.set_aspect(abs(val_lims[1] - val_lims[0])/abs(val_lims[1] - val_lims[0]))
           
def draw_diagonal(ax, plot_args={'color':'black','linestyle':':','alpha':0.5}):
    
    xlims,ylims = ax_lims(ax)
    dx = 0.01*(xlims[1] - xlims[0])
    dy = 0.01*(ylims[1] - ylims[0])
    
    xvals = np.linspace(xlims[0]-dx, xlims[1]+dx, 100)
    yvals = np.linspace(ylims[0]-dy, ylims[1]+dy, 100)
    
    ax.plot(xvals, yvals, **plot_args)
    
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
def add_grid(ax, grid_args = {'linestyle':':','alpha':0.2, 'color':'black'}):
    
    ax.grid('on', **grid_args)
    
def format_titer_str(titer_str):
    
    coef = int(np.ceil(np.log10(float(titer_str))))

    if coef>0:
        return int(titer_str)
    else:
        return np.round(titer_str,-coef+1)
    
def add_log_titer_ticklabels(ax, step=1):
    
    xlims,ylims = ax_lims(ax)
    
    x0 = np.ceil(xlims[0])
    x1 = np.floor(xlims[1])
    
    y0 = np.ceil(ylims[0])
    y1 = np.floor(ylims[1])
    
    xticks = np.arange(x0, x1+step, step)
    yticks = np.arange(y0, y1+step, step)
    
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    
    xticklabels = [format_titer_str(10*2**x) for x in xticks]
    yticklabels = [format_titer_str(10*2**x) for x in yticks]
    
    ax.set_xticklabels(xticklabels,rotation=90)
    ax.set_yticklabels(yticklabels)
    
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    
    
def convert_color_to_hex(color):
    
    color = color.replace('grey','gray')
    
    if color[0]=='#':
        return color
    else:
        return colors.color_name_to_rgb[color].hex_format()    
