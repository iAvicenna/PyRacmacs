#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 18:39:59 2021

@author: iavicenna
"""
import numpy as np
import matplotlib.patches as ptchs
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.cm import get_cmap
from matplotlib import colors as mplcolors


def get_text_dim(text, fig=None, ax=None, fontfamily = rcParams['font.family'], 
                 fontsize=rcParams['font.size'], rotation=0, relative_to = 'fig',
                 return_proportion=True):
    
    close_plot = False
    
    if fig is None and ax is None:
        fig,ax = plt.subplots(1,1, figsize=(10,10))
        close_plot = True
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    t = ax.text(.5*(xlim[0] + xlim[1]), .5*(ylim[0] + ylim[1]), text, 
                fontfamily=fontfamily, fontsize=fontsize, rotation=rotation)
    
    fw,fh = fig.get_size_inches()*fig.dpi
        
    bb = t.get_window_extent(renderer = fig.canvas.get_renderer())
    
    t.remove()
    
    if relative_to == 'fig':
        
        
        height = float(bb.y1 - bb.y0)
        width = float(bb.x1 - bb.x0)
        
        if close_plot:
            plt.close(fig)
        
        if not return_proportion:
            return (height/fig.dpi, width/fig.dpi)
        else:
            return (height/fh, width/fw)
        
        
    elif relative_to == 'ax':
        transf = ax.transData.inverted()
        bb_datacoords = bb.transformed(transf)
        
        height = float(bb_datacoords.y1 - bb_datacoords.y0)
        width = float(bb_datacoords.x1 - bb_datacoords.x0)
        
        if close_plot:
            plt.close(fig)
        
        return (height,width)

def generate_color_categories(categories, return_hex=False):
    """Generate n colors, useful for coloring different categories.
       author: https://www.andersle.no/
    """

    ncol = len(categories)
    if ncol <= 10:
        colors = [i for i in get_cmap('tab10').colors]
    elif 10 < ncol <= 20:
        colors = [i for i in get_cmap('tab20').colors]
    elif 25 < ncol <= 256:
        cmap = get_cmap(name='viridis')
        colors = cmap(np.linspace(0, 1, ncol))
    else:
        raise ValueError('Maximum 256 categories')
    color_map = {}
    for i, key in enumerate(categories):
        if not return_hex:
            color_map[key] = colors[i]
        else:
            color_map[key] = mplcolors.to_hex(colors[i])
        
        
    return color_map


def add_object(ax, object_type, position, options=None):
    
    if options is None:
            options = {}
    else:
        assert isinstance(options,dict)
    
    if object_type.upper() == 'VERTICAL STRIP':
        assert isinstance(position,(float,int))
        
        if 'width' not in options:
            width = 0.5
            options['width'] = 0.5
        else:
            width = options['width']
        if 'facecolor' not in options:
            options['facecolor'] = 'gray'
        if 'color' not in options:
            options['edgecolor'] = 'black'
        if 'alpha' not in options:
            options['alpha'] = 0.3
        if 'linestyle' not in options:
            options['linestyle'] = None
        if 'linewidth' not in options:
            options['linewidth'] = None
        
        
        options.pop('width')
        ax.axvspan(position-width, position+width,
                   **options)
    
    if object_type.upper() == 'TEXT BOX':
        
        assert 'text' in options
        text = options['text']
        fig = options['fig']

        x0 = position[0]
        x1 = position[1]
        y0 = position[2]
        y1 = position[3]
        
        if 'fontsize' not in options:
            options['fontsize']=12
            
        h,w = get_text_dim(text, fig, ax, fontsize=options['fontsize'], relative_to='ax')
        
        text_x = (x1 + x0 - w)/2
        text_y = (y1 + y0 - 0.7*h)/2 
         
        
        if 'facecolor' not in options:
            options['facecolor'] = 'gray'
        if 'alpha' not in options:
            options['alpha'] = 0.1
        if 'edgecolor' not in options:
            options['edgecolor'] = 'black'
        if 'linewidth' not in options:
            options['linewidth'] = 2
        if 'angle' not in options:
            options['angle'] = 0
        
        ax.text(x=text_x, y=text_y, s=text, fontsize=options['fontsize'], rotation=options['angle'])
        
        rect1 = ptchs.Rectangle((x0, y0), height=y1-y0, width=x1-x0,
                         
                         fill=True, facecolor=options['facecolor'], clip_on=False,
                         alpha=options['alpha'], edgecolor=options['edgecolor'],
                         linewidth=options['linewidth'])
  
        ax.add_patch(rect1)
        
        if options['edgecolor'] != 'black':
            
            bbox = get_ax_bbox(ax,fig)
            ratio = (bbox.y1-bbox.y0)/(bbox.x1-bbox.x0)
            
            yoffset=0.05
            xoffset = yoffset * ratio
            rect2 = ptchs.Rectangle((x0+xoffset, y0+yoffset), height=y1-y0-2*yoffset, width=x1-x0-2*xoffset,
                             fill='none', clip_on=False, facecolor='none',
                             alpha=1, edgecolor=options['edgecolor'],
                             linewidth=1.5*options['linewidth'])
        
            ax.add_patch(rect2)
        
    
    return ax

def get_ax_bbox(ax,fig):
    return ax.bbox.transformed(transform=fig.transFigure.inverted())


def fig_layout(ax_sizes, ax_rel_rects, htab=0.2, vtab=0.2):
    
    '''
    Given a matrix of axis sizes (the height and width of the axis including labels etc)
    and ax_rel_rects, the relative rectangle of the axis with respect to its height and width
    this function produces a grid figure from these. Inputs must be inches and must 
    be a 2d array of lists. The following would produce a layout with 2 axis on the first row.
    
    ax_sizes = [ [ [2.15, 2.51], [2.15, 2.51] ] ] => 1x2 array of lists
    ax_rel_rects = [ [ [ [0.51, 0.0, 2.0, 1.87], [0.51, 0.0, 2.0, 1.87] ] ] => 1x2 array of lists
    '''
    
    s1,s2 = ax_sizes.shape
    
    fig_height = 0
    fig_width = 0
    ax_rects = np.zeros((s1,s2),dtype=object)
    
    for i in range(s1):
        ax_widths = [x[1] for x in ax_sizes[i,:]]
        fig_width = max([fig_width, np.sum(ax_widths) + htab*s2])

    fig_width += htab

    for i in range(s2):
        ax_heights = [x[0] for x in ax_sizes[:,i]]
        fig_height = max([fig_height, np.sum(ax_heights) + vtab*s1])
        
    fig_height += vtab
    
    for i0 in range(s1):
        ax_widths = [x[1] for x in ax_sizes[i0,:]]
        
        for i1 in range(s2):
            ax_heights = [x[0] for x in ax_sizes[:,i1]]
            
            ax_x = np.sum(ax_widths[:i1]) + htab*(i1+1)
            ax_y = fig_height - (np.sum(ax_heights[:i0+1]) + vtab*(i0+1))
            
            ax_abs_x = ax_rel_rects[i0,i1][0]
            ax_abs_y = ax_rel_rects[i0,i1][1]
            ax_abs_width = ax_rel_rects[i0,i1][2]
            ax_abs_height = ax_rel_rects[i0,i1][3]
            
            ax_rects[i0,i1] = [(ax_x + ax_abs_x)/fig_width, (ax_y + ax_abs_y)/fig_height, ax_abs_width/fig_width, ax_abs_height/fig_height]
            
    
    return (fig_width, fig_height), ax_rects
            
    
def calculate_axis_loc(xtick_labels, ytick_labels, xlabel, ylabel, title, legend,
                      fontfamily = rcParams['font.family'], xrotation=90,
                      fontsize_dict= {'tick_labels':12, 'axis_label':14, 'title':14, 'legend':12},
                      xscale=1, yscale=1, show_x=True, show_y=True):
    
    '''
    This function produces ax height and width and its relative rectangle
    given information about its x and y axis objects, title and legend. ax height and
    width encompasses everything including labels, legends etc where as ax_rel_rect
    gives the relative position of the plotting box. Everything is calculated in inches.
    '''
    
    show_x = int(show_x)
    show_y = int(show_y)
    
    temp_fig, temp_ax = plt.subplots(1,1,figsize=(5,5))

    xtick_label_sizes = [get_text_dim(x,temp_fig, temp_ax, rotation=xrotation, 
                                      fontfamily=fontfamily, return_proportion=False,
                                      fontsize=fontsize_dict['tick_labels']) for x in xtick_labels]
    
    xtick_label_widths = [max(2*x[1],0.7) for x in xtick_label_sizes]
    
    xtick_label_heights = [1.5*x[0] for x in xtick_label_sizes]
    
    ytick_label_sizes = [get_text_dim(x, temp_fig, temp_ax,  fontfamily=fontfamily, return_proportion=False,
                                      fontsize=fontsize_dict['tick_labels']) for x in ytick_labels]
    ytick_label_widths = [2*x[1] for x in ytick_label_sizes]
    
    ytick_label_heights = [1.5*x[0] for x in ytick_label_sizes]
    
    legend_sizes = [get_text_dim(x,temp_fig, temp_ax, fontfamily=fontfamily, return_proportion=False,
                                      fontsize=fontsize_dict['legend'])
                         for x in legend]
    legend_widths = [1.5*x[1] for x in legend_sizes]
    
        
    
    title_height,title_width = get_text_dim(title, temp_fig, temp_ax, 
                                            fontsize=fontsize_dict['title'], return_proportion=False,)
    
    xlabel_height, xlabel_width = get_text_dim(xlabel, temp_fig, temp_ax, 
                                               fontsize=fontsize_dict['axis_label'], 
                                               rotation=xrotation, return_proportion=False,)
    
    ylabel_height, ylabel_width = get_text_dim(ylabel, temp_fig, temp_ax, 
                                               fontsize=fontsize_dict['axis_label'], 
                                               rotation=90, return_proportion=False,)
    
    title_height *= 2
     
    ax_width = xscale*sum(xtick_label_widths) + (max(ytick_label_widths) + ylabel_width)*show_y + max(legend_widths)
    
    ax_height = yscale*sum(ytick_label_heights) + max(xtick_label_heights)*show_x + title_height + xlabel_height

    ax_rel_rect = [(max(ytick_label_widths)+ylabel_width)*show_y, max(xtick_label_heights)*show_x, 
                   ax_width - show_y*(max(ytick_label_widths)+ylabel_width)-max(legend_widths), 
                   ax_height - show_x*max(xtick_label_heights) - title_height  ]

    plt.close(temp_fig)
    
    return [ax_height, ax_width], ax_rel_rect




aa_colors = [[],[]]

aa_colors[0] ={
    'G':"#FFA50080",
    'A':"#FFA50080",
    'S':"#FFA50080",
    'T':"#FFA50080",
    'C':"#00FF0080",
    'V':"#00FF0080",
    'I':"#00FF0080",
    'L':"#00FF0080",
    'P':"#00FF0080",
    'F':"#00FF0080",
    'Y':"#00FF0080",
    'M':"#00FF0080",
    'W':"#00FF0080",
    'N':"#FF00FF80",
    'Q':"#FF00FF80",
    'H':"#FF00FF80",
    'D':"#FF000080",
    'E':"#FF000080",
    'K':"#0000FF80",
    'R':"#0000FF80",
    '*':"#7F7F7F80",
    'X':"#7F7F7F80",
    '-':"#7F7F7F80",
    '?':"#7F7F7F80",
    '#':"#7F7F7F80",
    ' ':"#7F7F7F80",
    }


aa_colors[1] ={
    'H':'blue',
    'K':'blue',
    'R':'blue',
    'D':'red',
    'E':'red',
    'S':'green',
    'T':'green',
    'N':'green',
    'Q':'green',
    'A':'magenta',
    'V':'magenta',
    'L':'magenta',
    'I':'magenta',
    'M':'magenta',
    'F':'white',
    'W':'white',
    'Y':'white',
    'P':'brown',
    'G':'brown',
    'C':'orange',
    'B':'gray',
    'Z':'gray',
    '*':'gray',
    'X':'gray',
    '-':'gray',
    '?':'gray',
    '#':'gray',
    ' ':'gray'
    }


list_aa = list('ARNDBCEQZGHILKMFPSTWYV')
extra_chars = ['X','*','-','?','#',' ']