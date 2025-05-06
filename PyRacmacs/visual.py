#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:29:36 2022

@author: Sina Tureli
"""

import PyRacmacs as pr
import tempfile
import webbrowser
import matplotlib.pyplot as plt
import numpy as np
from . import Racmacs
from .structures import RacViewerOptions

def view(racmap: pr.RacMap,
         export_path = None,
         selfcontained = True,
         display=True,
         optimization_number=0,
         options = None
         ):

    if options is None:
        options = {}
    if isinstance(options, dict):
        options = RacViewerOptions(**options)

    if export_path is None:
        with tempfile.NamedTemporaryFile(suffix='.html') as temp_fp:
            export_path = temp_fp.name
    Racmacs.export_viewer(racmap._acmap_R,
                          export_path,
                          optimization_number=optimization_number+1,
                          options = options.options_R
                          )

    if display:
        webbrowser.open(export_path)


def plot(racmap: pr.RacMap,
         optimization_number=0,
         show_ag=True,
         show_sr=True,
         scale=1,
         make_square=False,
         alphas=None,
         linewidths=None,
         sizes=None,
         padding=0.015,
         plot_orders=None,
         grid_props=None,
         stress=None,
         map_lims=None
         ):

  if stress is None:
      stress = ''
  else:
      stress = f'S:{stress:.2f}'

  if grid_props is None:
    grid_props = {}

  grid_props = dict({'color':'grey', 'alpha':0.15, 'linestyle':'-'},
                    **grid_props)

  if map_lims is None:
    map_lims = _get_map_limits(racmap, show_ag, show_sr)
  else:
    assert 'x' in map_lims and 'y' in map_lims


  if make_square:
    lims = [min(map_lims['x'][0], map_lims['y'][0]) - 1 ,
            max(map_lims['x'][1], map_lims['y'][1]) + 1]
    map_lims['x'] = lims
    map_lims['y'] = lims


  if show_sr and show_ag:
    coordinates = racmap.coordinates
  elif show_sr and not show_ag:
    coordinates = racmap.sr_coordinates
  elif not show_sr and show_ag:
    coordinates = racmap.ag_coordinates
  else:
    raise ValueError("Can't have false for both show_sr and show_ag")

  s1,s2 = coordinates.shape
  fill_colors = racmap.ag_fills + racmap.sr_fills
  shapes = racmap.ag_shapes + racmap.sr_shapes
  outlines = racmap.ag_outlines + racmap.sr_outlines
  shapes = ['o' if shape.lower()=='circle' else 's' for shape in shapes]
  fill_colors = [color.replace('transparent',"None") for color in fill_colors]

  if alphas is None:
    alphas = [0.5 for _ in range(s1)]

  if fill_colors is None:
    fill_colors = ['grey' for _ in range(s1)]

  if shapes is None:
    shapes = ['o' for _ in range(s1)]

  if outlines is None:
    outlines = ['black' for _ in range(s1)]

  if linewidths is None:
    linewidths = [2 for _ in range(s1)]

  if sizes is None:
    sizes = [1000 for _ in range(s1)]



  if plot_orders is None:
      zorders = [i for i in range(s1+s2)]
  else:
      zorders = [s1+s2-po for po in plot_orders]

  ratio = (map_lims['x'][1] - map_lims['x'][0])/(map_lims['y'][1] - map_lims['y'][0])

  fig = plt.figure(figsize=(5*scale, 5*scale/ratio))
  ax = fig.add_axes([padding, padding*ratio, 1-2*padding, 1-2*padding*ratio])
  ax.set_aspect("equal", anchor="C")

  for i in range(s1):

      x = _project(coordinates[i,0], map_lims['x'])
      y = _project(coordinates[i,1], map_lims['y'])

      if x in map_lims['x'] or y in map_lims['y']:
          marker = _determine_marker(x, y, map_lims['x'], map_lims['y'])
          alpha = 0.2
      else:
          marker = shapes[i]
          alpha = alphas[i]

      ax.scatter(x,y,
                 alpha=alpha, linewidths=linewidths[i],
                 color=fill_colors[i], marker=marker,
                 edgecolors=outlines[i], s=sizes[i],
                 zorder=zorders[i], clip_on=False)


  ax.set_xlim(*map_lims['x'])
  ax.set_ylim(*map_lims['y'])
  _format_2d_plot(ax, grid_props)

  #ax.text(map_lims['x'][0]+1, map_lims['y'][0]+1, stress, fontsize=12)

  return fig,ax



def _project(val, lims):

    if val<lims[0]:
        val = lims[0]
    elif val>lims[1]:
        val = lims[1]

    return val

def _determine_marker(xval, yval, xlims, ylims):

    if xval == xlims[0]:
        marker = '<'
    elif xval == xlims[1]:
        marker = '>'
    elif yval == ylims[0]:
        marker = 'v'
    elif yval == ylims[1]:
        marker = '^'
    return marker


def _get_map_limits(racmap, show_ag, show_sr):

    if show_ag and show_sr:
        coordinates = racmap.coordinates
    elif show_ag:
        coordinates = racmap.ag_coordinates
    elif show_sr:
        coordinates = racmap.sr_coordinates
    else:
        raise ValueError('show_antigens and show_sera are both False')


    xlim = [np.nanmin(coordinates[:,0])-1, np.nanmax(coordinates[:,0])+1]
    ylim = [np.nanmin(coordinates[:,1])-1, np.nanmax(coordinates[:,1])+1]

    if racmap.number_of_dimensions==3:
        zlim = [np.nanmin(coordinates[:,2])-1, np.nanmax(coordinates[:,2])+1]
    else:
        zlim = [-1, 1]

    lims = {}

    lims['x'] = xlim
    lims['y'] = ylim
    lims['z'] = zlim

    return lims


def _format_2d_plot(ax, grid_props):

    xlims = ax.get_xlim()
    ylims = ax.get_ylim()

    ax.set_xticks(np.arange(xlims[0], xlims[1]+1,1))
    ax.set_yticks(np.arange(ylims[0], ylims[1]+1,1))

    ax.set_xticklabels(['' for _ in np.arange(xlims[0], xlims[1]+1, 1)])
    ax.set_yticklabels(['' for _ in np.arange(ylims[0], ylims[1]+1, 1)])

    ax.grid('on', **grid_props)

    _remove_tick_lines(ax)


def _remove_tick_lines(ax):

    for tick in ax.xaxis.get_major_ticks():
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)

    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)
