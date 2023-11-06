#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 22:37:11 2021

@author: iavicenna
"""

import numpy as np
import pandas as pd
import plotly as plotly

import matplotlib.pyplot as plt
import PyRacmacs as pr

from plotly.subplots import make_subplots as _make_subplots

from matplotlib.cm import get_cmap
from matplotlib import colors as mplcolors, colormaps

cm = colormaps["tab20"]
cm = [cm(i) for i in range(0,20,2)] + [cm(i) for i in range(1,20,2)]


__all__ = ['interactive_titer_plot_from_racmap', 'box_plot', # titer plots
           'bootstrap_bar', 'bootstrap_agcolored_map', 'bootstrap_agvol_vs_mean_sr_dist' # bootstrap related plots
           ]

try:
  import seaborn as sb
except ModuleNotFoundError:
  __all__.remove("box_plot")


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


## TITER PLOTTING #############################################################

def box_plot(racmap, groups, group_colors=None, plot_type='sera',
             invert_distances=True, xlabel=None, swarm=True, ax=None, display=False,
             ylims=None, xscale=1):

    '''
    groups is the group names for the antigens.
    group_colors is a dict mapping each group name to its color
    '''

    assert isinstance(plot_type,str),'Average type must be a string'
    plot_type = plot_type.upper()
    assert plot_type in ['SERA','ANTIGENS'], 'Average type must be either sera or antigens'
    sb.set_context("paper", font_scale=1.5)


    data = racmap.log_titers
    ylabel = 'Titer in log2 scale \n'

    if plot_type == 'ANTIGENS':
        data = data.T
        if xlabel is None:
            xlabel='Antigen Group'
        index_name = 'Serum'

        if group_colors is None:
           group_colors  = generate_color_categories(set(racmap.sr_names))

    elif xlabel is None:
        xlabel = 'Serum Group'
        index_name = 'Antigen'

        if group_colors is None:
           group_colors  = generate_color_categories(set(racmap.ag_names))


    s1,s2 = data.shape

    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=(xscale*2*len(groups)*len(group_colors)/16,5))

    index = list(racmap.ag_names)
    assert s2==len(groups)
    long_data = pd.DataFrame(np.zeros((s1*s2,3)),columns=['group','titers',index_name], dtype=object)
    counter = 0
    for i in range(s1):
        for j in range(s2):
            long_data.iloc[counter,0] = groups[j]
            long_data.iloc[counter,1] = data[i,j]
            long_data.iloc[counter,2] = index[i]
            counter += 1

    ax = sb.boxplot(x='group', y='titers', hue=index_name, data=long_data, palette=group_colors,ax=ax)

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    if ylims is None:
        ylims = ax.get_ylim()
    else:
        ax.set_ylim(ylims)

    yticks = range(int(ylims[0]), int(ylims[1])+1, 1)
    ytick_labels = [int(10*2**x) for x in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)

    plt.tight_layout()
    fig = plt.gcf()

    if not display:
        plt.close(fig)

    return fig,ax


def interactive_titer_plot_from_racmap(racmap, serum_group_ids, y_max=None, antigen_group_ids = None,
                                       line_widths=None, antigen_colors=None, line_styles=None,
                                       xlabel=None, ylabel=None, markers=None, lg_axis_width=250,
                                       xscale=1, scale=1, plot_type='line', cover_thresholded=None,
                                       threshold=10, ingroup_order=None, order_ingroups=False,
                                       revert_ingroup_order=None, invert_distances=False,
                                       transpose_grid=False, marker_sizes=None
                                       ):


    x_names = np.array(racmap.sr_names)
    set_serum_group_ids = list(set(serum_group_ids))
    I = np.argsort([serum_group_ids.index(x) for x in set_serum_group_ids])
    set_serum_group_ids = [set_serum_group_ids[x] for x in I]
    lessthans = racmap.lessthans.to_numpy()
    serum_II = []
    invert_yticklabels=False


    for serum_group_id in set_serum_group_ids:
        serum_II.append([ind for ind,x in enumerate(serum_group_ids) if x==serum_group_id])

    if antigen_group_ids is None:
        antigen_group_ids = racmap.ag_names

    if antigen_colors is None:
        antigen_colors = {x:y for x,y in zip(racmap.ag_names, racmap.ag_fills)}

    set_antigen_group_ids = list(set(antigen_group_ids))
    I = np.argsort([antigen_group_ids.index(x) for x in set_antigen_group_ids])
    set_antigen_group_ids = [set_antigen_group_ids[x] for x in I]
    antigen_II = []

    for antigen_group_id in set_antigen_group_ids:
        antigen_II.append([ind for ind,x in enumerate(antigen_group_ids) if x==antigen_group_id])


    if ylabel is None:
        ylabel = 'Titer in log2 scale'

    if cover_thresholded is None:
        cover_thresholded=True

    data = racmap.log_titers

    s1,s2 = data.shape

    if revert_ingroup_order is None:
        revert_ingroup_order=False


    width_ratios = [serum_group_ids.count(x)/len(serum_group_ids) for x in set_serum_group_ids]

    if order_ingroups:
        if ingroup_order is None:

            ingroup_order = np.array(list(range(len(serum_group_ids))))

            for serum_group_id in set_serum_group_ids:
                J = [ind for ind,x in enumerate(serum_group_ids) if x==serum_group_id]

                vals = np.nanmean(racmap.log_titers[:,J],axis=0)

                if any(not np.isnan(x) for x in vals):
                    K = np.argsort(vals).flatten()

                    if revert_ingroup_order:
                        K = K[::-1]

                    data[:,J] = data[:,[J[x] for x in K]]
                    lessthans[:,J] = lessthans[:,[J[x] for x in K]]
                    ingroup_order[J] = ingroup_order[[J[x] for x in K]]
                    x_names[J] = x_names[[J[x] for x in K]]
        else:
            data[:,:] = data[:,ingroup_order]
            lessthans[:,:] = lessthans[:,ingroup_order]
            x_names = x_names[ingroup_order]

    x_names = list(x_names)

    antigen_colors = [antigen_colors[x] for x in set_antigen_group_ids]

    data[lessthans] -= 0.5

    if plot_type == 'line':

        fig = _plotly_line_subplots(data=data, serum_group_indices=serum_II, serum_group_names=set_serum_group_ids,
                                    antigen_group_indices = antigen_II, antigen_group_names = set_antigen_group_ids,
                                    line_names=set_antigen_group_ids, x_names = x_names, y_max=y_max,
                                    width_ratios=width_ratios, xscale=xscale, log_scale=True,
                                    xlabel=xlabel, ylabel=ylabel, line_widths=line_widths, markers=markers,
                                    line_colors=antigen_colors, line_styles=line_styles,lg_axis_width=lg_axis_width,
                                    scale=scale, cover_thresholded=cover_thresholded, threshold=threshold,
                                    invert_yticklabels=invert_yticklabels, marker_sizes=marker_sizes)
    elif plot_type == 'box':

        fig = _plotly_box_subplots(data=data, serum_group_ids=serum_group_ids, serum_group_names=set_serum_group_ids,
                                  antigen_group_ids = antigen_group_ids, antigen_group_names = set_antigen_group_ids,
                                  x_names = set_antigen_group_ids, log_scale=True,
                                  ylabel=ylabel, xlabel=xlabel, box_colors=antigen_colors,
                                  xscale=xscale, cover_thresholded=cover_thresholded, threshold=threshold, invert_yticklabels=invert_yticklabels)

    return fig, ingroup_order


def _plotly_box_subplots(data, serum_group_ids, serum_group_names, x_names, antigen_group_ids, antigen_group_names,
                        size=400, xscale=1, yscale=1, log_scale=True, y_max=None, y_min=None, marker_size=None,
                        xlabel=None, ylabel=None, box_colors=None, cover_thresholded=False, threshold=10,
                        hide_boxes=False, add_lines=False, hide_xtick_labels=False, m=[40,0,10,0],
                        trace_line_colors=None,invert_yticklabels=False):

    assert set(serum_group_ids) == set(serum_group_names), 'As a set serum group ids must be identical to serum group names'
    assert set(antigen_group_ids) == set(antigen_group_names), 'As a set antigen group ids must be identical to antigen group names'

    data = data.copy()
    s1,s2 = data.shape
    s3 = len(antigen_group_names)

    if hide_xtick_labels:
        x_names = len(x_names)*['']

    if add_lines:
        jitter=0
    else:
        jitter=0.1

    assert s1==len(antigen_group_ids),'Number of antigen group ids should be equal to number of rows of data'
    assert s2==len(serum_group_ids),'Number of serum group ids should be equal to number of columns of data'

    fig_height = size*yscale
    fig_width = xscale*s3*len(serum_group_names)*size/5

    ratio = np.sqrt(fig_height*fig_width/len(serum_group_names))/360
    hs = ratio*10/fig_width

    box_ratio = ratio*8/len(x_names)

    if y_max is None:
        y_max = int(np.ceil(np.nanmax(data)))+1

    if y_min is None:
        y_min = int(np.floor(min([np.log2(threshold/10), np.nanmin(data)])))

    if box_colors is None:
        box_colors = generate_color_categories(range(s3),return_hex=True)
    else:
        assert len(box_colors)==s3, f'number of box colors ({len(box_colors)}) should equal number of groups ({s3})'

    if trace_line_colors is None:
        trace_line_colors = ['rgba(0,0,0,0.5)']*s2
    else:
        assert len(trace_line_colors) == s2, 'number of trace line colors should be equal to number of sera'

    fig = _make_subplots(rows=1, cols=len(serum_group_names),
                         horizontal_spacing=hs,
                         subplot_titles=[f'<b>{x}</b>' for x in serum_group_names])

    fig.update_layout(
        autosize=False,
        width=fig_width,
        height=fig_height,
        margin=dict(
            t=m[0]*ratio,
            b=m[1]*ratio,
            r=m[2]*ratio,
            l=m[3]*ratio,
        ),
    )

    if marker_size is None:
        marker_size = 12*ratio

    for j, group_id1 in enumerate(serum_group_names):

        I = [ind for ind,x in enumerate(serum_group_ids) if x==group_id1]

        if j==0:
            ytickvals = list(range(y_min,y_max))
            if log_scale:
                yticktext = [str(np.round(10*2**x,0)) for x in ytickvals]
            else:
                yticktext = ytickvals

            xaxis_title = xlabel
            yaxis_title = ylabel

        else:
            ytickvals  = list(range(y_min,y_max))
            yticktext = ['']*len(ytickvals)
            xaxis_title = ''
            yaxis_title = ''

        if invert_yticklabels:
            yticktext = yticktext[::-1]

        if add_lines:
            y_line_data = np.zeros((len(I),s3))*np.nan
            x_line_data = np.zeros((len(I),s3))*np.nan

        for k, group_id2 in enumerate(antigen_group_names):

            K = [ind for ind,x in enumerate(antigen_group_ids) if x==group_id2]
            y_data = np.reshape(data[np.ix_(K,I)],(len(I)*len(K),))

            if add_lines:
                y_line_data[:,k] = np.nanmean(np.reshape(data[np.ix_(K,I)], (len(K),len(I))),axis=0).T
                x_line_data[:,k] = k+1

            if hide_boxes:
                fillcolor='rgba(0,0,0,0)'
                marker = dict(size = marker_size, color=box_colors[k],
                          line=dict(
                    color='black',
                    width=1
                )
                )
                line = dict(color='rgba(0,0,0,0)', width=1*box_ratio)

            else:
                fillcolor=box_colors[k]
                marker = dict(size = marker_size)
                line = dict(color='black', width=1*box_ratio)

            fig.add_trace(plotly.graph_objects.Box(
                y=y_data,
                x=[k+1]*len(K)*len(I),
                name=x_names[k],
                fillcolor=fillcolor,
                showlegend=False,
                boxpoints='all',
                pointpos=0,
                jitter=jitter,
                line=line,
                marker=marker

            ),row=1, col=j+1)

        if add_lines:
            for k in range(len(I)):
                line = dict(color=trace_line_colors[I[k]], width=2*ratio)

                fig.add_trace(plotly.graph_objects.Scatter(x=x_line_data[k,:].flatten(), y=y_line_data[k,:].flatten(),
                                     showlegend=False, mode='lines',line=line,
                                     ), row=1, col=j+1)


        fig.update_xaxes(tickangle=-90)

        fig.update_layout(
            {
            f'yaxis{j+1}':{'range': [y_min, y_max], 'tickvals':ytickvals, 'ticktext':yticktext,
                           'autorange':False, 'gridcolor':'#f4eeee','mirror':True, 'title':yaxis_title},
            f'xaxis{j+1}':{'range':[0.5, s3+1],'gridcolor':'#f4eeee','mirror':True, 'title':xaxis_title,
                           'tickvals':list(range(1,s3+1)), 'ticktext':x_names},

            },
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(
                #size=18,
                color="Black"
            )
        )

        if cover_thresholded:
            x0 = -1
            x1 = s3+2

            fig.add_trace(plotly.graph_objects.Scatter(x=[x0,x0,x1,x1], y=[y_min,np.log2(threshold/10),np.log2(threshold/10),y_min,y_min],
                                                    fillcolor='rgba(128,128,128,0.75)',
                                                    showlegend=False,
                                                    fill="toself",
                                                    mode="lines",
                                                    line=dict(color='rgba(128,128,128,0)')), row=1, col=j+1)

    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

    return fig


def _plotly_line_subplots(data, serum_group_indices, serum_group_names, width_ratios,
                         antigen_group_indices, antigen_group_names, line_names, x_names,
                         size=400, xscale=1, log_scale=True, y_max=None, y_min=None,
                         xlabel=None, ylabel=None, line_widths=None, markers=None,
                         line_colors=None, line_styles=None, lg_axis_width=250, scale=1,
                         cover_thresholded=False, threshold=10, invert_yticklabels=False,
                         marker_sizes=None):

    s1,s2 = data.shape
    s3 = len(antigen_group_indices)

    size = size*scale
    T = serum_group_indices[-1][-1]

    if y_max is None:
        y_max = int(np.ceil(np.nanmax(data)))+1

    if y_min is None:
        y_min = int(np.floor(min([np.log2(threshold/10), np.nanmin(data)])))

    if y_min==y_max:
        print(f'Warning max value {y_max} is equal to min_value {y_min}. Lowering y_min by 2')
        y_min -= 2


    fig_height = (y_max - y_min)/10*size + 100*max([len(x) for x in x_names])/11
    fig_width = T/7*xscale*size + lg_axis_width

    hs = 10/fig_width

    widths = [fig_width*x  for x in width_ratios]
    widths[0] += 100

    widths = [lg_axis_width] + widths
    width_ratios = [x/sum(widths) for x in widths]

    serum_group_names = [''] + serum_group_names


    if line_widths is None:
        line_widths = [2]*s3
    else:
        assert len(line_widths) == s3, 'number of line widths should be same as number of antigen groups'
        assert all([isinstance(x,(int,float)) for x in line_widths]), 'line widths should be a list of numeric entries'

    if line_styles is None:
        line_styles = ['solid']*s3
    else:
        assert len(line_styles) == s3, 'number of line widths should be same as number of antigen groups'
        assert all([x in ['dot','dash','solid','dashdot']  for x in line_styles]), 'line styles can be only dash, dashdot, dot, solid'

    if markers is None:
        markers = [0]*s3
    else:
        assert len(markers) == s3, 'number of markers should be same as number of antigen groups'
        assert all([isinstance(x,(str)) for x in markers]), 'markers should be a list of integers'

    if marker_sizes is None:
        marker_sizes = [10]*s3
    else:
        assert len(marker_sizes) == s3, 'number of marker sizes should be same as number of antigen groups'

    if line_colors is None:
        line_colors = generate_color_categories(range(s3), return_hex=True)
    else:
        assert len(line_colors) == s3, 'number of line_colors should be same as number of antigen groups'

    fig = _make_subplots(rows=1, cols=len(serum_group_names),
                        column_widths=width_ratios, horizontal_spacing=hs,
                        subplot_titles=[f'<b>{x}</b>' for x in serum_group_names],
                        )

    fig.update_layout(
        autosize=False,
        width=fig_width,
        height=fig_height,
        margin=dict(
            t=40,
            r=10,
            l=0
        ),
        font=dict(
                #size=18,
                color="Black"
            )
    )

    fig.update_xaxes(tickangle=-45)


    for ind1,I in enumerate(serum_group_indices):

        s1,s2 = data.shape

        if ind1 == 0:
            showlegend=True
        else:
            showlegend=False

        subplot_x_names = [x_names[x] for x in I]
        x_vals = list(range(1,len(I)+1))

        if ind1==0:
            y_tickvals = list(range(y_min,y_max))
            if log_scale:
                y_ticktext = [str(np.round(10*2**x,0)) for x in y_tickvals]
            else:
                y_ticktext = y_tickvals

            xaxis_title = xlabel
            yaxis_title = ylabel
        else:
            y_tickvals  = list(range(y_min,y_max))
            y_ticktext = ['']*len(y_tickvals)
            xaxis_title = ''
            yaxis_title = ''

        if invert_yticklabels:
            y_tickvals = y_tickvals[::-1]

        for ind2, K in enumerate(antigen_group_indices):


            visible='legendonly'

            line = dict(color=line_colors[ind2], width=line_widths[ind2], dash=line_styles[ind2])


            y_data = np.nanmean(np.reshape(data[np.ix_(K,I)], (len(K),len(I))),axis=0)


            fig.add_trace(plotly.graph_objects.Scatter(name=line_names[ind2], x=x_vals, y=y_data,
                                     legendgroup=f'group{ind2}',showlegend=showlegend,
                                     marker={'size':marker_sizes[ind2],'line':{'width':2, 'color':'DarkSlateGrey'}},
                                     visible=visible, mode='lines+markers', marker_symbol=markers[ind2], line=line,
                                     ), row=1, col=ind1+2)



        fig.update_layout(
            {
            f'yaxis{ind1+2}':{'range': [y_min, y_max], 'tickvals':y_tickvals, 'ticktext':y_ticktext,
                           'autorange':False, 'gridcolor':'#f4eeee','mirror':True, 'title':yaxis_title},

            f'xaxis{ind1+2}':{'range': [0.5, len(I)+0.5], 'tickvals':x_vals, 'ticktext':subplot_x_names,
                             'autorange':False, 'gridcolor':'#f4eeee', 'mirror':True,
                             'title':xaxis_title, 'gridwidth':1, 'tickfont':dict(size=19)
                             },


            },
            plot_bgcolor='rgba(0,0,0,0)',

        )

        if cover_thresholded:
            x0 = -1
            x1 = len(I)+1
            fig.add_trace(plotly.graph_objects.Scatter(x=[x0,x0,x1,x1], y=[y_min,np.log2(threshold/10),np.log2(threshold/10),y_min,y_min],
                                                    fillcolor='rgba(128,128,128,0.75)',
                                                    showlegend=False,
                                                    fill="toself",
                                                    mode="lines",
                                                    line=dict(color='rgba(128,128,128,0)')), row=1, col=ind1+2)


    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0
        )   )

    return fig


## BOOTSTRAP RELATED PLOTS ####################################################


def _vol_to_radius(vol, ndims):

  if ndims==2:
    return vol**(1/2)/np.pi
  elif ndims==3:
    return (vol/(4/3*np.pi))**(1/3)


def bootstrap_bar(vol, ndims, input_max_y=None, transformation=None):

  '''
  transformation transforms the vol values before plotting.
  default is calculating radius in dimensions in ndims assuming the value
  is the volume of a spherical object. if you want to plot for instance
  maximum width of blobs, then suggested transformation would be
  lambda x: x/2.
  '''

  assert ndims in [2,3], f"ndims is {ndims} but can only be 2 or 3"

  if transformation is None:
    transformation = lambda x: _vol_to_radius(x, ndims)

  if input_max_y is not None:
    max_y = input_max_y
  else:
    max_y = -np.inf

  max_len = max([len(x) for x in vol])
  xbuffer = max_len*4.2/54 + 0.5

  height = 5 + xbuffer
  nbars = len(vol)
  width = 4+(nbars+1)*0.33

  x0 = 2/width
  y0 = xbuffer/height
  dx = (nbars+1)*0.33/width
  dy = 4.8/height

  fig = plt.figure(figsize=(width, height))

  ax = fig.add_axes([x0, y0, dx, dy])

  for indn,name in enumerate(vol):
    sizes = transformation(np.array(vol[name]))
    sizes = [size for size in sizes if size!=0]

    if input_max_y is None:
       max_y = max(max_y, np.sum(sizes))

    for indv,size in enumerate(sizes):
      ax.bar(indn, size, bottom=np.sum(sizes[:indv]), facecolor=cm[indv])

  ax.set_xticks(range(len(vol)))
  ax.set_xticklabels(vol, rotation=90)
  ax.set_xlim([-1, len(vol)])

  ax.set_ylim([0,max_y])
  step = max(int(max_y/10),1)
  ax.set_yticks(np.arange(0,max_y+step,step))
  ax.set_yticklabels(np.arange(0,max_y+step,step), fontsize=12)

  ax.set_ylabel("radius", fontsize=20)

  ax.grid(True, alpha=0.3)

  return fig,ax


def bootstrap_agcolored_map(vol, racmap: pr.RacMap, max_val=None, min_val=None):

  ndims = racmap.number_of_dimensions
  assert ndims in [2,3], f"map dimensions is {ndims} but can only be 2 or 3"

  sizes = []

  for indn,name in enumerate(vol):
    total_size = np.sum(_vol_to_radius(np.array(vol[name]), ndims))
    sizes.append(total_size)

  sizes = np.array(sizes)
  if max_val is None:
    max_val = np.max(sizes)
  if min_val is None:
    min_val = np.min(sizes)

  sizes[sizes>max_val] = max_val
  sizes[sizes<min_val] = min_val

  normalized_sizes = (sizes-min_val)/(max_val-min_val)

  colors = normalized_sizes[:,None]*np.array([255,0,0])[None,:] +\
    (1-normalized_sizes)[:,None]*np.array([0,0,255])[None,:]

  colors = [pr.utils.colors.rgb_to_hex(*x) for x in colors]

  racmap.ag_fills = colors[:racmap.num_antigens]

  return racmap


def bootstrap_agvol_vs_mean_sr_dist(vol, racmap: pr.RacMap):
  '''
  scatter plot of bootstrap blob volumes/areas vs
  the mean distance of the ag to the nearest ndims+1 sera.
  '''

  ndims = racmap.number_of_dimensions
  assert ndims in [2,3], f"map dimensions is {ndims} but can only be 2 or 3"

  sizes = []
  mean_sr_dists = []

  sr_coords = racmap.sr_coordinates

  for indn,name in enumerate(vol):
    total_size = np.sum(_vol_to_radius(np.array(vol[name]), ndims))
    sizes.append(total_size)

  I = (np.char.find(racmap.titer_table.values.astype(str),'<')==-1) &\
      (np.char.find(racmap.titer_table.values.astype(str),'>')==-1) &\
      (np.char.find(racmap.titer_table.values.astype(str),'*')==-1)


  for ind in range(racmap.num_antigens):

    ag_coord = racmap.ag_coordinates[ind,:]

    J = np.argwhere(I[ind,:]).flatten()

    dif = ag_coord[None,:] - sr_coords[J,:]
    dists = np.linalg.norm(dif, axis=1)
    K = np.argsort(dists)[:ndims+1]

    mean_sr_dists.append(np.mean(dists[K]))

  fig, ax = plt.subplots(1,1, figsize=(5, 5))

  ax.scatter(sizes, mean_sr_dists, c="black")
  ax.set_xlabel("radius")
  ax.set_ylabel("mean sr dist")
  ax.grid('on', alpha=0.2)

  return fig,ax
