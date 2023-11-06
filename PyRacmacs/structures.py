#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 22:10:09 2022

@author: Sina Tureli
"""

import rpy2
import rpy2.robjects as ro
import numpy as np
import pandas as pd
import warnings
import multiprocessing
from . import Racmacs

from .utils import conversion, regulations
from rpy2.robjects import pandas2ri
from rpy2.robjects import r

rNULL = ro.rinterface.NULL

STYLE_KEY_DICT = {
        'shown': '+',
        'aspect': 'a',
        'fill': 'F',
        'outline': 'O',
        'shape': 'S',
        'size': 's',
        'outline_width' : 'o',
        'rotation': 'r'
        }

OTHER_KEY_DICT = {
        'name' : 'N',
        'passage' : 'P',
        'group' : 'g',
        'id' : 'i',
        'sequence' : 'q'
    }

MAP_PROPS = ['pt_drawing_order']

STYLE_INV_DICT = {STYLE_KEY_DICT[key]:key for key in STYLE_KEY_DICT}
OTHER_INV_DICT = {OTHER_KEY_DICT[key]:key for key in OTHER_KEY_DICT}

check = regulations.check


class RacOptimizerOptions():

    def __init__(self, dim_annealing=False, method="L-BFGS",
                 maxit=3000, num_cores=max(1,multiprocessing.cpu_count() - 1),
                 report_progress=None, ignore_disconnected=False,
                 progress_bar_length=139
                 ):

      if report_progress is None:
          report_progress = rNULL

      self.dim_annealing = dim_annealing
      self.method = method
      self.maxit = maxit
      self.num_cores = num_cores
      self.report_progress = report_progress
      self.ignore_disconnected = ignore_disconnected
      self.progress_bar_length = progress_bar_length

    @property
    def options_R(self):
      return\
        Racmacs.RacOptimizer_options(dim_annealing=self.dim_annealing,
                                     method=self.method, maxit=self.maxit,
                                     num_cores=self.num_cores,
                                     report_progress=self.report_progress,
                                     ignore_disconnected=self.ignore_disconnected,
                                     progress_bar_length=self.progress_bar_length
                                     )


class RacMergeOptions():

    def __init__(self, sd_limit=1,
                 dilution_stepsize=1, method=None):

      "method should be conservative or likelihood"

      if method is None:
        method = rNULL



      self.dilution_stepsize = dilution_stepsize
      self.sd_limit = sd_limit
      self.method = method

    @property
    def options_R(self):
      _options_R=\
          Racmacs.RacMerge_options(sd_limit=self.sd_limit,
                                   dilution_stepsize=self.dilution_stepsize,
                                   method=self.method
                                   )

      #second option which is sth called function seems redundant and
      #probably is a bug. when that is supplied it complains that function
      #is an unused argument

      names = [_options_R.names[i] for i in [0,1,3]]
      vals =   [_options_R[i] for i in [0,1,3]]

      return ro.vectors.ListVector(zip(names, vals))



class RacViewerOptions():

    def __init__(self, opacity=None, viewer_controls="hidden",
                 grid_display="static", grid_col="#cfcfcf",
                 background_col = "#ffffff",
                 show_names = False, show_errorlines=False,
                 show_connectionlines=False, show_titers=False,
                 xlim=None, ylim=None, translation=None, rotation=None,
                 zoom=None):

        if xlim is None:
          xlim = rNULL
        if ylim is None:
          ylim = rNULL
        if opacity is None:
          opacity = ro.NA_Logical
        if zoom is None:
          zoom = rNULL
        if translation is None:
          translation = ro.FloatVector([0, 0 , 0])
        if rotation is None:
          rotation = ro.FloatVector([0, 0 ,0])


        self.viewer_controls = viewer_controls
        self.grid_display = grid_display
        self.grid_col = grid_col
        self.show_names = show_names
        self.show_errorlines = show_errorlines
        self.show_connectionlines = show_connectionlines
        self.show_titers = show_titers
        self.background_col = background_col
        self.xlim = xlim
        self.ylim = ylim
        self.opacity = opacity
        self.zoom = zoom
        self.translation = translation
        self.rotation = rotation


    @property
    def options_R(self):
        return\
            Racmacs.RacViewer_options(self.opacity, self.viewer_controls,
                                      self.grid_display, self.grid_col,
                                      self.background_col, self.show_names,
                                      self.show_errorlines,
                                      self.show_connectionlines,
                                      self.show_titers, self.xlim, self.ylim,
                                      self.translation, self.rotation, self.zoom)


class RacMap():
    '''
    The main attribute of this class is the _acmap_R object. This
    is a an object from rpy2.robject and therefore can directly
    be passed to functions imported from R via rpy2. All the other
    properties listed here allow accessing and modifying attributes
    of this object via the imported Racmacs package functionalities.
    '''

    def from_ace(read_path, optimization_number=None, sort_optimizations=True,
                 align_optimizations=True):

      from .io import read_racmap

      return read_racmap(read_path, optimization_number, sort_optimizations,
                         align_optimizations)


    def __init__(self, acmap_R=None, titer_table=None):

        if acmap_R is not None:
            self._acmap_R = acmap_R
        elif titer_table is not None:
            with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
                titer_table_r = ro.conversion.py2rpy(titer_table)
            self._acmap_R = Racmacs.acmap(titer_table=titer_table_r,
                                          ag_names = ro.StrVector(titer_table.index),
                                          sr_names = ro.StrVector(titer_table.columns))


    def get_bootstrap_data(self, optimization_number=0):
        return self.get_optimization(optimization_number)['bootstrap']

    def get_blobs(self):
        return _extract_blobs_from_R(self, "bootstrap")

    def get_triangulation_blobs(self, optimization_number=0):


        ag_tblobs, sr_tblobs = _extract_blobs_from_R(self, "triangulation",
                                                     optimization_number+1)

        return ag_tblobs, sr_tblobs

    def get_optimization(self, optimization_number):
        return conversion.recurse_r_tree(Racmacs.getOptimization(self._acmap_R, optimization_number+1))


    def get_optimization_coordinates(self, optimization_number):

        optimization = self.get_optimization(optimization_number)
        ag_coords = optimization['ag_base_coords']
        sr_coords = optimization['sr_base_coords']

        return np.concatenate([ag_coords, sr_coords])


    def get_red_chi2(self, std=1, optimization_number=0):
        return self.get_stress(optimization_number)/(std**2*self.df)


    def get_stress(self, optimization_number=0):
        return self.get_optimization(optimization_number)['stress'][0]


    def keep_optimizations(self, optimization_numbers):

        if isinstance(optimization_numbers,int):
            optimization_numbers = [optimization_numbers+1]
        else:
            optimization_numbers = [x+1 for x in optimization_numbers]

        new_map_R =  Racmacs.keepOptimizations(self._acmap_R,
                                               ro.IntVector(optimization_numbers))

        return RacMap(new_map_R)


    def order_antigens(self, ag_names):

        new_map_R = Racmacs.orderAntigens(self._acmap_R, ag_names)

        return RacMap(new_map_R)


    def order_sera(self, sr_names):

        new_map_R = Racmacs.orderSera(self._acmap_R, sr_names)

        return RacMap(new_map_R)


    def add_optimization(self, ag_coordinates, sr_coordinates, number_of_dimensions,
                         minimum_column_basis='none', fixed_column_bases = None,
                         ag_reactivity_adjustments = None
                         ):

        if sr_coordinates is None:
            sr_coordinates = np.nan*np.zeros((self.num_sera,
                                              number_of_dimensions))

        ag_coordinates_r = ro.numpy2ri.numpy2rpy(ag_coordinates)
        sr_coordinates_r = ro.numpy2ri.numpy2rpy(sr_coordinates)

        if fixed_column_bases is None:
            fixed_column_bases = rNULL
        else:
            fixed_column_bases = ro.FloatVector(fixed_column_bases)

        if ag_reactivity_adjustments is None:
            ag_reactivity_adjustments = rNULL
        else:
            ag_reactivity_adjustments = ro.FloatVector(ag_reactivity_adjustments)


        _acmap_R = Racmacs.addOptimization(self._acmap_R,  ag_coordinates_r,
                                            sr_coordinates_r, number_of_dimensions,
                                            minimum_column_basis, fixed_column_bases,
                                            ag_reactivity_adjustments)
        new_map = RacMap(_acmap_R)

        return new_map


    def recalculate_stress(self, optimization=0):

        return float(Racmacs.recalculateStress(self._acmap_R, optimization+1)[0])


    @property
    def df(self):
        Npoints = self.num_antigens + self.num_sera
        Nobs = self.number_of_measured_titers
        d = self.number_of_dimensions

        return Nobs - Npoints*d +  d + d*(d-1)/2

    @property
    def T(self, styles=None):
        new_map = RacMap(titer_table = self.titer_table.T)

        if styles is None:
            styles = ['outlines','fills','shapes','sizes', 'groups']

        for style in styles:
            new_map.__setattr__(f'sr_{style}',
                             self.__getattribute__(f'ag_{style}'))

            new_map.__setattr__(f'ag_{style}',
                             self.__getattribute__(f'sr_{style}'))

        return new_map

    @property
    def number_of_dimensions(self):
        return int(Racmacs.mapDimensions(self._acmap_R)[0])

    @property
    def min_col_basis(self):
        return

    @property
    def num_sera(self):
        return int(Racmacs.numSera(self._acmap_R)[0])

    @property
    def num_antigens(self):
        return int(Racmacs.numAntigens(self._acmap_R)[0])

    @property
    def num_optimizations(self):
        return int(Racmacs.numOptimizations(self._acmap_R)[0])

    @property
    def number_of_measured_titers(self):
        return len(self.measured_titer_locations[0])

    @property
    def number_of_well_coordinated_titers(self):
        return len(self.well_coordinated_titer_locations[0])

    @property
    def titer_table(self):
        titer_table_R = Racmacs.titerTable(self._acmap_R)

        row_names = conversion.return_list(titer_table_R.dimnames[0])
        col_names = conversion.return_list(titer_table_R.dimnames[1])
        titer_table = conversion.return_list(titer_table_R)


        s1 = len(row_names)
        s2 = len(col_names)

        titer_table = np.reshape(np.array(titer_table),(s2,s1)).T
        titer_table = pd.DataFrame(titer_table,index=row_names, columns=col_names)


        return titer_table

    @property
    def log_titer_table(self):

        num_titers = self.titer_table.applymap(conversion.titer_str_to_num).values
        log_titers = np.log2(num_titers/10)

        log_titer_table = pd.DataFrame(log_titers, index = self.ag_names,
                                       columns=self.sr_names)

        return log_titer_table

    @property
    def log_titers(self):

        return self.log_titer_table.values


    @property
    def below_LOD_titer_locations(self):
        '''
        indices of titers with < (LOD:limit of detection)
        '''
        indices = np.argwhere(np.char.find(self.titer_table.values.astype(str),'<')>-1)

        return tuple(np.array(indices).T.tolist())

    @property
    def above_LOD_titer_locations(self):
        '''
        indices of titers with > (LOD:limit of detection)
        '''
        indices = np.argwhere(np.char.find(self.titer_table.values.astype(str),'>')>-1)

        return tuple(np.array(indices).T.tolist())

    @property
    def detectable_titer_locations(self):
        '''
        titers which do not have < or > in them
        '''

        indices = np.argwhere((np.char.find(self.titer_table.values.astype(str),'<')==-1) &
                              (np.char.find(self.titer_table.values.astype(str),'>')==-1)
                              )

        return tuple(np.array(indices).T.tolist())

    @property
    def undetectable_titer_locations(self):
        '''
        titers which do have < or > in them
        '''

        indices = np.argwhere((np.char.find(self.titer_table.values.astype(str),'<')>-1) |
                              (np.char.find(self.titer_table.values.astype(str),'>')>-1)
                              )

        return tuple(np.array(indices).T.tolist())


    @property
    def measured_titer_locations(self):
        '''
        indices of titers which are measured i.e not *
        '''

        indices = np.argwhere(np.char.find(self.titer_table.values.astype(str),'*')==-1)

        return tuple(np.array(indices).T.tolist())

    @property
    def well_coordinated_serum_indices(self):

        '''
        sera for which number of measured titers > map dimensions
        '''

        s1 = self.shape[1]
        coordinated_serum_indices = []
        ndims = self.number_of_dimensions
        titers = self.titer_table.values.astype(str)

        for i in range(s1):
            if np.count_nonzero((np.char.find(titers[:,i],'*')==-1) &
                                (np.char.find(titers[:,i],'<')==-1)
                                )>ndims:
                coordinated_serum_indices.append(i)

        return coordinated_serum_indices


    @property
    def well_coordinated_antigen_indices(self):

        '''
        sera for which numer of measured titers > map dimensions
        '''

        s1 = self.shape[0]
        coordinated_antigen_indices = []
        ndims = self.number_of_dimensions
        titers = self.titer_table.values.astype(str)

        for i in range(s1):
            if np.count_nonzero((np.char.find(titers[i,:],'*')==-1) &
                                (np.char.find(titers[i,:],'<')==-1)
                                )>ndims:
                coordinated_antigen_indices.append(i)

        return coordinated_antigen_indices


    @property
    def well_coordinated_titer_locations(self):
        '''
        indices of titers for which both the antigen and serum have
        > map dimension number of measured titers
        '''

        well_coordinated_antigen_indices = self.well_coordinated_antigen_indices
        well_coordinated_serum_indices = self.well_coordinated_serum_indices

        measured_titer_locations = self.measured_titer_locations
        well_coordinated_titer_locations = [[],[]]

        for titer_location0,titer_location1 in (zip(measured_titer_locations[0],
                                                      measured_titer_locations[1])):

            if (titer_location0 in well_coordinated_antigen_indices and
                titer_location1 in well_coordinated_serum_indices):

                well_coordinated_titer_locations[0].append(titer_location0)
                well_coordinated_titer_locations[1].append(titer_location1)

        return well_coordinated_titer_locations


    @property
    def unmeasured_titer_locations(self):
        '''
        indices of unmeasured titers i.e those which are *
        '''

        indices = np.argwhere(np.char.find(self.titer_table.values.astype(str),'*')>-1)

        return tuple(np.array(indices).T.tolist())


    @property
    def is_measured(self):
      '''
      a dataframe reflecting whether a titer has * in or not
      '''

      return pd.DataFrame(np.core.defchararray.find(self.titer_table.values.astype(str),'*')!=-1,
                          index=self.ag_names, columns=self.sr_names)


    @property
    def is_below_LD(self):
        '''
        a dataframe reflecting whether a titer has < in or not
        '''

        return pd.DataFrame(np.core.defchararray.find(self.titer_table.values.astype(str),'<')!=-1,
                            index=self.ag_names, columns=self.sr_names)


    @property
    def lower_LD(self):
        '''
        return lower limit of detection by looking at values with < and returning
        minimum amongst them. if none found returns -np.inf
        '''

        undetected_titer_locations = self.below_LOD_titer_locations

        undetected_titers = list(set([int(x[1:]) for x in
                                      self.titer_table.values[undetected_titer_locations].flatten()]))

        if len(undetected_titers) == 0:
            return -np.inf
        elif len(undetected_titers)>1:
            warnings.warn(f'Multiple lower LOD titers: {sorted(undetected_titers)}, returning maximum.')

        return max(undetected_titers)


    @property
    def is_above_LD(self):
        '''
        a dataframe reflecting whether a titer has > in or not
        '''

        return pd.DataFrame(np.core.defchararray.find(self.titer_table.values.astype(str),'>')!=-1,
                            index=self.ag_names, columns=self.sr_names)


    @property
    def upper_LD(self):
        '''
        return lower limit of detection by looking at values with > and returning
        maximum amongst them. if none found returns np.inf
        '''

        undetected_titer_locations = self.above_LOD_titer_locations

        undetected_titers = list(set([int(x[1:]) for x in
                                      self.titer_table.values[undetected_titer_locations].flatten()]))

        if len(undetected_titers) == 0:
            return np.inf
        elif len(undetected_titers)>1:
            warnings.warn(f'Multiple upper LOD titers: {sorted(undetected_titers)}, returning maximum.')

        return max(undetected_titers)


    @property
    def sr_group_levels(self):
        if isinstance(Racmacs.srGroups(self._acmap_R),rpy2.rinterface_lib.sexp.NULLType):
            return []
        else:
            sr_group_levels = conversion.return_list(Racmacs.srGroups(self._acmap_R).levels)
            sr_groups = self.sr_groups
            I = np.argsort([sr_groups.index(x) for x in sr_group_levels])
            return [sr_group_levels[i] for i in I]

    @property
    def ag_group_levels(self):
        if isinstance(Racmacs.agGroups(self._acmap_R),rpy2.rinterface_lib.sexp.NULLType):
            return []
        else:
            ag_group_levels = conversion.return_list(Racmacs.agGroups(self._acmap_R).levels)
            ag_groups = self.ag_groups
            I = np.argsort([ag_groups.index(x) for x in ag_group_levels])
            return [ag_group_levels[i] for i in I]


    @property
    def shape(self):
        return (self.num_antigens,self.num_sera)


    @property
    def ag_group_names(self):

        ag_names = self.ag_names
        ag_groups = self.ag_groups

        return [f'{x}_{y}' for x,y in zip(ag_names, ag_groups)]

    @property
    def sr_group_names(self):

        sr_names = self.sr_names
        sr_groups = self.sr_groups

        return [f'{x}_{y}' for x,y in zip(sr_names, sr_groups)]

    @property
    def column_bases(self):
        return conversion.return_list(Racmacs.colBases(self._acmap_R))


    def map_distances(self, optimization_number=0):
        return np.array(Racmacs.mapDistances(self._acmap_R, optimization_number+1))


    def predicted_log_titers(self, optimization_number=0):

        map_distances = self.map_distances(optimization_number=optimization_number+1)
        s1,s2 = map_distances.shape
        column_bases = np.tile(np.reshape(self.column_bases,(s2,1)),(1,s1)).T

        return column_bases - map_distances


    def predicted_titers(self, optimization_number=0):

        predicted_log_titers = self.predicted_log_titers(optimization_number=optimization_number+1).flatten()

        predicted_titers = np.array([str(10*2**x) if not np.isnan(x) else str(np.nan)
                                for x in predicted_log_titers])

        predicted_titers = np.reshape(predicted_titers,self.titer_table.shape)

        return predicted_titers


    def sr_indices_of_sr_group(self, sr_group):

        sr_group_ids = self.sr_group_ids
        check(1,self,sr_group)
        return [ind for ind,x in enumerate(sr_group_ids) if x == sr_group]


    def subset_map(self, serum_indices=None, antigen_indices=None):

        if serum_indices is not None:
            serum_indices = [x+1 for x in serum_indices]
        else:
            serum_indices = list(range(1,self.num_sera+1))

        if antigen_indices is not None:
            antigen_indices = [x+1 for x in antigen_indices]
        else:
            antigen_indices =  list(range(1,self.num_antigens+1))

        antigen_indices = sorted(antigen_indices)
        serum_indices = sorted(serum_indices)

        antigen_indices = ro.IntVector(antigen_indices)
        serum_indices = ro.IntVector(serum_indices)
        subset_acmap_R = Racmacs.subsetMap(self._acmap_R, sera=serum_indices,
                                           antigens=antigen_indices)

        return RacMap(subset_acmap_R)


    def subset_map_bygroup(self, serum_groups=None, antigen_groups=None):

        if serum_groups is None:
            serum_indices =list(range(1,self.num_antigens+1))
        else:
            serum_indices = [ind+1 for ind,name in enumerate(self.sr_groups)
                             if name in serum_groups]

        if antigen_groups is None:
            antigen_indices =list(range(1,self.num_antigens+1))
        else:
            antigen_indices = [ind+1 for ind,name in enumerate(self.ag_groups)
                               if name in antigen_groups]

        antigen_indices = ro.IntVector(antigen_indices)
        serum_indices = ro.IntVector(serum_indices)

        subset_acmap_R = Racmacs.subsetMap(self._acmap_R, sera=serum_indices,
                                           antigens=antigen_indices)

        return RacMap(subset_acmap_R)


    def _get_sr_names(self):
        return conversion.return_list(Racmacs.srNames(self._acmap_R))

    def _get_ag_names(self):
        return conversion.return_list(Racmacs.agNames(self._acmap_R))


    def _get_sr_ids(self):
        return conversion.return_list(Racmacs.srIDs(self._acmap_R))


    def _get_ag_ids(self):
        return conversion.return_list(Racmacs.agIDs(self._acmap_R))


    def _get_ag_coordinates(self, optimization_number=0):
        return np.array(Racmacs.agCoords(self._acmap_R, optimization_number+1))


    def _get_sr_coordinates(self, optimization_number=0):
        return np.array(Racmacs.srCoords(self._acmap_R, optimization_number+1))


    def _get_coordinates(self, optimization_number=0):

        ag_coordinates = self._get_ag_coordinates()
        sr_coordinates = self._get_sr_coordinates()

        return np.concatenate([ag_coordinates, sr_coordinates])


    def _get_map_transformation(self, optimization_number=0):

        return np.array(Racmacs.mapTransformation(self._acmap_R,
                                                  optimization_number+1))


    def _get_sr_groups(self):

        sr_groups_factor = Racmacs.srGroups(self._acmap_R)

        return [sr_groups_factor.levels[x-1] for x in sr_groups_factor]


    def _get_ag_groups(self):

        ag_groups_factor = Racmacs.agGroups(self._acmap_R)

        return [ag_groups_factor.levels[x-1] for x in ag_groups_factor]


    def _get_dilution_stepsize(self):
        return Racmacs.dilutionStepsize(self._acmap_R)[0]


    def _get_sr_outlines(self):
        return conversion.return_list(Racmacs.srOutline(self._acmap_R))


    def _get_ag_outlines(self):
        return conversion.return_list(Racmacs.agOutline(self._acmap_R))


    def _get_sr_shapes(self):
        return conversion.return_list(Racmacs.srShape(self._acmap_R))


    def _get_ag_shapes(self):
        return conversion.return_list(Racmacs.agShape(self._acmap_R))


    def _get_sr_fills(self):
        return conversion.return_list(Racmacs.srFill(self._acmap_R))


    def _get_ag_fills(self):
        return conversion.return_list(Racmacs.agFill(self._acmap_R))


    def _get_sr_sizes(self):
        return conversion.return_list(Racmacs.srSize(self._acmap_R))


    def _get_ag_sizes(self):
        return conversion.return_list(Racmacs.agSize(self._acmap_R))


    def _get_ag_reactivity_adjustments(self):
        return conversion.return_list(Racmacs.agReactivityAdjustments(self._acmap_R))


    def _get_min_column_basis(self):
        return str(Racmacs.minColBasis(self._acmap_R))


    def _get_ag_sequences(self):
      seqs = np.array(Racmacs.agSequences(self._acmap_R))
      size = seqs.size
      nrows = self.num_antigens
      ncols = int(size/nrows)

      return np.reshape(seqs, (ncols, nrows))

    def _get_sr_sequences(self):
      seqs = np.array(Racmacs.srSequences(self._acmap_R))
      size = seqs.size
      nrows = self.num_sera
      ncols = int(size/nrows)

      return np.reshape(seqs, (nrows, ncols))

    def _set_sr_names(self, val):
        set_method = r("`srNames<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_names(self, val):
        set_method = r("`agNames<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_ids(self, val):
        set_method = r("`agIDs<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_sr_ids(self, val):
        set_method = r("`srIDs<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_sr_groups(self, val):
        set_method = r("`srGroups<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_groups(self, val):
        set_method = r("`agGroups<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_dilution_stepsize(self, val):
        set_method = r("`dilutionStepsize<-`")
        self._acmap_R = set_method(self._acmap_R, float(val))


    def _set_ag_reactivity_adjustments(self, val):
        set_method = r("`agReactivityAdjustments<-`")
        self._acmap_R = set_method(self._acmap_R, ro.FloatVector(val))


    def _set_sr_outlines(self, val):
        set_method = r("`srOutline<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_outlines(self, val):
        set_method = r("`agOutline<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_sr_shapes(self, val):
        set_method = r("`srShape<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_shapes(self, val):
        set_method = r("`agShape<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_sr_fills(self, val):
        set_method = r("`srFill<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_ag_fills(self, val):
        set_method = r("`agFill<-`")
        self._acmap_R = set_method(self._acmap_R, ro.StrVector(val))


    def _set_sr_sizes(self, val):
        set_method = r("`srSize<-`")
        self._acmap_R = set_method(self._acmap_R, ro.FloatVector(val))


    def _set_ag_sizes(self, val):
        set_method = r("`agSize<-`")
        self._acmap_R = set_method(self._acmap_R, ro.FloatVector(val))


    def _set_ag_coordinates(self, val, optimization_number=0):
        set_method = r("`agCoords<-`")

        ag_coordinates_r = ro.numpy2ri.numpy2rpy(val)
        self._acmap_R = set_method(self._acmap_R, optimization_number+1,
                                   ag_coordinates_r)


    def _set_sr_coordinates(self, val, optimization_number=0):
        set_method = r("`srCoords<-`")

        sr_coordinates_r = ro.numpy2ri.numpy2rpy(val)
        self._acmap_R = set_method(self._acmap_R, optimization_number+1,
                                   sr_coordinates_r)

    def _set_min_column_basis(self, val):

        set_method = r("`minColBasis<-`")
        self._acmap_R = set_method(self._acmap_R, val)


    def _set_coordinates(self, val, optimization_number=0):

        ag_coordinates = val[:self.num_antigens,:]
        sr_coordinates = val[self.num_antigens:self.num_antigens+self.num_sera,:]

        self._set_ag_coordinates(ag_coordinates, optimization_number)
        self._set_sr_coordinates(sr_coordinates, optimization_number)


    def _set_map_transformation(self, val, optimization_number=0):

        set_method = r("`mapTransformation<-`")
        transformation = ro.numpy2ri.numpy2rpy(val)
        self._acmap_R = set_method(self._acmap_R, optimization_number+1,
                                   transformation)

    def _set_ag_sequences(self, sequences):

        set_method = r("`agSequences<-`")
        sequences = ro.numpy2ri.numpy2rpy(sequences)
        self._acmap_R = set_method(self._acmap_R, sequences)


    def _set_sr_sequences(self, sequences):

        set_method = r("`srSequences<-`")

        sequences = ro.numpy2ri.numpy2rpy(sequences)

        self._acmap_R = set_method(self._acmap_R, sequences)


    ag_sequences = property(_get_ag_sequences, _set_ag_sequences)
    sr_sequences = property(_get_sr_sequences, _set_sr_sequences)
    ag_names = property(_get_ag_names, _set_ag_names)
    sr_names = property(_get_sr_names, _set_sr_names)
    dilution_stepsize = property(_get_dilution_stepsize, _set_dilution_stepsize)
    ag_ids = property(_get_ag_ids, _set_ag_ids)
    sr_ids = property(_get_sr_ids, _set_sr_ids)
    sr_groups = property(_get_sr_groups, _set_sr_groups)
    ag_groups = property(_get_ag_groups, _set_ag_groups)
    dilution_stepsize = property(_get_dilution_stepsize, _set_dilution_stepsize)
    ag_reactivity_adjustments = property(_get_ag_reactivity_adjustments, _set_ag_reactivity_adjustments)
    sr_outlines = property(_get_sr_outlines, _set_sr_outlines)
    ag_outlines = property(_get_ag_outlines, _set_ag_outlines)
    sr_shapes = property(_get_sr_shapes, _set_sr_shapes)
    ag_shapes = property(_get_ag_shapes, _set_ag_shapes)
    sr_fills = property(_get_sr_fills, _set_sr_fills)
    ag_fills = property(_get_ag_fills, _set_ag_fills)
    sr_sizes = property(_get_sr_sizes, _set_sr_sizes)
    ag_sizes = property(_get_ag_sizes, _set_ag_sizes)
    coordinates = property(_get_coordinates, _set_coordinates)
    ag_coordinates = property(_get_ag_coordinates, _set_ag_coordinates)
    sr_coordinates = property(_get_sr_coordinates, _set_sr_coordinates)
    map_transformation = property(_get_map_transformation, _set_map_transformation)



def _extract_blobs_from_R(map_with_blobs, blob_type, optimization_number=0):

    ag_names = map_with_blobs.ag_names
    sr_names = map_with_blobs.sr_names

    ag_blobs = {}
    sr_blobs = {}


    if map_with_blobs.number_of_dimensions==3:

      keys = ["vertices", "faces", "normals"]
      dtypes = [float, int, float]

      for ag in ag_names:
        ag_blobs[ag] = {}

        if blob_type == 'bootstrap':
          ag_blob = Racmacs.agBootstrapBlob(map_with_blobs._acmap_R, ag)
        elif blob_type == 'triangulation':
          ag_blob = Racmacs.agTriangulationBlob(map_with_blobs._acmap_R, ag,
                                                optimization_number+1)
        else:
          raise ValueError("blob_type can only be bootstrap or triangulation "
                           "but it was {blob_type}")

        for i in range(3):
          if ag_blob == rNULL or len(ag_blob)==0:
            ag_blobs[ag][keys[i]] = np.array([[]])
            continue

          ag_blobs[ag][keys[i]] = np.array(Racmacs.agBootstrapBlob(map_with_blobs._acmap_R,ag)[0][i],
                                           dtype=dtypes[i])

      for sr in sr_names:
        sr_blobs[sr] = {}

        if blob_type == 'bootstrap':
          sr_blob = Racmacs.srBootstrapBlob(map_with_blobs._acmap_R, sr)
        elif blob_type == 'triangulation':
          sr_blob = Racmacs.srTriangulationBlob(map_with_blobs._acmap_R, sr,
                                                optimization_number+1)
        else:
          raise ValueError("blob_type can only be bootstrap or triangulation "
                           "but it was {blob_type}")

        for i in range(3):

          if sr_blob == rNULL or len(sr_blob)==0:
            sr_blobs[sr][keys[i]] = np.array([[]])
            continue

          sr_blobs[sr][keys[i]] = np.array(Racmacs.srBootstrapBlob(map_with_blobs._acmap_R,sr)[0][i],
                                           dtype=dtypes[i])

    else:

      for ag in ag_names:
        ag_blobs[ag] = []

        if blob_type == "bootstrap":
          ag_blob_list = Racmacs.agBootstrapBlob(map_with_blobs._acmap_R,ag)
        elif blob_type == "triangulation":
          ag_blob_list = Racmacs.agTriangulationBlob(map_with_blobs._acmap_R,ag)

        if ag_blob_list==rNULL or len(ag_blob_list)==0:
          continue

        for ag_blob in ag_blob_list:
          x = np.array(list(ag_blob.items())[1][1])
          y = np.array(list(ag_blob.items())[2][1])

          ag_blobs[ag].append(np.array([x,y,np.zeros(x.shape)]).T)


      for sr in sr_names:
        sr_blobs[sr] = []

        if blob_type == "bootstrap":
          sr_blob_list = Racmacs.srBootstrapBlob(map_with_blobs._acmap_R, sr)
        elif blob_type == "triangulation":
          sr_blob_list = Racmacs.srTriangulationBlob(map_with_blobs._acmap_R, sr)

        if sr_blob_list==rNULL or len(sr_blob_list)==0:
          continue

        for sr_blob in sr_blob_list:
          x = np.array(list(sr_blob.items())[1][1])
          y = np.array(list(sr_blob.items())[2][1])

          sr_blobs[sr].append(np.array([x,y,np.zeros(x.shape)]).T)

    return ag_blobs, sr_blobs
