#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:29:36 2022

@author: Sina Tureli
"""

import PyRacmacs as pr
import tempfile
import webbrowser
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
