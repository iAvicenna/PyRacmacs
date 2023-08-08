#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 23:08:28 2022

@author: Sina Tureli
"""
import sys
import time
import rpy2.robjects as ro
import rpy2.rinterface_lib.callbacks
from .structures import Racmacs, RacMap

rNULL = ro.rinterface.NULL

class RProgressBar(list):

    def __init__(self, length=130):
        super(RProgressBar, self).__init__(['-' for _ in range(length)])
        self.stack = []

    def update(self,n):
        if '-' in self:
            i = self.index('-')
            self[i:i+n] = ['=' for _ in range(n)]

    def is_complete(self):

        return all(x=='=' for x in self)

    def reset(self):
        super(RProgressBar, self).__init__(['-' for _ in range(len(self))])
        self.stack = []

    def print(self):

        sys.stdout.write('\r'+''.join(self)+' ')
        sys.stdout.flush()


def capture_r_output(display_progress: bool, is_optimization: bool):
    """
    Will cause all the output that normally goes to the R console,
    to end up instead in a python list.

    racmacs prints progress bar for operations like bootstrap, dimensionality
    testing etc differently than optimization. so is_optimization (True/False)
    is used to tailor the function to these needs

    modified from a script by the author: xApple
    link: https://stackoverflow.com/questions/38793688/how-do-i-silence-output-to-stdout-in-rpy2

    """

    stdout = []
    stderr = []

    if not is_optimization:
      progress_bar = RProgressBar(80)
    else:
      progress_bar = RProgressBar(0)

    def add_to_stdout(line):

        if display_progress:
            print(line)
        stdout.append(line)

    def add_to_stderr_nonopt(line):

        if display_progress and line!='':

            if all(x in ['-'] for x in line)\
              and line.count('-') != len(progress_bar):

              progress_bar.__init__(line.count('-'))

            if not all(x in ['\n','\r','-','=',''] for x in line):
                if progress_bar is not None and\
                  progress_bar.is_complete():

                    print('')
                    progress_bar.reset()

                print(line, end='')

            elif all(x in ['=','-'] for x in line.strip('\n')) and\
                     len(line.strip('\n'))>0:

                progress_bar\
                  .update(n=line.count('=')-
                          progress_bar.count('='))


                progress_bar.print()


        stderr.append(line)

    def add_to_stderr_opt(line):

        if display_progress:

            if not all(x in ['\n','\r','-','='] for x in line):
              end = ''
              if line[-1]!='\n': end = '\n'
              print(line, end=end)
            elif line == '-':
              progress_bar.__init__(len(progress_bar)+1)
            elif line=='=':
              progress_bar.stack.append('=')
            elif (line=='\r' or
                  (len(progress_bar.stack)>0 and line not in ['=','\r'])):

              progress_bar.update(len(progress_bar.stack))
              progress_bar.print()

              if not all(x in ['\n','\r','-','='] for x in line):
                end = ''
                if line[-1]!='\n': end = '\n'
                print(line, end=end)

              progress_bar.reset()


        stderr.append(line)

    # Set the call backs #
    rpy2.rinterface_lib.callbacks.consolewrite_print = add_to_stdout

    if is_optimization:
      rpy2.rinterface_lib.callbacks.consolewrite_warnerror = add_to_stderr_opt
    else:
      rpy2.rinterface_lib.callbacks.consolewrite_warnerror = add_to_stderr_nonopt

    return stdout,stderr


def write_racmap(racmap: RacMap, save_path: str, compress=False,
                 round_titers=False, pretty=None):

    if pretty is None:
        pretty = not compress

    Racmacs.save_acmap(racmap._acmap_R, save_path, compress, pretty, round_titers)


def read_racmap(read_path: str, optimization_number=None, sort_optimizations=False,
                align_optimizations=False):

    if optimization_number is None:
        optimization_number = rNULL
    else:
        optimization_number += 1 # python starts from 0 R from 1

    acmap_R = Racmacs.read_acmap(read_path, optimization_number, sort_optimizations,
                                 align_optimizations)

    racmap = RacMap(acmap_R)

    return racmap


class PyProgressBar(list):

    def __init__(self, total_iterations, length=100, step=0.01):
        super(PyProgressBar, self).__init__(['-' for _ in range(length)])
        self.step = step
        self.t0 = time.time()
        self.proportion_complete = 0
        self.total_iterations = total_iterations
        self.has_completed = False


    @property
    def is_complete(self):

        return '-' not in self

    def grow(self, n=None):

        if not self.has_completed:
            previous_proportion_complete = self.proportion_complete

            if n is None: n=1

            self.proportion_complete += n/self.total_iterations

            if '-' in self:
                n = int(self.proportion_complete*len(self))
                self[:n+1] = ['=' for _ in range(n+1)]

            if int(previous_proportion_complete/self.step)<int(self.proportion_complete/self.step):
                self.print()


    def reset(self):
        super(PyProgressBar, self).__init__(['-' for _ in range(len(self))])

    def print(self):

        sys.stdout.write('\r'+''.join(self))
        sys.stdout.flush()
        if self.is_complete:
            print(f'\n{time.time()-self.t0:.2f} seconds.')
            self.has_completed = True



# type_check(expected)
