#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    cleaning.py
#
#    Utility functions to clean up unnecessary files after
#    running a protocol.
#
#    Copyright (C) 2020 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.


# standard library
import os
# RosettaDDGPrediction
from . import util



def clean_folders(process, stepname, wd, options, level):
    """Clean up folders after having run a protocol step.

    NB: 'process' is passed to ensure it is finished before
    the cleaning starts (because of the Dask task graph).
    """

    # if the step is 'relax' (cartddg protocols)
    if stepname == "relax":
        cleanfunc = clean_relax
    # if the step is 'cartesian' (cartddg protocols)
    elif stepname == "cartesian":
        cleanfunc = clean_cartesian
    # if the step is 'flexddg' (flexddg protocols)
    elif stepname == "flexddg":
        cleanfunc = clean_flexddg
    # call the appropriate cleaning function
    cleanfunc(wd = wd, \
              options = options, \
              level = level)


def clean_relax(wd, options, level):
    """Clean the folder where the 'relax' step was run."""

    # no cleaning procedure has been implemented so far
    if level is not None:
        errstr = "No cleaning procedure for 'relax' has " \
                 "been implemented yet."
        raise NotImplementedError(errstr)    


def clean_cartesian(wd, options, level):
    """Clean the folders where the 'cartesian' step was run."""

    # simply return if no cleaning was requested
    if level is None:
        return

    # if ony removal of PDB files has been requested
    elif level == "pdb":
        for item in os.listdir(wd):
            # remove all PDB files
            if item.endswith(".pdb"):
                os.remove(item)            


def clean_flexddg(wd, options, level):
    """Clean the folders where the 'flexddg' step was run."""

    # simply return if no cleaning was requested
    if level is None:
        return
    
    # if some cleaning has been requested
    elif level in ["structdbfile", "all"]:

        # get the RosettaScript options 
        rscriptoptions = \
            options[util.get_option_key(options = options, \
                                        option = "scriptvars")]
        # get the struct DB file name
        structdbname = \
            rscriptoptions[\
                util.get_option_key(options = rscriptoptions, \
                                    option = "structdbfile")]
        
        # get the struct DB file path
        structdbfile = os.path.join(wd, structdbname)
        # remove the struct DB file
        os.remove(structdbfile)

        # if the deepest level of cleaning has been requested
        if level == "all":
            for item in os.listdir(wd):
                # remove also PDB files and scorefiles
                if item.endswith(".pdb") or item.endswith(".sc"):
                    os.remove(item)
