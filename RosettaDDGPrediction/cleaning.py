#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    cleaning.py
#
#    Utility functions to clean up unnecessary files after
#    running a protocol.
#
#    Copyright (C) 2022 Valentina Sora 
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


# Standard library
import os
# RosettaDDGPrediction
from . import util



def clean_folders(process,
                  step_name,
                  wd,
                  options,
                  level):
    """Clean up folders after having run a protocol step.

    NB: 'process' is passed to ensure it is finished before
    the cleaning starts (because of the Dask task graph).
    """

    # If the step is 'relax' (cartddg protocols)
    if step_name in ("relax", "relax2020"):
        clean_func = clean_relax
    
    # If the step is 'cartesian' (cartddg protocols)
    elif step_name == "cartesian":
        clean_func = clean_cartesian
    
    # If the step is 'flexddg' (flexddg protocols)
    elif step_name == "flexddg":
        clean_func = clean_flexddg
    
    # Call the appropriate cleaning function
    clean_func(wd = wd,
               options = options,
               level = level)


def clean_relax(wd,
                options,
                level):
    """Clean the folder where the 'relax' step was run.
    """

    # No cleaning procedure has been implemented so far
    if level is not None:
        errstr = "No cleaning procedure for 'relax' has " \
                 "been implemented yet."
        raise NotImplementedError(errstr)    


def clean_cartesian(wd,
                    options,
                    level):
    """Clean the folders where the 'cartesian' step was run.
    """

    # Simply return if no cleaning was requested
    if level is None:
        return

    # If the removal of PDB files has been requested
    elif level == "pdb":
        for item in os.listdir(wd):
            # Remove all PDB files
            if item.endswith(".pdb"):
                os.remove(item)            


def clean_flexddg(wd,
                  options,
                  level):
    """Clean the folders where the 'flexddg' step was run.
    """

    # Simply return if no cleaning was requested
    if level is None:
        return
    
    # If some cleaning has been requested
    elif level in ["structdbfile", "all"]:

        # Get the RosettaScript options 
        r_script_options = \
            options[util.get_option_key(options = options,
                                        option = "scriptvars")]
        # Get the struct DB file name
        struct_db_name = \
            r_script_options[\
                util.get_option_key(options = r_script_options,
                                    option = "structdbfile")]
        
        # Get the struct DB file path
        struct_db_file = os.path.join(wd, struct_db_name)
        
        # Remove the struct DB file
        os.remove(struct_db_file)

        # If the deepest level of cleaning has been requested
        if level == "all":
            for item in os.listdir(wd):
                # Remove also PDB files and scorefiles
                if item.endswith(".pdb") or item.endswith(".sc"):
                    os.remove(item)
