#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    rosetta_ddg_check_run.py
#
#    Check that the Rosetta calculations have terminated without
#    errors.
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
import argparse
import itertools
import logging as log
import os
import os.path
import sys
# RosettaDDGProtocols
from .defaults import (
    CONFIG_RUN_DIR,
    MUTINFO_COLS,
    )
from . import util



def main():



    ######################### ARGUMENT PARSER #########################



    description = \
        "\nCheck that the Rosetta calculations have terminated " \
        "without errors.\n"

    # Create the argument parser
    parser = argparse.ArgumentParser(description = description)

    # Add arguments
    cr_help = f"Configuration file of the protocol that was " \
              f"run. If it is a name without extension, it is " \
              f"assumed to be the name of a YAML file in " \
              f"{CONFIG_RUN_DIR}."
    parser.add_argument("-cr", "--configfile-run",
                        type = str,
                        required = True,
                        help = cr_help)

    d_help = "Directory where the protocol was run. " \
             "Default is the current working directory."
    parser.add_argument("-d", "--running-dir",
                        type = str,
                        default = os.getcwd(),
                        help = d_help)

    mf_help = "File with info about the mutations (it " \
              "is created when running). Not needed if the " \
              "calculation to be checked did not perform " \
              "any mutations."
    parser.add_argument("-mf", "--mutinfofile",
                        type = str,
                        default = None,
                        help = mf_help)

    # Parse the arguments
    args = parser.parse_args()

    # Configuration files
    config_file_run = args.configfile_run
    
    # Directories
    run_dir = util.get_abspath(args.running_dir)
    
    # Others
    mutinfo_file = util.get_abspath(args.mutinfofile)



    ############################## LOGGING ############################



    # Basic logging configuration
    log.basicConfig(level = log.INFO)



    ########################## CONFIGURATION ##########################



    # Try to get the configuration of the run from the
    # corresponding configuration file
    try:
        
        config_run = util.get_config_run(config_file_run)
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = f"Could not parse the configuration file " \
                 f"{config_file_run}: {e}"
        log.error(errstr)
        sys.exit(errstr)



    ############################## CHECK ##############################



    # Get the family of the protocol that was run
    family = config_run["family"]


    # If a mutinfo file has been passed
    if mutinfo_file is not None:
        
        # Get info about the mutations
        try:
            
            mutinfo = util.get_mutinfo(mutinfo_file = mutinfo_file)
        
        # If something went wrong, report it and exit
        except Exception as e:
            
            errstr = \
                f"Could not load mutations' info from " \
                f"{mutinfo_file}: {e}"
            log.error(errstr)
            sys.exit(errstr)

        # Get the mutations' directories
        dirs_mutations = \
            [os.path.join(run_dir, d) for d \
             in mutinfo[MUTINFO_COLS["dir_name"]]]

        # If it is a cartddg protocol
        if family == "cartddg":

            # The directories where the calculations have been
            # performed are simply the directories containing
            # the mutations
            dirs_paths = dirs_mutations

        # If it is a flexddg protocol
        if family == "flexddg":

            # Get the number of structures generated
            n_struct = config_run["mutations"]["nstruct"]
            
            # Format the structure names as strings
            struct_nums = \
                [str(num) for num in range(1, n_struct + 1)]

            # The directories where the calculations have been
            # performed are sub-directories (one per structure
            # generated) of the mutations' directories
            dirs_paths = \
                list(itertools.chain(\
                    *[[os.path.join(d, sn) for sn in struct_nums] \
                    for d in dirs_mutations]))
    
    # If no mutinfo file has been passed
    else:

        # The only directory to be checked is the running directory
        dirs_paths = [run_dir]

    # Try to check the runs performed
    try:

        crashed_runs_paths = util.check_rosetta_run(dirs_paths)
        
    # If something went wrong, report it and exit
    except Exception as e:
            
        errstr = \
            f"Could not check the runs performed in {run_dir}: {e}"
        log.error(errstr)
        sys.exit(errstr)


    # If no run has crashed
    if not crashed_runs_paths:

        # Log it and report it to the standard output
        logstr = "No crashed run has been found."
        log.info(logstr)
        sys.stdout.write(logstr + "\n")

    # Otherwise
    else:

        # For each run that has crashed
        for crashed_run_path in crashed_runs_paths:

            # Log it and report it to the standard output
            logstr = \
                f"The run in the following directory " \
                f"reported a crash: {crashed_run_path}."
            log.warning(logstr)
            sys.stdout.write(logstr + "\n")
