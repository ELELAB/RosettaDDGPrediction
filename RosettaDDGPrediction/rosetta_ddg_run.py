#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    rosetta_ddg_run.py
#
#    Run Rosetta protocols for the prediction of the ΔΔG of
#    stability upon mutation of a monomeric protein or the
#    ΔΔG of binding upon mutation of a protein complex.
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
import logging as log
import os
import os.path
import sys
# Third-party packages
import dask
from distributed import (
    Client,
    fire_and_forget,
    LocalCluster
)
import yaml
# RosettaDDGProtocols
from . import cleaning
from .defaults import (
    CONFIG_RUN_DIR,
    CONFIG_SETTINGS_DIR,
    MUT_DIR_NAME,
    MUT_DIR_PATH,
    ROSETTA_PROTOCOLS
)
from . import pythonsteps
from . import util



def main():



    ######################### ARGUMENT PARSER #########################



    description = \
        "\nRun Rosetta protocols for the prediction of the ΔΔG of " \
        "stability upon mutation of a monomeric protein or the " \
        "ΔΔG of binding upon mutation of a protein complex.\n"

    # Create the argument parser
    parser = argparse.ArgumentParser(description = description)
    general_args = \
        parser.add_argument_group("General arguments")
    sat_args = \
        parser.add_argument_group("Saturation-related arguments")


    #---------------------------- General ----------------------------#
 

    p_help = "PDB file of the wild-type structure."
    general_args.add_argument("-p", "--pdbfile",
                              type = str,
                              required = True,
                              help = p_help)

    cr_help = f"Configuration file of the protocol to be " \
              f"run. If it is a name without extension, it is " \
              f"assumed to be the name of a YAML file in " \
              f"{CONFIG_RUN_DIR}."
    general_args.add_argument("-cr", "--configfile-run",
                             type = str,
                             required = True,
                             help = cr_help)

    cs_help = \
        f"Configuration file containing settings to be used for " \
        f"the run. If it is a name without extension, it is assumed " \
        f"to be the name of a YAML file in {CONFIG_SETTINGS_DIR}. "
    general_args.add_argument("-cs", "--configfile-settings",
                              type = str,
                              required = True,
                              help = cs_help)

    r_help = "Path to the Rosetta installation directory."
    general_args.add_argument("-r", "--rosettapath",
                              type = str,
                              required = True,
                              help = r_help)

    d_help = \
        "Directory where the protocol will be run. " \
        "Default is the current working directory."
    general_args.add_argument("-d", "--rundir",
                              type = str,
                              default = os.getcwd(),
                              help = d_help)

    l_help = \
        "File containing the list of selected mutations " \
        "(or positions to perform saturation mutagenesis)."
    general_args.add_argument("-l", "--listfile",
                              type = str,
                              default = None,
                              help = l_help)

    n_help = \
        "Number of processes to be started in parallel. " \
        "Default is one process (no parallelization)."
    general_args.add_argument("-n", "--nproc",
                              type = int,
                              default = 1,
                              help = n_help)


    #-------------------- Saturation mutagenesis ---------------------#


    saturation_help = \
        "Perform saturation mutagenesis on selected positions."
    sat_args.add_argument("--saturation",
                          action = "store_true",
                          help = saturation_help)

    reslistfile_help = \
        "File containing the list of residue types to be included " \
        "in the saturation mutagenesis. It is used only if " \
        "--saturation is provided."
    sat_args.add_argument("--reslistfile",
                          type = str,
                          default = None,
                          help = reslistfile_help)

    # Collect the arguments
    args = parser.parse_args()
    
    # Files
    pdb_file = util.get_abspath(args.pdbfile)
    rosetta_path = util.get_abspath(args.rosettapath)
    run_dir = util.get_abspath(args.rundir)
    list_file = util.get_abspath(args.listfile)
    res_list_file = util.get_abspath(args.reslistfile)
    
    # Configuration files
    config_file_run = args.configfile_run
    config_file_settings = args.configfile_settings
    
    # Others
    n_proc = args.nproc
    saturation = args.saturation



    ############################## LOGGING ############################



    # Basic logging configuration
    log.basicConfig(level = log.INFO)



    ############################# SATURATION ##########################



    # Ensure that if the saturation mutagenesis was
    # requested, a reslistfile has been passed
    if saturation and res_list_file is None:
        
        errstr = \
            "You requested a saturation mutagenesis scan but did " \
            "not provide a file containing residue types to be " \
            "included in the scan."
        log.error(errstr)
        sys.exit(errstr)



    ############################## CLIENT #############################



    # Try to get the run settings from the default YAML file
    try:
        
        settings = util.get_config_settings(config_file_settings)
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = \
            f"Could not parse the configuration file " \
            f"{config_file_settings}: {e}"
        log.error(errstr)
        sys.exit(errstr)


    # Create the local cluster
    cluster = LocalCluster(n_workers = n_proc,
                           **settings["localcluster"])
    
    # Open the client from the cluster
    client = Client(cluster)



    ############################# OPTIONS #############################


    
    # Try to get the protocol options from the YAML file
    try:
        
        options = util.get_config_run(config_file = config_file_run)
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = \
            f"Could not parse the configuration " \
            f"file {config_file_run}: {e}"
        log.error(errstr)
        sys.exit(errstr)

    # Get the protocol steps
    steps = options["steps"]

    # Get the protocol family
    family = options["family"]

    # Get the full Rosetta path to the executables (path to
    # where Rosetta is installed + path to where Rosetta
    # executables are usually stored)
    exec_path = os.path.join(rosetta_path,
                             settings["rosetta"]["execpath"])

    # Get the suffix the Rosetta executables of interest should
    # have (it differs according to whether they support MPI,
    # the compiler used, etc.)
    exec_suffix = settings["rosetta"]["execsuffix"]
    


    ########################### INPUT FILES ###########################

 

    # Try to load the PDB file
    try:
        
        # Set the current PDB file to the PDB passed by the user
        # (after checking it)
        curr_pdb_file = util.check_pdb_file(pdb_file = pdb_file,
                                            **options["pdb"])
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = f"Could not load the PDB file {pdb_file}: {e}"
        log.error(errstr)
        sys.exit(errstr)

    # If the file with the list of mutations was passed
    if list_file:
        
        # Try to retrieve the options to be used for
        # handling mutations
        try:
            
            mut_options = options["mutations"]

        # If something went wrong, report it and exit
        except KeyError:
            
            errstr = f"No options to handle mutations found in " \
                     f"the configuration file {config_file_run}."
            log.error(errstr)
            sys.exit(errstr)
        
        # Try to generate the list of mutations
        try:
            
            mutations, mutations_original = \
                util.get_mutations(\
                    list_file = list_file,
                    res_list_file = res_list_file,
                    pdb_file = curr_pdb_file,
                    res_numbering = mut_options["resnumbering"],
                    extra = mut_options["extra"],
                    n_struct = mut_options["nstruct"])
        
        # If something went wrong, report it and exit
        except Exception as e:
            
            errstr = f"Could not generate the list of mutations: {e}"
            log.error(errstr)
            sys.exit(errstr)
    
    # If no list of mutations was passed
    else:

        # Inform the user that no mutations will be performed
        logstr = \
            "No list of mutations was passed. Therefore, no " \
            "mutations will be performed."
        log.info(logstr)
        
        # No mutations will be performed
        mutations, mutations_original = [], [] 



    ############################### RUN ###############################


    
    # Create an empty list to keep track of the running futures
    futures = []

    # Create an empty list to keep track of all the directories
    # where the previous step was run
    prev_rundirs = []

    # For each step of the protocol that has to be run
    for step_name, step in steps.items():

        # If there are still pending futures from the previous step
        if futures:
            
            # Gather them
            client.gather(futures)
            
            # Clear the list of pending futures
            futures = []            
        
        # Try to Get the step options
        step_opts = step["options"]

        # Get the step features (hard-coded, they define what
        # the step is internally)
        step_features = ROSETTA_PROTOCOLS[family][step_name]

        # If the step is run via Rosetta commands
        if step_features["run_by"] == "rosetta":

            # Get the step cleaning level
            clean_level = step["cleanlevel"]

            # Get the role of the step
            role = step_features["role"]

            # Get the step working directory
            step_wd = step["wd"]

            # If the specified directory was "."
            if step_wd == ".":
                
                # Run in the current directory without generating
                # a sub-directory for the current step
                step_wd = run_dir
            
            else:
                
                # Set a new directory
                step_wd = os.path.join(run_dir, step_wd)


            # If it is a processing step
            if role == "processing":

                # If it is a relax step
                if step_name in ("relax", "relax2020"):
                
                    # Run the step
                    process = \
                        client.submit(util.run_relax,
                                      step_features = step_features,
                                      exec_path = exec_path,
                                      exec_suffix = exec_suffix,
                                      step = step,
                                      step_wd = step_wd,
                                      step_opts = step_opts,
                                      curr_pdb_file = curr_pdb_file,
                                      settings = settings,
                                      n_proc = n_proc)

                    # Append the process to the list of futures so that
                    # it gets gathered before the next step
                    futures.append(process)

                    # Submit also the post-run cleaning (we can fire and
                    # forget about this one since no other task depends
                    # on the cleaning)
                    fire_and_forget(\
                        client.submit(cleaning.clean_folders,
                                      step_name = step_name,
                                      wd = step_wd,
                                      options = step_opts,
                                      level = clean_level,
                                      wait_on = [process]))

                # Add the directory where the step was run to
                # the list of directories
                prev_rundirs.append(step_wd)

            # If it is a ΔΔG prediction step
            elif role == "ddg":

                # Write out the file mapping the directory
                # names to the mutations
                util.write_mutinfo_file(\
                    mutations_original = mutations_original,
                    out_dir = step_wd,
                    mutinfo_file = mut_options["mutinfofile"])

                # Log the order in which the mutations will be
                # performed
                _mut_list = \
                    [m[MUT_DIR_NAME] for m in mutations_original]
                mut_list = \
                    [m for m in list(dict.fromkeys(_mut_list))]
                logstr = f"The following mutations will be " \
                         f"performed:\n{', '.join(mut_list)}."
                log.info(logstr)

                # For each mutation
                for mut, mut_orig in zip(mutations, mutations_original):

                    # Set the path to the mutation directory
                    mut_wd = \
                        os.path.join(step_wd, mut_orig[MUT_DIR_PATH])              
                            
                    # If the step is cartesian ΔΔG calculation
                    if step_name in ("cartesian", "cartesian2020"):
                        
                        # Run the step
                        process = \
                            client.submit(\
                                util.run_cartesian,
                                step_features = step_features,
                                exec_path = exec_path,
                                exec_suffix = exec_suffix,
                                mut = mut,
                                mut_wd = mut_wd,
                                step = step,
                                step_opts = step_opts,
                                curr_pdb_file = curr_pdb_file,
                                settings = settings)

                    # If the step is Flex ddG ΔΔG calculation
                    elif step_name == "flexddg":

                        # Run the step
                        process = \
                            client.submit(\
                                util.run_flexddg,
                                step_features = step_features,
                                exec_path = exec_path,
                                exec_suffix = exec_suffix,
                                mut = mut,
                                mut_wd = mut_wd,
                                step = step,
                                step_opts = step_opts,
                                curr_pdb_file = curr_pdb_file,
                                settings = settings)                   

                    # Append the process to the list of futures so that
                    # it gets gathered before the next step
                    futures.append(process)

                    # Submit also the post-run cleaning
                    fire_and_forget(\
                        client.submit(cleaning.clean_folders,
                                      step_name = step_name,
                                      wd = mut_wd,
                                      options = step_opts,
                                      level = clean_level,
                                      wait_on = [process]))


        # If the step is run by Python
        elif step_features["run_by"] == "python":

            # If the step is structure selection
            if step_name == "structure_selection":

                # Run the step and return the path to the PDB
                # file to be used in the next step
                curr_pdb_file = \
                    run_structure_selection(step_opts = step_opts,
                                            curr_pdb_file = curr_pdb_file,
                                            prev_opts = prev_opts)


        # Store the previous step options
        prev_opts = step_opts

        # Store the previous step working directory
        prev_wd = step_wd


    # Gather the futures pending after running all steps
    client.gather(futures)


if __name__ == "__main__":
    main()