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
from distributed import Client, LocalCluster
import yaml
# RosettaDDGProtocols
from . import cleaning
from .defaults import (
    CONFIG_RUN_DIR,
    CONFIG_SETTINGS_DIR,
    MUT_DIR_NAME,
    MUT_DIR_PATH,
    ROSETTA_OPTIONS,
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
        
        errstr = f"Could not parse the configuration file " \
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
        
        options = client.submit(util.get_config_run,
                                config_file = config_file_run).result()
    
    # If something went wrong, report it and exit
    except Exception as e:
        
        errstr = f"Could not parse the configuration " \
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
        curr_pdb_file = client.submit(util.check_pdb_file,
                                      pdb_file = pdb_file,
                                      **options["pdb"]).result()
    
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
                client.submit(\
                    util.get_mutations,
                        list_file = list_file,
                        res_list_file = res_list_file,
                        pdb_file = curr_pdb_file,
                        res_numbering = mut_options["resnumbering"],
                        extra = mut_options["extra"],
                        n_struct = mut_options["nstruct"]).result()
        
        # If something went wrong, report it and exit
        except Exception as e:
            
            errstr = f"Could not generate the list of mutations: {e}"
            log.error(errstr)
            sys.exit(errstr)
    
    # If no list of mutations was passed
    else:
        
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
        
        # Get the step options
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

            # Get the executable needed to run the step
            executable = step_features["executable"]

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

            # Try to get the Rosetta executable
            try:
                
                executable = \
                    client.submit(util.get_rosetta_executable,
                                  exec_name = executable,
                                  exec_path = exec_path,
                                  exec_suffix = exec_suffix)
            
            # If something went wrong, report it and exit
            except Exception as e:
                
                errstr = f"Could not get Rosetta executable " \
                         f"'{executable}' from {exec_path}."
                log.error(errstr)
                sys.exit(errstr)

            # If it is a processing step
            if role == "processing":
                
                # Get and update the step options
                step_opts = client.submit(util.update_options,
                                          options = step_opts,
                                          pdb_file = curr_pdb_file)
                
                # Set the flags file and Rosetta output
                flagsfile = os.path.join(step_wd, step["flagsfile"])
                output = os.path.join(step_wd, step["output"])
                
                # Write the flagsfile
                flagsfile = client.submit(util.write_flagsfile,
                                          options = step_opts,
                                          flagsfile = flagsfile)
                
                # Submit the calculation
                # NB: all processes will be used as MPI processes
                # if MPI is available, since there is no need for
                # parallelization over the mutations (it is a
                # processing step).
                try:
                    
                    process = \
                        client.submit(\
                            util.run_rosetta,
                            executable = executable,
                            flagsfile = flagsfile,
                            output = output,
                            wd = step_wd,
                            use_mpi = settings["mpi"]["usempi"],
                            mpi_exec = settings["mpi"]["mpiexec"],
                            mpi_args = settings["mpi"]["mpiargs"],
                            mpi_n_proc = n_proc)
               
                # If something went wrong, report it and exit
                except Exception as e:
                    
                    errstr = f"'{step_name}' run in {step_wd} exited " \
                             f"with code {process['returncode']}. " \
                             f"The following exception occurred: {e}"
                    log.errstr(errstr)
                    sys.exit(errstr)

                # Submit also the post-run cleaning
                futures.append(\
                    client.submit(cleaning.clean_folders,
                                  process = process,
                                  step_name = step_name,
                                  wd = step_wd,
                                  options = step_opts,
                                  level = clean_level))

                # Add the directory where the step was run to
                # the list of directories
                prev_rundirs.append(step_wd)

            # If it is a ΔΔG prediction step
            elif role == "ddg":

                # Write out the file mapping the directory
                # names to the mutations
                futures.append(\
                    client.submit(\
                        util.write_mutinfo_file,
                        mutations_original = mutations_original,
                        out_dir = step_wd,
                        mutinfo_file = mut_options["mutinfofile"]))

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
                    
                    # Set the flags file and Rosetta output
                    flagsfile = os.path.join(mut_wd, step["flagsfile"])
                    output = os.path.join(mut_wd, step["output"])          

                    # If the step is cartesian ΔΔG calculation
                    if step_name == "cartesian":
                        
                        # Get the keyword used to specify the mutfile
                        # in the configuration file
                        mutfile_key = \
                            client.submit(util.get_option_key,
                                          options = step_opts,
                                          option = "mutfile").result()
                        
                        # Set the path to the mutfile that will be
                        # written
                        mutfile = \
                            os.path.join(mut_wd, step_opts[mutfile_key])
                        
                        # Write the mutfile
                        client.submit(util.write_mutfile,
                                      mut = mut,
                                      mutfile = mutfile).result()

                        # Update the options for the current mutation
                        # (add input PDB file and mutation-specific
                        # options)
                        opts_mut = \
                            client.submit(util.update_options,
                                          options = step_opts,
                                          pdb_file = curr_pdb_file)

                    # If the step is Flex ddG ΔΔG calculation
                    elif step_name == "flexddg":
                        
                        # Get the keyword used to specify the Rosetta
                        # script variables in the configuration file
                        script_vars_key = \
                            client.submit(\
                                util.get_option_key,
                                options = step_opts,
                                option = "script_vars").result()
                        
                        # Get the keyword used to specify the resfile
                        # in the configuration file
                        resfile_key = \
                            client.submit(\
                                util.get_option_key,
                                options = step_opts[script_vars_key],
                                option = "resfile").result()
                        
                        # Set the path to the resfile that will be
                        # written
                        resfile = \
                            os.path.join(\
                                mut_wd,
                                step_opts[script_vars_key][resfile_key])
                        
                        # Write the resfile
                        client.submit(util.write_resfile,
                                      mut = mut,
                                      resfile = resfile).result()

                        # Update the options for the current mutation
                        # (add input PDB file and mutation-specific
                        # options)
                        opts_mut = \
                            client.submit(util.update_options,
                                          options = step_opts,
                                          pdb_file = curr_pdb_file,
                                          mut = mut)                                

                    # Write the flagsfile
                    flagsfile = client.submit(util.write_flagsfile,
                                              options = opts_mut,
                                              flagsfile = flagsfile)

                    # Submit the calculation
                    # NB: only one MPI process will be used if MPI is
                    # available since there is no gain in using MPI
                    # with cartesian_ddg or the Flex ddG procedure.
                    # The parallelization is done over the mutations.
                    try:
                        
                        process = \
                            client.submit(\
                                util.run_rosetta,
                                executable = executable,
                                flagsfile = flagsfile,
                                output = output,
                                wd = mut_wd,
                                use_mpi = settings["mpi"]["usempi"],
                                mpi_exec = settings["mpi"]["mpiexec"],
                                mpi_args = settings["mpi"]["mpiargs"],
                                mpi_n_proc = 1)
                    
                    # If something went wrong, report it and exit
                    except Exception as e:
                        
                        errstr = \
                            f"'{step_name}' run in {mut_wd} exited " \
                            f"with code {process['returncode']}. " \
                            f"The following exception occurred: {e}"
                        
                        # Log the error and exit
                        log.errstr(errstr)
                        sys.exit(errstr)

                    # If the step was a flexddg step
                    if step_name == "flexddg":

                        # Get the options regarding the extraction
                        # of structures from the database file
                        extract_options = step["extract_structures"]

                        # The extraction of the structures is not
                        # a 'step' per se because it needs to be
                        # run in every folder where the flexddg
                        # step was run and before the unnecessary
                        # output files are deleted, since one
                        # of these files is needed to extract
                        # the structures

                        # Get whether the user has requested the
                        # extraction
                        if extract_options["extract"]:

                            # Check whether some steps of the extraction
                            # have already been performed
                            is_struct_extracted, is_struct_renamed = \
                                client.submit(\
                                    util.check_structures_extraction,
                                    wd = mut_wd).result()

                            # Get the RosettaScripts options specified
                            # in the flexddg protocol
                            r_script_options = \
                                step_opts[script_vars_key]

                            # Initialize the variable storing the
                            # extraction process to None
                            # (it is needed since it passed to the
                            # rename function to be sure the renaming
                            # occurs only after the extraction,
                            # but it can happen that only the 
                            # extraction has been performed but not
                            # the renaming, and we want only the
                            # renaming to be performed)
                            process_extract = None

                            # If the structures have not been
                            # extracted yet
                            if not is_struct_extracted:

                                # Get the name of the Rosetta 
                                # executable responsible
                                # for the extraction
                                exec_extract_name = \
                                    step_features["executable_extract"]

                                # Set the flags file
                                flagsfile_extract = \
                                    os.path.join(\
                                        mut_wd,
                                        extract_options["flagsfile"])

                                # Set the Rosetta output
                                output_extract = \
                                    os.path.join(\
                                        mut_wd,
                                        extract_options["output"])

                                # Try to get the Rosetta executable
                                # The 'wait_on' option does not do
                                # anything apart from telling Dask
                                # that this task should only be
                                # executed if the flexddg run
                                # completed successfully
                                try:
                                    
                                    executable_extract = \
                                        client.submit(\
                                            util.get_rosetta_executable,
                                            exec_name = exec_extract_name,
                                            exec_path = exec_path,
                                            exec_suffix = exec_suffix,
                                            wait_on = [process])

                            
                                # If something went wrong, report it
                                except Exception as e:
                                    
                                    errstr = \
                                        f"Could not get Rosetta " \
                                        f"executable " \
                                        f"'{exec_extract_name}' " \
                                        f"from {exec_path}."
                                    log.error(errstr)

                                # Get the options to be passed to the
                                # Rosetta executable
                                opts_extract = \
                                    extract_options["options"]

                                # Get the key used to define the database
                                # file containing the structures in the
                                # configuration file
                                # The 'wait_on' option does not do
                                # anything apart from telling Dask
                                # that this task should only be
                                # executed if the flexddg run
                                # completed successfully
                                db_name_key = \
                                    client.submit(\
                                        util.get_option_key,
                                        options = r_script_options,
                                        option = "struct_db_file",
                                        wait_on = [process]).result()

                                # Get the path to the database file
                                db_file = \
                                    os.path.join(\
                                        mut_wd,
                                        r_script_options[db_name_key])

                                # Get the key that will be used to define
                                # the database file when given as input
                                # for the extraction
                                db_file_key = \
                                    sorted(ROSETTA_OPTIONS["db_name"],
                                           key = len,
                                           reverse = True)[0]

                                # Update the options with the database file
                                opts_extract.update(\
                                    {db_file_key : db_file})

                                # Write the flagsfile
                                # The 'wait_on' option does not do
                                # anything apart from telling Dask
                                # that this task should only be
                                # executed if the flexddg run
                                # completed successfully
                                flagsfile_extract = \
                                    client.submit(\
                                        util.write_flagsfile,
                                        options = opts_extract,
                                        flagsfile = flagsfile_extract,
                                        wait_on = [process])

                                # Submit the extraction
                                try:
                                    
                                    process_extract = \
                                        client.submit(\
                                            util.run_rosetta,
                                            executable = executable_extract,
                                            flagsfile = flagsfile_extract,
                                            output = output_extract,
                                            wd = mut_wd,
                                            use_mpi = settings["mpi"]["usempi"],
                                            mpi_exec = settings["mpi"]["mpiexec"],
                                            mpi_args = settings["mpi"]["mpiargs"],
                                            mpi_n_proc = 1)
                                
                                # If something went wrong, report it
                                except Exception as e:
                                    
                                    errstr = \
                                        f"The extraction of the " \
                                        f"structures from the file " \
                                        f"'{db_file}' in {mut_wd} " \
                                        f"exited with code " \
                                        f"{process['returncode']}. " \
                                        f"The following exception " \
                                        f"occurred: {e}"
                                    
                                    # Log the error
                                    log.errstr(errstr)

                                # Turn on the flag indicating that the
                                # structures have now been extracted
                                is_struct_extracted = True

                            # If the structures have been extracted
                            # but not renamed
                            if (is_struct_extracted) and \
                            (not is_struct_renamed):

                                # Try to rename the structures
                                # The 'wait_on' option does not do
                                # anything apart from telling Dask
                                # that this task should only be
                                # executed if the extraction 
                                # completed successfully
                                try:
                                        
                                    process_rename = \
                                        client.submit(\
                                            util.rename_structures_flexddg,
                                            path = mut_wd,
                                            r_script_options = r_script_options,
                                            wait_on = [process_extract])
                                    
                                # If something went wrong, report it
                                except Exception as e:
                                        
                                    errstr = \
                                        f"Could not rename the  " \
                                        f"structures extracted from " \
                                        f"the file '{db_file}' in " \
                                        f"{mut_wd}." \
                                        f"The following exception " \
                                        f"occurred: {e}"
                                        
                                    # Log the error
                                    log.errstr(errstr)


                    # Submit also the post-run cleaning
                    futures.append(\
                        client.submit(cleaning.clean_folders,
                                      process = process,
                                      step_name = step_name,
                                      wd = mut_wd,
                                      options = opts_mut,
                                      level = clean_level))

        # If the step is run by Python
        elif step_features["run_by"] == "python":

            # If the step is structure selection
            if step_name == "structure_selection":

                # Get the input file, the file type and the
                # selection criterion
                in_file, in_file_type, select = \
                    step_opts["infile"], step_opts["infiletype"], \
                    step_opts["select"]

                # Get the path to the input file (assuming the
                # input file is specified as either a file name
                # or a relative path starting from the working
                # directory)
                in_file = os.path.join(run_dir, in_file)
                
                # Try to select the structure
                try:
                    
                    # Run the selection
                    struct = \
                        client.submit(pythonsteps.select_structure,
                                      in_file = in_file,
                                      in_file_type = in_file_type,
                                      select = select)
                
                # If something went wrong, report it and exit
                except Exception as e:
                    
                    errstr = f"Could not not perform the " \
                             f"structure selection: {e}"
                    log.error(errstr)
                    sys.exit(errstr)
                
                # Get the name of the current PDB file
                pdb_file = os.path.basename(curr_pdb_file)

                # Try to get the name of the selected PDB file by
                # using the options from the previous (Rosetta) step
                # to retrieve the correct PDB file
                try:
                    
                    curr_pdb_file_name = \
                        client.submit(util.get_out_pdb_name,
                                      options = prev_opts,
                                      pdb_file = pdb_file, 
                                      struct = struct).result()
                
                # If something went wrong, report the error and exit
                except Exception as e:
                    
                    errstr = f"Could not get the name of the PDB " \
                             f"of the selected structure: {e}"
                    log.error(errstr)
                    sys.exit(errstr)

                # Get the path to the PDB file
                curr_pdb_file = \
                    os.path.join(prev_wd, curr_pdb_file_name)

        # Store the previous step options
        prev_opts = step_opts

        # Store the previous step working directory
        prev_wd = step_wd

    # Gather the futures pending after running all steps
    client.gather(futures)


if __name__ == "__main__":
    main()