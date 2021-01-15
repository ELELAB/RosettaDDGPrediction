#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    rosetta_ddg_run.py
#
#    Run Rosetta protocols for the prediction of the ΔΔG of
#    stability upon mutation of a monomeric protein or the
#    ΔΔG of binding upon mutation of a protein complex.
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
import argparse
import logging as log
import os
import os.path
import sys
# third-party packages
import dask
from distributed import Client, LocalCluster
import yaml
# RosettaDDGProtocols
from . import cleaning
from .defaults import (
    CONFIGRUNDIR,
    CONFIGSETTINGSDIR,
    MUTDIRPATH,
    ROSETTAPROTOCOLS
)
from . import pythonsteps
from . import util



def main():



    ######################### ARGUMENT PARSER #########################



    description = \
        "\nRun Rosetta protocols for the prediction of the ΔΔG of " \
        "stability upon mutation of a monomeric protein or the " \
        "ΔΔG of binding upon mutation of a protein complex.\n"

    # create the argument parser
    parser = argparse.ArgumentParser(description = description)
    generalargs = parser.add_argument_group("General arguments")
    satargs = parser.add_argument_group("Saturation-related arguments")


    #---------------------------- General ----------------------------#
 

    p_help = "PDB file of the wild-type structure."
    generalargs.add_argument("-p", "--pdbfile", \
                             type = str, \
                             required = True, \
                             help = p_help)

    cr_help = f"Configuration file of the protocol to be " \
              f"run. If it is a name without extension, it is " \
              f"assumed to be the name of a YAML file in " \
              f"{CONFIGRUNDIR}."
    generalargs.add_argument("-cr", "--configfile-run", \
                             type = str, \
                             required = True, \
                             help = cr_help)

    cs_help = \
        f"Configuration file containing settings to be used for " \
        f"the run. If it is a name without extension, it is assumed " \
        f"to be the name of a YAML file in {CONFIGSETTINGSDIR}. "
    generalargs.add_argument("-cs", "--configfile-settings", \
                             type = str, \
                             required = True, \
                             help = cs_help)

    r_help = "Path to the Rosetta installation directory."
    generalargs.add_argument("-r", "--rosettapath", \
                             type = str, \
                             required = True, \
                             help = r_help)

    d_help = \
        "Directory where the protocol will be run. " \
        "Default is the current working directory."
    generalargs.add_argument("-d", "--rundir", \
                             type = str, \
                             default = os.getcwd(), \
                             help = d_help)

    l_help = \
        "File containing the list of selected mutations " \
        "(or positions to perform saturation mutagenesis)."
    generalargs.add_argument("-l", "--listfile", \
                             type = str, \
                             default = None, \
                             help = l_help)

    n_help = \
        "Number of processes to be started in parallel. " \
        "Default is one process (no parallelization)."
    generalargs.add_argument("-n", "--nproc", \
                             type = int, \
                             default = 1, \
                             help = n_help)


    #-------------------- Saturation mutagenesis ---------------------#


    saturation_help = \
        "Perform saturation mutagenesis on selected positions."
    satargs.add_argument("--saturation", \
                         action = "store_true", \
                         help = saturation_help)

    reslistfile_help = \
        "File containing the list of residue types to be included " \
        "in the saturation mutagenesis It is used only if " \
        "--saturation is provided."
    satargs.add_argument("--reslistfile", \
                         type = str, \
                         default = None, \
                         help = reslistfile_help)

    # collect the arguments
    args = parser.parse_args()
    # files
    pdbfile = util.get_abspath(args.pdbfile)
    rosettapath = util.get_abspath(args.rosettapath)
    rundir = util.get_abspath(args.rundir)
    listfile = util.get_abspath(args.listfile)
    reslistfile = util.get_abspath(args.reslistfile)
    # configuration files
    configfilerun = args.configfile_run
    configfilesettings = args.configfile_settings
    # others
    nproc = args.nproc
    saturation = args.saturation

    # ensure that if the saturation mutagenesis was
    # requested, a reslistfile has been passed
    if saturation and reslistfile is None:
        errstr = "You requested a saturation mutagenesis " \
                 "scan but did not provide a reslistfile."
        sys.exit(errstr)

    # get the name of the configuration file for running
    # the protocol
    configrunname = \
        os.path.basename(configfilerun).rstrip(".yaml")
    # get the name of the configuration file for the
    # running options
    configsettingsname = \
        os.path.basename(configfilesettings).rstrip(".yaml")

    # if the configuration file is a name without extension
    if configfilerun == configrunname:
        # assume it is a configuration file in the directory
        # storing configuration files for running protocols
        configfilerun = os.path.join(CONFIGRUNDIR, \
                                     configrunname + ".yaml")
    # otherwise assume it is a file name/file path
    else:
        configfilerun = util.get_abspath(configfilerun)

    # if the configuration file is a name without extension
    if configfilesettings == configsettingsname:
        # assume it is a configuration file in the directory
        # storing configuration files for run settings
        configfilesettings = os.path.join(CONFIGSETTINGSDIR, \
                                          configsettingsname + ".yaml")
    # otherwise assume it is a file name/file path
    else:
        configfilesettings = util.get_abspath(configfilesettings)



    ############################## CLIENT #############################



    # try to get the run settings from the default YAML file
    try:
        settings = yaml.safe_load(open(configfilesettings, "r"))
    # if something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not parse the configuration file " \
                 f"{configfilesettings}: {e}"
        log.error(errstr)
        sys.exit(errstr)

    # create the local cluster
    cluster = LocalCluster(n_workers = nproc, \
                           **settings["localcluster"])
    
    # open the client from the cluster
    client = Client(cluster)



    ############################# OPTIONS #############################


    
    # try to get the protocol options from the YAML file
    try:
        options = client.submit(util.get_config_run, \
                                configfile = configfilerun).result()
    # if something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not parse the configuration " \
                 f"file {configfilerun}: {e}"
        log.error(errstr)
        sys.exit(errstr)

    
    # get the protocol steps
    steps = options["steps"]

    # get the protocol family
    family = options["family"]

    # get the full Rosetta path to the executables (path to
    # where Rosetta is installed + path to where Rosetta
    # executables are usually stored)
    execpath = os.path.join(rosettapath, \
                            settings["rosetta"]["execpath"])

    # get the suffix the Rosetta executables of interest should
    # have (it differs according to whether they support MPI, to
    # the compiler used, etc.)
    execsuffix = settings["rosetta"]["execsuffix"]
    


    ########################### INPUT FILES ###########################

 

    # try to load the PDB file
    try:
        # set the current PDB file to the PDB passed by the user
        # (after checking it)
        currpdbfile = client.submit(util.check_pdbfile, \
                                    pdbfile = pdbfile, \
                                    **options["pdb"]).result()
    # if something went wrong, report it and exit
    except Exception as e:
        errstr = f"Could not load the PDB file {pdbfile}: {e}"
        log.error(errstr)
        sys.exit(errstr)
    
    
    # if the file with the list of mutations was passed
    if listfile:
        
        # get the mutations options
        mutoptions = options["mutations"]
        
        # try to generate the list of mutations
        try:
            mutations = \
                client.submit(\
                    util.get_mutations, \
                        listfile = listfile, \
                        reslistfile = reslistfile, \
                        pdbfile = currpdbfile, \
                        resnumbering = mutoptions["resnumbering"], \
                        extra = mutoptions["extra"], \
                        nstruct = mutoptions["nstruct"])
        
        # if something went wrong, report it and exit
        except Exception as e:
            errstr = f"Could not generate the list of mutations: {e}"
            log.error(errstr)
            sys.exit(errstr)
    
    # otherwise
    else:
        # no mutations will be performed
        mutations = [] 



    ############################### RUN ###############################


    
    # create an empty list to keep track of the running futures
    futures = [] 

    
    # for each step of the protocol that has to be run
    for stepname, step in steps.items():

        # if there are still pending futures from the previous step
        if futures:
            # gather them
            client.gather(futures)
            # clear the list of pending futures
            futures = []            
        
        # get the step options
        stepopts = step["options"]

        # get the step features
        stepfeatures = ROSETTAPROTOCOLS[family][stepname]

        # get the step cleaning level
        cleanlevel = step["cleanlevel"]


        # if the step is run via Rosetta commands
        if stepfeatures["runby"] == "rosetta":

            # get the role of the step
            role = stepfeatures["role"]

            # get the executable needed to run the step
            executable = stepfeatures["executable"]

            # if the specified directory was "."
            if step["wd"] == ".":
                # run in the current directory without generating
                # a sub-directory for the current step
                stepwd = rundir
            else:
                # set a new directory
                stepwd = os.path.join(rundir, step["wd"])

            # try to get the Rosetta executable
            try:
                executable = \
                    client.submit(util.get_rosetta_executable, \
                                  execname = executable, \
                                  execpath = execpath, \
                                  execsuffix = execsuffix)
            
            # if something went wrong, report it and exit
            except Exception as e:
                errstr = f"Could not get Rosetta executable " \
                         f"'{executable}' from {execpath}."
                log.error(errstr)
                sys.exit(errstr)

        
            # if it is a processing step
            if role == "processing":
                
                # get and update the step options
                stepopts = client.submit(util.update_options, \
                                         options = stepopts, \
                                         pdbfile = currpdbfile) 
                
                # set the flags file and Rosetta output
                flagsfile = os.path.join(stepwd, step["flagsfile"])
                output = os.path.join(stepwd, step["output"])
                
                # write the flagsfile
                flagsfile = client.submit(util.write_flagsfile, \
                                          options = stepopts, \
                                          flagsfile = flagsfile)
                
                # submit the calculation
                # NB: all processes will be used as MPI processes
                # if MPI is available, since there is no need for
                # parallelization over the mutations (it is a
                # processing step).
                try:
                    process = client.submit(util.run_rosetta, \
                                            executable = executable, \
                                            flagsfile = flagsfile, \
                                            output = output, \
                                            mpinproc = nproc, \
                                            wd = stepwd, \
                                            **settings["mpi"])
               
                # if something went wrong, report it and exit
                except Exception as e:
                    errstr = f"'{stepname}' run in {stepwd} exited " \
                             f"with code {process["returncode"]}. " \
                             f"The following exception occurred: {e}"
                    log.errstr(errstr)
                    sys.exit(errstr)

                # submit also the post-run cleaning
                futures.append(cleaning.clean_folders, \
                               process = process, \
                               step = stepname, \
                               wd = stepwd, \
                               options = stepopts, \
                               level = cleanlevel)

          
            # if it is a ΔΔG prediction step
            elif role == "ddg":

                # write out the file mapping the directory
                # names to the mutations
                futures.append(\
                    client.submit(util.write_mutinfofile, \
                                  mutations = mutations, \
                                  outdir = stepwd, \
                                  mutinfofile = mutoptions["mutinfofile"]))

                # log the order in which the mutations will be performed
                logstr = f"The following mutations will be " \
                         f"performed:\n{'\n'.join(mutations)}."

                # for each mutation
                for mut in mutations:
                    
                    # set the path to the mutation directory
                    mutwd = os.path.join(stepwd, mut[MUTDIRPATH])              
                    
                    # set the flags file and Rosetta output
                    flagsfile = os.path.join(mutwd, step["flagsfile"])
                    output = os.path.join(mutwd, step["output"])          

                    
                    # if the step is cartesian ΔΔG calculation
                    if stepname == "cartesian":
                        
                        # get the keyword used to specify the mutfile
                        # in the configuration file
                        mutfilekey = \
                            client.submit(util.get_option_key, \
                                          options = stepopts, \
                                          option = "mutfile").result()
                        
                        # set the path to the mutfile that will be written
                        mutfile = os.path.join(mutwd, stepopts[mutfilekey])
                        
                        # write the mutfile
                        client.submit(util.write_mutfile, \
                                      mut = mut, \
                                      mutfile = mutfile).result()

                        # update the options for the current mutation
                        # (add input PDB file and mutation-specific options)
                        optsmut = client.submit(util.update_options, \
                                                options = stepopts, \
                                                pdbfile = currpdbfile)


                    # if the step is Flex ddG ΔΔG calculation
                    elif stepname == "flexddg":
                        
                        # get the keyword used to specify the Rosetta
                        # script variables in the configuration file
                        scriptvarskey = \
                            client.submit(util.get_option_key, \
                                          options = stepopts, \
                                          option = "scriptvars").result()
                        
                        # get the keyword used to specify the resfile
                        # in the configuration file
                        resfilekey = \
                            client.submit(util.get_option_key, \
                                          options = stepopts[scriptvarskey], \
                                          option = "resfile").result()
                        
                        # set the path to the resfile that will be written
                        resfile = \
                            os.path.join(mutwd, \
                                         stepopts[scriptvarskey][resfilekey])
                        
                        # write the resfile
                        client.submit(util.write_resfile, \
                                      mut = mut, \
                                      resfile = resfile).result()

                        # update the options for the current mutation
                        # (add input PDB file and mutation-specific options)
                        optsmut = client.submit(util.update_options, \
                                                options = stepopts, \
                                                pdbfile = currpdbfile, \
                                                mut = mut)                                

                    
                    # write the flagsfile
                    flagsfile = client.submit(util.write_flagsfile, \
                                              options = optsmut, \
                                              flagsfile = flagsfile)


                    # submit the calculation
                    # NB: only one MPI process will be used if MPI is
                    # available since there is no gain in using MPI
                    # with cartesian_ddg or the Flex ddG procedure.
                    # The parallelization is done over the mutations.
                    try:
                        process = client.submit(util.run_rosetta, \
                                                executable = executable, \
                                                flagsfile = flagsfile, \
                                                output = output, \
                                                mpinproc = 1, \
                                                wd = mutwd, \
                                                **settings["mpi"])
                    
                    # if something went wrong, report it and exit
                    except Exception as e:
                        errstr = f"'{stepname}' run in {stepwd} exited " \
                                 f"with code {process["returncode"]}. " \
                                 f"The following exception occurred: {e}"
                        # log the error and exit
                        log.errstr(errstr)
                        sys.exit(errstr)

                    # submit also the post-run cleaning
                    futures.append(cleaning.clean_folders, \
                                   process = process, \
                                   step = stepname, \
                                   wd = mutwd, \
                                   options = optsmut, \
                                   level = cleanlevel)

      
        # if the step is run by Python
        elif stepfeatures["runby"] == "python":

            # if the step is structure selection
            if stepname == "structure_selection":

                # get the input file, the file type and the
                # selection criterion
                infile, infiletype, select = \
                    stepopts["infile"], stepopts["infiletype"], \
                    stepopts["select"]

                # get the path to the input file (assuming the
                # input file is specified as either a file name
                # or a relative path starting from the working
                # directory)
                infile = os.path.join(rundir, infile)
                
                # try to select the structure
                try:
                    # run the selection
                    struct = \
                        client.submit(pythonsteps.select_structure, \
                                      infile = infile, \
                                      infiletype = infiletype, \
                                      select = select)
                
                # if something went wrong, report it and exit
                except Exception as e:
                    errstr = f"Could not not perform the " \
                             f"structure selection: {e}"
                    log.error(errstr)
                    sys.exit(errstr)
                
                # get the name of the current PDB file
                pdbfile = os.path.basename(currpdbfile)

                # try to get the name of the selected PDB file by
                # using the options from the previous (Rosetta) step
                # to retrieve the correct PDB file
                try:
                    currpdbfilename = \
                        client.submit(util.get_outpdbname, \
                                      options = prevopts, \
                                      pdbfile = pdbfile, 
                                      struct = struct).result()
                # if something went wrong, report the error and exit
                except Exception as e:
                    errstr = f"Could not get the name of the PDB " \
                             f"of the selected structure: {e}"
                    log.error(errstr)
                    sys.exit(errstr)

                # get the path to the PDB file
                currpdbfile = os.path.join(prevwd, currpdbfilename)


        # store the previous step options
        prevopts = stepopts

        # store the previous step working directory
        prevwd = stepwd
  

    # gather the futures pending after running all steps
    client.gather(futures)


if __name__ == "__main__":
    main()