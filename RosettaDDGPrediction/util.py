#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    General utility functions.
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
import copy
import itertools
import logging as log
import operator
import os
import os.path
import re
import subprocess
# Third-party packages
import Bio.PDB as PDB
import matplotlib.font_manager as fm
import pandas as pd
import yaml
# RosettaDDGProtocols
from .dask_patches import reset_worker_logger
from .defaults import (
    CHAIN,
    CHAIN_SEP,
    COMP_SEP,
    CONFIG_AGGR_DIR,
    CONFIG_PLOT_DIR,
    CONFIG_RUN_DIR,
    CONFIG_SETTINGS_DIR,
    DIR_MUT_SEP,
    FLEXDDG_STATES,
    MULTI_MUT_SEP,
    MUT,
    MUT_DIR_NAME,
    MUT_DIR_PATH,
    MUTINFO_COLS,
    MUTR,
    MUT_SEP, 
    NUMR,
    POSR,
    ROSETTA_CRASH_LOG,
    ROSETTA_OPTIONS,
    ROSETTA_PROTOCOLS,
    ROSETTA_SCRIPTS_DIR,
    STRUCT,
    STRUCT_EXTRACTED_PATTERN,
    WTR
)


# Get the module logger
logger = log.getLogger(__name__)



########################### ROSETTA-RELATED ###########################



def get_rosetta_executable(exec_name,
                           exec_path,
                           exec_suffix,
                           **kwargs):
    """Get the path to a Rosetta executable from the Rosetta
    installation directory.
    """

    # Get all the executables available
    all_execs = os.listdir(exec_path)
    
    # If no suffix was provided, assume the executable name does not
    # have any suffix
    exec_comp_name = \
        exec_name + exec_suffix if exec_suffix is not None \
        else exec_name

    # Filter function to get the Rosetta executables
    filt = lambda x: x == exec_comp_name
    
    # Get all Rosetta executables satifying the criteria of
    # the filter function (you should get only one executable) 
    execs = list(filter(filt, all_execs))
    
    # If no executable was found
    if len(execs) == 0:
        # Raise an exception
        errstr = \
            f"No executable found for {exec_comp_name} in {exec_path}."
        raise ValueError(errstr)
    
    # If multiple conflicting executables were found
    elif len(execs) > 1:
        # Raise an exception
        errstr = \
            f"Multiple executables found for {exec_comp_name} in " \
            f"{exec_path}."
        raise ValueError(errstr)

    # Return the path to the executable
    return os.path.join(exec_path, execs[0])


def run_rosetta(executable,
                flagsfile,
                output,
                wd,
                use_mpi,
                mpi_exec,
                mpi_args,
                mpi_n_proc,
                **kwargs):
    """Run Rosetta.
    """


    # Make sure that the working directory exists. If not, create it.
    os.makedirs(wd, exist_ok = True)
    
    # Set an empty prefix (before the Rosetta command line)
    # to run with MPI (remains empty if you do not run with MPI)
    mpi_prefix = []
    
    # If MPI is requested
    if use_mpi:
        # Update the MPI prefix
        mpi_prefix.extend([mpi_exec, "-n", str(mpi_n_proc)] + mpi_args)
    
    # Set the arguments for the command line
    args = mpi_prefix + [executable, "@", flagsfile]
    
    # Launch the process
    popen = subprocess.Popen(args,
                             stdout = open(output, "w"),
                             stderr = subprocess.STDOUT,
                             cwd = wd)
    
    # Wait for the process to complete
    popen.wait()
    
    # Return the Popen attributes of interest (cannot
    # return Popen itself since it is not serializable and
    # we are launching the process with Dask)
    return {"args" : popen.args,
            "stdin" : popen.stdin,
            "stdout" : popen.stdout,
            "stderr" : popen.stderr,
            "pid" : popen.pid,
            "returncode" : popen.returncode}


def check_rosetta_run(dirs_paths, **kwargs):
    """Check that a Rosetta calculation has exited without errors.
    """

    # Create an empty list to store the paths to the runs that
    # have crashed
    crashed_runs_paths = []

    # For each run to be checked
    for dir_path in dirs_paths:

        # If Rosetta wrote a crash log in the directory
        # where the run was performed
        if ROSETTA_CRASH_LOG in os.listdir(dir_path):

            # Append the path to the directory to the final list  
            crashed_runs_paths.append(dir_path)

    # Return the list of paths to the crashed runs
    return crashed_runs_paths


def parse_scorefile_text(scorefile, **kwargs):
    """Parse a Rosetta scorefile in text format to get the 
    structure numbers associated with the corresponding scores.
    """

    with open(scorefile, "r") as f:
        
        # Initialize the parsing flag to False
        parse = False
        
        # Initialize the structures' count to 1
        nstruct = 1
        
        # Initialize an empty list to store the scores
        scores = []
        
        # For each line in the file
        for l in f:
            
            if l.startswith("SCORE: total_score"):
                # Start parsing the file
                parse = True
                # Ignore the current line
                continue
            
            if parse:
                
                # Ignore empty elements in the string
                l = [item for item in l.split(" ") if item != ""]
                
                # The score is the second element
                scores.append((nstruct, float(l[1])))
                
                # Update the structure number (each line
                # corresponds to a structure)
                nstruct += 1
        
        # Return the scores associated to the corresponding
        # structures
        return scores


def write_flagsfile(options,
                    flagsfile,
                    **kwargs):
    """Write a flags file with the Rosetta options
    that will be used for the run.
    """

    # Make sure that the specified path exists.
    # If not, create it.
    file_path, file_name = os.path.split(flagsfile)
    os.makedirs(file_path, exist_ok = True)

    with open(flagsfile, "w") as f:
        
        # Format for the options (option and corresponding
        # value are whitespace-separated)
        opt_fmt = "{:s} {:s}\n"
        
        # Format for the variable substitutions in a RosettaScript
        # (variable and corresponding substitution are separated
        # by an equal sign)
        var_fmt = "{:s}={:s}"
        
        # For each option, option value pair
        for key, val in options.items():
            
            # If the option value is a dictionary, there
            # are some variables associated to the option
            if type(val) == dict:
                # Format the variables and their substitutions
                val = " ".join(\
                    [var_fmt.format(k,v) for k, v in val.items()])
            
            # Write the option to the flags file
            f.write(opt_fmt.format(key, val))
        
        # Return the path to the flags file
        return os.path.abspath(flagsfile)


def write_mutfile(mut,
                  mutfile,
                  **kwargs):
    """Write a mutfile containing the mutation performed
    in the current run (if any).
    """

    # Make sure that the specified path exists.
    # If not, create it.
    file_path, file_name = os.path.split(mutfile)
    os.makedirs(file_path, exist_ok = True)

    with open(mutfile, "w") as out:
        
        # Get the mutation
        mut = mut[MUT]
        
        # Get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR)
        
        # Set the header line
        out.write(f"total {len(mut)}\n{len(mut)}")
        
        # Write the mutations
        for smut in mut:
            
            # Get the values corresponding to the mutation
            # attributes
            chain, wtr, numr, mutr = operator.itemgetter(*keys)(smut)
            
            # Write the mutation
            out.write(f"\n{wtr} {numr} {mutr}")
        
        # Return the path to the mutfile
        return os.path.abspath(mutfile)


def write_resfile(mut,
                  resfile,
                  **kwargs):
    """Write a resfile containing the mutation performed
    in the current run (if any).
    """

    # Make sure that the specified path exists.
    # If not, create it.
    filepath, filename = os.path.split(resfile)
    os.makedirs(filepath, exist_ok = True)

    with open(resfile, "w") as out:
        
        # Get the mutation
        mut = mut[MUT]
        
        # Get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR)
        
        # Set the header and start line
        out.write("NATAA\nstart")         
        
        # For each single mutation (each single mutation is
        # in the form ("A", "C151Y"))
        for smut in mut:
            chain, wtr, numr, mutr = operator.itemgetter(*keys)(smut)
            out.write(f"\n{numr} {chain} PIKAA {mutr}")
        
        # Return the path to the resfile
        return os.path.abspath(resfile)



########################### CONFIG/OPTIONS ############################



def _convert_rosetta_option_to_string(option, value):
    """Convert a single Rosetta option and the associated value
    to a string representing the option/value pair.
    """
    
    # If the value is a boolean
    if isinstance(value, bool):
        # Rosetta uses 'true' and 'false' all lowercase
        return str(value).lower()
    
    # If the value is a string, integer or float
    elif isinstance(value, (int, str, float)):
        # Convert it to a string
        return str(value)
    
    # Otherwise, raise an exception
    else:
        raise ValueError(f"Unrecognized type {type(value)} for " \
                         f"value '{value}' of option '{option}'.")


def _recursive_traverse(data,
                        actions,
                        keys = None,
                        func = None):
    """Recursively traverse a dictionary performing actions on
    its items. It is used to traverse and modify the dictionary
    of options retrieved when parsing a YAML configuration file.

    Actions than can be performed are:
    
    - pop_empty : removal of keys associated to `None` 
                  values.
    - substitute : substitution of values of specific
                   keys with the result of a function taking
                   as inputs the key and the value.
    - substitute_dict : substitution of an entire dictionary with
                        a the return value of a function taking 
                        as arguments the items in the dictionary.
    """

    # If data is a dictionary
    if isinstance(data, dict):
        
        # Keys of items on which the actions will be
        # performed. If no keys are passed, all keys
        # in the dictionary will be considered.
        sel_keys = keys if keys else data.keys()
        
        # For each key, value pair
        for k, v in list(data.items()):

            # If value is None
            if v is None:
                
                # If removal of None values has been requested
                if "pop_empty" in actions:
                    
                    # Remove the key from the dictionary
                    data.pop(k)
                    continue

            # If value is a dictionary
            elif isinstance(v, dict):
                
                # If the susbtitution concerns the entire
                # dictionary
                if "substitute_dict" in actions:
                    
                    # If the key is in the selected keys
                    if k in sel_keys:
                        # Substitute the value with a function
                        # of the current value, with the function
                        # taking as inputs the key, values pairs
                        # in the dictionary
                        data[k] = func(**v)
                
                # Recursively check the sub-dictionaries
                # of the current dictionary
                _recursive_traverse(data = v,
                                    actions = actions,
                                    keys = keys,
                                    func = func)
        
            # If value is something else than a dictionary
            else:
                
                # If a substitution of the current value
                # has been requested
                if "substitute" in actions:
                    
                    # If the key is in the list of selected keys
                    if k in sel_keys:
                        
                        # Substitute the value with a function
                        # of it, assuming the function takes
                        # the key and the value as arguments
                        data[k] = func(k, v)

        # Return the dictionary
        return data


def get_option_key(options,
                   option,
                   **kwargs):
    """Given a dictionary of Rosetta options and the name of
    a particular option as defined in ROSETTA_OPTIONS, get which
    one of the possible alterative keys to define that option
    has been used in the given dictionary.
    """

    # Reset the worker's logger so that log messsages reach
    # the output
    logger = reset_worker_logger()

    # Get the values of each of the command line arguments that
    # can define the option of interest in the options dictionary
    # (it would be an empty string if not found)
    keys = set(ROSETTA_OPTIONS[option]) & set(options.keys())
    
    # If the user passed multiple equivalent command line
    # arguments for the option (i.e. both -in:file:s and -s), 
    # the longest one will be used
    list_keys = sorted(list(keys), key = len, reverse = True)
    
    # If the key exists
    if list_keys:
        if len(list_keys) > 1:
            warnstr = \
                f"Multiple instances of equivalent options " \
                f"found ({', '.join(list_keys)}). The more " \
                f"verbose one will be used ({list_keys[0]})."
            logger.warning(warnstr)
        
        # Return the longest option
        return list_keys[0]
    
    # Return None otherwise
    else:
        return None


def update_options(options,
                   pdb_file,
                   mut = None,
                   **kwargs):
    """Update a dictionary of Rosetta options with the input
    PDB file and possibly replace specific placeholders with
    the corresponding attribute of a mutation.
    """

    # Create a brand new dictionary to store the updated options
    new_options = copy.deepcopy(options)
    
    # Get one of the option keys used to define the
    # input PDB file (the longest one available by default)
    in_pdb_file_opt = sorted(ROSETTA_OPTIONS["in_pdb_file"],
                             key = len,
                             reverse = True)[0]

    # Update the new dictionary with the input PDB file option
    new_options[in_pdb_file_opt] = pdb_file

    # Possible RosettaScript in the protocol
    rosetta_script_opt = get_option_key(options, "protocol")
    rosetta_script = \
        options[rosetta_script_opt] if rosetta_script_opt else None

    # If a RosettaScript exists
    if rosetta_script:

        # If the RosettaScript is just a XML file name
        if os.path.basename(rosetta_script) == rosetta_script:
            # Assume it is a file in the default directory
            # containing RosettaScripts
            rosetta_script = os.path.join(ROSETTA_SCRIPTS_DIR,
                                          rosetta_script)
        
        # Otherwise, take the absolute path of the RosettaScript
        else:
            rosetta_script = os.path.abspath(rosetta_script)
        
        # Add the option with the path to the RosettaScript option
        new_options[rosetta_script_opt] = rosetta_script

        # If a mutation was passed (i.e. the protocol defined
        # includes a mutagenesis step)
        if mut:
            
            # Get the key in the option dictionary defining the
            # variable substitutions to be performed inside the
            # Rosetta script
            script_vars_opt = get_option_key(options, "script_vars")

            # Create a dictionary of the variable substitutions
            # with values and variables swapped (the keyword of the 
            # mutation attribute would be a value)
            val2var = \
                {v : k for k, v in options[script_vars_opt].items()}
            
            # For each attribute of the mutation
            for mut_key, mut_val in mut.items():
                
                # Get the name of the variable that will store the
                # mutation attribute, if present
                script_var = val2var.get(mut_key, None)
                
                # If the variable exists
                if script_var:
                    
                    # Susbtitute its current value (the mutation
                    # attribute keyword) with the mutation attribute
                    # value
                    new_options[script_vars_opt][script_var] = mut_val
    
    # Return the updated dictionary of options
    return new_options


def get_out_pdb_name(options,
                     pdb_file,
                     struct = None,
                     **kwargs):
    """Given a dictionary of Rosetta options, the name of
    an input PDB file that was the starting structure for the
    generation of an ensemble and possibly the number
    of the structure of interest in that ensemble, get the
    correct name given by Rosetta to that structure.
    """

    # Get the PDB file name
    pdb_name = pdb_file.rstrip(".pdb")
    
    # Get the option defining the PDB prefix, if any
    prefix_opt = get_option_key(options, "out_prefix")
    
    # Get the defined prefix or set it to an empty string
    # if not present
    prefix = options[prefix_opt] if prefix_opt else ""
    
    # Get the option defining the PDB suffix, if any
    suffix_opt = get_option_key(options, "out_suffix")
    
    # Get the defined prefix or set it to an empty string
    # if not present
    suffix = options[suffix_opt] if suffix_opt else ""
    
    # If a structure number was specified
    if struct:
        return f"{prefix}{pdb_name}{suffix}_{struct:04d}.pdb"
    
    # If no structure number was specified
    else:
        return f"{prefix}{pdb_name}{suffix}.pdb"


def _get_config_run_version_1(config):
    """Get the configuration from version 1 YAML 
    configuration files.
    """

    # Get the protocol family
    family = config["family"]

    # If the protocol family is not recognized, raise
    # an error
    if not family in ROSETTA_PROTOCOLS.keys():
        errstr = f"Unrecognized protocol family {family}."
        raise ValueError(errstr)

    # For each step
    for step_name, step in config["steps"].items():
        
        # If the step name is not recognized, raise an error
        if not step_name in ROSETTA_PROTOCOLS[family].keys():
            errstr = \
                f"Unrecognized step name {step_name} " \
                f"for protocol family {family}."
            raise ValueError(errstr)

        # Get the step's fixed settings
        step_settings = ROSETTA_PROTOCOLS[family][step_name]
        
        # Only consider Rosetta steps
        if step_settings["run_by"] == "rosetta":
            
            # Create a copy of the configuration
            step_opts = dict(step["options"])
            
            # Recursively remove all options that map
            # to None and convert the other ones to strings
            new_step_opts = \
                _recursive_traverse(\
                    data = step_opts,
                    actions = ["pop_empty", "substitute"],
                    func = _convert_rosetta_option_to_string)
            
            # Update the dictionary of options with the new
            # options for the step
            config["steps"][step_name]["options"] = new_step_opts

            # If the extraction of the structures is part of
            # the step
            if "extract_structures" in list(step.keys()):
                
                # Create a copy of the configuration
                extract_opts = \
                    dict(step["extract_structures"]["options"])

                # Recursively remove all options that map
                # to None and convert the other ones to strings
                new_extract_opts = \
                    _recursive_traverse(\
                        data = extract_opts,
                        actions = ["pop_empty", "substitute"],
                        func = _convert_rosetta_option_to_string)
                
                # Update the dictionary of options with the new
                # options for the step
                config["steps"][step_name]["extract_structures"]["options"] = \
                    new_extract_opts

    # Return the configuration
    return config


def get_config_run(config_file, **kwargs):
    """Get the configuration for running the protocol.
    """

    # Get the name of the configuration file for running
    # the protocol
    config_file_name = os.path.basename(config_file).rstrip(".yaml")

    # If the configuration file is a name without extension
    if config_file == config_file_name:
        
        # Assume it is a configuration file in the directory
        # storing configuration files for running protocols
        config_file = os.path.join(CONFIG_RUN_DIR,
                                   config_file_name + ".yaml")
    
    # Otherwise, assume it is a file name/file path
    else:
        config_file = get_abspath(config_file)

    # Load the configuration from the file
    config = yaml.safe_load(open(config_file, "r"))

    # Check the version of the configuration file
    if config["version"] == 1:
        
        # Return the configuration written in version 1 format
        return _get_config_run_version_1(config = config)
    
    # Only version 1 supported so far
    else:
        errstr = \
            "Only version 1 configuration files " \
            "are supported for now."
        raise ValueError(errstr)


def get_config_settings(config_file, **kwargs):
    """Get the configuration for the run's settings.
    """

    # Get the name of the configuration file for the
    # run's settings
    config_file_name = os.path.basename(config_file).rstrip(".yaml")

    # If the configuration file is a name without extension
    if config_file == config_file_name:
        
        # Assume it is a configuration file in the directory
        # storing configuration files for run settings
        config_file = os.path.join(CONFIG_SETTINGS_DIR,
                                   config_file_name + ".yaml")
    
    # Otherwise, assume it is a file name/file path
    else:
        config_file = get_abspath(config_file)

    # Return the configuration
    return yaml.safe_load(open(config_file, "r"))


def get_config_aggregate(config_file, **kwargs):
    """Get the configuration for data aggregation.
    """
    
    # Get the name of the configuration file for data
    # aggregation
    config_file_name = os.path.basename(config_file).rstrip(".yaml")

    # If the configuration file is a name without extension
    if config_file == config_file_name:
        
        # Assume it is a file in the directory where
        # configuration files for data aggregation are stored
        config_file = os.path.join(CONFIG_AGGR_DIR,
                                   config_file_name + ".yaml")
    
    # Otherwise, assume it is a file name/file path
    else:
        config_file = get_abspath(config_file)

    # Return the configuration
    return yaml.safe_load(open(config_file, "r"))


def get_config_plot_version_1(config):
    """Parse the configuration file for plotting,
    version 1.
    """
    
    # Create a copy of the configuration
    new_config = dict(config)
    
    # Substitute the font properties definitions
    # with the corresponding FontProperties instances
    _recursive_traverse(data = new_config,
                        actions = ["substitute_dict"],
                        func = fm.FontProperties,
                        keys = {"fontproperties"})
    
    # Return the configuration
    return new_config


def get_config_plot(config_file, **kwargs):
    """Get the plotting configuration.
    """

    # Get the name of the configuration file for plotting
    config_file_name = os.path.basename(config_file).rstrip(".yaml")

    # If the configuration file is a name without extension
    if config_file == config_file_name:
        
        # Assume it is a file in the directory where
        # configuration files for plotting are stored
        config_file = os.path.join(CONFIG_PLOT_DIR,
                                   config_file_name + ".yaml")
    
    # Otherwise, assume it is a file name/file path
    else:
        config_file = get_abspath(config_file)
    
    # Load the configuration from the file
    config = yaml.safe_load(open(config_file, "r"))

    # Check the version of the configuration file
    if config["version"] == 1:
        # Return the configuration written in version 1 format
        return _get_config_plot_version_1(config = config)
    
    # Only version 1 is supported so far
    else:
        errstr = \
            "Only version 1 configuration files " \
            "are supported for now."
        raise ValueError(errstr)



########################## MUTATIONS-RELATED ##########################



def get_res_list(res_list_file):
    """Parse the file containing the list of residue types
    to be used for saturation mutagenesis scans.
    """
            
    with open(res_list_file, "r") as f:
        
        # Return all lines that are not empty
        return [l.rstrip("\n") for l in f \
                if not re.match(r"^\s*$", l)]


def _get_mut_list(list_file):
    """Parse the file containing the list of 
    positions/mutations.
    """

    # Reset the worker's logger so that log messsages reach
    # the output
    logger = reset_worker_logger()

    with open(list_file, "r") as f:
        
        # Initialize an empty list to store the mutations
        mut_list = []
        
        # For each line in the file
        for line in f:
            
            # Ignore empty lines
            if re.match(r"^\s*$", line):
                continue
            
            # Separate the mutations from extra data associated
            # to them     
            data, *extra_data = line.rstrip("\n").split(" ")
            
            # Create an empty list to store the mutation
            mut = []
            
            # For each single mutation potentially part of
            # a multiple mutation
            for single_mut in data.split(MUT_SEP):
                mut.append(tuple(single_mut.split(COMP_SEP)))
            
            # Convert the list to a tuple
            mut = tuple(mut)
            
            # Convert each piece of extra data to a string
            extra_data = [str(item) for item in extra_data]

            # Pack the mutation data into a tuple
            mut_data = tuple([(mut), *extra_data])
            
            # If the mutation has been already found in the file
            if mut_data in mut_list:

                _mut = line.rstrip('\n').split(' ')
                
                # Warn the user that the mutation will be performed
                # only once
                warnstr = \
                    f"Mutation '{_mut}' has been defined more than " \
                    f"once in {list_file}. It will be performed " \
                    f"only once."
                
                logger.warning(warnstr)

            # Otherwise, append the mutation data to the list
            else:
                mut_list.append(mut_data)

        # Return the list of mutations
        return mut_list



def _get_saturation_mut_list(pos_list, res_list):
    """Generate a mutation list for saturation mutagenesis.
    
    Every position specified as (chain, wild_type_residue,
    position) will be treated as a position to perform
    saturation mutagenesis on, while every mutation specified
    as (chain, wild_type_residue, position, mutated_residue)
    will be treated as a single mutation.
    
    If you define multiple positions to be simultaneously
    mutated (i.e., on the same line of the mutations' list file),
    all possible combinations of mutations of those positions
    will be performed (Cartesian product).
    
    Example
    -------
    If in the mutations' list file you have:

    A.R.10,A.R.11,A.F.52.A

    and in the residues' list file you have:

    A
    C

    The following mutations will be performed:

    A.R.10.A and A.R.11.A and F.52.A
    A.R.10.A and A.R.11.C and F.52.A
    A.R.10.C and A.R.11.A and F.52.A
    A.R.10.C and A.R.11.C and F.52.A
    """

    # Create an empty list to store the list of mutations
    # when saturation mutagenesis is involved
    sat_mut_list = []
    
    # For each combination of mutations/positions
    for pos_data, *extra_data in pos_list:

        # If there is only one mutation/postion
        if len(pos_data) == 1:

            # If it is a mutation (chain, wild-type residue,
            # position, and mutated residue specified)
            if len(pos_data[0]) == 4:

                # Just add the mutation to the list (it will
                # not be part of the saturation mutagenesis)
                sat_mut_list.append((pos_data, *extra_data))
            
            # If it is a position (chain, wild-type residue,
            # position, and mutated residue specified)
            elif len(pos_data[0]) == 3:
                
                # Get the chain, wild-type residue, and position
                chain, wtr, numr = pos_data[0]

                # Add all possible mutations for the current
                # position to the list
                sat_mut_list.extend(\
                    [(((chain, wtr, numr, res),), *extra_data) \
                     for res in res_list])

        # If there are multiple mutations/positions
        else:

            # For each position, generate all possible mutations.
            # For each mutation, simply add it.
            single_muts = \
                [[(*i, r) if len(i) == 3 else tuple(i) \
                 for r in res_list] for i in pos_data]

            # Compute the Cartesian product between the different
            # possibilities to obtain the final list of mutations
            sat_mut_list.extend(\
                ([(i, *extra_data) for i in \
                  list(itertools.product(*single_muts))]))

    # Return the list of mutations
    return sat_mut_list


def _convert_to_pose_numbering(mut_list, pdb_file):
    """Retun a copy of the mutation list with residue numbers
    changed to the Rosetta pose numbering.
    """

    # Create the PDB parser
    parser = PDB.PDBParser()
    
    # Get the structure fron the PDB file
    structure = parser.get_structure("structure", pdb_file) 
    
    # Create an empty dictionary to store the mapping
    # between the PDB numbering and the pose numbering
    pdbnum2posenum = {}
    
    # Set a chain offset to start at the correct index when
    # changing chain
    chain_offset = 0
    
    # For each chain
    for chain in structure[0]:
        
        # For each residue in the chain
        for ix, res in enumerate(chain, start = 1):
            
            # Unpack the residue full ID
            struc, mod, chainid, (het, resn, icode) = res.get_full_id()
            
            # Map the PDB chain and residue number to the numbering
            # that the residue will have in the mutfile (starts
            # at 1); convert the PDB residue number to a string
            pdbnum2posenum[(chainid, str(resn))] = \
                str(ix + chain_offset)
        
        # Update the chain offset
        chain_offset += ix
    
    # Generate the mutation list but with residue numbers changed
    # according to the mutfile convention
    pose_mut_list = []

    # For each mutation and associated extra data
    for mut, *extra_data in mut_list:

        # Get the mutation with the pose numbering 
        new_mut = \
            [(chain, wtr, pdbnum2posenum[(chain, numr)], mutr) \
              for (chain, wtr, numr, mutr) in mut]

        # Append it to the final list together with its associated
        # extra data
        pose_mut_list.append(tuple([tuple(new_mut), *extra_data]))
    
    # Return the new mutations' list
    return pose_mut_list


def _generate_mutation_dir_path(mut_dict):
    """Given a mutation, generate the corresponding directory
    name (since the Rosetta-compatible name of the mutation may
    not be usable as-it-is as a directory name).
    """

    # Get the keys of the attributes of the mutation
    keys = (CHAIN, WTR, NUMR, MUTR)
    
    # Initialize a list for the formatted mutation names
    fmt_muts = []      
    
    # Each single mutation is in the form ("A", "C151Y")
    for single_mut in mut_dict[MUT]:
        
        # Get the mutation attributes
        chain, wtr, numr, mutr = \
            operator.itemgetter(*keys)(single_mut)
        
        # Strip Rosetta identifiers of non canonical residues
        # from the name of the residues
        wtr = wtr.strip("X[]")
        mutr = mutr.strip("X[]")

        # If the chain has a chain ID      
        if chain != "_":
            fmt_muts.append(f"{chain}{CHAIN_SEP}{wtr}{numr}{mutr}")

        # If the chain has no ID       
        else:
            fmt_muts.append(f"{wtr}{numr}{mutr}")
    
    # Get the string representing the mutation directory name
    dir_name = DIR_MUT_SEP.join(fmt_muts)
    
    # If the run generates multiple structures
    if STRUCT in mut_dict.keys():
        # Add another level to the path, so that each
        # structure will have a dedicated sub-folder
        dir_path = os.path.join(dir_name, mut_dict[STRUCT])
    
    # If the run generates only one structure
    else:
        # There will be no sub-folders for the structures
        dir_path = dir_name
    
    # Return the directory path and name
    return dir_path, dir_name


def get_mutations(list_file,
                  res_list_file,
                  pdb_file,
                  res_numbering,
                  extra,
                  n_struct,
                  **kwargs):
    """Get the list of mutations to be performed.
    """

    # Check the pose numbering argument
    if not res_numbering in ("pose", "pdb"):
        errstr = \
            f"Residue numbering must be either 'pose' or " \
            f"'pdb', but '{res_numbering}' was passed."
        raise ValueError(errstr)


    #---------------------- Saturation or not? -----------------------#


    # Get the mutations/positions list
    mut_list = _get_mut_list(list_file)

    # If a list of residue types has been passed, assume it is a
    # saturation mutagenesis scan
    if res_list_file:

        # Get the list of residue types
        res_list = get_res_list(res_list_file)
        
        # Treat the mutations' list as a list of positions
        # and generate the new list of mutations
        mut_list = _get_saturation_mut_list(mut_list, res_list)

    # Generate a copy of the original mutation list (in case it
    # gets changed because of a different residue numbering)
    mut_list_original = copy.deepcopy(mut_list)


    #----------------------- Residue numbering -----------------------#


    # If the protocol requires the residue numbering to follow
    # the Rosetta pose numbering convention
    if res_numbering == "pose":
        
        # Convert the numbering from PDB to pose numbering
        mut_list = _convert_to_pose_numbering(mut_list, pdb_file)

    # Otherwise, assume it follows the PDB convention (requires no
    # conversion)


    #------------------------ List generation ------------------------#


    # Create an empty list to store the mutations
    mutations = []
    
    # For each mutations' list (original and final)
    for mutl in [mut_list, mut_list_original]:

        # Initialize an empty list to store the single mutations
        # defined in each mutation
        muts = []

        # ((("A","C","151","Y"), ("A","S","154","N")), *extra_data)
        for mut, *extra_data in mutl:

            # Create a dictionary where the mutation's attributes
            # will be stored
            mut_dict = {MUT : []}
            
            # For each mutation attribute
            for chain, wtr, numr, mutr in mut:
                
                # If the wild-type residue is noncanonical
                if len(wtr) > 1:
                    wtr = f"X[{wtr}]"
                
                # If the mutant residue is noncanonical
                if len(mutr) > 1:
                    mutr = f"X[{mutr}]"             
                
                # Update the mutation dictionary
                mut_dict[MUT].append({CHAIN : chain,
                                      WTR : wtr,
                                      NUMR : numr,
                                      MUTR : mutr})
            
            # If there are extra data associated with the mutation      
            if extra:
                
                # Add them with the corresponding key
                for key, data in zip(extra, extra_data):
                    mut_dict[key] = data
            
            # If multiple structures will be generatred for each mutation
            if n_struct:
                
                # Add the mutation as many times as the structures
                # requested, with a new key identifying the structure
                # number
                for struct in range(n_struct):
                    
                    # Create a fresh copy of the dictionary to be modified
                    mut_dict_copy = dict(mut_dict)
                    
                    # Update the copy
                    mut_dict_copy.update({STRUCT : str(struct + 1)})
                    
                    # Append the copy to the list of mutations
                    muts.append(mut_dict_copy)

            # If only one structure will be generated
            else:

                # Simply append the mutation dictionary
                muts.append(mut_dict)

        # For each mutation dictionary
        for mut_dict in muts:
            
            # Add the mutation directory name and directory path
            mut_dir_path, mut_dir_name = \
                _generate_mutation_dir_path(mut_dict)
            mut_dict[MUT_DIR_PATH] = mut_dir_path
            mut_dict[MUT_DIR_NAME] = mut_dir_name

        # Append the current list of mutations to the final list
        mutations.append(muts)
   
    # Return the list of mutations
    return mutations


def write_mutinfo_file(mutations_original,
                       out_dir,
                       mutinfo_file,
                       **kwargs):
    """Write a comma-separated file containing the names of the
    directories containing the results for all mutations, the
    names of the mutations as written in the mutations list file
    and the labels to be used for the mutations when aggregating
    data and plotting.
    """

    # Reset the worker's logger so that log messsages reach
    # the output
    logger = reset_worker_logger()

    # Get the name of the mutinfo file without the extension
    mutinfo_filename, mutinfo_ext = os.path.splitext(mutinfo_file)

    # If the output directory has already been created
    if os.path.exists(out_dir):

        # Check if a mutinfo file with the name specified in
        # the config file is already present (or others
        # already numbered)        
        mutinfo_files = \
            [f for f in os.listdir(out_dir) if f == mutinfo_file \
             or re.match(f"{mutinfo_filename}[0-9]{mutinfo_ext}", f)]

        # If (an)other mutinfo file(s) was (were) found
        if mutinfo_files:

            # Get the file numbers
            num_files = \
                [f.lstrip(mutinfo_filename).rstrip(mutinfo_ext) \
                 for f in mutinfo_files]

            # Sort the files in ascending order
            files_sorted = \
                sorted(list(zip(num_files, mutinfo_files)),
                       key = lambda x: x[0])

            # Set the name of the new mutinfo file, given the files
            # already found (name of the original mutinfo file name +
            # higherst file number found+1 + original mutinfo file
            # extension
            mutinfo_file = \
                f"{mutinfo_filename}{int(files_sorted[-1][0])+1}" \
                f"{mutinfo_ext}"

            # Warn the user about thw existing mutinfo files
            # and the new one that will be written
            warnstr = \
                f"The following mutinfo files have been found in " \
                f"{out_dir}: " \
                f"{', '.join([i[1] for i in files_sorted])}. A new " \
                f"file named {mutinfo_file} will be written for " \
                f"this run."
            logger.warning(warnstr)

    # Make sure that the directory exists. If not, create it.
    os.makedirs(out_dir, exist_ok = True)

    # Set the path to the output file
    mutinfo_file_path = os.path.join(out_dir, mutinfo_file)

    with open(mutinfo_file_path, "w") as out:
        
        # Keep track of the mutations already written
        # to the mutinfo file
        muts_written = set()

        # Get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR)
        
        # For each mutation in the original list
        for mut_orig in mutations_original:

            # Generate empty lists to store the single
            # mutations (to be joined later) composing each
            # mutation to be performed
            mutation_def = []
            mut_label_def = []
            pos_label_def = []

            # Hashable name for the mutation
            mut_hash = tuple([tuple(i.items()) for i in mut_orig[MUT]])

            # For each single mutation
            for i, single_mut_orig in enumerate(mut_orig[MUT]):

                # Get the attributes of the original mutation
                chain_orig, wtr_orig, numr_orig, mutr_orig = \
                    operator.itemgetter(*keys)(single_mut_orig)

                # Strip Rosetta identifiers of non canonical
                # residues from the name of the residues (so
                # that it is consistent with what it is written
                # in the mutations' list file)
                wtr_orig = wtr_orig.strip("X[]")
                mutr_orig = mutr_orig.strip("X[]")

                # Compose the mutation name            
                mutation_def.append(\
                    f"{chain_orig}{COMP_SEP}{wtr_orig}{COMP_SEP}" \
                    f"{numr_orig}{COMP_SEP}{mutr_orig}")
                    
                # Compose the mutation label (by default
                # without the chain ID)
                mut_label_def.append(\
                    f"{wtr_orig}{numr_orig}{mutr_orig}")
                    
                # Compose the position label (by default
                # without the chain ID)
                pos_label_def.append(\
                    f"{wtr_orig}{numr_orig}")
                
                # If the mutation has not been written yet
                if i == len(mut_orig[MUT])-1 \
                and mut_hash not in muts_written:
                
                    # Write out the directory name, the mutation
                    # and the label that will be used for the
                    # mutation in the aggregation/plot
                    out.write(f"{MULTI_MUT_SEP.join(mutation_def)},"
                              f"{mut_orig[MUT_DIR_NAME]},"
                              f"{DIR_MUT_SEP.join(mut_label_def)},"
                              f"{DIR_MUT_SEP.join(pos_label_def)}\n")

                    muts_written.add(mut_hash)
                


def get_mutinfo(mutinfo_file, **kwargs):
    """Create a data frame from a comma-separated file containing 
    the names of the directories containing the results for all 
    mutations, the names of the mutations as written in the 
    mutations list file and the labels to be used for the 
    mutations when aggregating data and plotting.
    """

    with open(mutinfo_file, "r") as f:
        
        # Create an empty list to store directory 
        # names and mutations
        mutinfo = []
        
        # For each line in the file
        for l in f:
            
            # Ignore empty lines
            if re.match(r"^\s*$", l):
                continue

            # Get mutation, directory name and mutation name
            mut_name, dir_name, mut_label, pos_label = \
                l.rstrip("\n").split(",")
            
            # Update the dictionary corresponding to
            # the current mutation
            mutinfo.append(\
                {MUTINFO_COLS["mut_name"] : mut_name,
                 MUTINFO_COLS["dir_name"] : dir_name,
                 MUTINFO_COLS["mut_label"] : mut_label,
                 MUTINFO_COLS["pos_label"] : pos_label})

        # Create a data frame from the list and return it
        return pd.DataFrame(mutinfo)



################################ OTHERS ###############################



def check_pdb_file(pdb_file,
                   allow_multi_chains,
                   allow_no_chain_ids,
                   **kwargs):
    """Check a PDB file before passing it to Rosetta.
    """
    
    # Create the PDB parser
    parser = PDB.PDBParser()
    
    # Try to get the structure
    try:
        structure = parser.get_structure("structure", pdb_file)
    
    # In case the PDB file was not found
    except FileNotFoundError:
        
        # Raise an error
        errstr = f"PDB file {pdb_file} not found."
        raise FileNotFoundError(errstr)
    
    # In case something went wrong in accessing/opening the file
    except IOError:
        
        # Raise an error
        errstr = f"Could not open PDB file {pdb_file}."
        raise IOError(errstr)
    
    # In case something else went wrong
    except Exception as e:
        
        # Re-raise the error
        raise Exception(e)

    # If the PDB file contains more than one model
    if len(structure) > 1: 
        
        # Raise an exception since multi-model structures
        # are not allowed
        errstr = "Multi-model PDB files are not allowed."
        raise ValueError(errstr)
    
    # Get the number of chains in the structure
    num_chains = len(structure[0])
    
    # Check if multiple chains are allowed
    if not allow_multi_chains and num_chains > 1:
        
        # Raise an exception if the structure contains multiple
        # chains but multiple chains are not allowed
        errstr = "Multi-chains structures are not allowed."
        raise ValueError(errstr)
    
    # Check if the absence of chain IDs is allowed
    if not allow_no_chain_ids:
        
        # For each chain in the first model of the structure (it
        # is already guaranteed that there is only one model)
        for chain in structure[0]:
            
            # Raise an exception if there are no chain IDs but
            # chain IDs were mandatory
            if not chain.get_id():
                errstr = "All chains must have a chain ID."
                raise ValueError(errstr)
    
    # Return the PDB file
    return pdb_file


def get_abspath(path, **kwargs):
    """Given a path, return its absolute path. Return
    None if the path given is None.
    """

    return os.path.abspath(path) if path is not None else path


def get_items(d,
              keys,
              default = None,
              **kwargs):
    """Similar to operator.itemgetter but defaults to
    a specific value if the key is not found (instead of
    throwing an exception).
    """

    return [d.get(k, default) for k in keys]


def check_structures_extraction(wd):
    """Check whether the extraction of the structures or their
    renaming has already happened. 
    """

    # Get the files in the directory
    files_in_wd = os.listdir(wd)

    # Get the structures extracted from the database file
    struct_extracted = \
        list(filter(lambda x: re.match(STRUCT_EXTRACTED_PATTERN, x),
                    files_in_wd))

    # If a structure whose name begins with the given prefix has been
    # found, the extraction has been performed (all structures
    # in the folder are extracted at once)
    is_struct_extracted = True if len(struct_extracted) > 0 else False

    # Get the structures renamed
    struct_renamed = \
        list(filter(lambda x: x.startswith(FLEXDDG_STATES) \
                    and x.endswith(".pdb"),
                    files_in_wd))

    # If a structure whose name begins with any of the prefixes
    # has been found, the renaming has been performed (all structures
    # in the folder are renamed at once)
    is_struct_renamed = True if len(struct_renamed) > 0 else False

    # If the flag for the structures' renaming is on
    if is_struct_renamed:

        # If also the flag for structures' extraction is on at
        # this point (which means some structures that have
        # been extracted but not renamed are left for some reason)
        if is_struct_extracted:

            # Turn off the is_struct_extracted flag
            is_struct_renamed = False

        # Otherwise (no leftover structure to rename)
        else:

            # Turn on the flag for the structures' extractuon
            # (since none of the structures has their original
            # name after renaming, the is_struct_extracted flag
            # would be now False even if the
            # structures have been in fact extracted)
            is_struct_extracted = True

    # Return the two flags reporting the status of the extraction
    return is_struct_extracted, is_struct_renamed


def rename_structures_flexddg(path,
                              r_script_options,
                              **kwargs):
    """Rename the structures extracted from the db3 file
    produced by the flexddg protocols.
    """

    # Reset the worker's logger so that log messsages reach
    # the output
    logger = reset_worker_logger()

    # Get the number of backrub trials
    backrub_n_trials = \
        r_script_options[get_option_key(options = r_script_options,
                                        option = "backrub_n_trials")]
                            
    # Get the backrub trajectory stride
    backrub_traj_stride = \
        r_script_options[get_option_key(options = r_script_options,
                                       option = "backrub_traj_stride")]

    # Compute the trajectory stride
    traj_stride = int(backrub_n_trials) // int(backrub_traj_stride)

    # Create an empty list to store the structures to be renamed
    # (for logging purposes)
    struct_to_rename = []

    # Create an empty list to store the structures that have been
    # renamed (for logging purposes)
    struct_renamed = []

    # For each file in the specified path
    for f in os.listdir(path):

        # Split the name of the file into path leading up to
        # the file and file base name
        f_path, f_basename = os.path.split(os.path.join(path, f))

        # Get the match between the pattern that we look for
        # the name of the structure
        struct_match = re.match(STRUCT_EXTRACTED_PATTERN, f_basename)
        
        # If the current file is a structure to be renamed
        if struct_match:

            # Get the structure ID
            struct_id = int(struct_match.group(1))

            # Set the new name for the structure
            new_fname = \
                "%s_%05d.pdb" % \
                    (FLEXDDG_STATES[(struct_id-1) % len(FLEXDDG_STATES)],
                    (((struct_id-1) // len(FLEXDDG_STATES)) + 1) \
                    * traj_stride)

            # Rename the structure
            os.rename(os.path.join(f_path, f_basename), \
                      os.path.join(f_path, new_fname))

            # Add the old name to the list of structures to
            # be renamed and the new name in the list of
            # structures renamed
            struct_to_rename.append(f_basename)
            struct_renamed.append(new_fname)

    # Inform the user about which structures were found and renamed
    infostr_to_rename = \
        f"The following structures (extracted from the protocol's " \
        f"output database) were found and will be renamed:" \
        f"\n{', '.join(struct_to_rename)}."
    logger.info(infostr_to_rename)

    infostr_renamed = \
        f"The following structures have been renamed: " \
        f"\n{', '.join(struct_renamed)}."
    logger.info(infostr_renamed)