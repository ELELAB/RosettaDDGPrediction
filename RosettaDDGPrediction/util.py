#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    General utility functions.
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
import copy
import operator
import os
import os.path
import re
import subprocess
# third-party packages
import Bio.PDB as PDB
import pandas as pd
import yaml
# RosettaDDGProtocols
from .dask_patches import reset_worker_logger
from .defaults import (
    CHAIN,
    COMPSEP,
    DIRCHAINSEP,
    DIRMUTSEP,
    MUT,
    MUTDIRNAME,
    MUTDIRPATH,
    MUTR,
    MUTSEP,
    NOMUTR, 
    NUMR,
    ROSETTADFCOLS,
    ROSETTAOPTIONS,
    ROSETTAPROTOCOLS,
    ROSETTASCRIPTSDIR,
    STRUCT,
    WTR
)



########################### ROSETTA-RELATED ###########################



def get_rosetta_executable(execname, execpath, execsuffix):
    """Get the path to a Rosetta executable from the Rosetta
    installation directory.
    """

    # reset the worker logger
    logger = reset_worker_logger()

    # get all the executables available
    allexecs = os.listdir(execpath)
    
    # if no suffix was provided, assume the executable name does not
    # have any suffix
    execcompname = execname + execsuffix if execsuffix is not None \
                   else execname

    # filter function to get the Rosetta executables
    filt = lambda x: x == execcompname
    
    # get all Rosetta executables satifying the criteria of
    # the filter function (you should get only one executable) 
    execs = list(filter(filt, allexecs))
    
    # if no executable was found
    if len(execs) == 0:
        # raise an exception
        raise ValueError(f"No executable found for {execcompname}.")
    # if multiple conflicting executables were found
    elif len(execs) > 1:
        # raise an exception
        raise ValueError(f"Multiple executables found for " \
                         f"{execcompname}.")

    # return the path to the executable
    return os.path.join(execpath, execs[0])


def run_rosetta(executable, \
                flagsfile, \
                output, \
                wd, \
                usempi, \
                mpiexec, \
                mpiargs, \
                mpinproc):
    """Run Rosetta."""

    # make sure that the directory exists. If not, create it.
    os.makedirs(wd, exist_ok = True)
    
    # set an empty prefix (before the Rosetta command line)
    # to run with MPI (remains empty if you do not run with MPI)
    mpiprefix = []
    # if MPI is requested
    if usempi:
        # update the MPI prefix
        mpiprefix.extend([mpiexec, "-n", str(mpinproc)] + mpiargs)
    # set the arguments for the command line
    args = mpiprefix + [executable, "@", flagsfile]
    
    # launch the process
    popen = subprocess.Popen(args, \
                             stdout = open(output, "w"), \
                             stderr = subprocess.STDOUT, \
                             cwd = wd)
    
    # wait for the process to complete
    popen.wait()
    
    # return the Popen attributes of interest (cannot
    # return Popen itself since it is not serializable and
    # we are launching the process with Dask)
    return {"args" : popen.args, \
            "stdin" : popen.stdin, \
            "stdout" : popen.stdout, \
            "stderr" : popen.stderr, \
            "pid" : popen.pid, \
            "returncode" : popen.returncode}


def parse_scorefile_text(scorefile):
    """Parse a Rosetta scorefile in text format to get the 
    structure numbers associated with the corresponding scores.
    """

    with open(scorefile, "r") as f:
        # initialize the parsing flag to False
        parse = False
        # initialize the structures count to 1
        nstruct = 1
        # initialize an empty list to store the scores
        scores = []
        for l in f:
            if l.startswith("SCORE: total_score"):
                # start parsing the file
                parse = True
                continue
            if parse:
                # ignore empty elements in the string
                l = [item for item in l.split(" ") if item != ""]
                # the score is the second element
                scores.append((nstruct, float(l[1])))
                # update the structure number (each line
                # corresponds to a structure)
                nstruct += 1
        
        # return the scores
        return scores


def write_flagsfile(options, flagsfile):
    """Write a flags file with the Rosetta options
    that will be used for the run.
    """

    # make sure that the specified path exists.
    # If not, create it.
    filepath, filename = os.path.split(flagsfile)
    os.makedirs(filepath, exist_ok = True)

    with open(flagsfile, "w") as f:
        # format for the options (option and corresponding
        # value are whitespace-separated)
        optfmt = "{:s} {:s}\n"
        # format for the variable substitutions in a RosettaScript
        # (variable and corresponding substitution are separated
        # by an equal sign)
        varfmt = "{:s}={:s}"
        # for each option, option value
        for key, val in options.items():
            # if the option value is a dictionary, there
            # are some variables
            if type(val) == dict:
                # format the variables and their substitutions
                val = " ".join(\
                    [varfmt.format(k,v) for k, v in val.items()])
            # write the option to the flags file
            f.write(optfmt.format(key, val))
        # return the path to the flags file
        return os.path.abspath(flagsfile)


def write_mutfile(mut, mutfile):
    """Write a mutfile containing the mutation performed
    in the current run (if any).
    """

    # make sure that the specified path exists.
    # If not, create it.
    filepath, filename = os.path.split(mutfile)
    os.makedirs(filepath, exist_ok = True)

    with open(mutfile, "w") as out:
        # get the mutation
        mut = mut[MUT]
        # get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR)
        # set the header line
        out.write(f"total {len(mut)}\n{len(mut)}")
        # write the mutations
        for smut in mut:
            # get the values corresponding to the mutation
            # attributes
            chain, wtr, numr, mutr = operator.itemgetter(*keys)(smut)
            # write the mutation
            out.write(f"\n{wtr} {numr} {mutr}")
        # return the path to the mutfile
        return os.path.abspath(mutfile)


def write_resfile(mut, resfile):
    """Write a resfile containing the mutation performed
    in the current run (if any).
    """

    # make sure that the specified path exists.
    # If not, create it.
    filepath, filename = os.path.split(resfile)
    os.makedirs(filepath, exist_ok = True)

    with open(resfile, "w") as out:
        # get the mutation
        mut = mut[MUT]
        # get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR)
        # set the header and start line
        out.write("NATAA\nstart")         
        # each single mutation is in the form ("A", "C151Y")
        for smut in mut:
            chain, wtr, numr, mutr = operator.itemgetter(*keys)(smut)
            out.write(f"\n{numr} {chain} PIKAA {mutr}")
        # return the path to the resfile
        return os.path.abspath(resfile)



########################### OPTIONS-RELATED ###########################



def convert_rosetta_option_to_string(option, value):
    """Convert a single Rosetta option to a string."""
    
    # if the value is a boolean
    if isinstance(value, bool):
        # Rosetta uses 'true' and 'false' all lowercase
        return str(value).lower()
    
    # if the value is a string, integer or float
    elif isinstance(value, (int, str, float)):
        # convert it to a string
        return str(value)
    
    # otherwise, raise an exception
    else:
        raise ValueError(f"Unrecognized type {type(value)} for " \
                         f"value '{value}' of option '{option}'.")


def recursive_traverse(data, \
                       actions, \
                       keys = None, \
                       func = None):
    """Recursively traverse a dictionary performing actions on
    its items. It is used to traverse and modify the dictionary
    of options retrieved when parsing a YAML configuration file.

    Actions than can be performed are:
    
    - pop_empty : removal of keys associated to `None` 
                  values.
    - substitute : substitution of values of specific
                   keys with a function of those values.
    - substitute_dict : substitution of an entire dictionary with
                        a the return value of a function taking 
                        as arguments the items in the dictionary.
    """

    # if data is a dictionary
    if isinstance(data, dict):
        
        # keys of items on which the actions will be
        # performed. If no keys are passed, all keys
        # in the dictionary will be considered.
        selkeys = keys if keys else data.keys()
        
        # for each key, value pair
        for k, v in list(data.items()):

            # if value is None
            if v is None:
                # if removal of None values has been requested
                if "pop_empty" in actions:
                    # remove the key from the dictionary
                    data.pop(k)
                    continue

            # if value is a dictionary
            elif isinstance(v, dict):
                # if the susbtitution concerns the entire
                # dictionary
                if "substitute_dict" in actions:
                    # if the key is in the selected keys
                    if k in selkeys:
                        # substitute the value with a function
                        # of the current value, with the function
                        # taking as inputs the key, values pairs
                        # in the dictionary
                        data[k] = func(**v)
                
                # recursively check the sub-dictionaries
                # of the current dictionary
                recursive_traverse(data = v, \
                                   actions = actions, \
                                   keys = keys, \
                                   func = func)
        
            # if value is something else
            else:
                # if substitution of the current value
                # has been requested
                if "substitute" in actions:
                    # if the key is in the list of selected keys
                    if k in selkeys:
                        # substitute the value with a function
                        # of it, assuming the function takes
                        # the key and the value as arguments
                        data[k] = func(k,v)

        # return the dictionary
        return data


def get_config_run_version_1(config):
    """Get the configuration from version 1 YAML 
    configuration files.
    """

    # get the protocol family
    family = config["family"]

    # if the protocol family is not recognized, raise
    # an error
    if not family in ROSETTAPROTOCOLS.keys():
        errstr = f"Unrecognized protocol family {family}."
        raise ValueError(errstr)

    # for each step
    for stepname, step in config["steps"].items():
        
        # if the stepname is not recognized, raise an error
        if not stepname in ROSETTAPROTOCOLS[family].keys():
            errstr = f"Unrecognized step name {stepname} " \
                     f"for protocol family {family}."
            raise ValueError(errstr)

        # get the step fixed settings
        stepsettings = ROSETTAPROTOCOLS[family][stepname]
        
        # only consider Rosetta steps
        if stepsettings["runby"] == "rosetta":
            # create a copy of the configuration
            stepopts = dict(step["options"])
            # recursively remove all options that map
            # to None and convert the other ones to strings
            newstepopts = recursive_traverse(\
                            data = stepopts, \
                            actions = ["pop_empty", "substitute"], \
                            func = convert_rosetta_option_to_string)
            # update the dictionary of options with the new
            # options for the step
            config["steps"][stepname]["options"] = newstepopts

    # return the configuration
    return config


def get_config_run(configfile):
    """Get the configuration for running the protocol."""

    # load the configuration from the file
    config = yaml.safe_load(open(configfile, "r"))

    # check the version of the configuration file
    if config["version"] == 1:
        # return the configuration written in version 1 format
        return get_config_run_version_1(config = config)
    else:
        errstr = "Only version 1 configuration files " \
                 "are supported for now."
        raise ValueError(errstr)


def get_option_key(options, option):
    """Given a dictionary of Rosetta options and the name of
    a particular option as defined in ROSETTAOPTIONS, get which
    one of the possible alterative keys to define that option
    has been used in the given dictionary.
    """

    # get the values of each of the command line arguments that
    # can define the option of interest in the options dictionary
    # (it would be an empty string if not found)
    keys = set(ROSETTAOPTIONS[option]) & set(options.keys())
    
    # if the user passed multiple equivalent command line
    # arguments for the option (i.e. both -in:file:s and -s), 
    # the longest one will be used
    listkeys = sorted(list(keys), key = len, reverse = True)
    
    # if the key exists, return it
    if listkeys:
        return listkeys[0]
    # return None otherwise
    else:
        return None


def update_options(options, \
                   pdbfile, \
                   mut = None):
    """Update a dictionary of Rosetta options with the input
    PDB file and possibly replace specific placeholders with
    the corresponding attribute of a mutation.
    """

    # create a brand new dictionary to store the updated options
    newoptions = copy.deepcopy(options)
    
    # get one of the option keys used to define the
    # input PDB file (the longest one available by default)
    inpdbfileopt = sorted(ROSETTAOPTIONS["inpdbfile"], \
                          key = len, \
                          reverse = True)[0]
    
    # update the new dictionary with the input PDB file option
    newoptions[inpdbfileopt] = pdbfile

    # possible RosettaScript in the protocol
    rosettascriptopt = get_option_key(options, "protocol")
    rosettascript = options[rosettascriptopt] if rosettascriptopt \
                    else None

    # if a RosettaScript exists
    if rosettascript:
        # if the RosettaScript is just a XML file name
        if os.path.basename(rosettascript) == rosettascript:
            # assume it is a file in the default directory
            # containing RosettaScripts
            rosettascript = os.path.join(ROSETTASCRIPTSDIR, \
                                         rosettascript)
        # otherwise take the absolute path of the RosettaScript
        else:
            rosettascript = os.path.abspath(rosettascript)
        # add the option with the path to the RosettaScript option
        newoptions[rosettascriptopt] = rosettascript
    
        # if a mutation was passed
        if mut:
            # get the key in the option dictionary defining the
            # variable substitutions to be performed inside the
            # Rosetta script
            scriptvarsopt = get_option_key(options, "scriptvars")
            # create a dictionary of the variable substitutions
            # with values and variables swapped (the keyword of the 
            # mutation attribute would be a value)
            val2var = {v : k for k,v in options[scriptvarsopt].items()}
            # for each attribute of the mutation
            for mutkey, mutval in mut.items():
                # get the name of the variable that will store the
                # mutation attribute, if present
                scriptvar = val2var.get(mutkey, None)
                # if the variable exists
                if scriptvar:
                    # susbtitute its current value (the mutation
                    # attribute keyword) with the mutation attribute
                    # value
                    newoptions[scriptvarsopt][scriptvar] = mutval
    
    # return the updated dictionary of options
    return newoptions


def get_outpdbname(options, pdbfile, struct = None):
    """Given a dictionary of Rosetta options, the name of
    an input PDB file that was the starting structure for the
    generation of an ensemble and possibly the number
    of the structure of interest in that ensemble, get the
    correct name given by Rosetta to that structure.
    """

    # get the PDB file name
    pdbname = pdbfile.rstrip(".pdb")
    
    # get the option defining the PDB prefix, if any
    prefixopt = get_option_key(options, "outprefix")
    # get the defined prefix or set it to an empty string
    # if not present
    prefix = options[prefixopt] if prefixopt else ""
    
    # get the option defining the PDB suffix, if any
    suffixopt = get_option_key(options, "outsuffix")
    # get the defined prefix or set it to an empty string
    # if not present
    suffix = options[suffixopt] if suffixopt else ""
    
    # if a structure number was specified
    if struct:
        return f"{prefix}{pdbname}{suffix}_{struct:04d}.pdb"
    # if no structure number was specified
    else:
        return f"{prefix}{pdbname}{suffix}.pdb"



########################## MUTATIONS-RELATED ##########################



def get_mutlist(listfile):
    """Parse the file containing the list of 
    positions/mutations.
    """

    with open(listfile, "r") as f:
        # initialize an empty list to store the mutations
        mutlist = []
        for line in f:
            # ignore empty lines
            if re.match(r"^\s*$", line):
                continue
            # separate the mutations from extra data associated
            # to them     
            data, *extradata = line.rstrip("\n").split(" ")
            # create an empty list to store the mutation
            mut = []
            # for each single mutation potentially part of
            # a multiple mutation
            for singlemut in data.split(MUTSEP):
                mut.append(tuple(singlemut.split(COMPSEP)))
            # convert the list to a tuple
            mut = tuple(mut)
            # convert each piece of extra data to a string
            extradata = [str(item) for item in extradata]
            # append the tuple containing the mutation to
            # the mutation list together with all extra data
            mutlist.append(tuple([(mut), *extradata]))
        
        # return the list of mutations
        return mutlist


def get_reslist(reslistfile):
    """Parse the file containing the list of residue types
    to be used for saturation mutagenesis scans.
    """
            
    with open(reslistfile, "r") as f:
        # return all lines that are not empty
        return [l.rstrip("\n") for l in f \
                if not re.match(r"^\s*$", l)]


def get_saturation_mutlist(poslist, reslist):
    """Generate a mutation list for saturation mutagenesis from
    a list of residue positions and a list of residue types.
    """

    # create an empty list to store the saturation mutagenesis
    # list of mutations
    satmutlist = []
    # for each residue position
    for posdata, *extradata in poslist:
        # check that you have only single positions
        if len(posdata) > 1:
            errstr = "You cannot have multiple simultaneous " \
                     "mutations in saturation mutagenesis scans."
            raise ValueError(errstr)
        # for each position in the list
        for (chain, wtr, pos) in posdata:
            # add all possible mutations for the current
            # position to the list
            satmutlist.extend(\
                [(((chain, wtr, pos, res),), *extradata) \
                 for res in reslist])
    # return the list of mutations
    return satmutlist


def convert_to_pose_numbering(mutlist, pdbfile):
    """Retun a copy of the mutation list with residue numbers
    changed to the Rosetta pose numbering.
    """

    # create the PDB parser
    parser = PDB.PDBParser()
    
    # get the structure fron the PDB file
    structure = parser.get_structure("structure", pdbfile) 
    
    # create an empty dictionary to store the mapping
    pdbnum2posenum = {}
    
    # set a chain offset to start at the correct index when
    # changing chain
    chainoffset = 0
    
    # for each chain
    for chain in structure[0]:
        # for each residue in the chain
        for ix, res in enumerate(chain, start = 1):
            # unpack the residue full ID
            struc, mod, chainid, (het, resn, icode) = res.get_full_id()
            # map the PDB chain and residue number to the numbering
            # that the residue will have in the mutfile (starts
            # at 1); convert the PDB residue number to a string
            pdbnum2posenum[(chainid, str(resn))] = \
                str(ix + chainoffset)
        
        # update the chain offset
        chainoffset += ix
    
    # generate the mutation list but with residue numbers changed
    # according to the mutfile convention
    posemutlist = []
    for mut, *extradata in mutlist:
        newmut = \
            [(chain, wtr, pdbnum2posenum[(chain, numr)], mutr) \
              for (chain, wtr, numr, mutr) in mut]
        posemutlist.append(tuple([tuple(newmut), *extradata]))
    
    # return the new list
    return posemutlist


def generate_mutation_dirpath(mutdict):
    """Given a mutation, generate the corresponding directory
    name (since the Rosetta-compatible name of the mutation may
    not be usable as-it-is as a directory name).
    """

    # get the keys of the attributes of the mutation
    keys = (CHAIN, WTR, NUMR, MUTR)
    
    # initialize a list for the formatted mutation names
    fmtmuts = []      
    
    # each single mutation is in the form ("A", "C151Y")
    for smut in mutdict[MUT]:
        # get the mutation attributes
        chain, wtr, numr, mutr = operator.itemgetter(*keys)(smut)
        # strip Rosetta identifiers of NCAA from the name
        # of the residues
        wtr = wtr.strip("X[]")
        mutr = mutr.strip("X[]")
        if chain != "_":
            # chain with chain IDs
            fmtmuts.append(f"{chain}{DIRCHAINSEP}{wtr}{numr}{mutr}")
        else:
            # one chain with no chain ID
            fmtmuts.append(f"{wtr}{numr}{mutr}")
    
    # get the string representing the mutation directory name
    dirname = DIRMUTSEP.join(fmtmuts)
    
    # if the run generates multiple structures
    if STRUCT in mutdict.keys():
        # add another level to the path
        dirpath = os.path.join(dirname, mutdict[STRUCT])
    else:
        dirpath = dirname
    
    # return the directory path and name
    return dirpath, dirname


def get_mutations(listfile, \
                  reslistfile, \
                  resnumbering, \
                  extra, \
                  nstruct, \
                  pdbfile):
    """Get the list of mutations to be performed."""

    # check the pose numbering argument
    if not resnumbering in ("pose", "pdb"):
        errstr = f"Residue numbering must be either 'pose' or " \
                 f"'pdb', but '{resnumbering}' was passed."
        raise ValueError(errstr)


    #---------------------- Saturation or not? -----------------------#


    # get the mutations/positions list
    mutlist = get_mutlist(listfile)
    # if a list of residue types has been passed
    if reslistfile:
        reslist = get_reslist(reslistfile)
        # treat mutlist as a list of positions
        # and generate the new list of mutations
        mutlist = get_saturation_mutlist(mutlist, reslist)


    #----------------------- Residue numbering -----------------------#


    # if the protocol requires the residue numbering to follow
    # the Rosetta pose numbering convention
    if resnumbering == "pose":
        # convert the numbering from PDB to pose numbering
        mutlist = convert_to_pose_numbering(mutlist, pdbfile)

    # otherwise assume it follows the PDB convention (requires no
    # conversion)


    #------------------------ List generation ------------------------#


    # create an empty list to store the mutations
    mutations = []
    
    # ((("A","C","151","Y"), ("A","S","154","N")), *extradata)
    for mut, *extradata in mutlist:
        
        # create a dictionary where the mutation's attributes
        # will be stored
        mutdict = {MUT : []}
        
        # for each mutation attribute
        for chain, wtr, numr, mutr in mut:
            # if the wild-type residue is noncanonical
            if len(wtr) > 1:
                wtr = f"X[{wtr}]"
            # if the mutant residue is noncanonical
            if len(mutr) > 1:
                mutr = f"X[{mutr}]"             
            # update the mutation dictionary
            mutdict[MUT].append({CHAIN : chain, \
                                 WTR : wtr, \
                                 NUMR : numr, \
                                 MUTR : mutr})
        
        # if there are extra data associated with the mutation      
        if extra:
            # add them with the corresponding key
            for key, data in zip(extra, extradata):
                mutdict[key] = data
        
        # if multiple structures will be generatred for each mutation
        if nstruct:
            # add the mutation as many times as the structures
            # requested, with a new key identifying the structure
            # number
            for struct in range(nstruct):
                # create a fresh copy of the dictionary to be modified
                mutdict1 = dict(mutdict)
                # update the copy
                mutdict1.update({STRUCT : str(struct + 1)})
                # append the copy to the list of mutations
                mutations.append(mutdict1)

        # if only one structure will be generated
        else:
            mutations.append(mutdict)

    for mutdict in mutations:
        # add the mutation directory name and directory path
        mutdirpath, mutdirname = generate_mutation_dirpath(mutdict)
        mutdict[MUTDIRPATH] = mutdirpath
        mutdict[MUTDIRNAME] = mutdirname 
   
    # return the list of mutations
    return mutations


def write_dirnames2mutations(mutations, d2mfile):
    """Write a comma-separated file mapping the name of
    each directory containing data for a mutation to the
    mutation itself, represented as in the mutations'
    list file.
    """

    with open(d2mfile, "w") as out:
        # get the keys of the attributes of the mutation
        keys = (CHAIN, WTR, NUMR, MUTR, MUTDIRNAME)
        # for each mutation
        for mut in mutations:
            # for each single mutation in the mutation (only
            # one expected since the function is used only in
            # saturation mutagenesis scans)
            for smut in mut:
                # get the attributes of the mutation
                chain, wtr, numr, mutr, dirname = \
                    operator.itemgetter(*keys)(smut)
                # exit the loop
                break
            # write out the directory name and the mutation
            out.write(f"{dirname},{chain}{COMPSEP}"\
                      f"{wtr}{COMPSEP}{numr}{COMPSEP}{mutr}\n")


def get_dirnames2mutations(d2mfile):
    """Create a dataframe from a comma-separated file
    mapping the name of each directory containing data
    for a mutation to the mutation itself, represented
    as in the mutations' list file.
    """

    with open(filename, "r") as f:
        # create an empty list to store directory 
        # names and mutations
        dirnames2mutations = []
        
        for l in f:
            # ignore empty lines
            if re.match(r"^\s*$", l):
                continue

            # get directory name and mutation
            dirname, mut = l.rstrip("\n").split(",")
            # get the different attributes of the mutation
            chain, wtr, numr, mutr = tuple(mut.split(COMPSEP))

            # append a dictionary mapping each attribute
            # name to the attribute itself to the list
            dirnames2mutations.append(\
                {ROSETTADFCOLS["mutation"] : dirname, \
                 CHAIN : chain, \
                 WTR : wtr, \
                 NUMR: numr, \
                 MUTR : mutr})

        # create a dataframe from the list
        df = pd.DataFrame(dirnames2mutations)
        # sort mutations first by chain ID, then by residue number
        # and finally alphabetically by wild-type residue
        df = df.sort_values(by = [CHAIN, NUMR, WTR])
        # create a column storing the position name (chain ID,
        # residue number and wild-type residue but no mutant
        # residue)
        df[NOMUTR] = df[CHAIN] + DIRCHAINSEP + df[WTR] + df[NUMR]
        # drop the columns containing chain IDs, residue numbers
        # and wild-type residues (were kept only for sorting
        # purposes)
        df = df.drop([CHAIN, NUMR, WTR], axis = 1)
        # return the dataframe
        return df



################################ OTHERS ###############################



def check_pdbfile(pdbfile, \
                  allow_multi_chains, \
                  allow_no_chain_ids):
    """Check a PDB file before passing it to Rosetta."""
    
    # create the PDB parser
    parser = PDB.PDBParser()
    
    # try to get the structure
    try:
        structure = parser.get_structure("structure", pdbfile)
    # in case the PDB file was not found
    except FileNotFoundError:
        # log the error and propagate the exception
        errstr = f"PDB file {pdbfile} not found."
        raise FileNotFoundError(errstr)
    # in case something went wrong in accessing/opening the file
    except IOError:
        # log the error and propagate the exception
        errstr = f"Could not open PDB file {pdbfile}."
        raise IOError(errstr)
    # in case something else went wrong
    except Exception as e:
        # log the error and propagate the exception
        raise Exception(e)

    # if the PDB file contains more than one model
    if len(structure) > 1:
        # raise an exception since multi-model structures
        # are not allowed
        errstr = "Multi-model PDB files are not allowed."
        raise ValueError(errstr)
    
    # get the number of chains in the structure
    numchains = len(structure[0])
    
    # check if multiple chains are allowed
    if not allow_multi_chains and numchains > 1:
        # raise an exception if the structure contains multiple
        # chains but multiple chains are not allowed
        errstr = "Multi-chains structures are not allowed."
        raise ValueError(errstr)
    
    # check if the absence of chain IDs is allowed
    if not allow_no_chain_ids:
        # for each chain in the first model of the structure (it
        # is already guaranteed that there is only one model)
        for chain in structure[0]:
            # raise an exception if there are no chain IDs but
            # chain IDs were mandatory
            if not chain.get_id():
                errstr = "All chains must have a chain ID."
                raise ValueError(errstr)
    
    # return the PDB file
    return pdbfile


def get_abspath(path):
    """Given a path, return its absolute path. Return
    None if the path given is None.
    """

    return os.path.abspath(path) if path is not None else path


def get_items(d, keys, default = None):
    """Similar to operator.itemgetter but defaults to
    a specific value if the key is not found (instead of
    throwing an exception).
    """

    return [d.get(k, default) for k in keys]
