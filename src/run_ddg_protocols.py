#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    run_ddg_protocols.py
#
#    User-friendly Python wrapper to run the Rosetta-based protocols
#    developed by Park et al. [1]_ and Barlow et al. to predict the
#    ΔΔG of stability upon mutation of a monomeric protein and the
#    ΔΔG of binding upon mutation of a protein complex.
#
#    .. [1] Park, Hahnbeom, et al. "Simultaneous optimization of
#    biomolecular energy functions on features from small molecules
#    and macromolecules." Journal of chemical theory and computation 
#    12.12 (2016): 6201-6212.
#
#    .. [2] Barlow, Kyle A., et al. "Flex ddG: Rosetta Ensemble-Based
#    Estimation of Changes in Protein–Protein Binding Affinity upon 
#    Mutation." The Journal of Physical Chemistry B 122.21 (2018): 
#    5389-5399.
#
#    Copyright (C) 2020 Valentina Sora <sora.valentina1@gmail.com>
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

import argparse
import atexit
import functools
import itertools
import logging as log
import multiprocessing as mp
import os
import os.path
import re
import signal
import subprocess

import Bio.PDB as PDB

# Currently implemented protocols.
PROTOCOLS = \
    {\
    "cartddg_talaris2014" : \
        {"family" : "cartddg", \
         "execnames" : ["relax", "cartesian_ddg"], \
         "steps" : ["relax", "cartesian_ddg"], \
         "ddgsteps" : ["cartesian_ddg"], \
         "parallelsteps" : ["cartesian_ddg"]}, \
    "cartddg_ref2015" : \
        {"family" : "cartddg", \
         "execnames" : ["relax", "cartesian_ddg"], \
         "steps" : ["relax", "cartesian_ddg"], \
         "ddgsteps" : ["cartesian_ddg"], \
         "parallelsteps" : ["cartesian_ddg"]}, \
    "flexddg_talaris2014" : \
        {"family" : "flexddg", \
         "execnames" : ["rosetta_scripts"], \
         "steps" : [], \
         "ddgsteps" : [], \
         "parallelsteps" : []}, \
    "flexddg_ref2015" : \
        {"family" : "flexddg", \
         "execnames" : ["rosetta_scripts"], \
         "steps" : [], \
         "ddgsteps" : [], \
         "parallelsteps" : []}, \
    }

# signals that determines the killing of child processes
SIGNALS = [signal.SIGTERM, signal.SIGINT]

# relative path (from the Rosetta installation directory)
# to the database
ROSETTADATABASE = "main/database"

# relative path (from the Rosetta installation directory)
# to the executables
ROSETTAEXECUTABLES = "main/source/bin"


def terminate_child_processes(pids):
    """Terminate child processes by process ID."""
    for pid in pids:
        # terminate the process
        subprocess.terminate(pid)


def run_rosetta(executable, \
                flagfile = None, \
                options = None, \
                scriptvars = None, \
                rosettaoutput = "rosetta.out", \
                mpiexec = None, \
                mpinproc = None, \
                leaveattached = False):
    """Run a Rosetta executable.
    """

    # set prefix (before Rosetta command line) to run with MPI
    # add --quiet option to mpirun since it exits with errors
    # otherwise
    mpiprefix = []
    if mpiexec is not None:
        mpiprefix.extend([mpiexec, "-n", mpinproc, "--quiet"])
        if leaveattached:
            # leave session attached
            mpiprefix.append("--leave-session-attached")
    
    # start setting the command line by adding MPI prefix
    # and Rosetta executable
    commandline = mpiprefix + [executable]
    
    # add flag file to the command line
    if flagfile is not None:
        commandline.extend(["@", flagfile])
    
    # add Rosetta command line options
    if options is not None:
        for option, value in options:
            # check identity, not equality
            if value is True:
                value = "true"
            elif value is False:
                value = "false"
            elif isinstance(value, list):
                value = " ".join([str(item) for item in value])
            # convert also int and float to str
            commandline.extend([option, str(value)])
    
    # add variable substitutions for a Rosetta XML script
    if scriptvars is not None:
        # enable variable substitution in the script
        commandline.append("-parser:script_vars")
        for scriptvar, value in scriptvars:
            # check identity, not equality
            if value is True:
                value = "true"
            elif value is False:
                value = "false"
            elif isinstance(value, list):
                value = " ".join([str(item) for item in value])
            # convert also int and float to str
            commandline.append(scriptvar + "=" + value)

    # using args as list and not as a single string
    # because, as per the Python documentation about
    # subprocess, a single string should only be used
    # if running an executable without arguments
    # using Popen to be able to terminate the processes
    # once launched
    popen = subprocess.Popen(commandline, \
                             stdout = open(rosettaoutput, "w"), \
                             stderr = subprocess.STDOUT)

    # running directory is also returned since it would
    # not be easily retrievable otherwise
    return popen, os.getcwd()


def run_cartddg_relax(wd, \
                      runrosettafunc, \
                      flagfile, \
                      executable, \
                      database, \
                      inpdbfile, \
                      outprefix, \
                      scorefile, \
                      dirname, \
                      rosettaoutput, \
                      mpiexec, \
                      mpinproc, \
                      leaveattached):
    """Run the `relax` step of the cartesian_ddg-based protocols.
    """

    #------------------- Setup the relax directory -------------------#

    # create a directory where to run the relaxation and enter it
    relaxdir = os.path.join(wd, dirname)  
    os.mkdir(relaxdir)
    os.chdir(relaxdir)                      
    # set the path to the output file
    relax_out = os.path.join(relaxdir, rosettaoutput)
    scorefile = os.path.join(relaxdir, scorefile)

    #-------------------- Setup the command line ---------------------#

    # set options
    options = \
        [("-in:file:s", inpdbfile), \
         ("-in:path:database", database), \
         ("-out:prefix", outprefix), \
         ("-out:file:scorefile", scorefile)]

    #------------------------------ Run ------------------------------#

    popen, rundir = runrosettafunc(executable = executable, \
                                   flagfile = flagfile, \
                                   options = options, \
                                   rosettaoutput = rosettaoutput, \
                                   mpiexec = mpiexec, \
                                   mpinproc = mpinproc, \
                                   leaveattached = leaveattached)

    # get back to the working directory
    os.chdir(wd)
    # return the results of the run
    return popen, rundir
    

def run_cartddg_cartesian_ddg(mutdef, \
                              wd, \
                              runrosettafunc, \
                              executable, \
                              flagfile, \
                              database, \
                              inpdbfile, \
                              mutfile, \
                              ncaalistfile, \
                              ddgoutput, \
                              rosettaoutput, \
                              mpiexec, \
                              mpinproc, \
                              leaveattached, \
                              pdbnum2mutfilenum):
    """Run the `cartesian_ddg` step of the cartesian_ddg-based
    protocols.
    """

    #------------------ Setup the mutation directory -----------------#
    
    # create a directory for the current mutation (single or
    # multiple) and enter it
    mutdirname = "_".join(mutdef).replace("[", "").replace("]", "")
    mutpath = os.path.join(wd, mutdirname)
    os.makedirs(mutpath)
    os.chdir(mutpath)  

    #---------------------- Generate the mutfile ---------------------#

    # write the file containing the directives for the mutation
    with open(mutfile, "w") as out:
        # set the header and start line
        header = "total {:d}\n{:d}"
        out.write(header.format(len(mutdef), len(mutdef)))      
        # define the mutation string
        mutstring = "\n{:s} {:s} {:s}"
        # write the mutations
        for mut in mutdef:
            chain, wtr, numr, mutr = mutdef
            out.write(mutstring.format(wtr, numr, mutr))

    #-------------------- Setup the command line ---------------------#

    # set options
    options = \
        [("-in:file:s", inpdbfile), \
         ("-in:path:database", database), \
         ("-ddg:mut_file", mutfile), \
         ("-ddg:out", ddgoutput)]

    # additional residues to be considered while packing
    if ncaalistfile is not None:
        options.extend(\
            [("-packing:packer_palette:extra_base_type_file", \
              ncaalistfile)])

    #------------------------------ Run ------------------------------#

    popen, rundir = runrosettafunc(executable = executable, \
                                   flagfile = flagfile, \
                                   options = options, \
                                   rosettaoutput = rosettaoutput, \
                                   mpiexec = mpiexec, \
                                   mpinproc = mpinproc, \
                                   leaveattached = leaveattached)

    # go to the working directory
    os.chdir(wd)
    # return the results of the run
    return popen, rundir


def run_flexddg(mutdef, \
                wd, \
                runrosettafunc, \
                executable, \
                flagfile, \
                rosettascript, \
                database, \
                inpdbfile, \
                resfile, \
                ncaalistfile, \
                ddgdbfile, \
                structdbfile, \
                rosettaoutput, \
                mpiexec, \
                mpinproc, \
                leaveattached):
    """Run a Flex ddG-based protocol.
    """

    # unpack the mutation definition.
    # For single mutations it is in the form:
    #   (((("A", "C151Y"),), "A"), 1)
    # For multiple simultaneous mutations it is in the form:
    #   (((("A", "C151Y"), ("A", "S154N")), "A"), 1)
    (((muts), chaintomove), struct) = mutdef

    #------------------ Setup the mutation directory -----------------#
  
    # create a directory for the current mutation (single or
    # multiple) and enter it
    mutdirname = \
        "_".join(["-".join(mut).replace("[", "").replace("]", "") \
                    for mut in muts])
    mutpath = os.path.join(wd, mutdirname)

    #----------------- Setup the structure directory -----------------#
    
    structpath = os.path.join(mutpath, str(struct))
    os.makedirs(structpath, exist_ok = True)
    os.chdir(structpath)

    #---------------------- Generate the resfile ---------------------#

    # write the file containing the directives for the mutation
    # (it follows the Rosetta resfile conventions)
    resfile = os.path.abspath(os.path.join(mutpath, resfile))
    # do not overwrite the resfile if you already have it
    if not os.path.exists(resfile):
        with open(resfile, "w") as out:
            # set the header and start line
            out.write("{:s}\n{:s}".format("NATAA", "start"))      
            # define the mutation string
            fstring = "\n{:s} {:s} PIKAA {:s}"          
            # each singlemut in the form ("A", "C151Y")
            for singlemut in muts:
                chain, wtr, numr, mutr = singlemut
                out.write(fstring.format(numr, chain, mutr))

    #-------------------- Setup the command line ---------------------#

    # set options
    options = \
        [("-in:file:s", inpdbfile), \
         ("-in:path:database", database), \
         ("-parser:protocol", rosettascript)]

    # additional residues to be considered while packing
    if ncaalistfile is not None:
        options.extend(\
            [("-packing:packer_palette:extra_base_type_file", \
              ncaalistfile)])

    # set substitutions for the script variables
    scriptvars = \
        [("chaintomove", chaintomove), \
         ("resfile", resfile), \
         ("ddgdbfile", ddgdbfile), \
         ("structdbfile", structdbfile)]

    # set protocol-related options and variable substitutions
    if flagfile == "flexddg_talaris2014":
        scriptvars.extend([("sfweights", "talaris2014")]) 
    elif flagfile == "flexddg_ref2015":
        scriptvars.extend([("sfweights", "ref2015")])

    #------------------------------ Run ------------------------------#

    popen, rundir = runrosettafunc(executable = executable, \
                                   flagfile = flagfile, \
                                   options = options, \
                                   scriptvars = scriptvars, \
                                   rosettaoutput = rosettaoutput, \
                                   mpiexec = mpiexec, \
                                   mpinproc = mpinproc, \
                                   leaveattached = leaveattached)

    # go back to the working directory
    os.chdir(wd)
    # return the results of the run
    return popen, rundir


def get_abspath(path):
    """Get the absolute path of an object. Return `None`if `None`
    was passed.
    """
    # return None if path is None
    if path is None:
        return path
    # remove trailing slash in the path if present (gives problems
    # in getting the basename of the path)
    elif path[-1] == "/":
        path = path[:-1]
    
    return os.path.abspath(path)


def get_reslist(reslistfile):
    """Get a list of residue types for saturation mutagenesis.
    """

    with open(reslistfile, "r") as f:
        reslist = []
        for line in f:
            if not re.match(r"^\s*$", line):
                # ignore empty lines
                res = line.rstrip("\n")
                reslist.append(res)

        return reslist


def get_mutlist(listfile, reslist = None, nstruct = None):
    """Generate a list of mutations from a file containing either 
    a list of selected mutations or a list of residue positions 
    to perform saturation mutagenesis.
    In this case, a set of residue types to be used for the
    mutagenesis must be passed.
    There can be extra data other than the mutation on each line
    (e.g. the chain to be moved in Flex ddG protocols), all
    whitespace-separated.

    Example of file containing a list of selected mutations
    (with extra data):

    A-C151Y,A-S154N A
    A-S2A A

    Example of file containing a list of positions (without
    extra data):
    
    A-C151
    A-S154
    """
    
    with open(listfile, "r") as f:
        mutlist = []
        for line in f:
            # ignore empty lines
            if not re.match(r"^\s*$", line):
                data = line.rstrip("\n").split(" ")
                muts = data[0].split(",")
                extradata = data[1:] if len(data) > 1 else []
                chainsmuts = tuple([tuple(m.split("-")) for m in muts])
                if reslist is not None:
                    chainsmuts = \
                        [(chainsmuts[0][0], chainsmuts[0][1]+res) \
                         for res in reslist]

                tmplist = []
                for singlemut in chainsmuts:
                    # chain, mut = ("A", "C151Y")
                    chain, mut = singlemut
                    # find any definition including NCAA
                    # i.e. "X[DALA]134A" becomes ["", "DALA]134A"]
                    # i.e. "A134X[DALA]" becomes ["A134", "DALA]"]
                    # i.e. "X[DALA]134X[DPHE] becomes 
                    #      ["", "DALA]134", "DPHE]"]
                    mut = mut.split("X[")          
                    
                    # only canonical residues
                    if len(mut) == 1:
                        wtr, numr, mutr = re.split(r"(\d+)", mut[0])
                    
                    # either the wt or the mutated residue is a NCAA              
                    elif len(mut) == 2:
                        # i.e. ["", "DALA]134A"] or ["A134", "DALA]"]
                        if mut[0] == "":
                            # the wild-type residue is non-canonical
                            # i.e. ["", "DALA]134A"]
                            # mut = ["", "DALA134A"]
                            mut = mut[1].replace("]", "")
                            # wt, numr, mutr = ["DALA", "134", "A"]
                            wtr, numr, mutr = re.split(r"(\d+)", mut)
                            wtr = "X[" + wtr + "]"
                        else:
                            # the mutated residue is non-canonical
                            # i.e. ["A134", "DALA]"]
                            # wtr, numr, void = ["A", "134", ""]
                            wtr, numr, void = re.split(r"(\d+)", mut[0])
                            # mutres = "[DALA]"
                            mutr = "X[" +  mut[1]
                    
                    # both wt and mutated residue are NCAAs             
                    elif len(mut) == 3:
                        # i.e. ["", "DALA]134", "DPHE]"]
                        # wtandnum = "DALA134"
                        wtandnum = mut[1].replace("]", "") 
                        # wtr, numr, void = ["DALA", "134", ""]
                        wtr, numr, void = re.split(r"(\d+)", wtandnum)
                        # mutr = "[DPHE]"
                        mutr = "X[" + mut[2]
                        # wtr = "[DALA]"
                        wtr = "X[" + wtr + "]"                

                    if reslist is None:
                        tmplist.append((chain, wtr, numr, mutr))
                    else:
                        tmplist.append((((chain, wtr, numr, mutr),), \
                                         *extradata))
                
                if reslist is None:
                    mutlist.append((tuple(tmplist), *extradata))
                else:
                    mutlist.extend(tmplist)

        if nstruct is not None:
            structmutlist = []
            for mut in mutlist:
                structmutlist.extend(\
                    [(mut, struct) for struct in range(1, nstruct+1)])
            mutlist = structmutlist
        
        # return the list of mutations
        return mutlist


def get_pose_mutlist(pdbfile, \
                     mutlist):
    """Retun a copy of the mutation list with residue numbers changed
    according to the Rosetta pose numbering (used in mutfiles)."""
    
    # create the PDB parser
    parser = PDB.PDBParser()
    # parse the structure (take only the first model)
    structure = parser.get_structure("structure", pdbfile)[0]
    # create an empty dictionary to store the mapping
    pdbnum2mutfilenum = {}
    # set a chain offset to start at the correct index when
    # changing chain
    chainoffset = 0
    for chain in structure:
        for index, res in enumerate(chain, start = 1):
            # unpack the residue full ID
            struc, mod, chainid, (het, resn, icode) = res.get_full_id()
            # map the PDB chain and residue number to the numbering
            # that the residue will have in the mutfile (starts at 1)
            # convert the PDB residue number to a string
            pdbnum2mutfilenum[(chainid, str(resn))] = index + chainoffset
        # update the chain offset
        chainoffset += index
    # generate the mutation list but with residue numbers changed
    # according to the mutfile convention
    posemutlist = []
    for muts, *extradata in mutlist:
        newmuts = \
            [(chain, wtr, pdbnum2mutfilenum[(chain, numr)], mutr) \
             for (chain, wtr, numr, mutr) in muts]
        posemutlist.append(tuple([newmuts, *extradata]))

    return posemutlist


def get_ncaa(ncaalistfile):
    """Get a list of non-canonical residues to be considered as 
    allowed while performing the mutations.
    """
    with open(ncaalistfile, "r") as f:
        for line in f:
            if not re.match(r"^\s*$", line):
                # ignore empty line
                # all NCAA should be on one line, whitespace-
                # separated
                return line.rstrip("\n").split(" ")


def get_struct_lowscore(scorefile):
    """Get the number of the structure having the lowest REU score
    (and the associated score) from a Rosetta scorefile.
    """
    
    with open(scorefile, "r") as f:
        parse = False
        nstruct = 1
        scores = []
        for line in f:
            if line.startswith("SCORE: total_score"):
                # start parsing the file
                parse = True
                continue

            if parse:
                line = \
                    [item for item in line.split(" ") if item != ""]
                scores.append((nstruct, float(line[1])))              
                nstruct += 1

        # sort scores in ascending order and take the first
        # element of the sorted list. If two structures have the
        # same score, the one generated first will be reported.
        return sorted(scores, \
                      key = lambda x: x[1])[0]



if __name__ == "__main__":

    
    ####################### SET ARGUMENT PARSER #######################

    description = \
        "\nScript to run Rosetta protocols for the prediction of " \
        "the ΔΔG of stability upon mutation of a monomeric protein " \
        "and the ΔΔG of binding upon mutation of a protein complex.\n"

    # create the argument parser
    parser = argparse.ArgumentParser(description = description)
    generalargs = parser.add_argument_group("General arguments")
    logargs = parser.add_argument_group("Logging-related arguments")
    mpiargs = parser.add_argument_group("MPI-related arguments")
    satargs = parser.add_argument_group("Saturation-related arguments")


    #---------------------------- General ----------------------------#

    p_choices = list(PROTOCOLS.keys())
    p_helpstr = "Protocol to be run."
    generalargs.add_argument("-p", "--protocol", \
                             type = str, \
                             choices = p_choices, \
                             required = True, \
                             help = p_helpstr)

    s_helpstr = \
        "Step(s) of the protocol to be run. Default is to run all " \
        "steps (it expects all the corresponding flag files if there " \
        "are multiple steps)."
    generalargs.add_argument("-s", "--steps", \
                             type = str, \
                             required = False, \
                             default = "all", \
                             help = s_helpstr)
    
    i_helpstr = "Input PDB file of the wild-type structure."
    generalargs.add_argument("-i", "--inpdbfile", \
                             type = str, \
                             required = True, \
                             help = i_helpstr)

    f_helpstr = \
        "Flag file(s) for the protocol to be run. Flag files must be " \
        "named after the protocols' names. If you want to run " \
        "multiple steps of the same protocol in one run, pass " \
        "the corresponding flag files in the correct order."
    # do not add the choices keyword here since the files may not
    # be located in the current working directory
    generalargs.add_argument("-f", "--flagfiles", \
                             type = str, \
                             required = True, \
                             nargs = "+", \
                             help = f_helpstr)

    rosettapath_helpstr = \
        "Path to the Rosetta installation directory."
    generalargs.add_argument("-rp", "--rosettapath", \
                             type = str, \
                             required = True, \
                             help = rosettapath_helpstr)

    d_helpstr = \
        "Working directory. Default is the current working directory."
    generalargs.add_argument("-d", "--workdir", \
                             type = str, \
                             required = False, \
                             default = os.getcwd(), \
                             help = d_helpstr)

    l_helpstr = \
        "File containing the list of selected mutations " \
        "or positions to perform saturation mutagenesis."
    generalargs.add_argument("-l", "--listfile", \
                             type = str, \
                             required = False, \
                             default = None, \
                             help = l_helpstr)

    ncaalistfile_helpstr = \
        "File containing a list of non-canonical residue types " \
        "(full names). It is needed if you want to specify " \
        "mutations to non-canonical residues."
    generalargs.add_argument("-nc", "--ncaalistfile", \
                             type = str, \
                             required = False, \
                             default = None, \
                             help = ncaalistfile_helpstr)

    rosettascript_helpstr = "XML Rosetta script."
    generalargs.add_argument("-rs", "--rosettascript", \
                             dest = "rosettascript", \
                             type = str, \
                             required = False, \
                             default = None, \
                             help = rosettascript_helpstr)

    nproc_helpstr = \
        "Number of processes to be started in parallel. " \
        "Default is one process (no parallelization)."
    generalargs.add_argument("-n", "--nproc", \
                             dest = "nproc",
                             type = int, \
                             required = False, \
                             default = 1, \
                             help = nproc_helpstr)

    #---------------------------- Logging ----------------------------#

    logfile_default = "run_ddg_protocols.log"
    logfile_helpstr = \
        "Name of the log file. Default is {:s}.".format(logfile_default)
    logargs.add_argument("-lf", "--logfile", \
                         dest = "logfile", \
                         required = False, \
                         default = logfile_default, \
                         help = logfile_helpstr)

    v_helpstr = "Verbosity level: verbose."
    logargs.add_argument("-v", \
                         action = "store_true", \
                         dest = "v", \
                         help = v_helpstr)

    vv_helpstr = "Verbosity level: debug."
    logargs.add_argument("-vv", \
                         dest = "vv", \
                         action = "store_true", \
                         help = vv_helpstr)

    #------------------------------ MPI ------------------------------#

    mpiexec_helpstr = \
        "Executable to run Rosetta with MPI (e.g. mpiexec, mpirun)."
    mpiargs.add_argument("--mpiexec", \
                         dest = "mpiexec",
                         type = str, \
                         required = False, \
                         default = None, \
                         help = mpiexec_helpstr)

    leaveattached_helpstr = "Leave session attached when using MPI."
    mpiargs.add_argument("--leave-session-attached", \
                         action = "store_true", \
                         dest = "leaveattached",
                         help = leaveattached_helpstr)

    #-------------------- Saturation mutagenesis ---------------------#

    saturation_helpstr = \
        "Perform saturation mutagenesis on selected positions."
    satargs.add_argument("--saturation", \
                         dest = "saturation" , \
                         action = "store_true", \
                         help = saturation_helpstr)

    r_helpstr = \
        "File containing the list of residue types " \
        "(one-letter name) to include in the saturation " \
        "mutagenesis. It is used only if --saturation is provided."
    satargs.add_argument("--reslistfile", \
                         type = str, \
                         required = False, \
                         default = None, \
                         help = r_helpstr)


    # parse the arguments
    args = parser.parse_args()
    # general arguments
    protocol = args.protocol
    steps = args.steps
    inpdbfile = get_abspath(args.inpdbfile)
    flagfiles = get_abspath(args.flagfiles)
    rosettapath = get_abspath(args.rosettapath)
    workdir = get_abspath(args.workdir)
    listfile = get_abspath(args.listfile)
    ncaalistfile = get_abspath(args.ncaalistfile)
    rosettascript = get_abspath(args.rosettascript)
    nproc = args.nproc
    # logging-related arguments
    v = args.v
    vv = args.vv
    logfilename = get_abspath(args.logfile)
    # MPI-related arguments
    mpiexec = args.mpiexec
    leaveattached = args.leaveattached
    # saturation mutagenesis-related arguments
    saturation = args.saturation
    reslistfile = get_abspath(args.reslistfile)

    # get the full path to the Rosetta database
    rosettadatabase = os.path.join(rosettapath, ROSETTADATABASE)


    #################### SET EXECUTABLES AND STEPS ####################

    # get the full path to all Rosetta executables 
    allexecs = os.listdir(rosettapath, ROSETTAEXECUTABLES)
    # create an empty dictionary to store the Rosetta executables
    # necessary to run the protocol
    rosettaexecs = {}
    for execname in PROTOCOLS[protocol]["execnames"]:
        # filter function to get the right Rosetta executable
        filterf = lambda x: x.startswith(execname)
        if mpiexec is not None:
            # modified filter function if the user wants to run with
            # MPI (different set of Rosetta executables)
            filterf = lambda x: x.startswith(execname) and ".mpi." in x
        # save the first occurrence of the list of Rosetta executables
        # satisfying the criteria of the filter function (the list
        # should contain only one element)
        rosettaexecs[execname] = list(filter(filterf, allexecs))[0]

    # set the flag that signals if we are running ΔΔG prediction
    ddgprediction = \
        steps in PROTOCOLS[protocol]["ddgsteps"] or steps == "all"
    # set the flag that signals if we are running a step that allows
    # multiprocessing (Python parallelization, different from MPI)
    pythonparallel = \
        steps in PROTOCOLS[protocol]["parallelsteps"] or ddgprediction


    ######################### SET THE LOGGER ##########################

    # set the log level (default is only warnings and errors)
    loglevel = log.WARNING
    if v:
        loglevel = log.INFO
    if vv:
        # if both -v and -vv are passed, -vv is used
        loglevel = log.DEBUG
    
    # set the log file
    logfile = os.path.join(os.getcwd(), logfilename)
    # set log and date format
    logfmt = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    datefmt = "%d-%b-%y %H:%M:%S"
    # configure the logger
    log.basicConfig(filename = logfile, \
                    filemode = "w", \
                    level = loglevel, \
                    format = logfmt, \
                    datefmt = datefmt)


    ######################### SET THE WORKDIR #########################

    # go to the correct working directory
    os.chdir(workdir)
    # set a general string to store info about the runs
    runinfostr = \
        "Process in {:s} exited with code {:d}. Command line: {:s}"


    #################### GENERATE THE MUTATION LIST ###################

    # the list of mutations is needed only if predicting the ΔΔG
    if ddgprediction:
        # initialize reslist to None
        reslist = None
        # initialize nstruct to None
        nstruct = None
        # if saturation mutagenesis has been requested
        if saturation:
            # get the residue types set
            reslist = get_reslist(reslistfile = reslistfile)
            # log the residue types that will be used for the
            # saturation mutagenesis
            logstr = "Residue types included in the saturation " \
                     "mutagenesis: \n{:s}"
            log.info(logstr.format(", ".join(reslist)))
        # if a Flex ddG protocol has been requested
        if protocol in protocols["flexddg"]:
            # number of structures in the ensemble
            nstruct = 35
        # generate the list of mutations
        mutlist = get_mutlist(listfile = listfile, \
                              reslist = reslist, \
                              nstruct = nstruct)


    ####################### SET UP PROCESSES ##########################

    # initialize the variable keeping track of how many
    # Python processes are requested
    pythonnproc = None
    # initialize the variable keeping track of how many MPI processes
    # for a step including a single Rosetta run are requested
    mpinproc = None
    # initialize the variable keeping track of how many MPI processes
    # for a step including multiple Rosetta runs are requestes
    mpinprocs = None

    # if MPI parallelization requested
    if mpiexec is not None:
        # the number of MPI processes for steps that require only
        # one Rosetta run is always equal to the number of total
        # processes requested
        mpinproc = nproc
        # if a step parallelizing via multiprocessing is requested
        if pythonparallel:
            # if a step is the ΔΔG prediction
            if ddgprediction:
                nmuts = len(mutlist)
                # one MPI process per mutation, since it does not scale
                mpinprocs = [1 for mut in mutlist]
                if nmuts >= nproc:
                    # the Pool will start with as many Python processes
                    # as the number of processes given by the user
                    pythonnproc = nproc
                else:
                    logstr = \
                        "Setting more processes than mutations does " \
                        "not provide any performance improvement. " \
                        "Only {:d} processes will be launched for " \
                        "the ΔΔG prediction."
                    log.warning(logstr.format(nmuts))            
                    # the Pool will start with one process per mutation
                    pythonnproc = nmuts

        # no steps parallelizing via multiprocessing requested  
        else:
            # only one Python process (no parallelization via
            # multiprocessing, only MPI) since Rosetta is run 
            # only once
            pythonnproc = 1

    # if no MPI requested
    else:
        # if a step parallelizing via multiprocessing is requested
        if pythonparallel:
            # if a step is the ΔΔG prediction
            if ddgprediction:
                # no MPI processes also for steps performing the
                # ΔΔG prediction (multiple Rosetta runs)
                mpinprocs = [None for mut in mutlist]
                # as many Python processes as the number of processes
                # requested by the user
                pythonnproc = nproc

        # no steps parallelizing via multiprocessing requested 
        else:
            if nproc > 1:
                logstr = \
                    "Requesting multiple processes for step {:s} " \
                    "without using MPI does not provide any " \
                    "performance improvement. Only one process will " \
                    "be launched."
                log.warning(logstr.format(step))
            # only one Python process
            pythonnproc = 1

    # set the empty list that will contain the list of IDs of child
    # processed started
    pids = []


    ######################### CHECK PDB FILE ##########################

    parser = PDB.PDBParser()
    structure = parser.get_structure(inpdbfile)
    if len(structure) > 1:
        logstr = \
            "Multi-model PDB files are not supported. Exiting ..."
        log.error(logstr)
        exit(1)


    ######################## CARTDDG PROTOCOLS ########################

    if protocol["family"] == "cartddg":

        #--------------------------- relax ---------------------------#

        if step in ("relax", "all"):
            # relax executable
            relaxexec = rosettaexecs["relax"]
            # prefix for the names of the output structures
            outprefix = "relaxed_"
            # name of the scorefile
            scorefile = "score.sc"
            # name of the relaxation directory to be created
            relaxdir = "relax"
            # name of the Rosetta output
            rosettaoutput = "relax.out"
            # run the relaxation
            popen, rundir = \
                run_cartddg_relax(wd = workdir, \
                                  runrosettafunc = run_rosetta, \
                                  executable = relaxexec, \
                                  flagfile = flagfile, \
                                  rosettaoutput = rosettaoutput, \
                                  inpdbfile = inpdbfile, \
                                  database = rosettadatabase, \
                                  outprefix = outprefix, \
                                  scorefile = scorefile, \
                                  dirname = relaxdir, \
                                  mpiexec = mpiexec, \
                                  mpinproc = mpinproc, \
                                  leaveattached = leaveattached)
                
            # format the information about the run
            runinfostrfilled = \
                runinfostr.format(rundir, popen.returncode, \
                                " ".join(popen.args))

            # add the process ID to the list of child processes IDs
            pids.append(popen.pid)

            if popen.returncode == 0:
                # inform the user about the run
                log.info(runinfostrfilled)
            else:
                # warn the user since something went wrong
                log.error(runinfostrfilled)
                # exit the script (if something went wrong in the
                # relaxation it does not make sense to proceed)
                exit(1)

            # get the number of the lowest-score structure and
            # the corresponding score
            structlowscore, score = \
                get_struct_lowscore(os.path.join(relaxdir, scorefile))
            # report the structure having the lowest score
            logstr = \
                "Relaxed structure number {:d} is the one having " \
                "the lowest score ({:.3f})."
            log.info(logstr.format(structlowscore, score))                
            # go back to the top directory
            os.chdir(workdir)

        #--------------- prepare to run cartesian_ddg ----------------#

        if step == "all":
            # get the name of the structure with lowest score
            structlowscorename = \
                outprefix + str(structlowscorename) + ".pdb"
            # this will be the new input PDB file for cartesian_ddg
            inpdbfile = os.path.join(relaxdir, structlowscore)

        #----------------------- cartesian_ddg -----------------------#

        if step in ("cartesian_ddg", "all"):
            # cartesian_ddg executable
            cartddgexec = rosettaexecs["cartesian_ddg"]
            # name of the mutation file to be created for each
            # mutation
            mutfile = "mutation.mutfile"
            # name of the output file storing the scores
            ddgoutput = "ddg_scores.out"
            # name of the Rosetta output file
            rosettaoutput = "cartesian_ddg.out"
            # change residue numbers in the mutation list to abide
            # by the mutfile convention
            try:
                mutlist = get_pose_mutlist(inpdbfile, mutlist)
            except KeyError:
                logstr = \
                    "The chain IDs and/or residue numbers in the " \
                    "mutation file do not correspond to those in " \
                    "the PDB file."
                log.error(logstr)
                exit(1)
            # set values for some arguments of the run function
            partfunc = \
                functools.partial(\
                    run_cartddg_cartesian_ddg, \
                    wd = workdir, \
                    runrosettafunc = run_rosetta, \
                    executable = cartddgexec, \
                    flagfile = flagfile, \
                    database = rosettadatabase, \
                    inpdbfile = inpdbfile, \
                    mutfile = mutfile, \
                    ncaalistfile = ncaalistfile, \
                    ddgoutput = ddgoutput, \
                    rosettaoutput = rosettaoutput, \
                    mpiexec = mpiexec, \
                    leaveattached = leaveattached)
            # set arguments for the partial function
            partfuncargs = zip(mutlist, mpinprocs)


    ######################## FLEXDDG PROTOCOLS ########################

    elif protocol["family"] == "flexddg":     
        # check the XML Rosetta script
        if rosettascript is None:
            logstr = "You must provide the XML script describing " \
                     "the Flex ddG procedure. Exiting ..."
            log.error(logstr)
            exit(1)

        # check the structure for missing chain IDs (take the first
        # model, since it has already been ensured that the structure
        # contains only one model)
        for chain in structure[0]:
            # if the chain ID is an empty string
            if not chain.get_id():
                logstr = "Error in input PDB file: some chains do " \
                         "not have an ID. Flex ddG protocols cannot " \
                         "be run if the PDB has missing chain IDs. " \
                         "Exiting ..."
                log.error(logstr)
                exit(1)

        # rosetta_scripts executable
        rosettascriptsexec = rosettaexecs["rosetta_scripts"]
        # name of the mutation file (in resfile format) to be
        # created for each mutation
        resfile = "mutation.resfile"
        # name of the Rosetta output file
        rosettaoutput = "rosettascripts.out"
        # name of the database file containing on the ΔΔG scores
        ddgdbfile = "ddg.db3"
        # name of the database file containing info on the structures
        structdbfile = "struct.db3"
        # set values for some arguments of the run function
        partfunc = \
            functools.partial(\
                run_flexddg, \
                wd = workdir, \
                runrosettafunc = run_rosetta, \
                executable = rosettascriptsexec, \
                database = rosettadatabase, \
                rosettascript = rosettascript, \
                flagfile = protocol, \
                inpdbfile = inpdbfile, \
                resfile = resfile, \
                ncaalistfile = ncaalistfile, \
                ddgdbfile = ddgdbfile, \
                structdbfile = structdbfile, \
                rosettaoutput = rosettaoutput, \
                mpiexec = mpiexec, \
                leaveattached = leaveattached)
        # set arguments for the partial function
        partfuncargs = zip(mutlist, mpinprocs)


    ############################### RUN ###############################

    # run a step that allows multiprocessing (does not necessarily
    # mean that the user requested it)
    if pythonparallel:
        # launch parallel Python processes
        if pythonnproc > 1:
            # set "spawn" method to start processes (it will avoid
            # all issues that "fork" has)
            mp.set_start_method("spawn")
            # create the Pool of workers
            pool = mp.Pool(nproc)
            # set the mapping function (parallel, asyncronous)
            starmapfunc = pool.starmap_async
        # do not launch parallel Python processes
        else:
            # set the mapping function (serial)
            starmapfunc = itertools.starmap

        # empty list to store the return codes
        returncodes = []
        # run
        for popen, rundir in starmapfunc(partfunc, partfuncargs):
            # format the information about the run
            runinfostrfilled = \
                runinfostr.format(rundir, popen.returncode, \
                                  " ".join(popen.args))

            # add the process ID to the list of child processes IDs
            pids.append(popen.pid)
            # save the return code for later
            returncodes.append(popen.returncode)
            # check the return code of the run
            if popen.returncode == 0:
                # inform the user about the run
                log.info(runinfostrfilled)
            else:
                # warn the user since something went wrong
                log.error(runinfostrfilled)

        if not all(code == 0 for code in returncodes):
            # exit the script with code 1 if somewhere,
            # something went wrong
            exit(1)              

        # go back to the main working directory
        os.chdir(workdir)             
        # close and clean up the pool
        pool.close()
        pool.join()


    # register the function to be called when the program exits
    # (for termination signals handled by Python)
    atexit.register(terminate_child_processes, pids)
    # set to kill child processes if signals are received from the OS
    # (for termination signals not handled by Python)
    for sig in SIGNALS:
        signal.signal(sig, terminate_child_processes)

