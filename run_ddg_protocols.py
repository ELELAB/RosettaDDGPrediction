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


import subprocess
import multiprocessing as mp
import argparse
import os
import os.path
import logging as log
import functools
import re


def run_rosetta(executable, \
                options = None, \
                scriptvars = None, \
                rosettaoutput = "rosetta.out"):
    """Run a Rosetta executable with options and variable
    substitutions for a Rosetta XML script (if running
    `rosettascripts`).

    Parameters
    ----------
    options : `list`
    """

    commandline = [executable]

    # add Rosetta command line options
    if options is not None:
        for option, value in options:
            # check identity, not equality
            if value is True:
                value = "true"
            elif value is False:
                value = "false"
            elif isinstance(value, list):
                value = \
                    " ".join([str(item) for item in value])

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
                value = \
                    " ".join([str(item) for item in value])

            # convert also int and float to str
            commandline.append(scriptvar + "=" + value)

    # using args as list and not as a single string
    # because, as per the Python documentation about
    # subprocess, a single string should only be used
    # if running an executable without arguments
    completedprocess = subprocess.run(commandline, \
                                      stdout = open(rosettaoutput, "w"), \
                                      stderr = subprocess.STDOUT)

    # running directory is also returned since it would
    # not be easily retrievable otherwise
    return completedprocess, os.getcwd()


def run_cartddg_relax(wd, \
                      runrosettafunc, \
                      protocol, \
                      executable, \
                      database, \
                      pdbfile, \
                      outprefix, \
                      scorefile, \
                      dirname, \
                      rosettaoutput):
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

    # set options for the inputs
    inputoptions = \
        [("-in:file:s", pdbfile), \
         ("-in:ignore_unrecognized_res", True), \
         ("-in:file:fullatom", True), \
         ("-in:path:database", database)]

    # set options for the outputs
    outputoptions = \
        [("-out:prefix", outprefix), \
         ("-out:file:scorefile", scorefile), \
         ("-out:nstruct", 20)]

    # set run options
    runoptions = \
        [("-run:ignore_zero_occupancy", False)]

    # set packing-related options
    packingoptions = \
        [("-packing:use_input_sc", True)]

    # set relax-related options
    #
    # keep an eye on the default relax script and
    # number of iterations in later Rosetta releases
    relaxoptions = \
        [("-relax:constrain_relax_to_start_coords", True), \
         ("-relax:coord_constrain_sidechains", True), \
         ("-relax:min_type", "lbfgs_armijo_nonmonotone"), \
         ("-relax:cartesian", True), \
         ("-relax:script", "MonomerRelax2019"), \
         ("-relax:default_repeats", 5)]

    # set protocol-related options
    if protocol == "cartddg_talaris2014":
        protocoloptions = \
            [("-score:weights", "talaris2014_cart"), \
             ("-corrections:restore_talaris_behavior", True)] 
    elif protocol == "cartddg_ref2015":
        protocoloptions = [("-score:weights", "ref2015_cart")]

    options = \
        inputoptions + outputoptions + runoptions + \
        packingoptions + relaxoptions + protocoloptions

    #------------------------------ Run ------------------------------#

    cproc, rundir = runrosettafunc(executable = executable, \
                                   rosettaoutput = rosettaoutput, \
                                   options = options)

    # get back to the working directory
    os.chdir(wd)

    # return the results of the run
    return cproc, rundir


def run_cartddg_cartesian_ddg(mutdef, \
                              wd, \
                              runrosettafunc, \
                              executable, \
                              protocol, \
                              database, \
                              pdbfile, \
                              mutfile, \
                              ncaalistfile, \
                              ddgoutput, \
                              rosettaoutput):
    """Run the `cartesian_ddg`step of the cartesian_ddg-based
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
        
        for mut in mutdef:
            # find any definition including non-canonical residues
            # i.e. "X[DALA]134A" becomes ["", "DALA]134A"]
            # i.e. "A134X[DALA]" becomes ["A134", "DALA]"]
            # i.e. "X[DALA]134X[DPHE] becomes 
            #      ["", "DALA]134", "DPHE]"]
            mut = mut.split("X[")              
            # handle the different scenarios
            if len(mut) == 1:
                # only canonical residues, no splitting happened
                wtres, numres, mutres = re.split(r"(\d+)", mut[0])
                
            elif len(mut) == 2:
                # either the wild-type or the mutated residue
                # is non-canonical
                # i.e. ["", "DALA]134A"] or ["A134", "DALA]"]
                if mut[0] == "":
                    # the wild-type residue is non-canonical
                    # i.e. ["", "DALA]134A"]
                        
                    # mut = ["", "DALA134A"]
                    mut = mut[1].replace("]", "")
                    # wtres, numres, mutres = ["DALA", "134", "A"]
                    wtres, numres, mutres = re.split(r"(\d+)", mut)
                    # do not bother adding parentheses to wtres
                    # since it is not used in the resfile
                else:
                    # the mutated residue is non-canonical
                    # i.e. ["A134", "DALA]"]

                    # wtres, numres, voidstr = ["A", "134", ""]
                    wtres, numres, voidstr = re.split(r"(\d+)", mut[0])
                    # mutres = "[DALA]"
                    mutres = "X[" +  mut[1]
                
            elif len(mut) == 3:
                # both the wild-type and the mutated residues are
                # non-canonical
                # i.e. ["", "DALA]134", "DPHE]"]

                # wtandnum = "DALA134"
                wtandnum = mut[1].replace("]", "") 
                # wtres, numres, voidstr = ["DALA", "134", ""]
                wtres, numres, voidstr = re.split(r"(\d+)", wtandnum)
                # mutres = "[DPHE]"
                mutres = "X[" + mut[2]
                # do not bother adding parentheses to wtres since
                # it is not used in the resfile

            out.write(mutstring.format(wtres, numres, mutres))

    #-------------------- Setup the command line ---------------------#

    # set options for the inputs
    inputoptions = \
        [("-in:file:s", pdbfile), \
         ("-in:ignore_unrecognized_res", True), \
         ("-in:file:fullatom", True), \
         ("-in:path:database", database)]
    
    # set ddg-related options
    ddgoptions = \
        [("-ddg:cartesian", True), \
         ("-ddg:iterations", 3), \
         ("-ddg:bbnbrs", 1), \
         ("-ddg:dump_pdbs", True), \
         ("-ddg:mut_file", mutfile), \
         ("-ddg:out", ddgoutput)]

    # set packing-related options
    # additional residues to be considered while packing
    if ncaalistfile is not None:
        packingoptions = \
            [("-packing:packer_palette:extra_base_type_file", \
              ncaalistfile)]
    else:
        packingoptions = []

    # set scoring-related options
    scoringoptions = \
        [("-score:fa_max_dis", 9.0)]

    # set protocol-related options
    if protocol == "cartddg_talaris2014":
        protocoloptions = [("-score:weights", "talaris2014_cart"), \
                            ("-restore_talaris_behavior", True)]
    elif protocol == "cartddg_ref2015":
        protocoloptions = [("-score:weights", "ref2015_cart")]

    # all options
    options = \
        inputoptions + ddgoptions + scoringoptions + \
        packingoptions + protocoloptions

    #------------------------------ Run ------------------------------#

    cproc, rundir = runrosettafunc(executable = executable, \
                                   rosettaoutput = rosettaoutput, \
                                   options = options)

    # go to the working directory
    os.chdir(wd)

    # return the results of the run
    return cproc, rundir


def run_flexddg(mutdef, \
                wd, \
                runrosettafunc, \
                executable, \
                rosettascript, \
                database, \
                protocol, \
                pdbfile, \
                resfile, \
                ncaalistfile, \
                ddgdbfile, \
                structdbfile, \
                rosettaoutput):
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
    if not os.path.exists(resfile):
        # do not overwrite the resfile if you already have it
        with open(resfile, "w") as out:
            # set the header and start line
            out.write("{:s}\n{:s}".format("NATAA", "start"))      
            # define the mutation string
            fstring = "\n{:s} {:s} PIKAA {:s}"
            
            # each singlemut in the form ("A", "C151Y")
            for singlemut in muts:
                # chain, mut = ("A", "C151Y")
                chain, mut = singlemut
                # find any definition including non-canonical residues
                # i.e. "X[DALA]134A" becomes ["", "DALA]134A"]
                # i.e. "A134X[DALA]" becomes ["A134", "DALA]"]
                # i.e. "X[DALA]134X[DPHE] becomes 
                #      ["", "DALA]134", "DPHE]"]
                mut = mut.split("X[")
                
                # handle the different scenarios
                if len(mut) == 1:
                    # only canonical residues, no splitting happened
                    wtres, numres, mutres = re.split(r"(\d+)", mut[0])
                
                elif len(mut) == 2:
                    # either the wild-type or the mutated residue
                    # is non-canonical
                    # i.e. ["", "DALA]134A"] or ["A134", "DALA]"]
                    if mut[0] == "":
                        # the wild-type residue is non-canonical
                        # i.e. ["", "DALA]134A"]
                        
                        # mut = ["", "DALA134A"]
                        mut = mut[1].replace("]", "")
                        # wtres, numres, mutres = ["DALA", "134", "A"]
                        wtres, numres, mutres = re.split(r"(\d+)", mut)
                        # do not bother adding parentheses to wtres
                        # since it is not used in the resfile
                    else:
                        # the mutated residue is non-canonical
                        # i.e. ["A134", "DALA]"]

                        # wtres, numres, voidstr = ["A", "134", ""]
                        wtres, numres, voidstr = \
                            re.split(r"(\d+)", mut[0])
                        # mutres = "[DALA]"
                        mutres = "X[" +  mut[1]
                
                elif len(mut) == 3:
                    # both the wild-type and the mutated residues are
                    # non-canonical
                    # i.e. ["", "DALA]134", "DPHE]"]

                    # wtandnum = "DALA134"
                    wtandnum = mut[1].replace("]", "") 
                    # wtres, numres, voidstr = ["DALA", "134", ""]
                    wtres, numres, voidstr = \
                        re.split(r"(\d+)", wtandnum)
                    # mutres = "[DPHE]"
                    mutres = "X[" + mut[2]
                    # do not bother adding parentheses to wtres since
                    # it is not used in the resfile

                out.write(fstring.format(numres, chain, mutres))


    #-------------------- Setup the command line ---------------------#

    # set options for the inputs
    inputoptions = \
        [("-in:file:s", pdbfile), \
         ("-in:ignore_unrecognized_res", True), \
         ("-in:file:fullatom", True), \
         ("-in:path:database", database)]

    # set packing-related options
    packingoptions = \
        [("-packing:ex1", True), \
         ("-packing:ex2", True)]
    # additional residues to be considered while packing
    if ncaalistfile is not None:
        packingoptions.append(\
            ("-packing:packer_palette:extra_base_type_file", ncaalistfile))

    # set script parser-related options
    parseroptions = \
        [("-parser:protocol", rosettascript)]

    # set general substitutions for the script variables
    generalscriptvars = \
        [("chaintomove", chaintomove), \
         ("resfile", resfile), \
         ("ddgdbfile", ddgdbfile), \
         ("structdbfile", structdbfile)]

    # set protocol-related options and variable substitutions
    if protocol == "flexddg_talaris2014":
        protocolscriptvars = [("sfweights", "talaris2014")]
        protocoloptions = \
            [("-corrections:restore_talaris_behavior", True)]
    
    elif protocol == "flexddg_ref2015":
        protocolscriptvars = [("sfweights", "ref2015")]
        protocoloptions = []

    options = \
        inputoptions + packingoptions + parseroptions + protocoloptions

    scriptvars = generalscriptvars + protocolscriptvars

    #------------------------------ Run ------------------------------#

    cproc, rundir = runrosettafunc(executable = executable, \
                                   rosettaoutput = rosettaoutput, \
                                   options = options, \
                                   scriptvars = scriptvars)

    # go back to the working directory
    os.chdir(wd)

    # return the results of the run
    return cproc, rundir


def get_abspath(path):
    """Get the absolute path of an object. Return `None`if `None`
    was passed.
    """
    
    return os.path.abspath(path) if path is not None else path


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


def get_mutlist_cartddg(listfile, reslist = None):
    """Generate a list of mutations for the cartesian_ddg-based
    protocols from a file containing either a list of selected
    mutations or a list of residue positions to perform saturation
    mutagenesis. In this case, a list of residue types to be used
    for the mutagenesis must be passed.

    Example of file containing a list of selected mutations:

    C151Y,S154N
    S2A

    Example of file containing a list of positions:

    C151
    S154
    """
    with open(listfile, "r") as f:
        rawlist = []
        for line in f:
            if not re.match(r"^\s*$", line):
                # ignore empty lines
                items = tuple(line.rstrip("\n").split(","))
                rawlist.append(items)

        if reslist is None:
            # return the list of selected mutations
            return rawlist
        else:
            mutlist = []
            for item in rawlist:
                if len(item) > 1:
                    errstr = \
                        "Error in parsing {:s}. Each line " \
                        "must contain only one position."
                    raise ValueError(errstr.format(listfile))  

                # map each position to all residue types
                # included in the saturation mutagenesis
                muts = [tuple((item[0] + res,)) for res in reslist]
                mutlist.extend(muts)
            
            return mutlist


def get_mutlist_flexddg(listfile, reslist = None):
    """Generate a list of mutations for the Flex ddG-based protocols
    from a file containing either a list of selected mutations or a
    list of residue positions to perform saturation mutagenesis.
    In this case, a set of residue types to be used for the
    mutagenesis must be passed.

    Example of file containing a list of selected mutations:

    A-C151Y,A-S154N A
    A-S2A A

    Example of file containing a list of positions:
    
    A-C151 A
    A-S154 A
    """
    
    with open(listfile, "r") as f:
        rawlist = []
        for line in f:
            if not re.match(r"^\s*$", line):
                # ignore empty lines
                muts, chaintomove = line.rstrip("\n").split(" ")
                muts = muts.split(",")
                chainsmuts = \
                    tuple([tuple(mut.split("-")) for mut in muts])
                rawlist.append(tuple([chainsmuts, chaintomove]))

        if reslist is None:
            # return the list of selected mutations
            mutlist = rawlist
        else:
            mutlist = []
            # ((("A", "C151"),), "A")
            for item, chaintomove in rawlist:
                if len(item) > 1:
                    errstr = \
                        "Error in parsing {:s}. Each line " \
                        "must contain only one position."
                    raise ValueError(errstr.format(listfile))  
                
                # map each position to all residue types
                # included in the saturation mutagenesis
                # "A", "C151" = (("A", "C151"),)[0]
                chain, pos = item[0]
                # each mutation will be in the form:
                #   ((("A", "C151Y"),), "A"))
                # make sure that chain and position are together
                # in a tuple of length one, and not separated in
                # a tuple of length 2
                muts = \
                    [tuple(\
                        [tuple(\
                            [tuple([chain, pos + res]),]), \
                        chaintomove]) \
                     for res in reslist]
                mutlist.extend(muts)

        return mutlist


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
                return line.split(" ")


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
    
    ######################### ARGUMENT PARSER #########################

    description = \
        "\nScript to run Rosetta protocols for the prediction of " \
        "the ΔΔG of stability upon mutation of a monomeric protein " \
        "and the ΔΔG of binding upon mutation of a protein complex.\n"

    # create the argument parser
    parser = argparse.ArgumentParser(description = description)
    general_args = parser.add_argument_group("General arguments")


    #---------------------------- general ----------------------------#
    
    f_helpstr = "PDB file of the input wild-type structure."
    general_args.add_argument("-f", "--pdbfile", \
                              type = str, \
                              required = True, \
                              help = f_helpstr)

    p_helpstr = "Name of the protocol to be run."
    general_args.add_argument("-p", "--protocol", \
                              type = str, \
                              required = True, \
                              help = p_helpstr)

    d_helpstr = \
        "Working directory. Default is the current working directory."
    general_args.add_argument("-d", "--workdir", \
                              type = str, \
                              required = False, \
                              default = os.getcwd(), \
                              help = d_helpstr)

    rosettapath_helpstr = \
        "Path to Rosetta installation directory."
    general_args.add_argument("-r", "--rosettapath", \
                              type = str, \
                              metavar = "rosettapath", \
                              required = True, \
                              help = rosettapath_helpstr)

    l_helpstr = \
        "File containing the list of selected mutations " \
        "or positions to perform saturation mutagenesis"
    general_args.add_argument("-l", "--listfile", \
                              type = str, \
                              required = False, \
                              default = None, \
                              help = l_helpstr)

    step_helpstr = \
        "Which step of the protocol to run, since some protocols " \
        "consist of multiple steps that must be run separately."
    general_args.add_argument("-s", "--step", \
                              type = str, \
                              required = False, \
                              default = None, \
                              help = step_helpstr)

    nproc_helpstr = \
        "Number of processes to be started in parallel. " \
        "Default is one process."
    general_args.add_argument("-n", "--nproc", \
                              dest = "nproc",
                              type = int, \
                              required = False, \
                              default = 1, \
                              help = nproc_helpstr)

    logfile_default = "run_ddg_protocols.log"
    logfile_helpstr = \
        "Name of the log file (default: {:s}).".format(logfile_default)
    general_args.add_argument("-lf", "--logfile", \
                              dest = "logfile", \
                              required = False, \
                              default = logfile_default, \
                              help = logfile_helpstr)

    v_helpstr = "Verbosity level: verbose."
    general_args.add_argument("-v", \
                              action = "store_true", \
                              dest = "v", \
                              help = v_helpstr)

    vv_helpstr = "Verbosity level: debug."
    general_args.add_argument("-vv", \
                              dest = "vv", \
                              action = "store_true", \
                              help = vv_helpstr)

    rosettascript_helpstr = "XML Rosetta script."
    general_args.add_argument("--rosettascript", \
                              dest = "rosettascript", \
                              type = str, \
                              required = False, \
                              default = None, \
                              help = rosettascript_helpstr)


    #-------------------- saturation mutagenesis ---------------------#

    saturation_helpstr = \
        "Perform saturation mutagenesis on selected positions."
    general_args.add_argument("--saturation", \
                              dest = "saturation" , \
                              action = "store_true", \
                              help = saturation_helpstr)

    r_helpstr = \
        "File containing the list of residue types " \
        "(one-letter name) to include in the saturation " \
        "mutagenesis. It is used only if --saturation is provided."
    general_args.add_argument("--reslistfile", \
                              type = str, \
                              required = False, \
                              default = None, \
                              help = r_helpstr)


    #-------------------- Non-canonical residues ---------------------#
    
    ncaalistfile_helpstr = \
        "File containing a list of non-canonical residue types " \
        "(full names) to be considered. It is needed if you " \
        "want to specify mutations to non-canonical residues."
    general_args.add_argument("--ncaalistfile", \
                              type = str, \
                              required = False, \
                              default = None, \
                              help = ncaalistfile_helpstr)

    # parse the arguments
    args = parser.parse_args()
    # general arguments
    pdbfile = get_abspath(args.pdbfile)
    protocol = args.protocol
    rosettapath = get_abspath(args.rosettapath)
    workdir = get_abspath(args.workdir)
    step = args.step
    listfile = get_abspath(args.listfile)
    nproc = args.nproc
    v = args.v
    vv = args.vv
    logfilename = get_abspath(args.logfile)
    rosettascript = get_abspath(args.rosettascript)
    # saturation mutagenesis arguments
    saturation = args.saturation
    reslistfile = get_abspath(args.reslistfile)
    # NCAA-related arguments
    ncaalistfile = get_abspath(args.ncaalistfile)

    # path where Rosetta executables are
    rosettabin = os.path.join(rosettapath, "main/source/bin")

    # path where Rosetta database is
    rosettadatabase = os.path.join(rosettapath, "main/database")

    # executable names for currently implemented protocols
    execname2execpath = \
        {"relax" : None,  \
         "cartesian_ddg" : None, \
         "rosetta_scripts" : None}

    # find the path to each executable
    for name in execname2execpath.keys():
        for item in os.listdir(rosettabin):
            if item.startswith(name):
                execname2execpath[name] = \
                    os.path.join(rosettabin, item)
                break


    ############################ PROTOCOLS ############################

    # Currently implemented protocols. Protocols are grouped by aim
    # (predicting the ΔΔG of stability or binding) and then by family
    # (each family contains variants of the same protocol). Each
    # protocol is associated with the list of steps it is composed of
    # or None, in case it consists only of one step.
    protocols = \
        {"cartddg" : \
                {"cartddg_talaris2014" : ["relax", "cartesian_ddg"], \
                 "cartddg_ref2015" : ["relax", "cartesian_ddg"]}, \
        "flexddg": \
                {"flexddg_talaris2014" : None, \
                 "flexddg_ref2015" : None}, \
        }

    # which protocol steps are actually performing the ΔΔG prediction
    # (some of them may be preprocessing the structure, and need to
    # be run separately). If the protocol has only one step, that is
    # the ΔΔG prediction.
    ddgpredictionsteps = ["cartesian_ddg", None]

    # turn on/off the flag that signals that ΔΔG prediction has been
    # requested
    runddgprediction = \
        True if step in ddgpredictionsteps else False


    ########################## LOGGER SETUP ###########################

    # set the log level (default is only warnings and errors)
    loglevel = log.WARNING
    if v:
        loglevel = log.INFO
    if vv:
        # if both -v and -vv are passed, -vv is used
        loglevel = log.DEBUG
    
    logfile = os.path.join(os.getcwd(), logfilename)
    logfmt = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    datefmt = "%d-%b-%y %H:%M:%S"
    log.basicConfig(filename = logfile, \
                    filemode = "w", \
                    level = loglevel, \
                    format = logfmt, \
                    datefmt = datefmt)


    ########################## GENERAL SETUP ##########################

    # go to the correct working directory
    os.chdir(workdir)
    # set a general string to store info about the runs
    runinfostr = \
        "Process in {:s} exited with code {:d}. Command line: {:s}"


    ##################### MUTATION LIST GENERATION ####################

    # the list of mutations is needed only if calculating the ΔΔG
    if runddgprediction:
        # initialize reslist to None
        reslist = None
        # if saturation mutagenesis has been requested
        if saturation:
            # get the residue types set
            reslist = get_reslist(reslistfile = reslistfile)
            # log the residue types that will be used for the
            # saturation mutagenesis
            logstr = "Residue types included in the saturation " \
                     "mutagenesis: \n{:s}"
            log.info(logstr.format("\n".join(reslist)))

        if protocol in protocols["cartddg"]:
            # generate the list of single mutations
            mutlist = get_mutlist_cartddg(listfile, reslist)
            
        elif protocol in protocols["flexddg"]:
            # number fo structures in the ensemble
            nstruct = 35
            # generate the list of single mutations
            tmpmutlist = get_mutlist_flexddg(listfile, reslist)
            mutlist = []
            # repeat each mutation as many times as the number of
            # structures constituting the final ensemble (each
            # structure will be generated via a separate run)
            for mut in tmpmutlist:
                mutlist.extend(\
                    [(mut, struct) for struct in range(1, nstruct+1)])

        ### other elif statements if more protocols are implemented ###

    
    ######################## CARTDDG PROTOCOLS ########################

    if protocol in protocols["cartddg"]:

        #----------------------- Check inputs ------------------------#
            
        # ensure that the user defined which step to run
        if step is None:
            errstr = \
                "You must specify the step you want to run. Steps " \
                "are 'relax' and 'cartesian_ddg'."
            raise ValueError(errstr.format(protocol))

        #--------------------------- relax ---------------------------#

        if step == "relax":
            # prefix for the names of the output structures
            outprefix = "relaxed_"
            # name of the scorefile
            scorefile = "score.sc"
            # name of the relaxation directory to be created
            relaxdir = "relax"
            # name of the Rosetta output
            rosettaoutput = "relax.out"

            # run the relaxation
            cproc, rundir = \
                run_cartddg_relax(wd = workdir, \
                                  runrosettafunc = run_rosetta, \
                                  executable = execname2execpath["relax"], \
                                  protocol = protocol, \
                                  rosettaoutput = rosettaoutput, \
                                  pdbfile = pdbfile, \
                                  database = rosettadatabase, \
                                  outprefix = outprefix, \
                                  scorefile = scorefile, \
                                  dirname = relaxdir)
                
            # format the information about the run
            runinfostrfilled = \
                runinfostr.format(rundir, cproc.returncode, \
                                " ".join(cproc.args))

            if cproc.returncode == 0:
                # inform the user about the run
                log.info(runinfostrfilled)
            else:
                # warn the user since something went wrong
                log.error(runinfostrfilled)
                # exit the script (if something went wrong in)
                # the relaxation it does not make sense to proceed
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


        #----------------------- cartesian_ddg -----------------------#

        elif step == "cartesian_ddg":
            # name of the mutation file to be created for each
            # mutation
            mutfile = "mutation.mutfile"
            # name of the output file storing the scores
            ddgoutput = "ddg_scores.out"
            # name of the Rosetta output file
            rosettaoutput = "cartesian_ddg.out"

            partfunc = \
                functools.partial(\
                    run_cartddg_cartesian_ddg, \
                    wd = workdir, \
                    runrosettafunc = run_rosetta, \
                    executable = execname2execpath["cartesian_ddg"], \
                    protocol = protocol, \
                    database = rosettadatabase, \
                    pdbfile = pdbfile, \
                    mutfile = mutfile, \
                    ncaalistfile = ncaalistfile, \
                    ddgoutput = ddgoutput, \
                    rosettaoutput = rosettaoutput)


    ######################## FLEXDDG PROTOCOLS ########################

    elif protocol in protocols["flexddg"]:     
        # check the XML Rosetta script
        if rosettascript is None:
            errstr = \
                "You must provide the XML script describing " \
                "the procedure."
            raise ValueError(errstr)

        # name of the mutation file (in resfile format) to be
        # created for each mutation
        resfile = "mutation.resfile"
        # name of the Rosetta output file
        rosettaoutput = "rosettascripts.out"
        # name of the database file containing on the ΔΔG scores
        ddgdbfile = "ddg.db3"
        # name of the database file containing info on the structures
        structdbfile = "struct.db3"

        partfunc = \
            functools.partial(\
                run_flexddg, \
                wd = workdir, \
                runrosettafunc = run_rosetta, \
                executable = execname2execpath["rosetta_scripts"], \
                database = rosettadatabase, \
                rosettascript = rosettascript, \
                protocol = protocol, \
                pdbfile = pdbfile, \
                resfile = resfile, \
                ncaalistfile = ncaalistfile, \
                ddgdbfile = ddgdbfile, \
                structdbfile = structdbfile, \
                rosettaoutput = rosettaoutput)


    ########################## ΔΔG PREDICTION #########################

    # Protocol-agnostic way of sending jobs for ΔΔG predictions.
    #
    # It requires:
    # - the function running the protocol with all parameters already
    #   set ("partfunc" created with functools.partial or
    #   equivalent) except the mutation (whose format depends only on
    #   the function, since it will be parsed internally);
    # - that the function running the protocol returns a tuple
    #   containing a CompletedInstance of the process ("cproc") and
    #   the running directory ("rundir")

    if runddgprediction:
        # create the pool of workers
        pool = mp.Pool(nproc)
        # empty list to store the return codes
        returncodes = []
        # run the predictions in parallel keeping the order of the
        # inputs (with imap() the iterable is not split into chunks)
        for cproc, rundir in pool.imap(partfunc, mutlist):
            # format the information about the run
            runinfostrfilled = \
                runinfostr.format(rundir, cproc.returncode, \
                                " ".join(cproc.args))

            # save the return code for later
            returncodes.append(cproc.returncode)

            if cproc.returncode == 0:
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

