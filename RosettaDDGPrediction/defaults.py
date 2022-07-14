#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    defaults.py
#
#    Hard-coded defaults. Should only be changed
#    for development purposes.
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
import collections
import os.path
from pkg_resources import resource_filename, Requirement



############################ STEPS/OPTIONS ############################



# Recognized protocols
ROSETTA_PROTOCOLS = \
    {"flexddg" : \
        {"flexddg" : {"role" : "ddg",
                      "run_by" : "rosetta",
                      "executable" : "rosetta_scripts",
                      "executable_extract" : "score_jd2"}},
    "cartddg" : \
        {"relax" : {"role" : "processing",
                    "run_by" : "rosetta",
                    "executable" : "relax"},
         "structure_selection" : {"run_by" : "python"},
         "cartesian" : {"role" : "ddg",
                        "run_by" : "rosetta",
                        "executable" : "cartesian_ddg"}},
    "cartddg2020" : \
        {"relax2020" : {"role" : "processing",
                        "run_by" : "rosetta",
                        "executable" : "rosetta_scripts"},
         "structure_selection" : {"run_by" : "python"},
         "cartesian" : {"role" : "ddg",
                        "run_by" : "rosetta",
                        "executable" : "cartesian_ddg"}}}


# Possible ways of defining various Rosetta options that
# are needed to correctly retrieve files/parameters
ROSETTA_OPTIONS = \
    collections.defaultdict(str,
        {"ddg_out" : ("-ddg:out", "-out"),
         "in_pdb_file" : ("-in:file:s", "-s"),
         "mutfile" : ("-ddg:mut_file", "-mut_file"),
         "out_prefix" : ("-out:prefix", "-prefix"),
         "out_suffix" : ("-out:suffix", "-suffix"),
         "script_vars" : ("-parser:script_vars", "-script_vars"),
         "scorefile" : ("-out:file:scorefile", "-scorefile"),
         "scf_name" : ("-score:weights", "-weights", "scfname"),
         "backrub_n_trials" : ("backrubntrials",),
         "backrub_traj_stride" : ("backrubtrajstride",),
         "ddg_db_file" : ("ddgdbfile",),
         "protocol" : ("-parser:protocol", "-protocol"),
         "resfile" : ("resfile",),
         "struct_db_file" : ("structdbfile",),
         "db_name" : ("-inout:dbms:database_name",)})


# Default name for the Rosetta crash log
ROSETTA_CRASH_LOG = "ROSETTA_CRASH.log"



########################## DIRECTORIES/FILES ##########################



# Directory containing configuration files for running the protocols
CONFIG_RUN_DIR = \
    resource_filename(Requirement("RosettaDDGPrediction"),
                      "RosettaDDGPrediction/config_run")

# Directory containing configuration files for aggregating results
CONFIG_AGGR_DIR = \
    resource_filename(Requirement("RosettaDDGPrediction"),
                      "RosettaDDGPrediction/config_aggregate")

# Directory containing configuration files for plotting
CONFIG_PLOT_DIR = \
    resource_filename(Requirement("RosettaDDGPrediction"),
                      "RosettaDDGPrediction/config_plot")

# Directory containing configuration files for run settings
CONFIG_SETTINGS_DIR = \
    resource_filename(Requirement("RosettaDDGPrediction"),
                      "RosettaDDGPrediction/config_settings")

# Directory containing RosettaScripts used in some protocols
ROSETTA_SCRIPTS_DIR = \
    resource_filename(Requirement("RosettaDDGPrediction"),
                      "RosettaDDGPrediction/RosettaScripts")

# Default configuration file for data aggregation
CONFIG_AGGR_FILE = os.path.join(CONFIG_AGGR_DIR, "aggregate.yaml")



############################## MUTATIONS ##############################



# Keywords used to refer to specific mutation attributes

# Mutation (maps to as many dictionaries as the single mutations
# defined)
MUT = "_mut_"

# Structure number/ID (in case the same mutation must be performed
# multiple times on the same structure to obtain an ensemble)
STRUCT = "_struct_"

# Wild-type residue in a single mutation
WTR = "_wtr_"

# Mutant residue in a single mutation
MUTR = "_mutr_"

# Residue number in a single mutation
NUMR = "_numr_"

# Chain ID in a single mutation
CHAIN = "_chain_"

# All attributes that defines the position of the residue
# (chain ID, wild-type residue and residue number)
POSR = "_nomutr_"

# Mutation directory path
MUT_DIR_PATH = "_dirpath_"

# Mutation directory name
MUT_DIR_NAME = "_dirname_"

# Separator for mutations that are performed simultaneously    
MUT_SEP = ","

# Separator for the components of a mutation
# (e.g. C.L.53.E for mutating L53 on chain C to E)
COMP_SEP = "."

# Separator for chain ID and the rest of the mutation (wild-type
# residue, residue number and mutant residue) in the directory name
CHAIN_SEP = "-"

# Separator for single mutations in the directory name (in case of
# multiple simultaneous mutations)
DIR_MUT_SEP = "_"

# Separator for multiple mutations in the mutinfo file
MULTI_MUT_SEP = ":"



######################### STRUCTURE EXTRACTION ########################



# 'States' of the structures to be extracted from the database file
# generated by flexddg protocols (wild-type, mutant, or part of the
# backrub trajectory)
FLEXDDG_STATES = ("backrub", "wt", "mut")

# Pattern describing the extracted structures before renaming
STRUCT_EXTRACTED_PATTERN = r"(\d+)_0001.pdb"


############################# AGGREGATION #############################



# Column names in the mutations' info file
MUTINFO_COLS = {"mut_name" : "mut_name",
                "dir_name" : "dir_name",
                "mut_label" : "mut_label",
                "pos_label" : "pos_label"}

# Column names in Rosetta .db3 files and aggregated dataframes
ROSETTA_DF_COLS = {"b_steps" : "backrub_steps",
                   "energy_unit" : "energy_unit",
                   "mutation" : "mutation",
                   "mut_label" : "mutation_label",
                   "pos_label" : "position_label",
                   "name" : "name",
                   "n_struct" : "nstruct",
                   "scf_name" : "score_function_name",
                   "sc_type" : "score_type_name",
                   "sc_value" : "score_value",
                   "state" : "state",
                   "struct_id" : "struct_id",
                   "struct_num" : "struct_num",
                   "tot_score" : "total_score"}



############################### PLOTTING ##############################



# Valid plot types
PLOT_TYPES = ["contributions_barplot", "dg_swarmplot",
              "total_heatmap", "total_heatmap_saturation"]
