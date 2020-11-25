#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    aggregation.py
#
#    Utility functions for data aggregation.
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
import os.path
import sqlite3
# third-party packages
import pandas as pd
# RosettaDDGProtocols
from .defaults import ROSETTADFCOLS



def parse_output_cartddg(ddgout, \
                         listcontributions, \
                         scfname):
    """Parse the output file from cartddg protocols."""

    # column names
    scfnamecol = ROSETTADFCOLS["scfname"]
    statecol = ROSETTADFCOLS["state"]
    structnumcol = ROSETTADFCOLS["structnum"]
    totscorecol = ROSETTADFCOLS["totscore"]
    
    # state names  
    wt = "wt"
    mut = "mut"

    with open(ddgout, "r") as f:
        # create empty lists to store the wild-type and mutants ΔGs
        wtdgs = []
        mutdgs = []
        # for each line
        for line in f:
            if line.startswith("COMPLEX:"):
                # ddg_stability output
                line = \
                    [item for item in line.strip("\n").split(" ") \
                     if item != ""]
                # WT or MUT
                mutstatus = line[2]
                # total ΔG
                dg = float(line[3])
                # energy contributions
                contributions = [float(item) for item in line[5::2]]
                # store the data in a dictionary
                dgdata = dict(zip(listcontributions, contributions))
                # update the dictionary with the total score
                dgdata.update({totscorecol : dg})
                # append it to the right list according to it
                # referring to a wild-type or to a mutant structure    
                if mutstatus.startswith("WT"):
                    wtdgs.append(dgdata)
                elif mutstatus.startswith("MUT"):
                    mutdgs.append(dgdata)
        
        # the protocol must have run the same number
        # of rounds for both WT and MUT
        if len(wtdgs) != len(mutdgs):
            errstr = \
                f"The number of rounds run for the wild-type " \
                f"structure must be equal to those run for the " \
                f"mutant, while the file you provided ({ddgoutfile}) " \
                f"contains {len(wtdgs)} rounds for the wild-type " \
                f"and {len(mutdgs)} for the mutant. " \
                f"Please check your run."
            raise ValueError(errstr)
        
        # get the number of structures generated
        nstructs = [str(i) for i in range(1, len(wtdgs)+1)]
        
        # scores for wild-type structures
        dfwt = pd.DataFrame(wtdgs)
        dfwt[statecol] = wt
        dfwt[structnumcol] = nstructs
        dfwt[scfnamecol] = scfname
        
        # scores for mutant structures
        dfmut = pd.DataFrame(mutdgs)
        dfmut[statecol] = mut
        dfmut[structnumcol] = nstructs
        dfmut[scfnamecol] = scfname
        
        # concatenate the results
        return pd.concat([dfwt, dfmut]).reset_index(drop = True)


def aggregate_data_cartddg(df, \
                           listcontributions):
    """Aggregate data for cartddg protocols."""

    # column names
    statecol = ROSETTADFCOLS["state"]
    structnumcol = ROSETTADFCOLS["structnum"]
    totscorecol = ROSETTADFCOLS["totscore"]
    mutationcol = ROSETTADFCOLS["mutation"]

    # columns storing energy scores
    scorecols = listcontributions + [totscorecol]

    # create a dataframe with the wild-type ΔG scores
    dgwt = df.loc[df[statecol] == "wt"].reset_index(drop = True)

    # create a dataframe with the mutant ΔG scores
    dgmut = df.loc[df[statecol] == "mut"].reset_index(drop = True)

    # create a dataframe with the ΔΔG scores
    ddg = dgmut.copy()
    ddg[scorecols] = dgmut[scorecols] - dgwt[scorecols]
    ddg[statecol] = "ddg"

    # return all dataframes
    return (dgwt, dgmut, ddg)


def parse_output_flexddg(db3out, \
                         trajstride, \
                         structnum, \
                         scfname):
    """Parse the .db3 output from flexddg protocols."""

    # column names
    structidcol = ROSETTADFCOLS["structid"]
    namecol = ROSETTADFCOLS["name"]
    statecol = ROSETTADFCOLS["state"]
    scfnamecol = ROSETTADFCOLS["scfname"]
    bstepscol = ROSETTADFCOLS["bsteps"]
    sctypecol = ROSETTADFCOLS["sctype"]
    scvaluecol = ROSETTADFCOLS["scvalue"]
    structnumcol = ROSETTADFCOLS["structnum"]
    
    # open the connection
    db3out = os.path.abspath(db3out)
    connection = sqlite3.connect(db3out)
    connection.row_factory = sqlite3.Row
    
    # create the cursor
    cursor = connection.cursor()
    
    # selection string to get the total number of batches
    selstr = f'SELECT max(batch_id) from batches'
    nbatches = cursor.execute(selstr).fetchone()[0]
    # selection string for the structure scores
    sel = \
        "batches.name, structure_scores.struct_id," \
        "score_types.score_type_name, structure_scores.score_value," \
        "score_function_method_options.score_function_name"
    # selection string for the inner join on batches
    injoinbatches = \
        "batches.batch_id=structure_scores.batch_id"
    # selection string for the inner join on the scoring function
    injoinscfun = \
        "score_function_method_options.batch_id=batches.batch_id"
    # selection string for the inner join on the score type
    injoinsctype = \
        "score_types.batch_id=structure_scores.batch_id " \
        "AND score_types.score_type_id=structure_scores.score_type_id"
    
    # assemble the query string
    query = \
        f"SELECT {sel} from structure_scores\n" \
        f"INNER JOIN batches ON {injoinbatches}\n" \
        f"INNER JOIN score_function_method_options ON {injoinscfun}\n" \
        f"INNER JOIN score_types ON {injoinsctype}"
    
    # read the query into a dataframe
    df = pd.read_sql_query(query, connection)
    
    # function to renumber the structure IDs
    getnewid = lambda x: trajstride * (1 + (int(x - 1) // nbatches))
    # function to rename the structures
    getnewname = lambda x: x.replace("_dbreport", "") \
                 if x.endswith('_dbreport') else x
    
    # renumber structure IDs
    df[bstepscol] = df[structidcol].apply(getnewid)
    # rename structures
    df[statecol] = df[namecol].apply(getnewname)
    
    # convert into a pivot table
    df = df.pivot_table(index = [statecol, bstepscol, scfnamecol], \
                        columns = sctypecol, \
                        values = scvaluecol).reset_index()
    
    # add column for the structure number
    df[structnumcol] = structnum
    # add complete score function name to the
    # score function name column
    df[scfnamecol] = scfname
    # remove unnecessary column names
    df.columns.names = [None]
    
    # close the connection
    connection.close()

    # return the dataframe
    return df


def aggregate_data_flexddg(df, \
                           listcontributions):
    """Aggregate data for flexddg protocols."""

    # column names
    statecol = ROSETTADFCOLS["state"]
    scfnamecol = ROSETTADFCOLS["scfname"]
    bstepscol = ROSETTADFCOLS["bsteps"]
    structnumcol = ROSETTADFCOLS["structnum"]
    nstructcol = ROSETTADFCOLS["nstruct"]
    totscorecol = ROSETTADFCOLS["totscore"]

    # get the score function name
    scfname = df[scfnamecol].unique()[0]
    
    # columns storing energy scores
    scorecols = listcontributions + [totscorecol]
    
    # select only the rows corresponding to the maximum number
    # of backrub steps performed
    bstepscond = (df[bstepscol] == df[bstepscol].max())
    
    # generate selections on binding and mutation states
    ubwt = df.loc[(df[statecol] == "unbound_wt") & bstepscond]
    ubmut = df.loc[(df[statecol] == "unbound_mut") & bstepscond]
    bwt = df.loc[(df[statecol] == "bound_wt") & bstepscond]
    bmut = df.loc[(df[statecol] == "bound_mut") & bstepscond]
    
    # get the scores (reset index to make it uniform to be able
    # to perform the subtraction later)
    ubwtscores = ubwt[scorecols].reset_index(drop = True)
    ubmutscores = ubmut[scorecols].reset_index(drop = True)
    bwtscores = bwt[scorecols].reset_index(drop = True)
    bmutscores = bmut[scorecols].reset_index(drop = True)
    
    # get the list of structures from one of the dataframes
    structnums = list(ubwt[structnumcol])
    
    # create a dataframe with the wild-type ΔG scores
    dgwt = bwtscores - ubwtscores
    
    # create a dataframe with the mutant ΔG scores
    dgmut = bmutscores - ubmutscores
    
    # create a copy of the dataframe containing the mutant ΔG scores
    ddg = dgmut.copy()
    # subtract the wild-type ΔG scores to obtain the ΔΔG scores
    ddg[scorecols] = dgmut[scorecols] - dgwt[scorecols]
    
    # add a column to all dataframes with the state
    dgwt[statecol], dgmut[statecol], ddg[statecol] = "wt", "mut", "ddg"
    # add a column to all dataframes with the score function name
    dgwt[scfnamecol], dgmut[scfnamecol], ddg[scfnamecol] = \
        scfname, scfname, scfname
    # add a column to all dataframes with the structure numbers
    dgwt[structnumcol], dgmut[structnumcol], ddg[structnumcol] = \
        structnums, structnums, structnums

    # return the dataframes
    return (dgwt, dgmut, ddg)


def generate_output_dataframes(dgwt, \
                               dgmut, \
                               ddg, \
                               mutname, \
                               rescale, \
                               listcontributions, \
                               convfact):
    
    # columns names
    mutationcol = ROSETTADFCOLS["mutation"]
    statecol = ROSETTADFCOLS["state"]
    totscorecol = ROSETTADFCOLS["totscore"]
    energyunitcol = ROSETTADFCOLS["energyunit"]
    structnumcol = ROSETTADFCOLS["structnum"]
    scfnamecol = ROSETTADFCOLS["scfname"]

    # columns storing energy scores
    scorecols = [totscorecol] + listcontributions

    # columns of the structures dataframe
    structdfcols = [mutationcol, structnumcol, statecol, \
                    energyunitcol, scfnamecol, *scorecols]

    # columns of the aggregate dataframe
    aggrdfcols = [mutationcol, statecol, \
                  energyunitcol, scfnamecol, *scorecols]

    #------------------------ Structures data ------------------------#

    # concatenate the two dataframes
    structdf = pd.concat([dgwt, dgmut, ddg])
    # add the column with the mutation name as first column
    structdf.insert(0, mutationcol, mutname)

    #------------------------- Aggregate data ------------------------#

    # group by mutation and state
    groupby = [mutationcol, statecol, scfnamecol]
    aggrdf = \
        structdf.groupby(groupby).mean().reset_index()

    #----------------------------- Rescale ---------------------------#

    # default energy units are Rosetta Energy Units
    energyunit = "REUs"
    if rescale:
        # rescale the ΔΔG scores to kcal/mol.
        aggrdf[scorecols] = aggrdf[scorecols] * 1.0/convfact
        structdf[scorecols] = structdf[scorecols] * 1.0/convfact
        # change the energy unit
        energyunit = "kcal/mol"

    #------------------------ Add units column -----------------------#

    # add the energy units columns
    aggrdf[energyunitcol] = energyunit
    structdf[energyunitcol] = energyunit

    #-------------------------- Sort columns -------------------------#

    aggrdf = aggrdf[aggrdfcols]
    structdf = structdf[structdfcols]

    # return the structures data and the aggregate data
    return (aggrdf, structdf)