#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    aggregation.py
#
#    Utility functions for data aggregation.
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
import os.path
import sqlite3
# Third-party packages
import pandas as pd
# RosettaDDGProtocols
from .defaults import ROSETTA_DF_COLS



def parse_output_cartddg(ddg_out,
                         list_contributions,
                         scf_name):
    """Parse the output file from cartddg protocols and
    return a data frame.
    """

    # Get the column names
    scf_name_col = ROSETTA_DF_COLS["scf_name"]
    state_col = ROSETTA_DF_COLS["state"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]
    
    # State names  
    wt = "wt"
    mut = "mut"

    with open(ddg_out, "r") as f:
        
        # Create empty lists to store the wild-type and mutants ΔGs
        wt_dgs = []
        mut_dgs = []
        
        # For each line
        for line in f:
            
            # If the line contains useful information
            if line.startswith("COMPLEX:"):
                
                # Split the line and remove empty elements
                line = \
                    [item for item in line.strip("\n").split(" ") \
                     if item != ""]
                
                # WT or MUT
                mut_status = line[2]
                
                # Total ΔG
                dg = float(line[3])
                
                # Energy contributions
                contributions = [float(item) for item in line[5::2]]
                
                # Store the data in a dictionary
                dg_data = dict(zip(list_contributions, contributions))
                
                # Update the dictionary with the total score
                dg_data.update({tot_score_col : dg})
                
                # Append it to the right list according to it
                # referring to a wild-type or to a mutant structure    
                if mut_status.startswith("WT"):
                    wt_dgs.append(dg_data)
                elif mut_status.startswith("MUT"):
                    mut_dgs.append(dg_data)
        
        # The protocol must have run the same number
        # of rounds for both WT and MUT
        if len(wt_dgs) != len(mut_dgs):
            errstr = \
                f"The number of rounds run for the wild-type " \
                f"structure must be equal to those run for the " \
                f"mutant, while the file you provided ({ddg_out}) " \
                f"contains {len(wt_dgs)} rounds for the wild-type " \
                f"and {len(mut_dgs)} for the mutant. " \
                f"Please check your run."
            raise ValueError(errstr)
        
        # Get the number of structures generated
        n_structs = [str(i) for i in range(1, len(wt_dgs)+1)]
        
        # Scores for wild-type structures
        df_wt = pd.DataFrame(wt_dgs)
        df_wt[state_col] = wt
        df_wt[struct_num_col] = n_structs
        df_wt[scf_name_col] = scf_name
        
        # Scores for mutant structures
        df_mut = pd.DataFrame(mut_dgs)
        df_mut[state_col] = mut
        df_mut[struct_num_col] = n_structs
        df_mut[scf_name_col] = scf_name
        
        # Concatenate the results in a data frame
        return pd.concat([df_wt, df_mut]).reset_index(drop = True)


def aggregate_data_cartddg(df,
                           list_contributions):
    """Aggregate data for cartddg protocols.
    """

    # Get the column names
    state_col = ROSETTA_DF_COLS["state"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]
    mutation_col = ROSETTA_DF_COLS["mutation"]

    # Columns storing energy scores
    score_cols = list_contributions + [tot_score_col]

    # Create a data frame with the wild-type ΔG scores
    dg_wt = df.loc[df[state_col] == "wt"].reset_index(drop = True)

    # Create a data frame with the mutant ΔG scores
    dg_mut = df.loc[df[state_col] == "mut"].reset_index(drop = True)

    # Create a data frame with the ΔΔG scores
    ddg = dg_mut.copy()
    ddg[score_cols] = dg_mut[score_cols] - dg_wt[score_cols]
    ddg[state_col] = "ddg"

    # Return all dataframes
    return (dg_wt, dg_mut, ddg)


def parse_output_flexddg(db3_out,
                         traj_stride,
                         struct_num,
                         scf_name):
    """Parse the .db3 output from flexddg protocols
    and return a data frame.
    """

    # Get the column names
    struct_id_col = ROSETTA_DF_COLS["struct_id"]
    name_col = ROSETTA_DF_COLS["name"]
    state_col = ROSETTA_DF_COLS["state"]
    scf_name_col = ROSETTA_DF_COLS["scf_name"]
    b_steps_col = ROSETTA_DF_COLS["b_steps"]
    sc_type_col = ROSETTA_DF_COLS["sc_type"]
    sc_value_col = ROSETTA_DF_COLS["sc_value"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    
    # Open the connection
    db3_out = os.path.abspath(db3_out)
    connection = sqlite3.connect(db3_out)
    connection.row_factory = sqlite3.Row
    
    # Create the cursor
    cursor = connection.cursor()
    
    # Selection string to get the total number of batches
    sel_str = f'SELECT max(batch_id) from batches'
    n_batches = cursor.execute(sel_str).fetchone()[0]
    
    # Selection string for the structure scores
    sel = \
        "batches.name, structure_scores.struct_id," \
        "score_types.score_type_name, structure_scores.score_value," \
        "score_function_method_options.score_function_name"
    
    # Selection string for the inner join on batches
    injoin_batches = \
        "batches.batch_id=structure_scores.batch_id"
    
    # Selection string for the inner join on the scoring function
    injoin_scfunc = \
        "score_function_method_options.batch_id=batches.batch_id"
    
    # Selection string for the inner join on the score type
    injoin_sctype = \
        "score_types.batch_id=structure_scores.batch_id " \
        "AND score_types.score_type_id=structure_scores.score_type_id"
    
    # Assemble the query string
    query = \
        f"SELECT {sel} from structure_scores\n" \
        f"INNER JOIN batches ON {injoin_batches}\n" \
        f"INNER JOIN score_function_method_options ON {injoin_scfunc}\n" \
        f"INNER JOIN score_types ON {injoin_sctype}"
    
    # Try to read the query into a data frame
    try:
        df = pd.read_sql_query(query, connection)
    
    # If something went wrong, raise an error
    except Exception as e:
        errstr = f"Could not query the .db3 file ({db3_out}): {e}"
        raise IOError(errstr)
    
    # Function to renumber the structure IDs
    get_new_id = \
        lambda x: traj_stride * (1 + (int(x - 1) // n_batches))
    
    # Function to rename the structures
    get_new_name = \
        lambda x: x.replace("_dbreport", "") \
        if x.endswith('_dbreport') else x
    
    # Renumber the structure IDs
    df[b_steps_col] = df[struct_id_col].apply(get_new_id)
    
    # Rename the structures
    df[state_col] = df[name_col].apply(get_new_name)
    
    # Convert the dataframe into a pivot table
    df = \
        df.pivot_table(index = [state_col, b_steps_col, scf_name_col],
                       columns = sc_type_col,
                       values = sc_value_col).reset_index()
    
    # Add a column for the structure number
    df[struct_num_col] = struct_num
    
    # Add the complete score function name to the
    # score function name column
    df[scf_name_col] = scf_name
    
    # Remove unnecessary column names
    df.columns.names = [None]
    
    # Close the connection
    connection.close()

    # Return the data frame
    return df


def aggregate_data_flexddg(df,
                           list_contributions):
    """Aggregate data for flexddg protocols.
    """

    # Get the column names
    state_col = ROSETTA_DF_COLS["state"]
    scf_name_col = ROSETTA_DF_COLS["scf_name"]
    b_steps_col = ROSETTA_DF_COLS["b_steps"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    n_struct_col = ROSETTA_DF_COLS["n_struct"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]

    # Get the score function name
    scf_name = df[scf_name_col].unique()[0]
    
    # Columns storing energy scores
    score_cols = list_contributions + [tot_score_col]
    
    # Select only the rows corresponding to the maximum number
    # of backrub steps performed
    b_steps_cond = (df[b_steps_col] == df[b_steps_col].max())
    
    # Generate selections on binding and mutation states
    ub_wt = df.loc[(df[state_col] == "unbound_wt") & b_steps_cond]
    ub_mut = df.loc[(df[state_col] == "unbound_mut") & b_steps_cond]
    b_wt = df.loc[(df[state_col] == "bound_wt") & b_steps_cond]
    b_mut = df.loc[(df[state_col] == "bound_mut") & b_steps_cond]
    
    # Get the scores (reset the index to make it uniform to be able
    # to perform the subtraction between mutant ΔGs and
    # wild-type ΔGs later)
    ub_wt_scores = ub_wt[score_cols].reset_index(drop = True)
    ub_mut_scores = ub_mut[score_cols].reset_index(drop = True)
    b_wt_scores = b_wt[score_cols].reset_index(drop = True)
    b_mut_scores = b_mut[score_cols].reset_index(drop = True)
    
    # Get the list of structures from one of the data frames
    struct_nums = list(ub_wt[struct_num_col])
    
    # Create a data frame with the wild-type ΔG scores
    dg_wt = b_wt_scores - ub_wt_scores
    
    # Create a data frame with the mutant ΔG scores
    dg_mut = b_mut_scores - ub_mut_scores
    
    # Create a copy of the dataframe containing the mutant ΔG scores
    ddg = dg_mut.copy()
    
    # Subtract the wild-type ΔG scores from the mutant ΔG scores
    # to obtain the ΔΔG scores
    ddg[score_cols] = dg_mut[score_cols] - dg_wt[score_cols]
    
    # Add a column to all data frames with the state
    dg_wt[state_col], dg_mut[state_col], ddg[state_col] = "wt", "mut", "ddg"
    
    # Add a column to all data frames with the score function name
    dg_wt[scf_name_col], dg_mut[scf_name_col], ddg[scf_name_col] = \
        scf_name, scf_name, scf_name
    
    # Add a column to all data frames with the structure numbers
    dg_wt[struct_num_col], dg_mut[struct_num_col], ddg[struct_num_col] = \
        struct_nums, struct_nums, struct_nums

    # Return the data frames
    return (dg_wt, dg_mut, ddg)


def generate_output_dataframes(dg_wt,
                               dg_mut,
                               ddg,
                               mutation,
                               mut_label,
                               pos_label,
                               rescale,
                               list_contributions,
                               conv_fact,
                               family):
    
    # Get the columns names
    mutation_col = ROSETTA_DF_COLS["mutation"]
    mut_label_col = ROSETTA_DF_COLS["mut_label"]
    pos_label_col = ROSETTA_DF_COLS["pos_label"]
    state_col = ROSETTA_DF_COLS["state"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]
    energy_unit_col = ROSETTA_DF_COLS["energy_unit"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    scf_name_col = ROSETTA_DF_COLS["scf_name"]

    # Columns storing energy scores
    score_cols = [tot_score_col] + list_contributions

    # Columns of the structures dataframe
    struct_df_cols = \
        [mutation_col, mut_label_col, pos_label_col,
         struct_num_col, state_col, energy_unit_col,
         scf_name_col, *score_cols]

    # Columns of the aggregate dataframe
    aggr_df_cols = \
        [mutation_col, mut_label_col, pos_label_col, state_col,
         energy_unit_col, scf_name_col, *score_cols]


    #------------------------ Structures data ------------------------#


    # Concatenate the three dataframes
    struct_df = pd.concat([dg_wt, dg_mut, ddg])

    # Add the columns with the mutation name and the mutation labels
    struct_df.insert(0, mutation_col, mutation)
    struct_df.insert(1, mut_label_col, mut_label)
    struct_df.insert(2, pos_label_col, pos_label)


    #------------------------- Aggregate data ------------------------#


    # Get the column you want to group by
    group_by = \
        [mutation_col, mut_label_col, pos_label_col,
         state_col, scf_name_col]
    

    #
    # In some cases, 'struct_num_col' is parsed as string instead of integer.
    # So, this line was added to force the conversion of values to integer.
    #
    # Added by @alexandrefassio to fix issue #62.
    #
    struct_df[struct_num_col] = struct_df[struct_num_col].astype('int64')


    # If the protocol is the updated version of the cartddg protocol
    if family == "cartddg2020":

        # Sort the data frame so that the row corresponding to the
        # structure pair with the lowest energy score for the
        # mutant is the first one
        low_ddg_score_first = \
            struct_df.loc[struct_df[state_col] == "mut"].sort_values(\
                by = tot_score_col,
                axis = 0,
                ascending = True)

        # Get the number of the structure pair with lowest energy
        # score for the mutant
        struct_num_low = \
            low_ddg_score_first.iloc[0][struct_num_col]

        # The aggregate data frame will only contain data about
        # the structure pair having the lowest energy score for the
        # mutant
        aggr_df = \
            struct_df.loc[struct_df[struct_num_col] == struct_num_low]

        # Drop the column containing the structure number information
        aggr_df = aggr_df.drop(struct_num_col, axis = 1)

    # Otherwise
    elif family in ("cartddg", "flexddg"):
        
        # The aggregate data frame will contain data about the
        # average scores
        aggr_df = struct_df.groupby(group_by).mean().reset_index()


    #----------------------------- Rescale ---------------------------#


    # Default energy units are Rosetta Energy Units
    energy_unit = "REUs"
    
    # If rescaling has been requested
    if rescale:
        
        # Rescale the ΔΔG scores in REUs to kcal/mol
        aggr_df[score_cols] = aggr_df[score_cols] * 1.0/conv_fact
        struct_df[score_cols] = struct_df[score_cols] * 1.0/conv_fact
        
        # Change the energy unit
        energy_unit = "kcal/mol"


    #------------------------ Add units column -----------------------#


    # Add the energy units column to both dataframes
    aggr_df[energy_unit_col] = energy_unit
    struct_df[energy_unit_col] = energy_unit


    #-------------------------- Sort columns -------------------------#


    # Sort columns in both dataframes
    aggr_df = aggr_df[aggr_df_cols]
    struct_df = struct_df[struct_df_cols]

    # Return the structures data and the aggregate data
    return (aggr_df, struct_df)


def write_mutatex_df(dfs,
                     mutatex_file,
                     index,
                     family):

    # Get the name of the column containing the total scores
    # and the state
    tot_score_col = ROSETTA_DF_COLS["tot_score"]
    state_col = ROSETTA_DF_COLS["state"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]

    # Create an empty dictionary to store data that will be part
    # of the final data frame
    final_df_dict = {}

    # For each mutation, data frame
    for mutr, df in dfs.items():

        # If the protocol belongs to the new cartddg2020 family
        if family == "cartddg2020":

            # Sort the data frame so that the row corresponding to the
            # structure pair with the lowest energy score for the
            # mutant is the first one
            low_ddg_score_first = \
                df.loc[df[state_col] == "mut"].sort_values(\
                    by = tot_score_col,
                    axis = 0,
                    ascending = True)

            # Get the number of the structure pair with lowest energy
            # score for the mutant
            struct_num_low = \
                low_ddg_score_first.iloc[0][struct_num_col]

            # The ΔΔG score will be the one associated with the 
            # structure pair having the lowest energy score for
            # the mutant
            avg_val = \
                df.loc[(df[struct_num_col] == struct_num_low) & \
                       (df[state_col] == "ddg"), tot_score_col].item()

            # There will be only one value, which will also
            # correspond to the minimum and maximum value
            # (the standard deviation will be 0.0)
            std_val = 0.0
            min_val = avg_val
            max_val = avg_val

            # Set the furst line that will be printed to
            # the output file
            first_line = "# low\tstd\tmin\tmax\n"
        
        elif family in ("cartddg", "flexddg"):

            # Get only the ΔΔG scores
            ddg_scores = df[df[state_col] == "ddg"][tot_score_col]
            
            # Get the average, standard deviation, minimum and
            # maximum value
            avg_val = ddg_scores.mean()
            std_val = ddg_scores.std()
            min_val = ddg_scores.min()
            max_val = ddg_scores.max()

            # Set the furst line that will be printed to
            # the output file
            first_line = "# avg\tstd\tmin\tmax\n"

        # Append the average, standard deviation, minimum and
        # maximum value to the dictionary
        final_df_dict[mutr] = [avg_val, std_val, min_val, max_val]

    # Generate the final data frame, sorting the values according
    # to the mutated residues
    final_df = \
        pd.DataFrame.from_dict(final_df_dict,
                               orient = "index").reindex(index)

    # Open the output file
    out = open(mutatex_file, "a")

    # Write out the first line
    out.write(first_line)

    # Save the data frame to the output file
    final_df.to_csv(out,
                    index = None,
                    header = False,
                    sep = " ",
                    float_format = "%.5f")

    # Close the output file
    out.close()