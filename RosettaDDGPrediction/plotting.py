#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    plotting.py
#
#    Utility functions for plotting.
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
import operator
import re
# Third-party packages
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# RosettaDDGProtocols
from .defaults import (
    CHAIN,
    CHAIN_SEP,
    COMP_SEP,
    MUTR,
    NUMR,
    ROSETTA_DF_COLS, 
    WTR
)
from .util import ( 
     get_items,
     recursive_traverse
)



############################# DATA LOADING ############################



def load_aggregated_data(in_file,
                         saturation = False):
    """Load aggregated data from a CSV file.
    """

    # Columns containing data of interest in the aggregated
    # dataframe
    mutation_col = ROSETTA_DF_COLS["mutation"]

    # Load the dataframe
    df = pd.read_csv(in_file)

    # If the aggregated data are from a saturation mutagenesis
    # scan
    if saturation:
        
        # Create new columns 
        new_cols = \
            pd.DataFrame(\
                df[mutation_col].str.split(COMP_SEP, 3).tolist(),
                columns = [CHAIN, WTR, NUMR, MUTR])
        
        # Add the new columns
        df = pd.concat([df, new_cols], axis = 1)

    # Return the dataframe
    return df



############################# AESTHETICS ##############################



def generate_ticks_positions(values,
                             config):
    """Generate the positions that the ticks
    will have on a plot axis/colorbar/etc.
    """
    
    # Get the configurations
    int_type = config.get("type")
    rtn = config.get("round_to_nearest")
    top = config.get("top")
    bottom = config.get("bottom")
    steps = config.get("steps")
    spacing = config.get("spacing")
    ciz = config.get("center_in_zero")
    
    # If no rounding has been specified and the
    # interval is continuous
    if rtn is None and int_type == "continuous":
        
        # Default to rounding to the nearest 0.5
        rtn = 0.5

    # If the maximum of the ticks interval has not
    # been specified
    if top is None:
        
        if int_type == "discrete":
            # Default top value is the maximum
            # of the values provided
            top = int(max(values))
        elif int_type == "continuous":
            # Default top value is the rounded up 
            # maximum of the values
            top = np.ceil(max(values)*(1/rtn))/(1/rtn)

    # If the minimum of the ticks interval has not
    # been specified
    if bottom is None:
        
        if int_type == "discrete":
            # Default bottom value is the minimum
            # of the values provided
            bottom = int(min(values))
        elif int_type == "continuous":
            # Default bottom value is the rounded down 
            # minimum of the values
            bottom = np.floor(min(values)*(1/rtn))/(1/rtn)

    # If the number of steps the interval should have
    # has not been specified
    if steps is None:
        
        if int_type == "discrete":
            # Default number of steps is the 
            # number of values provided
            steps = 1
        elif int_type == "continuous":
            # Default is 5 steps
            steps = 5

    # If the interval spacing has not been specified
    if spacing is None:
        
        if int_type == "discrete":
            # Default spacing is the one between two steps,
            # rounded up
            spacing = \
                int(np.ceil(np.linspace(bottom,
                                        top,
                                        steps,
                                        retstep = True)[1]))
        elif int_type == "continuous":
            # Default spacing is the one between two steps
            spacing = np.linspace(bottom,
                                  top,
                                  steps,
                                  retstep = True)[1]

    # If the two extremes of the interval coincide
    if top == bottom:
        
        # Return only one value
        return np.array([bottom])

    # If the interval needs to be centered in zero
    if ciz:
        
        # Get the highest absolute value
        absval = \
            np.ceil(top) if top > bottom else np.floor(bottom)
        
        # Top and bottom will be opposite numbers with
        # absolute value equal to absval
        top, bottom = absval, -absval
        
        # Return an evenly spaced interval
        # between top and bottom values
        return np.linspace(bottom, top, steps)

    # Return the ticks interval
    return np.arange(bottom, top + spacing, spacing)


def generate_heatmap_annotations(df,
                                 config):
    """Generate the annotations to be plotted on 
    a heatmap (each cell is annotated with the
    corresponding value).
    """

    # If the configuration is empty
    if config == dict():
        # Return a tuple filled with None values
        return (None, None)

    # Get the configuration for the style of the annotations
    # and for the number of decimals to be kept in the
    # annotations
    annot = config.get("annot")
    n_decimals = config.get("ndecimals", 2)
    
    # If no annotation is requested, leave the dictionary
    # for the annotation properties empty
    annot_kws = {}

    # If annotations are requested
    if config.get("annot"):
        
        # Create a function to set the annotations 
        # to the desired precision
        annot_func = lambda x : np.around(x, n_decimals) 
        
        # Vectorize the function
        annot_transform = np.vectorize(annot_func)
        
        # Create annotations for all cells of the heatmap
        annot = annot_transform(df.values)
        
        # Get the style of the annotations
        annot_kws = config["style"]

    # Return annotations and style
    return (annot, annot_kws)


def generate_barplot_annotations(ax,
                                 config,
                                 sum_pos,
                                 total,
                                 yticks):
    """Generate annotations to be plotted over the
    bars of a bar plot.
    """

    # Get whether the annotations are requested or not
    annot = config.get("annot")
    
    # Get the number of decimals to be kept in the annotations
    n_decimals = config.get("ndecimals")
    
    # Get the style of the annotations
    style = config.get("style")

    # Initialize an empty list to store annotations
    annots = []

    # If the annotation is requested   
    if annot:
        
        # Get the spacing between two ticks on the y-axis
        spacing = yticks[1] - yticks[0]
        
        # For each bar
        for patch, p_sum, score in zip(ax.patches, sum_pos, total):
            
            # Get the position of the annotation on the x-axis
            # (centered on the bar)
            x = patch.get_x() + patch.get_width()/2.0
            
            # Get the position of the annotation on the y-axis
            # (slightly above the bar)
            y = p_sum + spacing/5
            
            # Round the annotation to the desired precision
            s = round(score, n_decimals)
            
            # Add the annotation over the corresponding bar
            ax.text(x, y, s, **style)


def generate_axhline(ax,
                     config,
                     df):
    """Generate a horizontal line passing through 0.0.
    """
    
    # Set the intercept of the line
    y = 0
    
    # Set the length of the line
    xmax = (len(df) - 0.25) / len(df)
    
    # Plot the line
    plt.axhline(y = y,
                xmax = xmax,
                **config)


def generate_colorbar(mappable,
                      ticks,
                      config):
    """Generate a colorbar associated to a mappable
    (for example, a heatmap).
    """
 
    # Plot the colorbar
    cbar = plt.colorbar(mappable, **config["colorbar"])

    # Get the colorbar orientation
    orient = config["colorbar"].get("orientation") if \
             config["colorbar"].get("orientation") is not None \
             else "vertical"

    # If there is an axis label (horizontal orientation)
    if config["label"].get("xlabel"):
        # Set the colorbar label  
        cbar.ax.set_xlabel(**config.get("label"))
    
    # If there is an axis label (vertical orientation)
    elif config["label"].get("ylabel"):
        # Set the colorbar label     
        cbar.ax.set_ylabel(**config.get("label"))
    
    # Set the colorbar ticks and ticks labels
    # setting ticks on cbar.ax raises a UserWarning,
    # but setting tick labels does not
    cbar.set_ticks(ticks)

    # Format the ticklabels (can behave weirdly if the position
    # of a tick is represented by a number which has no precise
    # binary representation)
    ticklabels = [f"{float(np.round(t,2)):g}" for t in ticks]
    
    # If the orientation of the colorbar is horizontal
    if orient == "horizontal":
        # Set the x-axis ticks
        cbar.ax.set_xticklabels(ticklabels,
                                **config.get("ticklabels"))

    # If the orientation of the colorbar is vertical
    elif orient == "vertical":
        # Set the y-axis ticks
        cbar.ax.set_yticklabels(ticklabels,
                                **config.get("ticklabels"))

    # Return the colorbar
    return cbar


def generate_legend(ax,
                    config):
    """Generate a legend for the current plot.
    """

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Draw the legend
    ax.legend(handles = handles,
              labels = labels,
              bbox_transform = plt.gcf().transFigure,
              **config)

    # Retutn the ax the legend has been plotted on
    return ax


def generate_mask_nancells(ax,
                           cells,
                           config):
    """Generate a mask to mark differently cells
    in a heatmap containing NaN values.
    """

    # For each NaN cell
    for y,x in cells:
        
        # Add a rectangulat patch over the cell
        ax.add_patch(mpatches.Rectangle(xy = (x,y), **config))

    # Return the ax the patches have been plotted on
    return ax


def set_axis(ax,
             axis,
             config,
             ticks = None,
             ticklabels = None):
    """Set up the x- or y-axis.
    """
    
    # If no tick positions were provided
    if ticks is None:
        
        if axis == "x":
            # Default to the tick locations already present
            ticks = plt.xticks()[0]
        
        elif axis == "y":
            # Default to the tick locations already present
            ticks = plt.yticks()[0]
    
    # If no tick labels were provided
    if ticklabels is None:
        
        # Default to the string representations
        # of the ticks' locations
        ticklabels = [str(t) for t in ticks]

    # Set the tick labels' configuration
    ticklabels_config = {}
    if config.get("ticklabels"):
        ticklabels_config = config["ticklabels"]
    
    # If it is the x-axis
    if axis == "x":    
        
        # If there is an axis label
        if config.get("label"):
            
            # Set the axis label
            ax.set_xlabel(**config["label"])        
        
        # Set the ticks
        ax.set_xticks(ticks = ticks)        
        
        # Set the tick labels
        ax.set_xticklabels(labels = ticklabels,
                           **ticklabels_config)
        
        # If tick positions were provided
        if ticks != []:      
            
            # Set the axis boundaries
            ax.spines["bottom"].set_bounds(ticks[0],
                                           ticks[-1])
    
    # If it is the y-axis
    elif axis == "y":        
        
        # If there is an axis label
        if config.get("label"):
            
            # Set the axis label
            ax.set_ylabel(**config["label"])        
        
        # Set the ticks
        ax.set_yticks(ticks = ticks)        
        
        # Set the tick labels
        ax.set_yticklabels(labels = ticklabels,
                           **ticklabels_config)
        
        # If tick positions were provided
        if ticks != []:       
            
            # Set the axis boundaries
            ax.spines["left"].set_bounds(ticks[0],
                                         ticks[-1])

    # If a configuration for the tick parameters was provided
    if config.get("tick_params"):
        
        # Apply the configuration to the ticks
        ax.tick_params(axis = axis,
                       **config["tick_params"])

    # Return the axis
    return ax



################################ PLOT #################################



def plot_total_heatmap(df,
                       config,
                       out_file,
                       out_config,
                       saturation = False):
    """Plot either:

    - a one-row heatmap where all mutations are displayed on the
      x-axis, each one displayed as a cell color-coded according to
      the total ΔΔG score corresponding to that mutation;

    - a 2D heatmap containing the total ΔΔG scores for positions on
      which a saturation mutagenesis scan was performed.
    """
    
    # Clear the figure
    plt.clf()

    # Close the current figure window
    plt.close()


    #----------------------- Data preprocessing ----------------------#


    # Get the names of the columns of the aggregate data frame
    # containing data of interest
    mutation_col = ROSETTA_DF_COLS["mutation"]
    pos_label_col = ROSETTA_DF_COLS["pos_label"]
    mut_label_col = ROSETTA_DF_COLS["mut_label"]
    state_col = ROSETTA_DF_COLS["state"] 
    tot_score_col = ROSETTA_DF_COLS["tot_score"]
    
    # If the data are from a saturation mutagenesis scan
    if saturation:
        
        # Take only the total ΔΔG values
        final_df = df[df[state_col] == "ddg"][[pos_label_col,
                                               MUTR, tot_score_col]]
        
        # Create a data frame where rows are positions (identified
        # by all mutation elements that are not the mutated residue,
        # such as chain, residue number and wild-type residue) and
        # columns are mutated residues
        final_df = final_df.pivot(index = pos_label_col,
                                  columns = MUTR,
                                  values = tot_score_col).transpose()
        
        # Y-axis tick labels will be the row names
        y_ticklabels = final_df.index.values.tolist()
        
        # Set y-axis ticks to None so that ticks are
        # automatically placed
        y_ticks = None
        
        # The x-axis tick labels will be the residue positions
        positions = df[pos_label_col].unique()
        
        # Reorder the data frame columns according to the positions
        final_df = final_df[positions]

    # If the data are not from a saturation mutagenesis scan
    else:
        
        # Keep only the label column and total score
        final_df = \
            df[df[state_col] == "ddg"][[mut_label_col, tot_score_col]]
        
        # Set the mutations as index, drop the column containing
        # them and transpose the dataframe so that the mutations
        # are on the x-axis
        final_df = final_df.set_index(mut_label_col).transpose()
        
        # The y-axis should have neither ticks nor tick labels,
        # since the heatmap will have only one row
        y_ticks, y_ticklabels = [], []

    # X-axis tick labels will be the column names
    x_ticklabels = final_df.columns.values.tolist()
    
    # Flatten the array so that we are dealing only with
    # a list of values
    values = final_df.values.flatten()
    
    # Drop NaN values, values shown on the y-axis will be 
    # the total ΔΔG values
    y_values = values[~np.isnan(values)]
    
    # Get the cells where the value is NaN (mutations for
    # which the ΔΔG is not available, i.e. if you have run
    # an incomplete saturation mutagenesis on some positions)
    nan_cells = np.argwhere(np.isnan(final_df.values))


    #------------------------- Configuration -------------------------#
    
    
    # Get the configuration for the heatmap, colorbar, NaN cells,
    # x- and y-axis
    h_config, c_config, nan_config, x_config, y_config = \
        get_items(config, ("heatmap", "colorbar", "nancells", 
                  "xaxis", "yaxis"), {})

    # Get the configuration for the heatmap grid and for the
    # annotations
    h_config_heat, h_config_annot = \
        get_items(h_config, ("heatmap", "annot"), {})

    # Get the configuration for the interval to be represented
    # on the colorbar
    c_config_int = \
        get_items(c_config, ("interval",), {})[0]


    #----------------------------- Plot ------------------------------#


    # Get the colorbar ticks positions
    c_ticks = generate_ticks_positions(values = y_values,
                                       config = c_config_int)
    
    # Get maximum and minimum value from the interval
    vmin, vmax = c_ticks[0], c_ticks[-1]

    # Generate the heatmap annotations
    annots = generate_heatmap_annotations(df = final_df,
                                          config = h_config_annot)

    # Generate the heatmap
    ax = sns.heatmap(data = final_df,
                     cbar = False,
                     annot = annots[0],
                     annot_kws = annots[1],
                     vmin = vmin,
                     vmax = vmax,
                     center = (vmax+vmin)/2,
                     **h_config_heat)
    
    # Add a mask to the NaN cells
    generate_mask_nancells(ax = ax,
                           cells = nan_cells,
                           config = nan_config)

    # Add the colorbar to the heatmap
    generate_colorbar(mappable = ax.get_children()[0],
                      ticks = c_ticks, \
                      config = c_config)

    # Set the x-axis
    set_axis(ax = ax,
             axis = "x",
             ticklabels = x_ticklabels,
             config = x_config)

    # Set the y-axis
    set_axis(ax = ax,
             axis = "y",
             ticks = y_ticks,
             ticklabels = y_ticklabels,
             config = y_config)

    # For top and right spine of the plot
    for spine in ["top", "right"]:
        
        # Hide it
        ax.spines[spine].set_visible(False)

    # Save the plot to the output file
    plt.savefig(out_file, **out_config)



def plot_dg_swarmplot(df,
                      config,
                      out_file,
                      out_config):
    """Plot a swarmplot showing, for each mutation, the ΔG score of
    each wild-type and mutant structure present in the ensemble of
    structures generated by the protocol.
    """
    
    # Clear the figure
    plt.clf()

    # Close the current figure window
    plt.close()


    #----------------------- Data preprocessing ----------------------#

    
    # Get the names of the columns of the data frame
    # containing data of interest
    mut_label_col = ROSETTA_DF_COLS["mut_label"]
    state_col = ROSETTA_DF_COLS["state"]
    struct_num_col = ROSETTA_DF_COLS["struct_num"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]

    # Take both wild type and mutant ΔG values, but
    # not ΔΔG values
    df = df.loc[df[state_col].isin(["wt", "mut"])]
    
    # Create a new data frame containing only the 
    # columns of interest
    df = \
        df[[mut_label_col, tot_score_col,
            state_col, struct_num_col]]

    # Get the x-axis tick labels
    x_ticklabels = df[mut_label_col].unique()
    
    # Get the x-axis tick positions
    x_ticks = range(len(x_ticklabels))
    
    # Numeric values of interest will be the ΔG values
    y_values = df[tot_score_col].values.flatten()
    
    # Set x, y and hue columns that will be used
    # to generate the swarmplot
    x, y, hue = mut_label_col, tot_score_col, state_col


    #------------------------- Configuration -------------------------#


    # Get the configuration for the swarmplot, the legend, the x- and
    # the y-axis
    s_config, l_config, x_config, y_config = \
        get_items(config, ("swarmplot", "legend", "xaxis", "yaxis"),
                  {})

    # Get the configuration for the interval to be represented
    # on the x-axis
    x_config_int = x_config.get("interval", {})

    # Get the configuration for the interval to be represented
    # on the y-axis
    y_config_int = y_config.get("interval", {})


    #----------------------------- Plot ------------------------------#


    # Get the positions of the ticks on the y-axis
    y_ticks = generate_ticks_positions(values = y_values,
                                       config = y_config_int)

    # Generate the swarmplot
    ax = sns.swarmplot(data = df,
                       x = x,
                       y = y,
                       hue = hue,
                       **s_config)

    # Add the legend
    generate_legend(ax = ax,
                    config = l_config)

    # Set the x-axis
    set_axis(ax = ax,
             axis = "x",
             config = x_config,
             ticks = x_ticks,
             ticklabels = x_ticklabels)

    # Set the y-axis
    set_axis(ax = ax,
             axis = "y",
             config = y_config,
             ticks = y_ticks)

    # For top and right spine of the plot
    for spine in ["top", "right"]:
        
        # Hide it
        ax.spines[spine].set_visible(False)

    # Save the plot to the output file
    plt.savefig(out_file, **out_config)



def plot_contributions_barplot(df,
                               config,
                               contributions,
                               out_file,
                               out_config):
    """Plot a bar plot with stacked bars representing the different
    energy contributions making up the total ΔΔG scores. 
    Positive contributions are stacked upon the y positive
    semiaxis while negative contributions are stacked upon
    the y negative semiaxis.
    """

    # Clear the figure
    plt.clf()

    # Close the current figure window
    plt.close()


    #----------------------- Data preprocessing ----------------------#

    
    # Get the names of the columns of the data frame
    # containing data of interest
    mut_label_col = ROSETTA_DF_COLS["mut_label"]
    pos_label_col = ROSETTA_DF_COLS["pos_label"]
    state_col = ROSETTA_DF_COLS["state"]
    tot_score_col = ROSETTA_DF_COLS["tot_score"]

    # Get the total ΔΔG values
    df_total = df[df[state_col] == "ddg"][tot_score_col]

    # Get the values for all energy contributions
    df_contr = df[df[state_col] == "ddg"][contributions]
    
    # Retrieve the sum of all positive contributions and of all 
    # negative contributions per each case (= row). Will be used to 
    # position the labels over the stacked bars and to set the y 
    # axis limits.
    sum_pos = []
    sum_neg = []
    for case_name, data in df_contr.iterrows():
        sum_pos.append(data[data>0].sum())
        sum_neg.append(data[data<0].sum())

    # Values represented on the y-axis should span the entire
    # range of values, from the lowest negative to the highest
    # positive
    y_values = sum_pos + sum_neg


    #------------------------- Configuration -------------------------#


    # Get the configuration for the bar plot, the legend, the line
    # marking y coordinate 0.0 and the x- and the y-axis
    b_config, l_config, axh_config, x_config, y_config = \
        get_items(config,
                  ("barplot", "legend", "axhline", "xaxis", "yaxis"),
                  {})

    # Get the configuration for the bars and for the annotations
    # displayed over the bars
    b_config_bars, b_config_annot = \
        get_items(b_config, ("bars", "annot"), {})

    # Get the configuration for the interval that will be
    # represented on the y-axis
    y_config_int = y_config.get("interval", {})


    #----------------------------- Plot ------------------------------#

    
    # Open the multi-page PDF
    with PdfPages(out_file) as pdf:

        # Get the unique positions
        positions = df[pos_label_col].unique().tolist()
        
        # Set the number of pages of the PDF equal to the
        # number of positions scanned
        num_pages = len(positions)

        # For each page
        for page in range(num_pages):

            # Create a new figure
            plt.figure()

            # Get a sub-data frame of the original one containing
            # data only for the current position of interest
            sub_df = df.loc[df[pos_label_col] == positions[page]]

            # The labels of the ticks on the x-axis will be the names
            # of the mutations
            x_ticklabels = sub_df[mut_label_col].unique()

            # Get the positions of the ticks on the x-axis
            x_ticks = range(len(x_ticklabels))

            # Get the energy contributions for the current position
            sub_df_contr = \
                sub_df[sub_df[state_col] == "ddg"][contributions]
            
            # Get the total scores for the current position
            sub_df_total = \
                sub_df[sub_df[state_col] == "ddg"][tot_score_col]
            
            # Get the sum of total positive contributions for the
            # mutations of the current position
            sub_sum_pos = []
            for case_name, data in sub_df_contr.iterrows():
                sub_sum_pos.append(data[data>0].sum())

            # Draw a stacked bar plot in which positive contributions
            # are stacked on the positive semiaxis and negative 
            # contributions are stacked on the negative semiaxis
            ax = sub_df_contr.plot(kind = "bar", **b_config_bars)

            # Get the positions of the ticks on the y-axis
            y_ticks = generate_ticks_positions(values = y_values,
                                               config = y_config_int)

            # Generate text annotations with the total ΔΔG scores
            # that will appear on the top of the bars
            generate_barplot_annotations(ax = ax,
                                         config = b_config_annot,
                                         sum_pos = sub_sum_pos,
                                         total = sub_df_total,
                                         yticks = y_ticks)

            # Generate an horizontal line that passes through
            # y coordinate 0.0
            generate_axhline(ax = ax,
                             config = axh_config,
                             df = sub_df_total)

            # Add the legend
            generate_legend(ax = ax,
                            config = l_config)

            # Set the x-axis
            set_axis(ax = ax,
                     axis = "x",
                     config = x_config,
                     ticks = x_ticks,
                     ticklabels = x_ticklabels)

            # Set the y-axis
            set_axis(ax = ax,
                     axis = "y",
                     config = y_config,
                     ticks = y_ticks)

            # For top and right spine of the plot
            for spine in ["top", "right"]:
                
                # Hide it
                ax.spines[spine].set_visible(False)

            # Save the figure to the PDF page
            pdf.savefig(**out_config)
            
            # Clear the figure
            plt.clf()
            
            # Close the current figure window
            plt.close()
