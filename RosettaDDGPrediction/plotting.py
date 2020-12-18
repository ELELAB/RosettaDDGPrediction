#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    plotting.py
#
#    Utility functions for plotting.
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
import collections
import operator
import re
# third-party packages
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
# RosettaDDGProtocols
from .defaults import (
    ROSETTADFCOLS, 
    WTR, 
    NUMR, 
    MUTR, 
    CHAIN, 
    NOMUTR
)
from .util import (
     recursive_traverse, 
     get_items, 
     get_dirnames2mutations
)



############################ CONFIGURATION ############################



def get_config_plot_version_1(config):
    """Parse the configuration file for plotting,
    version 1.
    """
    
    # create a copy of the configuration
    newconfig = dict(config)
    # substitute the font properties definitions
    # with the corresponding FontProperties instances
    recursive_traverse(data = newconfig, \
                       actions = ["substitute_dict"], \
                       func = fm.FontProperties, \
                       keys = {"fontproperties"})
    # return the configuration
    return newconfig


def get_config_plot(configfile):
    """Get the plotting configuration."""
    
    # load the configuration from the file
    config = yaml.safe_load(open(configfile, "r"))

    # check the version of the configuration file
    if config["version"] == 1:
        # return the configuration written in version 1 format
        return get_config_plot_version_1(config = config)
    else:
        raise ValueError("Only version 1 configuration files " \
                         "are supported for now.")



############################# AESTHETICS ##############################



def generate_ticks_positions(values, config):
    """Generate the positions that the ticks
    will have on a plot axis/colorbar/etc.
    """
    
    # get the configurations
    inttype = config.get("type")
    rtn = config.get("round_to_nearest")
    top = config.get("top")
    bottom = config.get("bottom")
    steps = config.get("steps")
    spacing = config.get("spacing")
    ciz = config.get("center_in_zero")
    
    # if no rounding has been specified and the
    # interval is continuous
    if not rtn and inttype == "continuous":
        # default to rounding to the nearest 0.5
        rtn = 0.5

    # if the maximum of the ticks interval has not
    # been specified
    if not top:
        if inttype == "discrete":
            # default top value is the maximum
            # of the values provided
            top = int(max(values))
        elif inttype == "continuous":
            # default top value is the rounded up 
            # maximum of the values
            top = np.ceil(max(values)*(1/rtn))/(1/rtn)

    # if the minimum of the ticks interval has not
    # been specified
    if not bottom:
        if inttype == "discrete":
            # default bottom value is the minimum
            # of the values provided
            bottom = int(min(values))
        elif inttype == "continuous":
            # default bottom value is the rounded down 
            # minimum of the values
            bottom = np.floor(min(values)*(1/rtn))/(1/rtn)

    # if the number of steps the interval should have
    # has not been specified
    if not steps:
        if inttype == "discrete":
            # default number of steps is the 
            # number of values provided
            steps = 1
        elif inttype == "continuous":
            # default is 5 steps
            steps = 5

    # if the interval spacing has not been specified
    if not spacing:
        if inttype == "discrete":
            # default spacing is the one between two steps,
            # rounded up
            spacing = int(np.ceil(np.linspace(bottom, \
                                              top, \
                                              steps, \
                                              retstep = True)[1]))
        elif inttype == "continuous":
            # default spacing is the one between two steps
            spacing = np.linspace(bottom, \
                                  top, \
                                  steps, \
                                  retstep = True)[1]

    # if the two extremes of the interval coincide
    if top == bottom:
        # return only one value
        return np.array([bottom])

    # if the interval needs to be centered in zero
    if ciz:
        # get the highest absolute value
        absval = np.ceil(top) if top > bottom \
                 else np.floor(bottom)
        # top and bottom will be opposite numbers with
        # absolute value equal to absval
        top, bottom = absval, -absval
        # return an evenly spaced interval
        # between top and bottom values
        return np.linspace(bottom, top, steps)

    # return the ticks interval
    return np.arange(bottom, top + spacing, spacing)


def generate_heatmap_annotations(df, config):
    """Generate the annotations to be plotted on 
    a heatmap (each cell is annotated with the
    corresponding value).
    """

    # if the configuration is empty
    if config == dict():
        # return a tuple filled with None values
        return (None, None)

    # get the configuration for the style of the annotations
    # and for the number of decimals to be kept in the
    # annotations
    annot = config.get("annot")
    ndecimals = config.get("ndecimals", 2)
    # if no annotation is requested, leave the dictionary
    # for the annotation properties empty
    annotkws = {}

    # if annotations are requested
    if config.get("annot"):
        # create a function to set the annotations 
        # to the desired precision
        annotfunc = lambda x : np.around(x, ndecimals) 
        # vectorize the function
        annottransform = np.vectorize(annotfunc)
        # create annotations for all cells of the heatmap
        annot = annottransform(df.values)
        # get the style of the annotations
        annotkws = config["style"]

    # return annotations and style
    return (annot, annotkws)


def generate_barplot_annotations(ax, config, sumpos, total, yticks):
    """Generate annotations to be plotted over the
    bars of a bar plot.
    """

    # get whether the annotations are requested or not
    annot = config.get("annot")
    # get the number of decimals to be kept in the annotations
    ndecimals = config.get("ndecimals")
    # get the style of the annotations
    style = config.get("style")

    # if the annotation is requested
    annots = []
    if annot:
        # get the spacing between two ticks on the y-axis
        spacing = yticks[1] - yticks[0]
        # for each bar
        for patch, psum, score in zip(ax.patches, sumpos, total):
            # get the position of the annotation on the x-axis
            # (centered on the bar)
            x = patch.get_x() + patch.get_width()/2.0
            # get the position of the annotation on the y-axis
            # (slightly above the bar)
            y = psum + spacing/5
            # round the annotation to the desired precision
            s = round(score, ndecimals)
            # add the annotation over the corresponding bar
            ax.text(x, y, s, **style)


def generate_axhline(ax, config, df):
    """Generate a horizontal line passing through 0.0."""
    
    # set the intercept of the line
    y = 0
    # set the length of the line
    xmax = (len(df) - 0.25) / len(df)
    # plot the line
    plt.axhline(y = y, \
                xmax = xmax, \
                **config)


def generate_colorbar(mappable, \
                      ticks, \
                      config):
    """Generate a colorbar associated to a mappable
    (for example, a heatmap).
    """
 
    # plot the colorbar
    cbar = plt.colorbar(mappable, **config["colorbar"])

    # if there is an axis label
    if config["label"].get("ylabel"):
        # set the colorbar label     
        cbar.ax.set_ylabel(**config.get("label"))
    
    # set the colorbar ticks and ticks labels
    # setting ticks on cbar.ax raises a UserWarning, but setting
    # tick labels does not
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(ticks, **config.get("ticklabels"))

    # return the colorbar
    return cbar


def generate_legend(ax, config):
    """Generate a legend for the current plot."""

    # get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    # draw the legend
    ax.legend(handles = handles, \
              labels = labels, \
              bbox_transform = plt.gcf().transFigure, \
              **config)

    # retutn the ax the legend has been plotted on
    return ax


def generate_mask_nancells(ax, cells, config):
    """Generate a mask to mark differently cells
    in a heatmap containing NaN values.
    """

    # for each NaN cell
    for y,x in cells:
        # add a rectangulat patch over the cell
        ax.add_patch(mpatches.Rectangle(xy = (x,y), **config))

    # return the ax the patches have been plotted on
    return ax


def set_axis(ax, \
             axis, \
             config, \
             ticks = None, \
             ticklabels = None):
    """Set up the x- or y-axis."""
    
    if ticks is None:
        if axis == "x":
            # default to the tick locations already present
            ticks = plt.xticks()[0]
        elif axis == "y":
            # default to the tick locations already present
            ticks = plt.yticks()[0]
    
    if ticklabels is None:
        # default to the string representations
        # of the ticks' locations
        ticklabels = [str(t) for t in ticks]

    # set the tick labels
    ticklabelsconfig = {}
    if config.get("ticklabels"):
        ticklabelsconfig = config["ticklabels"]
    
    # if it is the x-axis
    if axis == "x":    
        # if there is an axis label
        if config.get("label"):
            # set the axis label
            ax.set_xlabel(**config["label"])        
        # set the ticks
        ax.set_xticks(ticks = ticks)        
        # set the tick labels
        ax.set_xticklabels(labels = ticklabels, \
                           **ticklabelsconfig)
        if ticks != []:      
            # set the axis boundaries
            ax.spines["bottom"].set_bounds(ticks[0], \
                                           ticks[-1])
    
    # if it is the y-axis
    elif axis == "y":        
        # if there is an axis label
        if config.get("label"):
            # set the axis label
            ax.set_ylabel(**config["label"])        
        # set the ticks
        ax.set_yticks(ticks = ticks)        
        # set the tick labels
        ax.set_yticklabels(labels = ticklabels, \
                           **ticklabelsconfig)
        if ticks != []:       
            # set the axis boundaries
            ax.spines["left"].set_bounds(ticks[0], \
                                         ticks[-1])

    # if a configuration for the tick parameters was provided
    if config.get("tick_params"):
        # apply the configuration to the ticks
        ax.tick_params(axis = axis, \
                       **config["tick_params"])

    # return the axis
    return ax



################################ PLOT #################################



def plot_total_heatmap(df, \
                       config, \
                       saturation = False, \
                       d2mfile = None):
    """Plot either:

    - a one-row heatmap where all mutations are displayed on the
      x-axis, each one displayed as a cell color-coded according to
      the total ΔΔG score corresponding to that mutation;

    - a 2D heatmap containing the total ΔΔG scores for positions on
      which a saturation mutagenesis scan was performed.
    """
    
    # clear the figure
    plt.clf()


    #----------------------- Data preprocessing ----------------------#


    # get the names of the columns of the dataframe
    # containing data of interest
    mutationcol = ROSETTADFCOLS["mutation"]
    statecol = ROSETTADFCOLS["state"] 
    totscorecol = ROSETTADFCOLS["totscore"]

    # take only the total ΔΔG values
    finaldf = df[df[statecol] == "ddg"][[mutationcol, totscorecol]]
    
    # if the data are from a saturation mutagenesis scan
    if saturation:
        # get a dataframe containing both mutation names
        # and all mutation elements
        d2m = get_dirnames2mutations(d2mfile)
        # merge the two dataframes on the mutation names
        finaldf = pd.merge(d2m, finaldf, on = mutationcol)
        # drop the mutation column
        finaldf = finaldf.drop(mutationcol, axis = 1)
        # create a dataframe where rows are positions (identified
        # by all mutation elements that are not the mutated residue,
        # such as chain, residue number and wild-type residue) and
        # columns are mutated residues
        finaldf = finaldf.pivot(index = NOMUTR, \
                                columns = MUTR, \
                                values = totscorecol).transpose()
        # y-axis tick labels will be the row names
        yticklabels = finaldf.index.values.tolist()
        # set y-axis ticks to None so that ticks are
        # automatically placed
        yticks = None

    # if the data are not from a saturation mutagenesis scan
    else:
        # set the mutations as index, drop the column containin
        # them and transpose the dataframe so that the mutations
        # are on the x-axis
        finaldf = finaldf.set_index(mutationcol).transpose()
        # the y-axis should have neither ticks nor tick labels,
        # since the heatmap will have only one row
        yticks, yticklabels = [], []

    # x-axis tick labels will be the column names
    xticklabels = finaldf.columns.values.tolist()
    # flatten the array so that we are dealing only with
    # a list of values
    values = finaldf.values.flatten()
    # drop NaN values, values shown on the y-axis will be 
    # the total ΔΔG values
    yvalues = values[~np.isnan(values)]
    # get the cells where the value is NaN (mutations for
    # which the ΔΔG is not available, i.e. if you have run
    # an incomplete saturation mutagenesis on some positions)
    nancells = np.argwhere(np.isnan(finaldf.values))


    #------------------------- Configuration -------------------------#
    
    
    # get the configuration for the heatmap, colorbar, NaN cells,
    # x- and y-axis
    hconfig, cconfig, nanconfig, xconfig, yconfig = \
        get_items(config, ("heatmap", "colorbar", "nancells", 
                  "xaxis", "yaxis"), {})

    # get the configuration for the heatmap grid and for the
    # annotations
    hconfigheat, hconfigannot = \
        get_items(hconfig, ("heatmap", "annot"), {})

    # get the configuration for the interval to be represented
    # on the colorbar
    cconfigint = \
        get_items(cconfig, ("interval",), {})[0]


    #----------------------------- Plot ------------------------------#


    # get the colorbar ticks positions
    cticks = generate_ticks_positions(values = yvalues, \
                                      config = cconfigint)
    
    # get maximum and minimum value from the interval
    vmin, vmax = cticks[0], cticks[-1]

    # generate the heatmap annotations
    annots = generate_heatmap_annotations(df = finaldf, \
                                          config = hconfigannot)

    # generate the heatmap
    ax = sns.heatmap(data = finaldf, \
                     cbar = False, \
                     annot = annots[0], \
                     annot_kws = annots[1], \
                     vmin = vmin, \
                     vmax = vmax, \
                     center = (vmax+vmin)/2, \
                     **hconfigheat)
    
    # add a mask to the NaN cells
    generate_mask_nancells(ax = ax, \
                           cells = nancells, \
                           config = nanconfig)

    # add the colorbar to the heatmap
    generate_colorbar(mappable = ax.get_children()[0], \
                      ticks = cticks, \
                      config = cconfig)

    # set the x-axis
    set_axis(ax = ax, \
             axis = "x", \
             ticklabels = xticklabels, \
             config = xconfig)

    # set the y-axis
    set_axis(ax = ax, \
             axis = "y", \
             ticks = yticks, \
             ticklabels = yticklabels, \
             config = yconfig)

    # return the plot axis
    return ax



def plot_dg_swarmplot(df, config):
    """Plot a swarmplot showing, for each mutation, the ΔG score of
    each wild-type and mutant structure present in the ensemble of
    structures generated by the protocol.
    """
    
    # clear the figure
    plt.clf()


    #----------------------- Data preprocessing ----------------------#

    
    # get the names of the columns of the dataframe
    # containing data of interest
    structnumcol = ROSETTADFCOLS["structnum"]
    totscorecol = ROSETTADFCOLS["totscore"]
    mutationcol = ROSETTADFCOLS["mutation"]
    statecol = ROSETTADFCOLS["state"] 

    # take both wild type and mutant ΔG values, but
    # not ΔΔG values
    df = df.loc[df[statecol].isin(["wt", "mut"])]
    # create a new dataframe containing only the 
    # columns of interest
    df = df[[mutationcol, totscorecol, \
             statecol, structnumcol]]

    # get the x-axis tick labels
    xticklabels = df[mutationcol].unique()
    # get the x-axis tick positions
    xticks = range(len(xticklabels))
    
    # numeric values of interest will be the ΔG values
    yvalues = df[totscorecol].values.flatten()
    # set x, y and hue columns that will be used
    # to generate the swarmplot
    x, y, hue = mutationcol, totscorecol, statecol


    #------------------------- Configuration -------------------------#


    # get the configuration for the swarmplot, the legend, the x- and
    # the y-axis
    sconfig, lconfig, xconfig, yconfig = \
        get_items(config, ("swarmplot", "legend", "xaxis", "yaxis"), \
                  {})

    # get the configuration for the interval to be represented
    # on the x-axis
    xconfigint = xconfig.get("interval", {})

    # get the configuration for the interval to be represented
    # on the y-axis
    yconfigint = yconfig.get("interval", {})


    #----------------------------- Plot ------------------------------#


    # get the positions of the ticks on the y-axis
    yticks = generate_ticks_positions(values = yvalues, \
                                      config = yconfigint)

    # generate the swarmplot
    ax = sns.swarmplot(data = df, \
                       x = x, \
                       y = y, \
                       hue = hue, \
                       **sconfig)

    # add the legend
    generate_legend(ax = ax, \
                    config = lconfig)

    # set the x-axis
    set_axis(ax = ax, \
             axis = "x", \
             config = xconfig, \
             ticks = xticks, \
             ticklabels = xticklabels)

    # set the y-axis
    set_axis(ax = ax, \
             axis = "y", \
             config = yconfig, \
             ticks = yticks)

    # return the plot axis
    return ax


def plot_contributions_barplot(df, contributions, config):
    """Plot a bar plot with stacked bars representing the different
    energy contributions making up the total ΔΔG scores. 
    Positive contributions are stacked upon the y positive
    semiaxis while negative contributions are stacked upon
    the y negative semiaxis.
    """

    # clear the figure
    plt.clf()


    #----------------------- Data preprocessing ----------------------#

    
    # get the names of the columns of the dataframe
    # containing data of interest
    mutationcol = ROSETTADFCOLS["mutation"]
    statecol = ROSETTADFCOLS["state"]
    totscorecol = ROSETTADFCOLS["totscore"]

    # get the total ΔΔG values
    dftotal = df[df[statecol] == "ddg"][totscorecol]
    # get the values for all energy contributions
    dfcontr = df[df[statecol] == "ddg"][contributions]
    
    # retrieve the sum of all positive contributions and of all 
    # negative contributions per each case (= row). Will be used to 
    # position the labels over the stacked bars and to set the y 
    # axis limits.
    sumpos = []
    sumneg = []
    for casename, data in dfcontr.iterrows():
        sumpos.append(data[data>0].sum())
        sumneg.append(data[data<0].sum())

    # values represented on the y-axis should span the entire
    # range of values, from the lowest negative to the highest
    # positive
    yvalues = sumpos + sumneg

    # the labels of the ticks on the x-axis will be the names
    # of the mutations
    xticklabels = df[mutationcol].unique()

    # get the positions of the ticks on the x-axis
    xticks = range(len(xticklabels))


    #------------------------- Configuration -------------------------#


    # get the configuration for the bar plot, the legend, the line
    # marking y coordinate 0.0 and the x- and the y-axis
    bconfig, lconfig, axhconfig, xconfig, yconfig = \
        get_items(config, ("barplot", "legend", "axhline", \
                  "xaxis", "yaxis"), {})

    # get the configuration for the bars and for the annotations
    # displayed over the bars
    bconfigbars, bconfigannot = \
        get_items(bconfig, ("bars", "annot"), {})

    # get the configuration for the interval that will be
    # represented on the y-axis
    yconfigint = yconfig.get("interval", {})


    #----------------------------- Plot ------------------------------#


    # draw a stacked bar plot in which positive contributions are
    # stacked on the positive semiaxis and negative contributions are
    # stacked on the negative semiaxis
    ax = dfcontr.plot(kind = "bar", **bconfigbars)

    # get the positions of the ticks on the y-axis
    yticks = generate_ticks_positions(values = yvalues, \
                                      config = yconfigint)

    # generate text annotations with the total ΔΔG scores that will
    # appear on the top of the bars
    generate_barplot_annotations(ax = ax, \
                                 config = bconfigannot, \
                                 sumpos = sumpos, \
                                 total = dftotal, \
                                 yticks = yticks)

    # generate an horizontal line that passes through
    # y coordinate 0.0
    generate_axhline(ax = ax, \
                     config = axhconfig, \
                     df = dftotal)

    # add the legend
    generate_legend(ax = ax, \
                    config = lconfig)

    # set the x-axis
    set_axis(ax = ax, \
             axis = "x", \
             config = xconfig, \
             ticks = xticks, \
             ticklabels = xticklabels)

    # set the y-axis
    set_axis(ax = ax, \
             axis = "y", \
             config = yconfig, \
             ticks = yticks)

    # return the plot axis
    return ax