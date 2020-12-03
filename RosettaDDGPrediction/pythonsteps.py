#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pythonsteps.py
#
#    Protocol steps performed by Python.
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
import bisect
import random
import statistics
# RosettaDDGProtocols
from . import util


def select_structure(select, \
                     infile, \
                     infiletype):
    """Select a structure (= return its path) based on data stored
    in different Rosetta output files. The structres themselves are
    assumed to live in the same directory as the output file to be
    parsed.
    """

    #-------------------------- File parsing -------------------------#

    
    # if getting data about the structure from a scorefile
    if infiletype == "scorefile_text":
        # open and parse the file
        data = util.parse_scorefile_text(scorefile = infile)
    else:
        # no parsing for other file types has been implemented so far
        raise NotImplementedError

    
    #---------------------- Structure selection ----------------------#

    
    # select a random structure
    if select == "random":
        struct, sdata = random.choice(data)
    # select the first structure
    elif select == "first":
        struct, sdata = data[0]
    # select the last structure
    elif select == "last":
        struct, sdata = data[-1]
    # select the structure with the lowest associated value
    elif select == "lowest":
        struct, sdata = \
            sorted(data, key = lambda x: x[1], reverse = False)[0]
    # select the structure with the highest associated value    
    elif select == "highest":
        struct, sdata = \
            sorted(data, key = lambda x: x[1], reverse = True)[0]
    # select the structure whose value is closest to the mean value
    elif select == "closest_to_mean":
        # sort the data in ascending order
        sorteddata = \
            sorted(data, key = lambda x: x[1], reverse = False)
        structs, ssdata = zip(*sorteddata)
        # compute the mean
        avg = statistics.mean(ssdata)
        # get where the mean would be in the sorted list of data
        posavg = bisect.bisect_right(ssdata, avg)
        # assess if the value closest to the mean is the one on
        # the left or on the right of where the mean would be in
        # the sorted list of data
        lowercloser = \
             abs(ssdata[posavg-1]-avg) < abs(ssdata[posavg+1]-avg)
        posclosest = posavg-1 if lowercloser else posavg+1
        struct, sdata = structs[posclosest], ssdata[posclosest]
    # get the structure with median value
    elif select == "median":
        # sort the data in ascending order
        sorteddata = \
            sorted(data, key = lambda x: x[1], reverse = False)
        structs, ssdata = zip(*sorteddata)
        # compute the median
        median = statistics.median(ssdata)
        # get where the median would be in the sorted list of data
        posmedian = bisect.bisect_right(ssdata, median)
        # the structure with value equal to the median is the one
        # immediately preceding where the median would be in the
        # list of sorted data (if two values are equal, bisect
        # right returns an insertion point that is after any
        # existing instance of the value in the list)
        posclosest = posmedian-1
        struct, sdata = structs[posclosest], ssdata[posclosest]
    # raise an error if an unrecognized rule was provided
    else:
        raise ValueError(f"No rule with name {select} exists.")
    # return the structure
    return struct