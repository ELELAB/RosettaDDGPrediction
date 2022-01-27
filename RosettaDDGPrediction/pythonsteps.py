#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pythonsteps.py
#
#    Protocol steps performed by Python.
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
import bisect
import random
import statistics
# RosettaDDGProtocols
from . import util



def select_structure(select,
                     in_file,
                     in_file_type):
    """Select a structure (= return its path) based on data stored
    in different Rosetta output files. The structres themselves are
    assumed to live in the same directory as the output file to be
    parsed.
    """


    #-------------------------- File parsing -------------------------#

    
    # If data about the structure are taken from a scorefile in text
    # format
    if in_file_type == "scorefile_text":
        
        # Open and parse the file
        data = util.parse_scorefile_text(scorefile = in_file)
    
    else:
        # No parsing for other file types has been implemented so far
        errstr = \
            f"You requested the parsing of a file of type " \
            f"'{in_file_type}' in order to select the relaxed " \
            f"structure, but only parsing of scorefiles in text " \
            f"format has been implemented so far."
        raise NotImplementedError(errstr)

    
    #---------------------- Structure selection ----------------------#

    
    # Select a random structure
    if select == "random":
        struct, s_data = random.choice(data)
    
    # Select the first structure
    elif select == "first":
        struct, s_data = data[0]
    
    # Select the last structure
    elif select == "last":
        struct, s_data = data[-1]
    
    # Select the structure with the lowest associated value
    elif select == "lowest":
        struct, s_data = \
            sorted(data, key = lambda x: x[1], reverse = False)[0]
    
    # Select the structure with the highest associated value    
    elif select == "highest":
        struct, s_data = \
            sorted(data, key = lambda x: x[1], reverse = True)[0]
    
    # Select the structure whose value is closest to the mean value
    elif select == "closest_to_mean":
        
        # Sort the data in ascending order
        sorted_data = \
            sorted(data, key = lambda x: x[1], reverse = False) 
        structs, ss_data = zip(*sorted_data)
        
        # Compute the mean
        avg = statistics.mean(ss_data)
        
        # Get where the mean would be in the sorted list of data
        pos_avg = bisect.bisect_right(ss_data, avg)
        
        # Assess if the value closest to the mean is the one on
        # the left or on the right of where the mean would be in
        # the sorted list of data
        lower_closer = \
             abs(ss_data[pos_avg-1]-avg) < abs(ss_data[pos_avg+1]-avg)
        pos_closest = pos_avg-1 if lower_closer else pos_avg+1
        struct, s_data = structs[pos_closest], ss_data[pos_closest]
    
    # Get the structure with median value
    elif select == "median":
        
        # Sort the data in ascending order
        sorted_data = \
            sorted(data, key = lambda x: x[1], reverse = False)
        structs, ss_data = zip(*sorted_data)
        
        # Compute the median
        median = statistics.median(ss_data)
        
        # Get where the median would be in the sorted list of data
        pos_median = bisect.bisect_right(ss_data, median)
        
        # The structure with value equal to the median is the one
        # immediately preceding where the median would be in the
        # list of sorted data (if two values are equal, bisect
        # right returns an insertion point that is after any
        # existing instance of the value in the list)
        pos_closest = pos_median-1
        struct, s_data = structs[pos_closest], ss_data[pos_closest]
    
    # Raise an error if an unrecognized rule was provided
    else:
        raise ValueError(f"No rule with name {select} exists.")
    
    # Return the structure
    return struct