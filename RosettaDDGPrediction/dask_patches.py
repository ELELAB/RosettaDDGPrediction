#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    dask_patches.py
#
#    Patches to fix unexpected Dask's behaviors.
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
import logging
# Third-party packages
import dask
from distributed.utils import DequeHandler



# To address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger():
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # New level
    NEW_LEVEL = logging.INFO
    
    # Get the logger
    logger = logging.getLogger("distributed.worker")
    
    # Define the handlers to keep
    h_to_keep = \
        [h for h in logger.handlers if type(h).__name__ == \
         DequeHandler.__name__]
    
    # Remove all the handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    
    # Add the handlers to keep
    for h in h_to_keep:
        
        # Set the new level
        h.setLevel(NEW_LEVEL)
        
        # Add the handler to the logger
        logger.addHandler(h)
    
    # Reset the logger level to the new level
    logger.setLevel(NEW_LEVEL)
    
    # Return the new logger
    return logger