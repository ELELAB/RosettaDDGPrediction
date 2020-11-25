#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

import logging

import dask
from distributed.utils import DequeHandler

# to address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger():
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # new level
    NEWLEVEL = logging.INFO
    # get the logger
    logger = logging.getLogger("distributed.worker")
    # define the handlers to keep
    htokeep = [h for h in logger.handlers if type(h).__name__ == \
               DequeHandler.__name__]
    # remove all the handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    # add the handlers to keep
    for h in htokeep:
        # set the new level
        h.setLevel(NEWLEVEL)
        # add the handler to the logger
        logger.addHandler(h)
    # reset the logger level to the new level
    logger.setLevel(NEWLEVEL)
    # return the new logger
    return logger