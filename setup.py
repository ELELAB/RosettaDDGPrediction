#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "RosettaDDGProtocols"
url = "https://github.com/ELELAB/RosettaDDGProtocols"
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"
version = "0.0.1"
description = "Python wrapper of Rosetta-based protocols for ΔΔG calculation."
#package_data = {"RosettaDDGProtocols" : ["files/*", "config_run/*", "config_aggregate/*", "config_plot/*"]}
package_dir = {"RosettaDDGProtocols" : "RosettaDDGProtocols"}
packages = ["RosettaDDGProtocols"]
entry_points = \
    {"console_scripts" : \
        ["rosetta_ddg_run = RosettaDDGProtocols.rosetta_ddg_run:main", \
         "rosetta_ddg_aggregate = RosettaDDGProtocols.rosetta_ddg_aggregate:main", \
         "rosetta_ddg_plot = RosettaDDGProtocols.rosetta_ddg_plot:main"], \
    }
install_requires = ["biopython", \
                    "dask", \
                    "distributed", \
                    "matplotlib", \
                    "MDAnalysis", \
                    "numpy", \
                    "pandas", \
                    "pyyaml", \
                    "seaborn"]

setup(name = name, \
      url = url, \
      author = author, \
      version = version, \
      description = description, \
      include_package_data = True, \
      package_dir = package_dir, \
      packages = packages, \
      entry_points = entry_points, \
      install_requires = install_requires)