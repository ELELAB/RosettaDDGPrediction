#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "RosettaDDGPrediction"
url = "https://github.com/ELELAB/RosettaDDGPrediction"
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"
version = "0.0.1"
description = "Python wrapper of Rosetta-based protocols for ΔΔG calculation."
package_dir = {"RosettaDDGPrediction" : "RosettaDDGPrediction"}
packages = ["RosettaDDGPrediction"]
package_data = {'RosettaDDGPrediction' : ['config_aggregate/*',
                                          'config_plot/*',
                                          'config_run/*',
                                          'config_settings/*',
                                          'files/*',
                                          'RosettaScripts/*']}
entry_points = \
    {"console_scripts" : \
        ["rosetta_ddg_run = RosettaDDGPrediction.rosetta_ddg_run:main", \
         "rosetta_ddg_aggregate = RosettaDDGPrediction.rosetta_ddg_aggregate:main", \
         "rosetta_ddg_plot = RosettaDDGPrediction.rosetta_ddg_plot:main"], \
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
      package_data = package_data, \
      package_dir = package_dir, \
      packages = packages, \
      entry_points = entry_points, \
      install_requires = install_requires)

