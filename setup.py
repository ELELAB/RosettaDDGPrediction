#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    setup.py
#
#    RosettaDDGPrediction setup.
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
from setuptools import setup


# Setup options
name = "RosettaDDGPrediction"

url = "https://github.com/ELELAB/RosettaDDGPrediction"

author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"

version = "0.0.1"

description = \
    "Python wrapper of Rosetta-based protocols for ΔΔG calculation."

package_dir = {"RosettaDDGPrediction" : "RosettaDDGPrediction"}

packages = ["RosettaDDGPrediction"]

package_data = \
    {"RosettaDDGPrediction" : ["config_aggregate/*",
                               "config_plot/*",
                               "config_run/*",
                               "config_settings/*",
                               "files/*",
                               "RosettaScripts/*"]}

entry_points = \
    {"console_scripts" : \
        ["rosetta_ddg_run = RosettaDDGPrediction.rosetta_ddg_run:main",
         "rosetta_ddg_check_run = RosettaDDGPrediction.rosetta_ddg_check_run:main",
         "rosetta_ddg_aggregate = RosettaDDGPrediction.rosetta_ddg_aggregate:main",
         "rosetta_ddg_plot = RosettaDDGPrediction.rosetta_ddg_plot:main"],
    }

install_requires = ["biopython",
                    "dask",
                    "distributed",
                    "matplotlib",
                    "MDAnalysis",
                    "numpy",
                    "pandas",
                    "pyyaml",
                    "seaborn"]


# Run the setup
setup(name = name,
      url = url,
      author = author,
      version = version,
      description = description,
      include_package_data = True,
      package_data = package_data,
      package_dir = package_dir,
      packages = packages,
      entry_points = entry_points,
      install_requires = install_requires)
