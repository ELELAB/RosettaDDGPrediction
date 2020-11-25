# RosettaDDGPrediction

**Rosetta version**: 3.12

**Last updated**: November 25th, 2020

## Overview

RosettaDDGPrediction is a Python package to run Rosetta-based protocols for the prediction of the ΔΔG of stability upon mutation of a monomeric protein or the ΔΔG of binding upon mutation of a protein complex and analyze the results.

## Caveats

### Rosetta options and defaults

Rosetta defaults for command line options and the options themselves may (and some are likely to) change over time.

The specific "Rosetta version" and the "Last updated" timestamp specified above can be useful to address inconsistencies that may arise between the content of this file and more recent Rosetta versions.

We will try to keep the protocols up-to-date with respect to changes in Rosetta, adjusting them in case they affect the expected behaviour of our code. However, we advise every user not to take RosettaDDGPrediction as a "black box" but instead spending a bit of time understanding the workflow of the protocols and the rationale behind them, which is also useful to understand their applicability to specific cases of study and their inherent limitations.

### Working with post-translational modifications and non-canonical residues

**Important: none of the protcols implemented so far has been benchmarked for post-translational modifications or non-canonical residues, therefore the reliability of the calculations on such residues must be taken with a grain of salt.**

Unfortunately, none of the protocols included in the package so far is capable of handling mutations to some post-traslationally modified residues (such as phosphorylated ones) because of design choices in the Rosetta suite itself, which considers them *variants* of the canonical amino acids instead of proper non-canonical residues.

However, all protocols implemented so far are able to handle mutations to non-canonical residues, such as D-amino acids.

For `cartddg` protocols, there seems to be a bug in Rosetta that causes a crash when either a mutation to a non-canonical residue is specified after a mutation to a canonical residue (when it is a simultaneous mutation, see the "Input files" section below) or when multiple simultaneous mutations to non-canonical residues are requested. We posted the issue on the RosettaCommons forum, and the resulting thread can be seen by going to [this page](https://www.rosettacommons.org/node/10846).

## Requirements

The user must have Python v3.7 or higher installed, together with the Rosetta modelling suite (which can be found [here](https://www.rosettacommons.org/software/license-and-download)).

Required Python dependencies, if not already present, will be installed along with RosettaDDGPrediction.

## Installation

1. (Optional) We recommend installing this package inside a virtual environment. Please see [here](https://docs.python.org/3/tutorial/venv.html) how to create a Python virtual environment. Upon successful creation, activate the virtual environment before moving to the second step.

2. To install the package, download and unzip this folder, enter the folder and run the following command:

   `python3.7 setup.py install`

Upon successful installation, you should have three executable (`rosetta_ddg_run`, `rosetta_ddg_aggregate` and `rosetta_ddg_plot`) available to perform the various steps of data collection and analysis.

## Usage

Please refer to the user guide for how to use RosettaDDGPrediction.