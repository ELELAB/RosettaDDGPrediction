# run_ddg_protocols.py

**Rosetta version**: Rosetta 3.11

**Last updated**: 24/05/2020

## The run_ddg_protocols.py script

The script is a user-friendly Python wrapper around a series of Rosetta commands, whose combinations reproduce the procedures developed by Park et al. [^park2016] and Barlow et al. [^barlow2018] for the prediction of the ΔΔG of stability of a monomeric protein upon mutation and the prediction of the ΔΔG of binding of a protein complex upon mutation, respectively. 

Other Rosetta-based procedures for ΔΔG predictions will be incorporated if need arises.

These procedures are hereafter referred to as "protocols".

## Caveats

### Rosetta options and defaults

Rosetta defaults for command line options and the options themselves may (and some are likely to) change over time. 

The specific Rosetta version and the "last updated" timestamp specified above may be useful to address the inconsistencies that may arise between the content of the README and more recent Rosetta versions. 

We will try to keep the script up-to-date with respect to changes in Rosetta, adjusting it in case they affect the behaviour of our code. However, we advise every user not to take it as a "black box" but instead spending a bit of time understanding the workflow of the protocols and the rationale behind them, which is also useful to understand their applicability to specific cases of study and inherent limitations.

### Working with post-translational modifications and non-canonical residues

**Important: none of the protcols implemented so far has been benchmarked for post-translational modifications or non-canonical residues, therefore the reliability of the results must be taken with a grain of salt.**

Unfortunately, all protocols included in the script so far are not capable of handling mutations to some post-traslationally modified residues (such as phosphorylated ones) because of design choices in the Rosetta suite itself that considers them *variants* of the canonical amino acids instead of proper non-canonical residues. 

However, all protocols implemented so far are able to handle mutations to non-canonical residues, such as D-amino acids. See the "Options" section below for more details on how to work with them. 

For `cartddg` protocols, there seems to be a bug in Rosetta that causes a crash when either a mutation to a non-canonical residue is specified after a mutation to a canonical residue (when it is a simultaneous mutation, see the "Input files" section below) or when multiple simultaneous mutations to non-canonical residues are requested. We posted the issue on the RosettaCommons forum, and the resulting thread can be seen by going to [this page](https://www.rosettacommons.org/node/10846). 

## ΔΔG of stability

#### Overview

We implemented two protocols for the prediction of the ΔΔG of stability of a monomeric protein upon mutation, named `cartddg_ref2015` and `cartddg_talaris2014`, that are actually versions of the same protocol using different scoring functions. They differ in the Rosetta energy function used for the calculations.

Both protocols include two steps:

1. A preprocessing of the input PDB structure, that consists in a relaxation in Cartesian space of the input PDB structure, performed with the Rosetta `relax` application.

2. The actual ΔΔG prediction, including the selection of the best rotameric side-chain for the wild-type and the mutated residue and another relaxation in Cartesian space, performed with the Rosetta `cartesian_ddg` application. 

For a set of mutations to be performed on a single PDB structure, step 1 should only be run once in order to produce the relaxed structure. Then, this structure would serve as input for as man runs as needed to calculate predict the ΔΔG for the mutations of interest.

For this reason, our script separates the run of the first step from the run of the second. You can select whether you want to run the `relax` step or the `cartesian_ddg` step when launching the script.

## ΔΔG of binding

#### Overview

We incorporated two protocols for the prediction of the ΔΔG of binding of a protein complex upon mutation, named `flexddg_ref2015` and `flexddg_talaris2014`, both coming from the work of Barlow et al. [^barlow2018]. They differ in the Rosetta energy function used for the calculations.

These two protocols make use of a XML script implementing the procedure devised by Barlow et al. (`Flex_ddG.xml`).

Both protocols consist in only one step, therefore there is no need to specify a `step` option when running them (see the "Options" section below).

## Installation

The script does not require any installation. There is no need for the script to be placed inside the directory from where you launch it, since you can call it from another directory and specify via command line the directory where you want the output and log files to be stored (see "Options" below).

## Requirements

`python3.7` or above.

## Usage

### Command line

`run_ddg_protocols.py [-h] -f PDBFILE -p PROTOCOL [-d WORKDIR] -r ROSETTAPATH [-l LISTFILE] [-s STEP] [-n NPROC] [-lf LOGFILE] [-v] [-vv] [--rosettascript ROSETTASCRIPT] [--saturation] [--reslistfile RESLISTFILE] [--ncaalistfile NCAALISTFILE]`

### Options

| Option                | Meaning                                                      |
| --------------------- | ------------------------------------------------------------ |
| `-h`, `--help`        | Show the help message and exit.                              |
| `-f`, `--pdbfile`     | PDB file of the wild-type input structure.                   |
| `-p`, `--protocol`    | Name of the protocol to be run.                              |
| `-d`, `--workdir`     | Working directory. Default is the current working directory. |
| `-r`, `--rosettapath` | Path to the Rosetta installation directory.                  |
| `-l`, `--listfile`    | File containing the list of selected mutations or positions to perform saturation mutagenesis. |
| `-s`, `--step`        | Which step of the protocol to run, since some protocols consist of multiple steps that must be run separately. |
| `-n`, `--nproc`       | Number of processes to be start in parallel. Default is one process. |
| `-lf`, `--logfile`    | Name of the log file (default is `run_ddg_protocols.log`)    |
| `-v`                  | Verbosity level: verbose.                                    |
| `-vv`                 | Verbosity level: debug.                                      |
| `--rosettascript`     | XML Rosetta script (for protocols that require one).         |
| `--saturation`        | Perform saturation mutagenesis on selected positions.        |
| `--reslistfile`       | File containing the list of residue types (one-letter name) to include in the saturation mutagenesis. It is used only if `--saturation` is provided. |
| `--ncaalistfile`      | File containing the list of non-canonical residue types (full names) to be considered. It is needed if you want to specify mutations to non-canonical residues. |

### Input files

#### PDB file

There are no restrictions on the input PDB file, only the requirement that it must be Rosetta-compatible. That means for example that Rosetta can throw errors if the file does not respect some Rosetta-specific conventions. Please check the Rosetta documentation for further details.

#### (Mutations/positions) list file

The format of this file depends on two factors:

* the protocol that is run
* if the run is on selected mutations or if it is a saturation mutagenesis on selected positions

File formats for the protocols implemented so far are shown below.

Non-canonical residue types are specified as `X[###]` where `###` is replaced by the full name of the residue type in Rosetta (i.e. `DALA` for D-alanine).

##### `cartddg` protocols (`cartddg_ref2015` and `cartddg_talaris2014`)

WARNING: because of format requirements in the `cartesian_ddg` application, residues should be numbered as if they started from 1 and were numbered progressively in the PDB file, regardless of the real PDB numbering and distinct chains. For this reason, no chain specification is needed when selecting the residue.

###### Selected mutations

Newline-separated list of mutations. If multiple simultaneous mutations are requested, they should be on the same line, separated by a comma.

```
C151Y,S154N
S2X[DALA]
```

###### Saturation mutagenesis

Newline-separated list of positions.

```
C151
S2
```

##### `flexddg` protocols (`flexddg_ref2015` and `flexddg_talaris2014`)

WARNING: residue numbering and chain specification must correspond to those in the input PDB file. Furthermore, because of a requirement of the Flex ddG procedure itself for the calculation of the ΔΔG of binding, a chain ID for each mutation must be specified. This is the chain that will be "moved away" from its partner during ΔΔG calculation, since the final value of the ΔΔG of binding is obtained by scoring the chains in complex and then individually and then subtracting the scores of the individual chains from that of the complex. If your complex is made by more than two chains, you just have to specify which chain(s) make(s) up the interface you are interested in. For example, if you have a complex including chains A,B and C and you are interested in the ΔΔG of binding between chain A and B, you can either move away chains A and C or B and C according to what makes more sense given the topology of your complex. In this case, the chain specification would be "AC" or "BC".

In case one of your chains did not have an ID, the corresponding symbol would be an underscore ("_").

###### Selected mutations

Newline-separated list of mutations. If multiple simultaneous mutations are requested, they should be on the same line, separated by a comma. The "chain(s)-to-be-moved" specification is separated from the mutations by a whitespace.

```
A-C151Y,A-S154N A
A-S2X[DALA] A
```

###### Saturation mutagenesis

Newline-separated list of positions. The "chain(s)-to-be-moved" specification is separated from each position by a whitespace.

```
A-C151 A
A-S2 A
```

#### Residue list file

Newline-separated list of one-letter residue names needed to specify which residues to include in a saturation mutagenesis.

```
A
C
D
E
```

#### Non-canonical Amino Acids (NCAA) list file

Whitespace-separated list of full names of non-canonical residue types (see the "Working with post-translational modifications and non-canonical residues" section above).

```
DALA DPHE
```

# References

[^park2016]: Park, Hahnbeom, et al. "Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules." *Journal of chemical theory and computation* 12.12 (2016): 6201-6212.
[^barlow2018]: Barlow, Kyle A., et al. "Flex ddG: Rosetta Ensemble-Based Estimation of Changes in Protein–Protein Binding Affinity upon Mutation." *The Journal of Physical Chemistry B* 122.21 (2018): 5389-5399.