# Developer notes

**Rosetta version**: Rosetta 3.11

**Last updated**: 24/05/2020

## General

* An arbitrary definition of "family" of protocols is implemented, where protocols distinguished only by minor differences (i.e. the scoring function used) are grouped together. In this case, a single `run` function for the whole group can be implemented to avoid code redundancy, and differences in options used can be resolved inside the function.
* An arbitrary definition of protocol "step" is also implemented, since some protocols include some conceptually different stages that must be run independently (i.e. a step to preprocess the input structure before doing ΔΔG calculations, where the preprocessed structure may be used to run other ΔΔG calculations later). In these cases, the steps are implemented in separate `run` functions. The implementation is such that multiple runs of the script are necessary to run all the steps.

## `run` functions

* All `run` functions for Rosetta protocols/steps should return a `subprocess.CompletedProcess` instance instaed of a simple return code, so that it is possible to have access to all the attributes of the instance if needed. They should also return the directory where the function was run. This allows for an easy trace back and logging of what happened during single runs. The return code and command line associated to the process can be retrieved via the `returncode` and `args` attributes of the `CompletedProcess` instance.
* In case of procedures that are designed to be run multiple times within a single run of the associated protocol (i.e. the Flex ddG  being applied a number of times to generate a number of structures), the corresponding `run` function should run only one instance at a time.
* All `run` functions should accept only one mutation (single or mutiple) at a time. A mutation can be defined in custom formats (tuples, strings, etc.), and the parsing of such formats should happen within the function. For procedures that must be run multiple times to generate a number of structures (i.e. Flex ddG), the definition of the mutation should include some reference to the specific instance (i.e. the structure) associated with that specific run.

## Input files

* The `listfile` is intended as a list of selected mutations or positions to perform saturation mutagenesis. It can be in any format according to the protocol, as far as a corresponding `get_mutlist` function is implemented. The format of the returned objects should match the one accepted and parsed by the corresponding `run` function.

## Parallelization

* The `map()` method is used for the `Pool` when multiprocessing so that the results are returned in the same order as the input, so that the parallelization can happen either mutation-wise (each mutation run in parallel), position-wise (for saturation mutagenesis, each position is run in parallel and mutations of that positions are run sequentially) or structure-wise (for procedures generating several structures with independent runs) according to the order of the inputs.

## Rosetta-related

#### Relax script

The relax script describing how the relaxation procedure should be run is`MonomerRelax2019.txt`, which is also the default, hard-coded one for the FastRelax runs within the `cartesian_ddg` application.

The choice of using the `legacy` relax script was previously made in order to be consisent with the one in use at the time of development and benchmarking of the `cartddg` protocols [^park2016], but since it is only possible to use it in `relax` but not in `cartesian_ddg` (because the selection is hard-coded in the latter) we decided to use the default Rosetta script also in `relax` for consistency between the two steps.

# References

[^park2016]: Park, Hahnbeom, et al. "Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules." *Journal of chemical theory and computation* 12.12 (2016): 6201-6212.

