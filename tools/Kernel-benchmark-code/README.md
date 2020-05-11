# CP2K QM/MM kernel benchmark

## Description

The CP2K QM/MM kernel benchmark provided here runs the most expensive subroutines identified during profiling of the whole-application biomolecular QM/MM benchmark specified in this repository under `../whole-application-benchmark`, namely `qmmm_forces_with_gaussian_LG` (or `qmmm_forces_with_gaussian_LR` for non-periodic systems). This is by far the most costly single subroutine in the execution of the whole-application benchmark, taking almost 30% of overall runtime when running at a reasonable parallel scale. It is similar in form to the first part of the subroutine that computes the corresponding forces (`qmmm_elec_with_gaussian_LG`), which is also the next most expensive subroutine during execution of the whole-application benchmark. Any optimisations found to speed up the kernel benchmark are therefore likely to also yield improvement in the QM/MM energies subroutine and hence taken together in almost half of the overall execution time when running at a reasonably efficient parallel scale.

The kernel subroutine computes the force contributions due to the long-range electrostatic interaction between QM and MM regions. It does this by treating the QM-MM coupling using the method of multigrids and Gaussian Expansion of the QM/MM Electrostatic Potential (GEEP) described in [1] and taking into account periodic boundary conditions as described in [2]. 

During regular execution of CP2K running the whole-application benchmark each MPI rank - one per core - executes the kernel subroutine independently without any communication between processes taking place within the call. Moreover,  execution times are almost identical between ranks. The code provided for the kernel benchmark therefore runs the kernel subroutine in serial in a way that is fully representative of a single MPI rank in the parallel execution context. Optimisations within the kernel subroutine should benefit performance of each rank and thereby directly overall runtime if changes are integrated back into the full application.


### References

[1] An Efficient Real Space Multigrid QM/MM Electrostatic Coupling. Laino, T.; Mohamed, F.; Laio, A.; Parrinello, M.  J. Chem. Theory Comput. 2005, 1, 1176-1184. https://doi.org/10.1021/ct050123f  

[2] An Efficient Linear-Scaling Electrostatic Coupling for Treating Periodic Boundary Conditions in QM/MM Simulations. Laino, T.; Mohamed, F. Laio, A. and Parrinello, M. J. Chem. Theory Comput. 2006 2 (5), 1370-1378 https://doi.org/10.1021/ct6001169  

## Requirements

The benchmark kernel is entirely self contained - the kernel code itself (`qmmm_gpw_forces.f90`), a  wrapper (`kernel_benchmark.f90`), some dependencies adapted from the CP2K codebase, and input files corresponding to the whole-application benchmark also contained in this repository are all included (more information in `./src/README.md` and `./data/README.md`). All that is required to build the kernel benchmark is a Fortran 2008-compatible compiler (gfortran is recommended, as it is for CP2K as a whole). If a Fortran 2008-compatible compiler is not available, you will have to modify a few file IO statements in `./src/kernel_benchmark.f90` to not make use of the 'newunit' feature. The code runs serially so no MPI library is required. The kernel benchmark has negligible memory requirements (~20MB memory high watermark). 

## How to build

- In the Makefile in this directory, specify a Fortran compiler that supports Fortran 2008 (in particular the `newunit` IO feature), and specify an optimisation level by setting FCFLAGS (e.g. to `-O3`)

- Type `make`

- The resulting executable kernel_benchmark will now appear in this directory

- `make clean` can be executed in this directory to clear all object and module files as well as the `kernel_benchmark` executable

## How to run

The executable `kernel_benchmark` should be run without any input parameters. It expects the `data` subdirectory to be located in the working directory where it executes. The executable prints out times obtained using calls to `system_clock()` to time the initialisation and the time spent running the kernel. 

## Checking results and expected runtime

A reference output file `./outputs/Forces.out` is provided for each of the kernel benchmarks to check against the `Forces.out` produced by running the benchmark. This can be compared using the `diff` command. The reference output was produced by compiling with gfortran version 8.2.0 with `-g` and default optimisation level (`-O0`) and running the benchmark on a single core of a quad-core Intel i7-3820QM@2.7GHz on a macOS laptop, for which the kernel call takes just over 3 minutes. No difference in `Forces.out` was encountered when running on the same platform with `-O3`, for which the kernel call takes around 60 seconds.

## Adapting the kernel benchmark

One might like to run the kernel benchmark code with different inputs. This is facilitated by inclusion of the file `./tools/qmmm_gpw_forces.F`, which was used to write to file the data structures needed to run the kernel benchmark and which are included here under `./data`. For more information see `./tools/README.md`