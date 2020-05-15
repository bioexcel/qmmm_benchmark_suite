# QM/MM Benchmark Suite

This benchmark suite contains a number of QM/MM systems designed for benchmarking CP2K and the Grm+CP2K interface.

The benchmarks are organised 


## Systems

|Name      |QM atoms       |Total atoms |Functional  |Basis set       |MD run type  |Periodic?|
|----------|---------------|------------|------------|----------------|-------------|---------|
|MQAE      |34             |16,396      |BLYP	     |DZVP-MOLOPT-GTH |NVE          |Y        |
|ClC-19    |19             |150,925	    |BLYP	     |DZVP-MOLOPT-GTH |NVE	        |Y        |
|ClC-253   |253            |150,925	    |BLYP	     |DZVP-MOLOPT-GTH |NVE          |Y        |
|CBD_PHY   |68             |167,922	    |PBE         |DZVP-MOLOPT-GTH |NVE          |Y        |
|GFP_QM-77 |20, 32, 53, 77 |28,264      |BLYP        |DZVP-GTH-BLYP   |NVT          |N        |


### MQAE

The MQAE system is a solute-solvent system consisting of a N-(6-methoxyquinolyl)
acetoethyl ester in solution taken from [ref]. All 34 atoms of the ester are 
treated with QM whereas the remaining water atoms are treated with MM. This 
system is chosen to represent a system with a small QM and a small MM subsystem.
The MM parameters are rebuilt using the Amber forcefields, where the parameters 
for the organic molecule are created using the General Amber Force Field (GAFF) 
[REF] and the water molecules are modelled using the SPCE model [REF]. We 
selected the BLYP functional as the XC functional in keeping with the 
simulations performed in [ref]. An energy cut-off of 300 Ry for the plane waves
was found to be suitable.


### ClC-19 and ClC-253

The ClC-19 and ClC-253 systems consist of a (ClC-ec1) chloride ion channel
embedded in a lipid bilayer (PDB-ID: 1KPK), which is solvated in water. The
ClC-19 and ClC-253 systems contain 19 and 253 QM atoms respectively, and 
therefore these systems represent a small and large QM subsystem within a large
MM subsystem (150,925 atoms in total). These systems are taken from [ref] but 
adapted slightly to reduce the number of QM atoms (by removing waters treated 
with QM). The QM regions are modelled using the GPW method with the 
DZVP-MOLOPT-GTH basis set and the BLYP XC functional and the corresponding
pseudopotentials. An energy cut-off for the plane waves of 300 Ry was found to
be suitable. The Amber14 forcefield is used for the protein and lipid14 
forcefield is used for the lipid molecules [REF], and water molecules are 
treated using the TIP3P model [REF].

### CBD_PHY

This system is an example of contains a fluorescent protein, namely a 
phytochrome dimer solvated in water (adapted from PBD-ID: 4O0P). There are 68 QM
atoms in this system and 167,922 atoms in total, and therefore it represents a 
fairly large QM subsystem (68 QM atoms) combined within a large MM subsystem 
(167,922 atoms in total). The QM region is modelled using GPW method with 
DZVP-MOLOPT-GTH basis set and PBE XC functional. An energy cut-off for the 
plane waves of 400 Ry was found to be suitable. For the MM part the Amber03
forcefields is used for the protein and are used with the TIP3P model for the 
water molecules. 

### GFP_ScaleQM

The GFP_QM systems contains a green fluorescent protein (GFP) in solution
(PDB-ID: 1GFL). They represent a range of different QM subsystem sizes within a
small MM subsystem (28,264 atoms in total). For the QM atoms the DZVP-GTH-BLYP 
basis set is used along with the BLYP XC functional. An energy cut-off for the 
plane waves of 300 Ry was found to be suitable. For the MM part of the system, 
the Amber03 forcefield is used for the protein and TIP3P model for the water 
molecules. Unlike the other systems the QM subsystem is treated non-periodically
by using the Poisson solver for the electrostatics. This means that the 
non-periodic versions of the GEEP routines will be used.

## Test cases

### CP2K Whole Application MD

This test case consists of a short MD simulation for each system treated with 
QM/MM using CP2K. For each system 5 steps are performed with a time step of 1
fs. The benchmarks in this test case can each be run with the standard release
version of CP2K.

For all systems the QM/MM coupling is described with the Gaussian Expansion of
the Electrostatic Potential (GEEP) method, and any bonds between QM and MM atoms
are treated using the Generalized Hybrid Orbital (GHO) method. The treatment of
the QM atoms may vary between systems as the basis sets and exchange correlation
functionals have been chosen depending on the system in question.

For each system the CP2K input file (.inp), the initial input coordinates (.pdb)
and the Amber MM forcefield (.prmtop) are provided in the WholeApp-MD directory.
More details can be found in the README file for each system. 

### CP2K Kernel

This test case aims to benchmark the performance of the CP2K QM/MM kernel, 
qmmm_forces_with_gaussian_LG. This kernel is one of the routines involved in 
calculation of the Gaussian expansion of the electrostatic potential (GEEP) 
which is required in the calculation of the QM+MM contribution to the atomic
forces for periodic systems.

This benchmark was created by running the entire code in order to generate the
data that feeds into the kernel for a single process. This means that the data
will be dependent on the process it was generated on, and the total number of
processes used. This benchmark runs the kernel subroutine in serial in a way 
that is fully representative of a single MPI rank in the parallel execution 
context. The code and instructions to generate the relevant data can be found
in the /tools directory of the benchmark suite.

For each system, this required data is provided in the /data directory. 
The code to run the benchmark itself is provided in the /src directory,
and the required Makefile is found in the top level of the kernel benchmark. 
The /outputs directory contains the force output from running the kernel
benchmark and can be used to check the correctness of any changes made to the
kernel. 

When run the kernel benchmark proceeds by reading in the generated data, and 
then calling the kernel in question, qmmm_forces_with_gaussian_LG, to calculate
the forces. The forces are written out in order to check for correctness and the
run time for the kernel subroutine itself is reported using the system_clock().
This benchmark is designed to run using OpenMP and can be used to determine the 
speed up on a single process as a function of the number of threads used.


### Grm+CP2K interface

The test case benchmarks the Gromacs/CP2K interface for performing an MD 
simulation. Gromacs is the main driver for the interface with libcp2k called to 
calculate the QM/MM energies and forces each MD step. Input parameters are
passed through the cp2k.inp file, which is used to generate a force environment
within CP2K. This is maintained throughout the simulation and is used to improve 
the performance by using the density function from the previous step to 
extrapolate the wavefunction of the next step (using the always stable 
predictor corrector method - ASPC). At each step, updated atomic positions 
are passed from Gromacs through libcp2k and the atomic forces and the energy 
are then calculated and returned.

This benchmark requires the use of Gromacs built with the CP2K interface. To do
this, CP2K must first be built as a library. When building Gromacs, libcp2k.h 
needs to be included in the header search path and the path to libcp2k.a needs 
to be added to the library path.

The interface is run in almost the same way as Gromacs, with the required files
being the topology file (.top), the configuration file (.gro), the MD
paramemters file (grompp.mdp) and the index (.ndx) file. The major differences
are in the .mdp file where values specific to QM/MM simulation are provided.
The QM atoms group may be specified here, with the atom numbers given in the
index file. QM/MM parameters can also be supplied to CP2K in this file. It is
also possible to directly provide a CP2K input file (cp2k.inp) containing all
parameters for the QM/MM calculations by setting QMMMInput=INPUT in the .mdp 
file.

For each system, all the required files are provided in the GRM+CP2K directory.
To ensure consistency been the native CP2K and the interface the Amber 
forcefields used in CP2K benchmark have been converted into Gromacs format using
Parmed [ref]. The benchmark is set up to perform X MD steps, with a time step of
1 fs.


## References