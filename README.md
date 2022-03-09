# QM/MM Benchmark Suite

This repository contains QM/MM simulation benchmarks for a number of 
biomolecular systems designed to run using CP2K and the GROMACS/CP2K 
interface. 

The suite is organised with different QM/MM systems at the top level. 

For each biomolecular system the suite contains up to three different
benchmarks, corresponding to different execution scenarios:

1. CP2K QM/MM MD - `system-name/CP2K/WholeApp-MD`
2. CP2K QM/MM kernel - `system-name/CP2K/Kernel`
3. GROMACS+CP2K QM/MM MD - `system-name/GRM+CP2K`
 
In addition there is a `/tools` directory which contains the code required to run the CP2K kernel benchmark.

## Biomolecular systems

|Name      |Type                |QM atoms       |Total atoms |Functional      |Basis set       |MD run type  |Periodic?|
|----------|--------------------|---------------|------------|----------------|----------------|-------------|---------|
|MQAE      |solute-solvent      |34             |16,396      |BLYP	           |DZVP-MOLOPT-GTH |NVE          |Y        |
|ClC       |ion channel         |19, 253        |150,925     |BLYP	           |DZVP-MOLOPT-GTH |NVE	         |Y        |
|CBD_PHY   |phytochrome         |68             |167,922     |PBE             |DZVP-MOLOPT-GTH |NVE          |Y        |
|GFP_QM-77 |fluorescent protein |20, 32, 53, 77 |28,264      |BLYP            |DZVP-GTH-BLYP   |NVT          |N        |


### MQAE

The MQAE system is a solute-solvent system consisting of a N-(6-methoxyquinolyl)
acetoethyl ester in solution taken from [1]. All 34 atoms of the ester are 
treated with QM whereas the remaining water atoms are treated with MM. This 
system is chosen to represent a system with a small QM and a small MM subsystem.
The MM parameters are rebuilt using the Amber forcefields, where the parameters 
for the organic molecule are created using the General Amber Force Field (GAFF) 
[2] and the water molecules are modelled using the SPCE model [3]. We 
selected the BLYP functional as the XC functional in keeping with the 
simulations performed in [1]. An energy cut-off of 400 Ry for the plane waves
was found to be suitable.


### ClC

ClC consists of a (ClC-ec1) chloride ion channel
embedded in a lipid bilayer (PDB-ID: 1KPK), which is solvated in water. Two
variants are included for this system - ClC-19 and ClC-253 which differ only
in having respectively 19 and 253 atoms treated quantum mechanically, 
representing a small and large QM subsystem within a large MM subsystem
(150,925 atoms in total). These systems are taken from [1] but adapted 
slightly to reduce the number of QM atoms (by removing waters treated 
with QM). The QM regions are modelled using the GPW method with the 
DZVP-MOLOPT-GTH basis set and the BLYP XC functional and the corresponding
pseudopotentials. An energy cut-off for the plane waves of 300 Ry was found to
be suitable. The Amber14 forcefield [2] is used for the protein and lipid14 
forcefield is used for the lipid molecules [4], and water molecules are 
treated using the TIP3P model [5].

### CBD_PHY

This system contains a phytochrome dimer (PBD-ID: 4O0P) with a bound chromophore, 
solvated in water. There are 68 QM atoms in this system and 167,922 atoms in total, 
and therefore it represents a fairly large QM subsystem (68 QM atoms) combined 
within a large MM subsystem (167,922 atoms in total). The QM region is modelled 
using the GPW method with DZVP-MOLOPT-GTH basis set and PBE XC functional. An 
energy cut-off for the plane waves of 400 Ry was found to be suitable. For the MM 
part the Amber03 forcefield is used for the protein and water molecules are treated
using the TIP3P model. 

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

## Execution scenerios

### CP2K QM/MM MD

These benchmarks consist of a short MD simulation for each system treated with 
QM/MM executed using CP2K. For each system 5 steps are performed with a time step 
of 1 fs. The benchmarks can each be run with the standard release version of CP2K.

For all systems the QM/MM coupling is described with the Gaussian Expansion of
the Electrostatic Potential (GEEP) method [6][7], and any bonds between QM and MM atoms
are treated using the Generalized Hybrid Orbital (GHO) method. The treatment of
the QM atoms may vary between systems as the basis sets and exchange correlation
functionals have been chosen depending on the system in question.

For each system the CP2K input file (.inp), the initial input coordinates (.pdb)
and the Amber MM forcefield (.prmtop) are provided in the WholeApp-MD directory.
More details can be found in the README file for each system. 

### CP2K QM/MM kernel

These benchmark the performance of the CP2K QM/MM subroutine
qmmm_forces_with_gaussian_LG isolated from the CP2K codebase. This subroutine 
computes the contribution to the forces due to the long range part of the
QM/MM potential using the Gaussian Expansion of the Electrostatic Potential
(GEEP) with periodic boundary conditions. As described in previous
deliverable D1.2 including in [9], it and closely related routines are
large contributors to overall CP2K QM/MM simulation runtime and isolating
it in a kernel facilitates performance analysis and optimisation. As many 
of the QM/MM routines have a similar structure, optimisations made here 
can be tested before implementation into related routines. 

The version of qmmm_forces_with_gaussian_LG included in the benchmark
suite kernel code includes OpenMP threading optimisations developed by
BioExcel and reported in the CP2K section of this deliverable. 
 
The kernel benchmarks were created by running the corresponding whole
application benchmark using a modified version of CP2K that generates the 
data that feeds into the kernel on a single MPI rank. The data therefore 
depends on the rank it was generated by and the total number of processes
used. The kernel benchmarks run the subroutine in serial in a way that is
representative of its execution by a single MPI rank in the original 
parallel execution context. The CP2K code modifications and instructions
to regenerate the relevant data (or indeed generate new kernel benchmarks)
can be found in the /tools directory of the benchmark suite.
 
The data required to run each kernel benchmark is provided in the /data 
directory. The code to run the benchmark itself is provided in the /src 
directory, and the required Makefile is found in the top level of the
kernel benchmark. The /outputs directory contains the force output from
running the kernel benchmark and can be used to check the correctness of
any changes made to the kernel. 
 
When run, the kernel benchmark proceeds by reading in the generated data,
and then calls qmmm_forces_with_gaussian_LG to calculate the forces. The
forces are written out in order to check correctness, and the run time 
for the kernel subroutine itself is reported. 




### GROMACS+CP2K QM/MM MD

These benchmark the GROMACS/CP2K interface [8] for performing a QM/MM MD 
simulation. GROMACS is the main driver for the interface, with libcp2k called to 
calculate the QM/MM energies and forces each MD step. CP2K input parameters are
passed through the cp2k.inp file, which is used to generate a force environment
within CP2K. This is maintained throughout the simulation, which aids performance
by enabling the reuse of density function from the previous step to 
extrapolate the wavefunction of the next step (using the always stable 
predictor corrector method - ASPC). At each step, updated atomic positions 
are passed from GROMACS through calls to libcp2k and the atomic forces and the 
energy are then calculated and returned.

This benchmark requires the use of GROMACS built with the CP2K interface. To do
this, CP2K must first be built as a library. When building GROMACS, libcp2k.h 
needs to be included in the header search path and the path to libcp2k.a needs 
to be added to the library path. Details of how to do this can be found on the 
Gromacs website [9].

The interface is run in almost the same way as GROMACS, with the required files
being the topology file (.top), the configuration file (.gro), the MD
parameters file (.mdp) and the index (.ndx) file. The major differences
are in the .mdp file where values specific to QM/MM simulation are provided.
The QM atoms group may be specified here, with the atom numbers given in the
index file. QM/MM parameters can also be supplied to CP2K in this file. It is
also possible to directly provide a CP2K input file (cp2k.inp) containing all
parameters for the QM/MM calculations by setting QMMMInput=INPUT in the .mdp 
file.

For each system, all the required files are provided in the GRM+CP2K directory.
To ensure consistency between the CP2K standalone benchmark and its corresponding 
GROMACS/CP2K interface benchmark, the Amber forcefields used in the CP2K benchmark 
were converted into GROMACS format using ParmEd [10]. The benchmark is set up to 
perform 5 MD steps, with a time step of 1 fs.


## References

[1] Extreme Scalability of DFT-Based QM/MM MD Simulations Using MiMiC.
Viacheslav Bolnykh, Jógvan Magnus Haugaard Olsen, Simone Meloni, Martin P. Bircher, Emiliano Ippoliti, Paolo Carloni, and Ursula Rothlisberger.
Journal of Chemical Theory and Computation 2019 15 (10), 5601-5613,
https://doi.org/10.1021/acs.jctc.9b00424

[2] ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB.
James A. Maier, Carmenza Martinez, Koushik Kasavajhala, Lauren Wickstrom, Kevin E. Hauser, and Carlos Simmerling.
Journal of Chemical Theory and Computation 2015 11 (8), 3696-3713
https://doi.org/10.1021/acs.jctc.5b00255

[3] The missing term in effective pair potentials.
H. J. C. Berendsen, J. R. Grigera and T. P. Straatsma.
Journal Physical Chemistry 1987, 91, 24, 6269-6271.
https://doi.org/10.1021/j100308a038

[4] Lipid14: The Amber Lipid Force Field.
J. Dickson, Benjamin D. Madej, Åge A. Skjevik, Robin M. Betz, Knut Teigen, Ian R. Gould, and Ross C. Walker.
Journal of Chemical Theory and Computation 2014 10 (2), 865-879.
https://doi.org/10.1021/ct4010307

[5] Comparison of simple potential functions for simulating liquid water.
William L. Jorgensen, Jayaraman Chandrasekhar, and Jeffry D. Madura.
Journel Chemicak Physics 1983 79, 926.
https://doi.org/10.1063/1.445869

[6] An Efficient Real Space Multigrid QM/MM Electrostatic Coupling.
Laino, T. Mohamed, F. Laio, A. Parrinello.
M. J. Chem. Theory Comput. 2005, 1 1176-1184. 
https://doi.org/10.1021/ct050123f

[7] An Efficient Linear-Scaling Electrostatic Coupling for Treating Periodic Boundary Conditions in QM/MM Simulations.
Laino, T.; Mohamed, F. Laio, A. and Parrinello.
M. J. Chem. Theory Comput. 2006 2 (5), 1370-1378
https://doi.org/10.1021/ct6001169

[8] https://manual.gromacs.org/documentation/2022-beta1/reference-manual/special/qmmm.html

[9] https://manual.gromacs.org/documentation/2022-beta1/install-guide/index.html#building-with-cp2k-qm-mm-support

[10] https://github.com/ParmEd/ParmEd


