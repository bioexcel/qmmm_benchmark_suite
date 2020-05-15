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
acetoethyl ester in water.All 34 atoms of the ester are treated with QM whereas 
the remaining water atoms are treated with MM. This system is chosen to 
represent a system with a small QM and a small MM subsystem. The MQAE system is
taken from [ref], however the forcefield is rebuilt using the Amber12 
forcefields, where the ester forcefield is created using the general Amber 
forcefield and the waters are modelled using the SPCE water model. For the XC
functional the BLYP functional is selected in keeping with the simulations
performed in [ref]. An energy cut-off of 300 Ry for the plane waves was found to
be suitable by checking the energy convergence.


### ClC-19 and ClC-253

The ClC-19 and ClC-253 benchmarks consist of a chloride ion channel embedded in a
lipid bilayer, which is solvated in water. There are 150,925 atoms in total, 
representing a large MM subsystem. The ClC-19 and ClC-253 systems contain 19 
and 253 QM atoms respectively, and therefore these systems represent a small 
and large QM subsystem and a large MM subsystem. These systems are taken from 
[ref] but adapted slightly to reduce the number of QM atoms (by removing waters
treated with QM). The Amber12 forcefields are used, and the water is treated 
using the TIP3P model. For the QM atoms the GPW method is used with the
DZVP-MOLOPT-GTH basis set and the BLYP XC functional, and the corresponding 
pseudopotentials. An energy cut-off for the plane waves of 300 Ry was found
to be suitable.

### CBD_PHY

This system is an example of a fluorescent protein, namely a phytochrome dimer
solvated in water. There are 68 QM atoms in this system and 167,922 atoms in
total, and therefore it represents a fairly large QM subsystem combined with a
large MM subsystem. This system the PBE XC functional and the DZVP-MOLOPT-GTH 
basis set. The energy cut-off for the plane waves is set to 400 Ry. For the MM
part the Amber03 forcefields are used with the TIP3P water model. 

### QMScale

## Test cases

### CP2K Whole Application MD

Contains the input files to run an MD simulation in native CP2K.
For all systems the QM/MM coupling is described with the Gaussian expansion 
of the Electrostatic Potential (GEEP) method, and any bonds between QM and MM
atoms are treated using the generalized hybrid orbital method. The treatment of
the QM atoms may vary between systems as the basis sets and exchange correlation
functionals have been chosen depending on the system in question. 


### CP2K Kernel

This test case aims to benchmark the performance of the CP2K QM/MM kernel, 
qmmm_forces_with_gaussian_LG. This kernel is one of the routines involved in 
calculation of the Gaussian expansion of the electrostatic potential (GEEP) 
which is required in the calculation of the QM+MM contribution to the atomic
forces for periodic systems.

### Grm+CP2K interface



## References