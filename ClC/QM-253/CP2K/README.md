# ClC-253

* PBD ID: 1KPK
* 150,925 atoms
* 253 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short QM/MM simulation using CP2K. 
5 MD steps are performed with a time step of 1 fs. The following XC functional 
set ups are included:

* BLYP - using DVZP-MOLOPT-GTH
* B3LYP - using EMSL: 6-31Gxx
* PBE0 - using EMSL: 6-31Gxx




``ClC-253-cp2k.inp`` - The CP2K MD input file. Contains set up parameters for the MD run 
and QM parameters. 

``ClC-253-cp2k-wfn.inp`` - The CP2K input file to generate the inital SCF wavefunctions. 
Does a single energy point calculation using RUN_TYPE ENERGY.

``ClC.prmtop`` - Amber forcefield for MM atoms. The Amber12 forcefield and
the TIP3P water model are used.

``ClC.pdb`` - Atomic input coordinates.




### To Run: 

    mpi_exec -n $NPROCS cp2k.psmp -i ClC-253-cp2k.inp -o output.log

