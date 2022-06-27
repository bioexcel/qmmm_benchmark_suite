# ClC-19

* PBD ID: 1KPK
* 150,925 atoms
* 19 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short MD simulation of a QM/MM system using CP2K. 
5 MD steps are performed with a time step of 1 fs. The following XC functional
set ups are included:

* BLYP - using DVZP-MOLOPT-GTH
* B3LYP - using EMSL: 6-31Gxx, HFX_BASIS and ADMM with DZVP-MOLOPT-GTH
* PBE0 - using EMSL: 6-31Gxx


``ClC-19-cp2k.inp`` - The CP2K MD input file. Contains set up parameters for the MD run 
and QM parameters. 

``ClC-19-cp2k-wfn.inp`` - The CP2K input file to generate the inital SCF wavefunctions. Does
a single energy point calculation using RUN_TYPE ENERGY.

``ClC.prmtop`` - Amber forcefield for MM atoms. The Amber12 forcefield and
the TIP3P water model are used.

``ClC.pdb`` - Atomic input coordinates.




### To Run: 

Thw ``ClC.pdb`` and ``ClC.prmtop`` files ahould be in the same directory as the input file.

    mpi_exec -n $NPROCS cp2k.psmp -i ClC-19-cp2k.inp -o output.log
