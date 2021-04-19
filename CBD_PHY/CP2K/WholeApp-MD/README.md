# CBD_PHY

* PBD ID: 4O0P (adapted)
* 167,922 atoms
* 68 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short MD simulation of a QM/MM system using CP2K. 
5 MD steps are performed with a time step of 1 fs. There are PBE and PBE0 
directories containing the inputs for the respective functionals.


``CBD_PHY-cp2k.inp`` - The CP2K input file. Contains set up parameters for the MD run 
and QM parameters. The DZVP-MOLOPT-GTH basis set and the PBE/PBE0 XC functional are used.

``CBD_PHY-cp2k-wfn.inp`` - The CP2K input file for generating the initial SCF
wave function.

``CBD_PHY.prmtop`` - Amber forcefield for MM atoms. The Amber03 forcefield and
the TIP3P water model are used.

``CBD_PHY.pdb`` - Atomic input coordinates.




### To Run: 

    mpi_exec -n $NPROCS cp2k.psmp -i CBD_PHY-cp2k.inp -o output.log

On 24 processes (1 thread) the an MD steps takes ~190 seconds on ARCHER
