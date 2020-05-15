# GFP_ScaleQM - 20

* PBD ID: 1GFL
* 28,264 atoms
* 20 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short MD simulation of a QM/MM system using CP2K. 
5 MD steps are performed with a time step of 1 fs.

``qmmm-1.inp`` - The CP2K input file. Contains set up parameters for the MD run 
and QM parameters. The DZVP-GTH-BLYP basis set and the BLYP XC functional are used.

``gfp_new1.prmtop`` - Amber forcefield for MM atoms. The Amber03 forcefield and
the TIP3P water model are used.

``NPT-1.restart`` - CP2K restart file for input of equilibriated coordinates and velocities.



### To Run: 

    mpi_exec -n $NPROCS cp2k.psmp -i qmmm-1.inp -o output.log

On 24 processes (1 thread) an MD step takes X seconds on ARCHER