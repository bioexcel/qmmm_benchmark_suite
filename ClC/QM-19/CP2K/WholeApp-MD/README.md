# ClC-19

* PBD ID: 1KPK
* 150,925 atoms
* 19 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short MD simulation of a QM/MM system using CP2K. 
5 MD steps are performed with a time step of 1 fs.



``ClC-19-cp2k.inp`` - The CP2K input file. Contains set up parameters for the MD run 
and QM parameters. The DZVP-MOLOPT-GTH basis set and the BLYP XC functional are used.

``ClC.prmtop`` - Amber forcefield for MM atoms. The Amber12 forcefield and
the TIP3P water model are used.

``ClC.pdb`` - Atomic input coordinates.




### To Run: 

    mpi_exec -n $NPROCS cp2k.psmp -i ClC-19-cp2k.inp -o output.log

On 24 processes (1 thread) the benchmark takes ~57 seconds on ARCHER