# MQAE

* N-(6-methoxyquinolyl) acetoethyl ester in solution
* 16,396 atoms
* 34 QM atoms

## QM/MM CP2K Whole Application benchmark

### About

This benchmark consists of a short MD simulation of a QM/MM system using CP2K. 
5 MD steps are performed with a time step of 1 fs.

``mqae-cp2k.inp`` - The CP2K input file. Contains set up parameters for the MD run 
and QM parameters. The DZVP-MOLOPT-GTH basis set and the BLYP XC functional are used.

``mqae.prmtop`` - Amber forcefield for MM atoms. The Amber12 forcefield and
the SPCE water model are used.

``mqae.pdb`` - Input coordinates file.


### To Run: 

    mpi_exec -n $NPROCS cp2k.psmp -i mqae-cp2k.inp -o output.log

On 24 processes (1 thread) an MD step takes ~26 seconds on ARCHER