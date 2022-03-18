# MQAE

* N-(6-methoxyquinolyl) acetoethyl ester in solution
* 16,396 atoms
* 34 QM atoms

## QM/MM GMX_CP2K benchmark

### About

This benchmark consists of a short QM/MM simulation using the Gromacs/CP2K
interface. 
5 MD steps are performed with a time step of 1 fs. The following XC functional 
set ups are included:

* BLYP - using DVZP-GTH-BLYP


``mqae.top`` - The Gromacs topology file.

``mqae.ndx`` - The Gromacs index file.

``mqae.gro`` - The Gromacs coordinates and velocities file.

``mqae.mdp`` - The Gromacs MD parameter file.

``mqae.inp`` - The CP2K input file. Contains QM parameters. 

``mqae_cp2k.pdb`` - The pdb coordinates for CP2K. 


### To Run: 

    gmx grompp -f mqae.mdp -p mqae.top -c mqae.gro -n mqae.ndx -qmi mqae_cp2k.inp -o mqae.tpr
    gmx_mpi mdrun -s mqae.tpr
