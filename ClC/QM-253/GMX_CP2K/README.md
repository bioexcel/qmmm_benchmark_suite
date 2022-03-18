# ClC-253

* PBD ID: 1KPK
* 150,925 atoms
* 253 QM atoms

## QM/MM GMX_CP2K benchmark

### About

This benchmark consists of a short QM/MM simulation using the Gromacs/CP2K
interface. 
5 MD steps are performed with a time step of 1 fs. The following XC functional 
set ups are included:

* BLYP - using DVZP-MOLOPT-GTH


``ClC.top`` - The Gromacs topology file.

``ClC.ndx`` - The Gromacs index file.

``ClC.gro`` - The Gromacs coordinates file.

``ClC.mdp`` - The Gromacs Md parameter file.

``ClC_cp2k.inp`` - The CP2K input file. Contains QM parameters. 

``ClC_cp2k.pdb`` - The pdb coordinates for CP2K. 


### To Run: 

    gmx grompp -f ClC.mdp -p ClC.top -c ClC.gro -n ClC.ndx -qmi ClC_cp2k.inp -o ClC.tpr
    gmx_mpi mdrun -s ClC.tpr
