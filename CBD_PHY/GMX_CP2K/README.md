# CBD-PHY

* PBD ID: 4O0P (adapted)
* 167,922 atoms
* 68 QM atoms

## QM/MM GMX_CP2K benchmark

### About

This benchmark consists of a short QM/MM simulation using the Gromacs/CP2K
interface. 
5 MD steps are performed with a time step of 1 fs. The following XC functional 
set ups are included:

* PBE - using DVZP-MOLOPT-GTH


``CBD_PHY.top`` - The Gromacs topology file.

``CBD_PHY.ndx`` - The Gromacs index file.

``CBD_PHY.gro`` - The Gromacs coordinates file.

``CBD_PHY.mdp`` - The Gromacs MD parameter file.

``CBD_PHY_cp2k.inp`` - The CP2K input file. Contains QM parameters. 

``CBD_PHY_cp2k.pdb`` - The pdb coordinates for CP2K. 


### To Run: 

    gmx grompp -f CBD_PHY.mdp -p CBD_PHY.top -c CBD_PHY.gro -n CBD_PHY.ndx -qmi CBD_PHY_cp2k.inp -o CBD_PHY.tpr
    gmx_mpi mdrun -s CBD_PHY.tpr
