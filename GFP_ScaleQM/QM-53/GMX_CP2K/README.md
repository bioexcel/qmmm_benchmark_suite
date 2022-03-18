# GFP_ScaleQM - 53

* PBD ID: 1GFL 
* 28,264 atoms
* 53 QM atoms

## QM/MM GMX_CP2K benchmark

### About

This benchmark consists of a short QM/MM simulation using the Gromacs/CP2K
interface. 
5 MD steps are performed with a time step of 1 fs. The following XC functional 
set ups are included:

* BLYP - using DVZP-GTH-BLYP


``gfp_qm53.top`` - The Gromacs topology file.

``gfp_qm53.ndx`` - The Gromacs index file.

``gfp_qm53.vel.gro`` - The Gromacs coordinates and velocities file.

``gfp_qm53.mdp`` - The Gromacs MD parameter file.

``gfp_qm53.inp`` - The CP2K input file. Contains QM parameters. 

``gfp_qm53_cp2k.pdb`` - The pdb coordinates for CP2K. 


### To Run: 

    gmx grompp -f gfp_qm53.mdp -p gfp_qm53.top -c gfp_qm53.gro -n gfp_qm53.ndx -qmi gfp_qm53_cp2k.inp -o gfp_qm53.tpr
    gmx_mpi mdrun -s gfp_qm53.tpr
