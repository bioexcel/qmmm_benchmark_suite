# Conversion of CP2K (amber) inputs to Gromacs inputs with Parmed

This document describes how to convert CP2K input file with Amber
topologies into the GROMACS format for the GROMACS(+CP2K) interface.

## Install Parmed

```
wget https://github.com/ParmEd/ParmEd/archive/refs/tags/3.4.3.tar.gz
tar xvf 3.4.3.tar.gz
cd ParmEd-3.4.3
python setup.py install --user
```

## Using Parmed to convert amber topology to .top

```
import parmed as pmd
pdb = pmd.load_file("ClC.pdb")
pdb.save("ClC.inpcrd", format="rst7")
amber = pmd.load_file("ClC.prmtop", "ClC.inpcrd")
amber.save("ClC.top")
amber.save("ClC.gro")
```

You may want to adjust the box size at the end of the .gro file.

## Generate ndx file

```
gmx make_ndx -f ClC.gro -o ClC.ndx
'q': save and quit
```

Edit to add QMatoms e.g.
```
...

 [ QMatoms ]
2887 2888 2080 2079 2880 2883 2886 2072 2075 2078 2884 2885 2881 2882 2074 2073 2077 2076 2081
```

## Create mdp file

For example:
(pay attention to the QMMM options)

```
; md.mdp
integrator  = md    ; MD using leap-frog integrator
dt          = 0.001 ; 1fs time-step
nsteps      = 100   ; 100 fs simulation

; Set output frequency to each step
nstxout                  = 1 ; Coordinates to trr
nstlog                   = 1 ; Energies to md.log
nstcalcenergy            = 1 ; Energies to ener.edr
nstenergy                = 1 ; Energies to ener.edr

; Set cut-offs
rlist               = 0.2 ; NB-search cut-off
rcoulomb	    = 0.2 ; Short-range electrostatic cut-off
rvdw		    = 0.2 ; Short-range Van der Waals cut-off

;Temperature coupling options
tcoupl                   = v-rescale 
nsttcouple               = 1
tc-grps                  = System
tau-t                    = 0.1
ref-t                    = 300

; GENERATE VELOCITIES FOR STARTUP RUN
gen-vel                  = yes
gen-temp                 = 300 
gen-seed                 = -1

; CP2K QMMM parameters
qmmm-cp2k-active              = true   ; Activate QMMM MdModule
qmmm-cp2k-qmgroup             = QMatoms; Index group of QM atoms
qmmm-cp2k-qmmethod            = PBE    ; Method to use
qmmm-cp2k-qmcharge            = 1      ; Charge of QM system
qmmm-cp2k-qmmultiplicity      = 1      ; Multiplicity of QM system
```

## Running with Gromacs

```
export NAME=ClC
gmx_mpi grompp -f ${NAME}.mdp -p ${NAME}.top -c ${NAME}.gro -n ${NAME}.ndx  -o ${NAME}.tpr
srun gmx_mpi mdrun -ntomp 1 -s ${NAME}.tpr
```
