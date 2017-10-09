#!/usr/bin/env bash
#This code prepares folders and creates an initial index file of groups we like. LIPIDS, SOLVENTS, and SYSTEM, which is the LIPIDS and SOLVENTS combined.
source /pool/butter/gromacs5/bin/GMXRC.bash #loading gmx commands

mkdir -p framesRaw #frames stored here
mkdir -p framesCount #Binned python data here.
mkdir -p framesPublish #Put ready data in here. (graphs, txt)

# this leaves groups POP2, DOPC, W and WF. (0,1,2,3)
gmx make_ndx -f "../solution_eq_323.gro" -o "selectedAtoms.ndx"<<EOF
r DOPC POP2
name 19 LIPIDS
r W WF
name 20 SOLVENTS
del 0-18
r DOPC POP2 W WF
name 2 SYSTEM

q
EOF
# del 2
