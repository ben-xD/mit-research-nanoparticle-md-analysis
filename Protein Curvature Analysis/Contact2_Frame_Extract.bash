#!/usr/bin/env bash
#This code loops over every frame, and creates the gro file for each frame.
source /pool/butter/gromacs5/bin/GMXRC.bash #loading gmx commands

#BASH: XTC to GRO. Frame separation/ extraction
#for number in {200..400..2}; do #Loop over 100 frames. every 2 frames, from 200 to 400
for number in {200..400..2}; do
frame=$((${number}*300)) #first frame 200, would be 60000
upperframe=$(((${number}+1)*300))
# let frame=${number}*300 #declaring integer, then calculating using rpm's bc package
# Bc is an arbitrary precision numeric processing arithmetic language.
outputFrameName="framesRaw/frame${number}.gro"
gmx trjconv -f 	"../solution_md_323.xtc" -o "${outputFrameName}" \
-b $frame \
-e $upperframe \
-s "../solution_md_323.tpr" \
-n "selectedAtoms.ndx" <<EOF
2
EOF

#The code above only uses GROUP 2. we lose group 0 and 1 information.

#in time units (10000*0.03*frame number) So frame 100 is 30000
#frame 100 is 30 nanoseconds, 30000ps
done