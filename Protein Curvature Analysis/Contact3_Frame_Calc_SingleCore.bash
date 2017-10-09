#!/usr/bin/env bash
#Initial/ old code. Loops over all frame files, ran on one node.
source /pool/butter/gromacs5/bin/GMXRC.bash #loading gmx commands

#finds and deletes recovered overwritten files made by gromacs
find . -name "#binSelection.ndx.*" -exec rm {} \;
find . -name "#mindist.xvg.*" -exec rm {} \;
cd framesCount
rm countOnlyFrame*
rm *#numCountFrame*
rm *numCountFrame*
rm xyDataFrameTemp.txt
rm xyframescount.txt
cd ../
xyData="framesCount/xyframescount.txt"

#This script bins over small x and y areas over each frame. It then calculates the number of contacts between solvent and lipid molecules.

#Bin size? 0 to 400 in x, 0 to 400 in y
# let binCount=50
# let binSize=40/${binCount}
let binSize=1

#need to ensure matching of number of bins in below small section and later on binning
for ((xDim=0; xDim<=40; xDim=xDim+$binSize)) {
#Loop over each yvalue
for ((yDim=0; yDim<=40; yDim=yDim+$binSize)) {
echo "${xDim} ${yDim}" | cat >>${xyData}
}
}
#To select frames to calculate on, provide copy of frame in frameRaw folder.
for frame in framesRaw/*.gro; do
echo $frame
#column 3 for frame 1. col 4 for frame 2...
#Create text file with x, y values
#Run binning on GRO file
#loop over each x value
for ((xDim=1; xDim<=40; xDim=xDim+$binSize)) {
for ((yDim=1; yDim<=40; yDim=yDim+$binSize)) {
#1. Need to create a new binSelection every bin. between xmiddle, ymiddle, zmax
#http://manual.gromacs.org/documentation/5.1/onlinehelp/selections.html

#A constant position can be defined as [XX, YY, ZZ]
let XX1=$xDim #lower end
let XX2=$xDim+$binSize #higher end
let YY1=$yDim
let YY2=$yDim+$binSize
#gmx select to restrict region of interest. Creating an index file.
gmx select -s $frame -selrpos whole_res_com -on binSelection.ndx -select "x>$XX1 and x<=$XX2 and y>$YY1 and y<=$YY2" 

#[cutoff REAL] cutoff is the cutoff distance above which all numbers are not used. REAL is distance which is used as index
#Using selection keywords to select atoms in specific X and Y positions. Any Z is ok.
#Use parameter expansion to select part of string and use it in naming index string; originally $frame=framesRaw/frameXXX.gro
#We just want the number!
#framesRaw/frame200.gro
text=${frame%.*}
text2=${text#*/}
number=${text2#frame}

indexFileName="framesCount/numCountFrame${number}_X${xDim}_Y${yDim}.xvg"
xyDataFrameTemp="framesCount/xyDataFrameTemp.txt"

#Calculating contacts using mindist, for each bin
gmx mindist \
-f ${frame} \
-n binSelection.ndx \
-on ${indexFileName} \
-d 0.8 #Distance for contact counting. Units are in Angstroms? NO, its Nanometers
#Check now -on saves the data, then extract data below #2

#2. Load data from on file. Extract number of contacts only, put it into neater. Use parameter expansion to organise data neatly in rows and columns
#*********************NEED TO ADD EXTRACTION OF VALUE ONLY: adding to new column (x=1, cycle through y)
extractedText=`tail -1 ${indexFileName}`
#Remove first number before tab and tab (or 17 characters)
#This works only when the file doesn't exist? OTHERWISE no data is output :(
echo ${extractedText:17}|cat >>framesCount/countOnlyFrame${number}.txt

#*********************NEED TO CREATE DATA FILE of VALUES ABOVE
#concatenate the file with data each round: X,Y,#contactsF1, #contactsF2, contactsF3....
#make a new text temporary file, which then is pasted into previous file
# echo "${xDim} ${yDim}" | cat >>framesCount/.tempframe.txt
#pastes the temporary frame file into the working frame count file
# pr -mts" " "framesCount/.tempframe.txt" "framesCount/xyframescount.txt"
# paste -d" " "framesCount/.tempframe.txt" "framesCount/xyframescount.txt"
}
}
#the below gets looped for every frame. So every frame is a new column

#finds and deletes recovered overwritten files made by gromacs
find . -name "#binSelection.ndx.*" -exec rm {} \;
find . -name "#mindist.xvg.*" -exec rm {} \;

# paste ${xyData} ${xyDataFrameTemp}| cat >>${xyData}

done