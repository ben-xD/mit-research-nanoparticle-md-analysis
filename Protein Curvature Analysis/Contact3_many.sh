#!/bin/bash
#SBATCH -J Contact3
#SBATCH -o framesOutput/contact3.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p extended-cpu

source /pool/butter/gromacs5/bin/GMXRC.bash
#Bin size? 0 to 400 in x, 0 to 400 in y
# let binCount=40
#change binCount back to 200 or something! Just testing
# let binSize=40/${binCount}
let binSize=4
nmbr=fn
frameText="framesRaw/frame${nmbr}.gro"
echo $frameText

#To select frames to calculate on, provide copy of frameText in frameRaw folder.
# for frameText in framesRaw/*.gro; do
#column 3 for frameText 1. col 4 for frameText 2...
#Create text file with x, y values
#Run binning on GRO file
#loop over each x value
for ((xDim=0; xDim<=40; xDim=xDim+$binSize)) {
for ((yDim=0; yDim<=40; yDim=yDim+$binSize)) {
#1. Need to create a new binSelection every bin. between xmiddle, ymiddle, zmax
#http://manual.gromacs.org/documentation/5.1/onlinehelp/selections.html
#The select position units may be completely different...
let XX1=$xDim #lower end
let XX2=$xDim+$binSize #higher end
let YY1=$yDim
let YY2=$yDim+$binSize
#gmx select to restrict region of interest. Creating an index file.
#selectedAtoms.ndx has 3 groups, Solvent, Lipids, and SYSTEM (everything). 
gmx select -s $frameText -n selectedAtoms.ndx -selrpos whole_res_com -on binSelection.ndx -select "x>=$XX1 and x<=$XX2 and y>=$YY1 and y<=$YY2" 
#Removed AND group SOLVENT LIPIDS
#NEED TO CREATE TWO GROUPS, so mindist works on the difference between the two!
#Change selection reference position to something better

#[cutoff REAL] cutoff is the cutoff distance above which all numbers are not used. REAL is distance which is used as index
#Using selection keywords to select atoms in specific X and Y positions. Any Z is ok.
#Use parameter expansion to select part of string and use it in naming index string; originally $frameText=framesRaw/frameXXX.gro
#We just want the number!
#frameText = "framesRaw/frame200.gro"
text=${frameText%.*}
# echo $text
text2=${text#*/}
# echo $text2
number=${text2#frame}
# echo $number

indexFileName="framesCount/numCountFrame${number}_X${xDim}_Y${yDim}.xvg"
xyDataFrameTemp="framesCount/xyDataFrameTemp.txt"

#Calculating contacts using mindist, for each bin
#This function doesn't get to create an index file for all the bins! WHY: There is one group in the index only. The areas are too small, and for the small bins, only one can be found?
gmx mindist \
-f ${frameText} \
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
#pastes the temporary frameText file into the working frameText count file
# pr -mts" " "framesCount/.tempframe.txt" "framesCount/xyframescount.txt"
# paste -d" " "framesCount/.tempframe.txt" "framesCount/xyframescount.txt"
}
}

#finds and deletes recovered overwritten files made by gromacs
find . -name "#binSelection.ndx.*" -exec rm {} \;
find . -name "#mindist.xvg.*" -exec rm {} \;

# paste ${xyData} ${xyDataFrameTemp}| cat >>${xyData}

# done