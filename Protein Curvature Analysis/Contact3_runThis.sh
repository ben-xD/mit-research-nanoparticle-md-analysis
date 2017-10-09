#!/bin/bash
#This code will copy Contact3_many repeatedly, into Contact3_many_tmp.
#Contact3_many_tmp will then be sent to the cluster with sbatch.

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

let binSize=4

#need to ensure matching of number of bins in below small section and later on binning
for ((xDim=0; xDim<=40; xDim=xDim+$binSize)) {
#Loop over each yvalue
for ((yDim=0; yDim<=40; yDim=yDim+$binSize)) {
echo "${xDim} ${yDim}" | cat >>${xyData}
}
}

# for frame in framesRaw/*.gro; do
for frame in framesRaw/frame202.gro; do
echo $frame
text=${frame%.*}
echo $text
text2=${text#*/}
echo $text2
number=${text2#frame}
echo $number

cp Contact3_many.sh Contact3_many_tmp.sh
grep -l 'nmbr=' Contact3_many_tmp.sh | xargs sed -i "s/^nmbr\=.*/nmbr=$number/"
sbatch Contact3_many_tmp.sh
done
