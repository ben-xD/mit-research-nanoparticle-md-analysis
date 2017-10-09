# -*- coding: utf-8 -*- #character encoding
"""Created on Thurs Jul 20 10:18: 2017
@author: Adrian Butterworth"""
import numpy as np
import math
import matplotlib as mpl
mpl.use("agg"); mpl.rcParams["agg.path.chunksize"]=20000; mpl.rcParams["lines.linewidth"]=1
import matplotlib.pyplot as plt
#import MDAnalysis #to calc min distances

#load data
minTime, minDist = np.genfromtxt("min_dist_heads.xvg", skip_header=19, unpack=True, dtype="float_")
protTime, protDist, protHeight = np.genfromtxt("protein_r_pos.txt", unpack=True, dtype="float_")

#GRAPH 0: Protein Minimum Distance to Bilayer (distance against time)
plt.figure(1); plt.xlabel("Time / ps"); plt.ylabel("Distance / nm")
plt.scatter(minDist,minTime, c=protTime, label="Gradient_6nm_v1") #plt.legend(); 
plt.title("Protein Minimum Distance to Bilayer"); plt.legend()
plt.show(); plt.savefig("Graphs/G0_Protein_Minimum_Distance_to_Bilayer.png")

#calc desired cutoff amount to count as IN CONTACT
#buffer=4 #buffer. Check graph to see cutoff above average value.
# cutoff=np.mean(protDist)+buffer
# #Remove non-contact data points
# removal=np.where(protDist<=cutoff)
# protDistTouch=protDist[removal]
# protHeightTouch=protHeight[removal]
# protTimeTouch=protTime[removal]
cutoff=np.mean(minDist)+500000
#Remove non-contact data points
removal=np.where(protDist<=cutoff)
protDistTouch=protDist[removal]
protHeightTouch=protHeight[removal]
protTimeTouch=protTime[removal]

#Use contact data points to calculate graphs
#GRAPH 1: Scatter - Height vs Radial Distance (repeat this graph for all data: FullPlaneScatter)
plt.figure(2); plt.xlabel("Radial Distance / unit"); plt.ylabel("Height / unit")
plt.scatter(protDistTouch, protHeightTouch, c=protTimeTouch, label="Gradient_6nm_v1");
plt.title("Protein Minimum Distance to Bilayer")
plt.colorbar(); plt.legend(); 
plt.show(); plt.savefig("Graphs/G1_Protein_Minimum_Distance_to_Bilayer.png")

#GRAPH 2: Histogram - Frequency against radial distances
plt.figure(3); plt.xlabel("Radial Distance / unit"); plt.ylabel("Frequency")
plt.hist(protDistTouch, bins=15,label="Gradient_6nm_v1")
plt.title("Frequency against Radial Distance"); plt.legend(); limitX=150;
plt.xlim(0, limitX);
plt.show(); plt.savefig("Graphs/G2_Frequency against radial distance.png")

#GRAPH 3: PseudoEnergetics - Flip histogram and make more continuous. Fit? (First repeat previous graph)
plt.figure(4); 
plt.subplot(211); plt.hist(protDistTouch, protTimeTouch, label="Gradient_6nm_v1") #make this a subplot with Graph 2
plt.xlabel("Radial Distance / unit"); plt.ylabel("Frequency")
plt.title("Frequency against Radial Distance")
plt.xlim(0, limitX);
protDistTouchOrg=protDistTouch[:]; protDistTouchOrg.sort()
hist, bin_edges=np.histogram(protDistTouch, bins=15, density=1)
print(hist)
hist=np.append([0],hist)
print(hist)
#wanna plot frequency distance
width=np.amax(protDistTouch)/bin_edges;
x=0; protDistTouchN=[0]
for counter in protDistTouchOrg:
	x=x+width
	protDistTouchN=np.append(protDistTouchN,x); 
print("x size is " + str(np.size(x)))
print("hist size is "+str(np.size(hist)))
plt.subplot(212); plt.scatter(x, hist, label="Gradient_6nm_v1")
#plt.title("PseudoEnergetics") #title for bottom graph will overlap
plt.xlim(0, limitX);
plt.show(); plt.savefig("Graphs/G3_Frequency against radial distance AND PseudoEnergetics.png")

#GRAPH 4: Density Map - GRAPH 1 but Density Map


#GRAPH 5: Top Down Density Map - need to analyse gromacs data again


#GRAPH 6: Histogram - Number of contacts against radial distance