#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ProtPos2: Protein minimum distance calculations."""
"""Created on Thurs Jul 20 10:18: 2017 @author: Adrian Butterworth"""
import numpy as np
import math, os
from scipy.optimize import leastsq
import matplotlib as mpl
# from scipy.interpolate import spline
mpl.use("agg"); mpl.rcParams["agg.path.chunksize"]=40000; mpl.rcParams["lines.linewidth"]=1
import matplotlib.pyplot as plt
#import MDAnalysis #to calc min distances

#Ensuring graph folder is present.
if not os.path.exists("graphs/plotB/"):
    os.makedirs("graphs/plotB/")

#load data
minTime, minDist = np.genfromtxt("min_dist_heads.xvg", skip_header=19, unpack=True, dtype="float_")
protTime, protDist, protHeight = np.genfromtxt("protein_r_pos.txt", unpack=True, dtype="float_")

#GRAPH 0: Protein Minimum Distance to Bilayer (distance against time)
plt.figure(1); plt.xlabel("Time / ps"); plt.ylabel("Distance / nm")
plt.scatter(minTime,minDist, c=protTime, label="Gradient_6nm_v1") #plt.legend(); 
plt.title("Protein Minimum Distance to Bilayer"); plt.legend()
plt.show(); plt.savefig("graphs/plotB/G0_Protein_Minimum_Distance_to_Bilayer.png")

#calc desired cutoff amount to count as IN CONTACT
#buffer=4 #buffer. Check graph to see cutoff above average value.
# cutoff=np.mean(protDist)+buffer
# #Remove non-contact data points
# removal=np.where(protDist<=cutoff)
# protDistClose=protDist[removal]
# protHeightClose=protHeight[removal]
# protTimeClose=protTime[removal]
cutoff=np.percentile(np.sort(minDist), 15)+0.5
removal=np.where(minDist<=cutoff) #***change cutoff to higher to include ALL data.
#Remove non-contact data points
#removal=np.where(protDist<=cutoff)
protDistClose=protDist[removal]
protHeightClose=protHeight[removal]
protTimeClose=protTime[removal]
#Use contact data points to calculate graphs
#GRAPH 1.1: Scatter - Height vs Radial Distance (repeat this graph for all data: FullPlaneScatter)
plt.figure(2); plt.xlabel("Radial Distance / unit"); plt.ylabel("Height / unit")
plt.scatter(protDistClose, protHeightClose, c=protTimeClose, label="Gradient_6nm_v1");
plt.title("Protein Minimum Distance to Bilayer")
plt.colorbar(); plt.legend(); 
plt.show(); plt.savefig("graphs/plotB/G1_Protein_Minimum_Distance_to_Bilayer.png")
#GRAPH 1.2:
plt.figure(3); plt.xlabel("Radial Distance / unit"); plt.ylabel("Height / unit")
plt.scatter(protDist, protHeight, c=protTime, label="Gradient_6nm_v1");
plt.title("Protein Minimum Distance to Bilayer")
plt.colorbar();
plt.show(); plt.savefig("graphs/plotB/G1.2_Protein_Minimum_Distance_to_Bilayer_AllData.png")
#GRAPH 1.3: Density Map - GRAPH 1.1 but Density Map


#GRAPH 2: Histogram - Frequency against radial distances
binCount=100;
plt.figure(4); plt.xlabel("Radial Distance / unit"); plt.ylabel("Frequency")
g2counts,g2bins, g2bars=plt.hist(protDistClose, bins=binCount, label="Gradient_6nm_v1")
plt.title("Frequency against Radial Distance"); plt.legend(); limitX=150;
plt.xlim(0, limitX);
plt.show(); plt.savefig("graphs/plotB/G2_Frequency against radial distance.png")

#GRAPH 3: Energetics - Flip histogram and make more continuous. Fit? (First repeat previous graph)
#GRAPH 3.1: Repeat GRAPH 2 (Subplot 1)
plt.figure(5); 
plt.subplot(211); plt.hist(protDistClose, bins=binCount, label="Gradient_6nm_v1")
plt.xlabel("Radial Distance / unit"); plt.ylabel("Frequency")
plt.title("Frequency against Radial Distance")
plt.xlim(0, limitX);
#GRAPH 3.2: Energetics (Subplot 1)
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#array of bin means
plt.subplot(212); plt.scatter(binsMean, -g2counts, label="Gradient_6nm_v1")
#plt.title("Energetics") #title for bottom graph will overlap
plt.xlim(0, limitX);

energyP = np.polyfit(binsMean, -1*g2counts, 5)
binsMean=np.asarray(binsMean)
g2countsEst=energyP[0]*binsMean**5 + energyP[1]*binsMean**4 + energyP[2]*binsMean**3 + energyP[3]*binsMean**2 + energyP[4]*binsMean + energyP[5]
plt.plot(binsMean,g2countsEst)
plt.show(); plt.savefig("graphs/plotB/G3_Frequency_against_radial_distance_AND_Energetics.png")