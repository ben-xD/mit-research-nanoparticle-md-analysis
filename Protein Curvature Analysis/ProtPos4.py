#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protein Position Analysis and Graphs"""
"""Created on Thurs Jul 20 10:18: 2017 @author: Adrian Butterworth"""

import numpy as np
import numpy.ma as ma
import math, os
import scipy.interpolate as inter
from scipy.optimize import leastsq
import matplotlib as mpl
mpl.rcParams["agg.path.chunksize"]=40000;
import matplotlib.pyplot as plt

#PROBLEMS? get the prepared matplotlibrc file in folder
# Angstrom is given by "+r"$\AA$"

#Ensuring graph folder is present.
if not os.path.exists("graphs/ProtPos/"):
    os.makedirs("graphs/ProtPos/")

#load data
minTime, minDist = np.genfromtxt("min_dist_heads.xvg", skip_header=19, unpack=True, dtype="float_")
protTime, protDist, protHeight = np.genfromtxt("protein_r_pos.txt", unpack=True, dtype="float_")
tProt, xProt, yProt, zProt = np.genfromtxt("protein_xyz_pos.txt", unpack=True, dtype="float_")
#Use tProt, xProt, yProt for solvent <-> Tail contacts.


#GRAPH 0: Protein Minimum Distance to Bilayer (distance against time)
plt.figure(1); plt.xlabel("Time / ps"); plt.ylabel("Distance / nm")
plt.scatter(minTime,minDist, c=protTime, label="Gradient_6nm_v1") #plt.legend(); 
plt.title("Protein Minimum Distance to Bilayer"); plt.legend()
plt.show(); plt.savefig("graphs/ProtPos/G0_Protein_Minimum_Distance_to_Bilayer.png")

#calc desired cutoff amount to count as IN CONTACT
buffer=0.5 #buffer. Check graph to see cutoff above average value.
cutoff=np.percentile(np.sort(minDist), 15)+buffer
removal=(minDist<=cutoff)
# Remove non-contact data points
protDistClose=protDist[removal]
protHeightClose=protHeight[removal]
protTimeClose=protTime[removal]
tProt=tProt[removal]
xProt=xProt[removal]
yProt=yProt[removal]
zProt=zProt[removal]
#Below 3 lines are used to graph all data, including non-contact
# protDistClose=protDist
# protHeightClose=protHeight
# protTimeClose=protTime


#Use contact data points to calculate graphs
#GRAPH 1.1: Scatter - Height vs Radial Distance (repeat this graph for all data: FullPlaneScatter)
plt.figure(2); plt.xlabel("Radial Distance / "+r"$\AA$"); plt.ylabel("Height / "+r"$\AA$")
plt.scatter(protDistClose, protHeightClose, c=protTimeClose, label="Gradient_6nm_v1");
plt.title("Protein Minimum Distance to Bilayer (In Contact Only)")
plt.colorbar(); plt.legend(); 
plt.show(); plt.savefig("graphs/ProtPos/G1.1_Protein_Minimum_Distance_to_Bilayer.png")
#GRAPH 1.2:
plt.figure(3); plt.xlabel("Radial Distance / "+r"$\AA$"); plt.ylabel("Height / "+r"$\AA$")
plt.scatter(protDist, protHeight, c=protTime, label="Gradient_6nm_v1");
plt.title("Protein Minimum Distance to Bilayer")
plt.colorbar();
plt.show(); plt.savefig("graphs/ProtPos/G1.2_Protein_Minimum_Distance_to_Bilayer_AllData.png")
#GRAPH 1.3: Heat Map (Height vs Radial Distance)
binCount=100;
countHeat,xEdges,yEdges=np.histogram2d(protDistClose,protHeightClose,bins=100)
plt.figure(4); plt.xlabel("Radial Distance / "+r"$\AA$"); plt.ylabel("Height / "+r"$\AA$")
plt.pcolormesh(xEdges,yEdges,countHeat,cmap="Blues",alpha=0.8, vmin=countHeat.min(),vmax=countHeat.max())
# plt.xlim(0,xEdges.max());plt.ylim(0,yEdges.max());
plt.title("Height against Radial Distance - Heat Map (In Contact Only)")
plt.colorbar();
plt.show(); plt.savefig("graphs/ProtPos/G1.3_Frequency_Against_Position_Heat_Map.png")


#GRAPH 2: Histogram - Frequency against radial distances
plt.figure(5); plt.xlabel("Radial Distance / "+r"$\AA$"); plt.ylabel("Frequency")
g2counts,g2bins, g2bars=plt.hist(protDistClose, bins=binCount, label="Gradient_6nm_v1")
plt.title("Frequency against Radial Distance"); plt.legend(); limitX=150;
plt.xlim(0, limitX);
plt.show(); plt.savefig("graphs/ProtPos/G2_Frequency against radial distance.png")


#GRAPH 3: Energetics - Flip histogram and make more continuous. Fit? (First repeat previous graph)
#GRAPH 3.1: Repeat GRAPH 2 (Subplot 1)
plt.figure(6); 
plt.subplot(211); plt.hist(protDistClose, bins=binCount, label="Gradient_6nm_v1")
plt.ylabel("Frequency")
plt.title("Frequency against Radial Distance")
plt.xlim(0, limitX);
#GRAPH 3.2: Energetics (Subplot 1)
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#array of bin means
plt.subplot(212); plt.scatter(binsMean, -g2counts, label="Gradient_6nm_v1")
plt.xlabel("Radial Distance / "+r"$\AA$"); plt.ylabel("Energy / Arbitrary Units")
#plt.title("Energetics") #title for bottom graph will overlap
plt.xlim(0, limitX);

energyP = np.polyfit(binsMean, -1*g2counts, 5)
binsMean=np.asarray(binsMean)
g2countsEst=energyP[0]*binsMean**5 + energyP[1]*binsMean**4 + energyP[2]*binsMean**3 + energyP[3]*binsMean**2 + energyP[4]*binsMean + energyP[5]
plt.plot(binsMean,g2countsEst)
plt.show(); plt.savefig("graphs/ProtPos/G3_Frequency_against_radial_distance_AND_Energetics.png")


#GRAPH 4: First top down scatter plot
plt.figure(7);plt.xlabel("X Position / "+r"$\AA$"); plt.ylabel("Y Position / "+r"$\AA$")
plt.title("Protein Position over time");
#Why does this scatter plot just give a straight line?? Because of radial_traj_v2
plt.scatter(xProt,yProt, c=tProt, label="Gradient_6nm_v1");plt.colorbar()
#plt.xlim();plt.ylim();#can change limits to make it a circle
plt.show(); plt.savefig("graphs/ProtPos/G4_Top_Down_Protein.png")


#GRAPH 5: Birds Eye View Heat Map
plt.figure(8); #The figure module provides top level artist (figure), containing all plot elements
fig, ax = plt.subplots() #need to point to ax, to add circles later. This line creates a figure!
plt.xlabel("X Position / "+r"$\AA$"); plt.ylabel("Y Position / "+r"$\AA$")
Z, xBins, yBins = np.histogram2d(xProt,yProt, bins=150);#************************change to affect heatmap
#Plotting circles
circle1 = plt.Circle((200, 200), 60, alpha=0.5, label="Inner Crystal Wall", color='k', fill=False)#60 and 110 are radii
circle2 = plt.Circle((200, 200), 110, alpha=0.5, label="Outer Crystal Wall", color='k', fill=False)
ax.add_artist(circle1); ax.add_artist(circle2)
#calculating middle x and y value of each bin
xbinsMean=[0.5*(xBins[i]+xBins[i+1]) for i in range(len(Z))]
ybinsMean=[0.5*(yBins[i]+yBins[i+1]) for i in range(len(Z))]
#Interpolation for smoothing. Smoothing is not good. Then i broke the code.
# YBins, XBins=np.meshgrid(xbinsMean,ybinsMean) #NOTICE X and Y got SWAPPED
# f=inter.interp2d(XBins,YBins,Z, kind="cubic")  #new interpolation function
# f=inter.RectBivariateSpline(xBins,yBins,Z)
#calculate the interpolated values
# ZInter=f(xBins,yBins)
#Reassign new calculated values to Z
# Z=ZInter

cmap = plt.get_cmap("Reds")
cmap.set_under("white")

percentileSelect=np.percentile(Z,78) #calculates percentile. 78 is good
Z=ma.masked_array(Z, mask=Z<percentileSelect) #Want less than a certain percentile (less frequency). Theyre the masked ones.
plt.pcolormesh(xBins,yBins,Z.T,cmap=cmap,alpha=1) #3D map! pcolormesh! Does this even work correctly. Its flipped
plt.xlim(0,400); plt.ylim(0,400); plt.colorbar()
plt.title("Birds Eye Positions Heat Map");
plt.show(); plt.savefig("graphs/ProtPos/G5_2D_Heat_Map.png")


#GRAPH 6: Protein Movements over time (plot every few frames)
plotEvery=60#frames
sparseT=np.array([])
sparseX=np.array([]); sparseY=np.array([]); sparseZ=np.array([])
for counter,time in enumerate(tProt):
	if counter%plotEvery==0:
		sparseT=np.append(sparseT,time)
		sparseX=np.append(sparseX,xProt[counter])
		sparseY=np.append(sparseY,yProt[counter])
		sparseZ=np.append(sparseZ,zProt[counter])
#Why plot and scatter? plot gives the lines
plt.plot(sparseX, sparseY, alpha=0.2,c="k") #reduce opacity of black line between sparse data points
#Scatter gives coloured data points
plt.scatter(sparseX,sparseY, alpha=0.7, c=sparseT, marker=".", edgecolors='none'); plt.colorbar();
plt.title("Protein Movements over time (plot every %i frames)" % plotEvery);
plt.show(); plt.savefig("graphs/ProtPos/G6_Time_AND_G5.png")


#GRAPH 7: Save G6 separately
plt.figure()
fig2, ax2 = plt.subplots() #need to point to ax, to add circles later. This line creates a figure!
plt.xlabel("X Position / "+r"$\AA$"); 
plt.ylabel("Y Position / "+r"$\AA$")
plt.plot(sparseX, sparseY, alpha=0.2,c="k")
plt.scatter(sparseX,sparseY, alpha=1, c=sparseT, marker="o", edgecolors='none'); plt.colorbar();
circle1 = plt.Circle((200, 200), 60, alpha=0.5, label="Inner Crystal Wall", color='k', fill=False)#60 and 110 are radii
circle2 = plt.Circle((200, 200), 110, alpha=0.5, label="Outer Crystal Wall", color='k', fill=False)
ax2.add_artist(circle1)
ax2.add_artist(circle2)
plt.xlim(0, 400)
plt.ylim(0, 400)
plt.title("Protein Movements over time (plot every %i frames)" % plotEvery);
plt.show()
plt.savefig("graphs/ProtPos/G7_Time.png")