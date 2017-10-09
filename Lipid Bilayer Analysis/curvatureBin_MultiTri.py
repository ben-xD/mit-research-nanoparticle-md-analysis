#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Binning lipid bilayer to allow more consistent curvature calculations."""

"""Bins lipid bilayer positions (radial distances) and turns
 histogram data into scatter data points. 
 Then curvature calculations are performed at every data point. 
 Every calculation uses the data point and the 2 neighbours.
 For each data point, 3 circles are generated, from the (1st, 2nd, 3rd)
 or (1st 3rd 4th) or (1st 2nd 4th) data points"""
 
#Created on Thu Jul 27 07:33:05 2017 @author: Adrian Butterworth
#import built in modules then custom, thats the convention.
import math, os
from scipy import optimize, odr, stats #odr is orthogonal distance regression.
import numpy as np
import matplotlib as mpl #first comment specifies the encoding. 8-bit Unicode Transformation Format
mpl.use("agg") #make sure the backend is right to save a file to check after running script.
mpl.rcParams["agg.path.chunksize"]=20000 
mpl.rcParams["lines.linewidth"]=1 #rcParams: rc means configuration file or run commands parameters.  
#rcParams is a dictionary like variable
#backend is not case sensitive. Agg can save PNGs. Agg stands for antigrain geometry
import matplotlib.pyplot as plt #the backend is specifically for pyplot.
r, z = np.genfromtxt("heads_r_pos.txt", unpack=True,dtype="float_")

#Ensuring graph folder is present.
if not os.path.exists("graph_binned/"):
    os.makedirs("graph_binned/")
	
#Definition
def circum_circle_radius(a, b, c):
	#function from http://www.adamfranco.com/2012/12/05/curvature-py/
  # Circumcircle radius calculation from http://www.mathopenref.com/trianglecircumcircle.html
  """If you know the length (a,b,c) of the three sides of a triangle,
  the radius of its circumcircle is given by the formula:"""
  if a > 0 and b > 0 and c > 0:
    try:
      divider = math.sqrt(math.fabs((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)))
      return (a * b * c) / divider
    except ZeroDivisionError:
      return 10000
  else:
    return 10000
	
#need to make the average line through middle. Then select upper leaf only.
maxRadius=150 #won't fit to values with radial distance, r  > 150
selection = np.where(r<=maxRadius)#gives element numbers of positions smaller than 150
fit=np.polyfit(r[selection],z[selection],2) #polynomial fit. Returns an array of coefficients.
def split(x):
	y=fit[0]*x**2+fit[1]*x+fit[2]
	return y
upperLeaf = np.where(z > split(r)) #array of elements numbers of upper leaflet.
rUpperLeaf=r[upperLeaf]; zUpperLeaf=z[upperLeaf] #work on upper leaf only.

#Binning data, and plotting mean/ median heights
#if within walls, then calculate average for that wall
binCount=20;#tailor this to make it smooth. 50 or below shows good curvature
binNumber = np.arange(1,binCount+1,1)
plt.figure(1); g2counts,g2bins,g2bars=plt.hist(rUpperLeaf, bins=binCount, label="Gradient_6nm_v1")
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#gives the middle r value of the bin
#^takes average of the wall values of each bin.^
plt.xlabel("Radial Distance / A"); plt.ylabel("Frequency")
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graph_binned/G1_Histogram") #includes 3 half plane plots

meanBinHeights,binEdges,binInts=stats.binned_statistic(rUpperLeaf,zUpperLeaf, bins=binCount)#gives the mean or median
plt.figure(2); plt.subplot(211)
plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
# plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
# plt.savefig("graph_binned/G2_Binned_Mean_Heights_against_Binned_Positions")
#binsMean, meanBinHeights = zip(*sorted(zip(binsMean, meanBinHeights)))

#Calculate continuous curvature of 3 data points
curvature=np.array([]); curvature1=np.array([]); curvature2=np.array([]); curvature3=np.array([])
#loop over every data point (starting from element 1(2nd))
for counter, _ in enumerate(binsMean):
	if counter >= 1 and counter < len(binsMean)-1: #can't calculate for data points on the edge
		#first calculate triangle side lengths
		#First 3 data points
		aSide1=math.sqrt((binsMean[counter-1]-binsMean[counter-2])**2+(meanBinHeights[counter-1]-meanBinHeights[counter-2])**2)
		bSide1=math.sqrt((binsMean[counter]-binsMean[counter-1])**2+(meanBinHeights[counter]-meanBinHeights[counter-1])**2)
		cSide1=math.sqrt((binsMean[counter]-binsMean[counter-2])**2+(meanBinHeights[counter]-meanBinHeights[counter-2])**2)
		#First, second, fourth data points (-2, -1, 0>1)
		aSide2=math.sqrt((binsMean[counter-1]-binsMean[counter-2])**2+(meanBinHeights[counter-1]-meanBinHeights[counter-2])**2)
		bSide2=math.sqrt((binsMean[counter]-binsMean[counter-1])**2+(meanBinHeights[counter]-meanBinHeights[counter-1])**2)
		cSide2=math.sqrt((binsMean[counter]-binsMean[counter-2])**2+(meanBinHeights[counter]-meanBinHeights[counter-2])**2)
		#First, third, fourth data points (-2, -1->0, 0->1)
		aSide3=math.sqrt((binsMean[counter]-binsMean[counter-2])**2+(meanBinHeights[counter]-meanBinHeights[counter-2])**2)
		bSide3=math.sqrt((binsMean[counter+1]-binsMean[counter])**2+(meanBinHeights[counter+1]-meanBinHeights[counter])**2)
		cSide3=math.sqrt((binsMean[counter+1]-binsMean[counter-2])**2+(meanBinHeights[counter+1]-meanBinHeights[counter-2])**2)
		#calc curvature
		curvature1=np.append(curvature,circum_circle_radius(aSide1,bSide1,cSide1))
		curvature2=np.append(curvature,circum_circle_radius(aSide2,bSide2,cSide2))
		curvature3=np.append(curvature,circum_circle_radius(aSide3,bSide3,cSide3))
		#average them
		curvature=np.append(curvature, np.mean([curvature1,curvature2,curvature3]))

plt.figure(2); plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("Curvature / A-1")
#print "The max curvature is ", np.amax(curvature)
# print len(binsMean), "vs", len(binsMean[1:-1]), "vs", len(curvature)
plt.scatter(binsMean[1:-1], 1/curvature)#boundaries?
plt.ylim(0,0.0002)
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graph_binned/G3_Curvature_against_Binned_Positions")