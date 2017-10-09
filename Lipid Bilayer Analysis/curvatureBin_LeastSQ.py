#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Binning lipid bilayer to allow more consistent least-square method fitting of curvature."""

"""Bins lipid bilayer positions (radial distances) and turns
 histogram data into scatter data points. 
 Then curvature calculations are performed at every data point. 
 Every calculation uses the data point and the 2 neighbours."""
 
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

#Loading data
r, z = np.genfromtxt("heads_r_pos.txt", unpack=True,dtype="float_")
#Ensuring graph folder is present.
if not os.path.exists("graph_binned/leastSQ"):
    os.makedirs("graph_binned/leastSQ")
	
#Function definitions*********************************************
def calc_R(rc, zc): #Used for circle fitting later
	#distance of 2d points from center rc, zc. Returns a vector
	return np.sqrt((lsqR-rc)**2 + (lsqZ-zc)**2) #c means center.
def circum_circle_radius(a, b, c):
	"""Calculates radius of circumcircle for Triangle (3 point -> Circle) Fitting """
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
		"""calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc)"""
def f(c):
	#calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc)
	Ri = calc_R(*c) #what is the meaning of this. *c can be a set of arguments
	return Ri - Ri.mean()
	#End of Function definitions**************************************

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
binCount=50;#tailor this to make it smooth. 50 or below shows good curvature
binNumber = np.arange(1,binCount+1,1); plt.figure(1)
g2counts,g2bins,g2bars=plt.hist(rUpperLeaf, bins=binCount, label="Gradient_6nm_v1")
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#gives the middle r value of the bin
#^takes average of the wall values of each bin.^

meanBinHeights,binEdges,binInts=stats.binned_statistic(rUpperLeaf,zUpperLeaf, bins=binCount)#gives the mean or median
plt.figure(2)
plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graph_binned/leastSQ/G1_Heights_against_Positions")
binsMean, meanBinHeights = zip(*sorted(zip(binsMean, meanBinHeights)))

#LeastSQ fitting (#circleStepSize points)
circleStepSize=5
curvatureLSQ=np.array([]); lsqRBins=np.array([]);
print "Number of circles = ", binCount/circleStepSize
for count, value in enumerate(binsMean):
	if count%circleStepSize==0 & count < (len(binsMean)- (circleStepSize+1)):
		lsqR=binsMean[count:count+(circleStepSize+1)]
		lsqZ=meanBinHeights[count:count+(circleStepSize+1)]
		#initial guess for circle fit
		center_estimate = np.mean(lsqR), np.mean(lsqZ) #estimate initial
		#LSQ********************
		center, ier = optimize.leastsq(f, center_estimate) #send estimate, and function in. Sum of squares of equations minimised.
		rC_2, zC_2 = center #solutions
		Ri = calc_R(rC_2, zC_2)
		#End of LSQ*************
		curvatureLSQ=np.append(curvatureLSQ,Ri.mean())
		lsqRBins=np.append(lsqRBins,binsMean[count])
		
plt.figure(3); 
plt.subplot(211); plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
plt.title("Curvature of Lipid Bilayer (Custom v0)");
plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("Curvature / A-1")
# print len(lsqRBins), "vs", len(curvatureLSQ)
plt.scatter(lsqRBins, 1/curvatureLSQ) #boundaries?
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graph_binned/leastSQ/G2_Curvature_against_Positions_AND_G1")

print "The max curvature from Least Square Fitting is ", np.amax(1/curvatureLSQ)