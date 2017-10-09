#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Binning lipid bilayer to allow more consistent odr-fitting:curvature calculations."""

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
import matplotlib.pyplot as plt
#having PROBLEMS? make sure you have a copy of
#matplotlibrc file in the same folder.

#Loading data
r, z = np.genfromtxt("heads_r_pos.txt", unpack=True,dtype="float_")
#Ensuring graph folder is present.
if not os.path.exists("graphs/binned/odr"):
    os.makedirs("graphs/binned/odr")
	
#Function definitions*********************************************
def circle(beta, x): 
	"""Circle function defined implicitly"""
	return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2
def calc_R(rc, zc):
	"""Calculates radius (distance of fitted points from circle center). Used during ODR fitting."""
	#distance of 2d points from center rc, zc. Returns a vector
	return np.sqrt((odrR-rc)**2 + (odrZ-zc)**2) #c means center.
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
binCount=70;#tailor this to make it smooth. 50 or below shows good curvature, but not enough datapoints
#Too large bin count causes a lot of noise to be shown. Subsequent curvature calculations are erroneous
g2counts,g2bins,g2bars=plt.hist(rUpperLeaf, bins=binCount, label="Gradient_6nm_v1")
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#gives the middle r value of the bin
#^takes average of the wall values of each bin.^
meanBinHeights,binEdges,binInts=stats.binned_statistic(rUpperLeaf,zUpperLeaf, bins=binCount)#gives the mean or median

#ODR fitting (5 points)
circleStepSize=5
curvatureODR=np.array([]); odrRBins=np.array([]);
print "Number of circles = ", binCount/circleStepSize
for count, value in enumerate(binsMean):
	if count%circleStepSize==0 & count < (len(binsMean)- 6):
		odrR=binsMean[count:count+6]
		odrZ=meanBinHeights[count:count+6]
		#fit data to circle
		#initial guess for circle fit
		R_m = calc_R(np.mean(odrR), np.mean(odrZ)).mean()
		beta0 = [ np.mean(odrR), np.mean(odrZ), R_m]#initial parameters
		#ODR********************
		lsc_data  = odr.Data(np.ma.row_stack([odrR, odrZ]), y=1)#what does this line do? The data to fit.
		#row_stack stacks arrays in sequence vertically. 
		#scalar input into odr.Data for y= value indicates implicit model
		lsc_model = odr.Model(circle, implicit=True) #The Model class stores information about the function you wish to fit.
		lsc_odr   = odr.ODR(lsc_data, lsc_model, beta0) #The ODR class gathers all information and coordinates the running of the main fitting routine.
		lsc_out   = lsc_odr.run()
		rc2, zc2, R = lsc_out.beta #estimated parameter values
		Ri		  = calc_R(rc2,zc2)
		#End of ODR*************
		curvatureODR=np.append(curvatureODR,Ri.mean())
		odrRBins=np.append(odrRBins,binsMean[count])

plt.figure(1);
plt.subplot(211); plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
#plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Radius of Curvature / A")
#plt.scatter(odrRBins, curvatureODR)
# plt.xlim((50,150));
plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Curvature / A-1")
plt.scatter(odrRBins, 1/curvatureODR)
# plt.xlim((50,150));
# plt.ylim((0,0.1));
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.tight_layout(); plt.show()
plt.savefig("graphs/binned/odr/G1_ODR_curvature_against_position_AND_height")

plt.figure(2); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Curvature / A-1")
plt.scatter(odrRBins, 1/curvatureODR)
# plt.xlim((50,150));
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graphs/binned/odr/G2_ODR_curvature_against_position")

curvatureText="The max curvature from ODR Fitting is "+str(max(1/curvatureODR))
print(curvatureText)

#Saving data to text file
f=open("graphs/binned/odr/curvature.txt","w+")
f.write(curvatureText)
f.close()