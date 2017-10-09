#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Both, meaning ODR and Tri. Not MultiTri or LeastSQ.
Binning lipid bilayer to allow more consistent curvature calculations."""

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
import matplotlib.pyplot as plt #the backend is specifically for pyplot.
r, z = np.genfromtxt("heads_r_pos.txt", unpack=True,dtype="float_")
#PROBLEMS? make sure matplotlibrc is present in dir

#Ensuring graph folder is present.
if not os.path.exists("graphs/binned/mixed"):
    os.makedirs("graphs/binned/mixed")
	
#Function definitions*********************************************
#For ODR Define implicitly, circle function for fitting
def circle(beta, x): #beta are the parameters
	return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2
def calc_R(rc, zc): #Used for circle fitting later
	#distance of 2d points from center rc, zc. Returns a vector
	return np.sqrt((odrR-rc)**2 + (odrZ-zc)**2) #c means center.
#For Triangle (3 point -> Circle) Fitting 
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
binCount=50;#tailor this to make it smooth. 50 or below shows good curvature
binNumber = np.arange(1,binCount+1,1)
plt.figure(1); g2counts,g2bins,g2bars=plt.hist(rUpperLeaf, bins=binCount, label="Gradient_6nm_v1")
binsMean=[0.5*(g2bins[i]+g2bins[i+1]) for i in range(len(g2counts))]#gives the middle r value of the bin
#^takes average of the wall values of each bin.^
plt.xlabel("Radial Distance / A"); plt.ylabel("Frequency")
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
# plt.savefig("graphs/binned/mixed/G0_Histogram") #includes 3 half plane plots

meanBinHeights,binEdges,binInts=stats.binned_statistic(rUpperLeaf,zUpperLeaf, bins=binCount)#gives the mean or median
plt.figure(2)
plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graphs/binned/mixed/G1_Binned_Mean_Heights_against_Binned_Positions")
binsMean, meanBinHeights = zip(*sorted(zip(binsMean, meanBinHeights)))

#Triangle -> Circle Fitting (3 points)
curvatureTri=np.array([])
#loop over every data point (starting from element 1(2nd))
for counter, _ in enumerate(binsMean):
	if counter >= 1 and counter < len(binsMean): #can't calculate for data points on the edge
		#first calculate triangle side lengths
		aSide=math.sqrt((binsMean[counter-1]-binsMean[counter-2])**2+(meanBinHeights[counter-1]-meanBinHeights[counter-2])**2)
		bSide=math.sqrt((binsMean[counter]-binsMean[counter-1])**2+(meanBinHeights[counter]-meanBinHeights[counter-1])**2)
		cSide=math.sqrt((binsMean[counter]-binsMean[counter-2])**2+(meanBinHeights[counter]-meanBinHeights[counter-2])**2)
		#calc curvature
		curvatureTri=np.append(curvatureTri,circum_circle_radius(aSide,bSide,cSide))

#ODR fitting (5 points)
circleStepSize=3
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
plt.figure(3); 
plt.subplot(211); plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
plt.title("Curvature of Lipid Bilayer (Custom v0)");
plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("Curvature / A-1")
plt.scatter(binsMean[1:], 1/curvatureTri)#boundaries?
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.tight_layout();
plt.savefig("graphs/binned/mixed/G2_Tri_Curvature_against_Binned_Positions_AND_G1")

# plt.figure(4); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Radius of Curvature / A")
# plt.scatter(odrRBins, curvatureODR)
plt.figure(4); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Curvature / A-1")
plt.scatter(odrRBins, 1/curvatureODR)
# plt.xlim((50,150));
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graphs/binned/mixed/G3_ODR_Curvature_against_Binned_Positions")

plt.figure(5);
plt.subplot(211); plt.xlabel("Radial Distance / A"); plt.ylabel("Mean Bin Heights")
plt.scatter(binsMean,meanBinHeights)
#plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Radius of Curvature / A")
#plt.scatter(odrRBins, curvatureODR)
# plt.xlim((50,150));
plt.subplot(212); plt.xlabel("Radial Distance / A"); plt.ylabel("ODR Curvature / A-1")
plt.scatter(odrRBins, 1/curvatureODR)
# plt.xlim((50,150));
plt.ylim((0,0.1));
plt.tight_layout();
plt.title("Curvature of Lipid Bilayer (Custom v0)"); plt.show()
plt.savefig("graphs/binned/mixed/G4_ODR_Curvature_against_Binned_Positions_AND_G1")

curvatureTextODR="The max curvature from ODR Fitting is "+ str(max(1/curvatureODR))
curvatureTextTRI="The max curvature from Triangle -> Circle Fitting is "+ str(np.amax(1/curvatureTri))
print curvatureTextODR
print curvatureTextTRI

#Saving data to text file
f=open("graphs/binned/mixed/curvature.txt","w+")
f.write(curvatureTextODR+"\n"+curvatureTextTRI)
f.close()