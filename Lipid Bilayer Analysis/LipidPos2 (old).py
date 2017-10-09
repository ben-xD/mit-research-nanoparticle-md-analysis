# -*- coding: utf-8 -*-
"""Created on Tue Jul 18 09:08:20 2017
@author: Adrian Butterworth"""
import numpy as np
import math, os
# import scipy.optimize.curve_fit as curve_fit
from scipy.optimize import curve_fit
import matplotlib as mpl #first comment specifies the encoding. 8-bit Unicode Transformation Format
mpl.use("agg") #make sure the backend is right to save a file to check after running script.
mpl.rcParams["agg.path.chunksize"]=20000 
mpl.rcParams["lines.linewidth"]=1 #rcParams: rc means configuration file or run commands parameters.  
#rcParams is a dictionary like variable
#backend is not case sensitive. Agg can save PNGs. Agg stands for antigrain geometry
import matplotlib.pyplot as plt #the backend is specifically for pyplot.
r, z = np.genfromtxt("heads_r_pos_FullPlane.txt", unpack=True,dtype="float_")
#Need to remove lines, just have scatter points.
#Duplicate, and then merge into Half Plane:
rPos=abs(r) #make all r values positive.
#need to make the average line through middle. SPLIT THEN FIT
#SPLIT
MaxRadius=150 #won't fit to values with radius > 150
selection = np.where(rPos<=150)#gives element numbers of positions smaller than 150
fit=np.polyfit(rPos[selection],z[selection],2) #polynomial fit. Returns an array of coefficients.
def split(x):
	y=fit[0]*x**2+fit[1]*x+fit[2]
	return y
upperLeaf = np.where(z > split(r)) #binary array

#Defining sigmoid function to fit (use initial guess to improve match)
def sigmoid(x,x0,k,a,c): #the sigmoid curve maths. S(x)=1/(1+e^<-x>)
	y = a / (1 + np.exp(-k*(x-x0))) + c #1 doesnt need to be explicitly float
	return y #return y was missing. it would not provide anything back, being a useless NoneType!
	#that missing return brought a warning else where, so was difficult to find out!
#Fitting sigmoid upper leaf
initial_guess = [60,1e-2,100,130] #from Alexis
#for curve fit, need 1 row, many columns. or list
rCurve=np.ndarray.tolist(r[upperLeaf])
zCurve=np.ndarray.tolist(z[upperLeaf])
pOpt, pCov = curve_fit(sigmoid, rCurve, zCurve, p0=initial_guess)
#pOpt is parameters which minimise squared sum residuals.
#pCov is estimated covariance of popt

#PLOT: Full plane height vs radial distance
plt.scatter(r,z,label="Full Plane (FP)")
plt.xlabel("Radial Distance / A")
plt.ylabel("Height")
plt.title("Curvature of Lipid Bilayer (Custom v0)")
plt.legend()
plt.show()
if not os.path.exists("Graphs/"):
    os.makedirs("Graphs/")
plt.savefig("Graphs/Full_Plane_Separate")

split=np.vectorize(split)#to allow running over each element in NDarray
sigmoid=np.vectorize(sigmoid)
#PLOT: half plane scatter (values)
plt.figure()
plt.scatter(rPos,z,label="Half Plane (HP)")
#PLOT: middle curve (parameters):
minR, maxR = 0, max(r) #find max radius, then plot the functions to this value.
points = int(1e4) #no. of points
SplitR = map(lambda x: float(maxR-minR)*x/points, range(points+1))#map returns a list. later convert to NDarray
#Range(x) makes a list of 0 to (x-1). Lambda is an anon. function

splitRnp=np.array(SplitR)
SplitZ = split(splitRnp) #SplitZ = map(split(SplitR), SplitR)
plt.plot(SplitR, SplitZ, label="Midline Split (HP)")
#PLOT: upper leaf curve (parameters): Sigmoid
ZUpper=sigmoid(splitRnp, pOpt[0],pOpt[1],pOpt[2],pOpt[3])
#ZUpper=map(sigmoid(SplitR, pOpt[0],pOpt[1],pOpt[2],pOpt[3]), SplitR)
plt.plot(SplitR,ZUpper,label="Upper Leaf Fit (HP)")

plt.xlabel("Radial Distance / A")
plt.ylabel("Height")
plt.title("Curvature of Lipid Bilayer (Custom v0)")
plt.legend()
plt.show()
plt.savefig("Graphs/Half_Plane_SetSeparate") #includes 3 half plane plots

plt.figure() #refresh, to prevent hiding the first graph by the scatter.
plt.xlabel("Radial Distance / A"); plt.ylabel("Height");plt.title("Curvature of Lipid Bilayer (Custom v0)")
plt.scatter(r,z,label="Full Plane (FP)")
plt.plot(SplitR,ZUpper,label="Upper Leaf Fit (HP)")
plt.plot(SplitR, SplitZ, label="Midline Split (HP)")
plt.legend()
plt.show()
plt.savefig("Graphs/All_Graphs")

#Also need to output the fitting parameters. Otherwise the fit isn't very useful.

#--------------------------------------------------------------------------------------#
#Select region #can make the variables argument to this script later.
crystalWidth=50 #v1=5.0, v2=2.5, v3=1.5, v4=1.0
minSelect=60; maxSelect=minSelect+crystalWidth
zUpper=z[upperLeaf]; rUpper=r[upperLeaf]
selection=np.where((rUpper>=minSelect) & (rUpper<=maxSelect))
zSelect=zUpper[selection]; rSelect=rUpper[selection]
#Fit region. Define fit function, then curve_fit. (Cant just calculate curvature on discrete/ spread out scatter points)
#need to find a better fitting function!
def	sqrtFit(x,K,c):
	return K*x**0.5+c
init=[1,100]
pOpt, pCov=curve_fit(sqrtFit,rSelect,zSelect,p0=init); K=pOpt[0]
def curvature(x):
	return (0.5*K*x**(-0.5))/((1+(0.5*K*x**(-0.5))**2)**1.5)
curvatureArray=curvature(rSelect)
#plot curvature with radial distance. Another subplot for contact point/ presence
plt.figure(); plt.xlabel("Radial Distance / A"); plt.ylabel("Curvature")
plt.scatter(rSelect, curvatureArray)
plt.show(); plt.savefig("Graphs/Curvature_against Radial Distance")
print("Curvatures are " +str(curvatureArray))