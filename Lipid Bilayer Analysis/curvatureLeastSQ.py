# -*- coding: utf-8 -*- #first comment specifies the encoding. 8-bit Unicode Transformation
"""This code calculates curvature for one region in the data. Using Least Squares method.
It is an edited version of planePlot, which generates and saves the data. This script generates but just uses the data."""
import time,sys,math,random,os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import optimize, odr
from scipy.optimize import curve_fit
import MDAnalysis
from MDAnalysis.coordinates.XTC import XTCReader
#having PROBLEMS? make sure you have a copy of
#matplotlibrc file in the same folder.

#*****Importing data
GRO="../solution_eq.gro"#topology: which atoms, connections, constant properties (partial charge,..)
XTC="../solution_md.xtc"#trajectory: changing coordinate info
u=MDAnalysis.Universe(GRO,XTC)
print u.atoms #Universe.atoms is a AtomGroup and can be thought of as list consisting of Atom objects. 
print "Number of Frames = ", len(u.trajectory) ##frames or time step
#The Atom is the elementary and fundamental object in MDAnalysis.

#*****Selecting atom groups. Details provided by ASE documentation******
#ASE = Atomic Simulation Environment
PO4 = u.select_atoms("name PO4")
heads = u.select_atoms("name PO4 NC3")
u_heads = heads.atoms[0:len(heads)/2] 
l_heads = heads.atoms[len(heads)/2:len(heads)] 
tails = u.select_atoms("name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B") 
crystal = u.select_atoms("resname SI")
cryst_com = crystal.center_of_mass()
print cryst_com
#********************************************************************
#Creating empty arrays
count=np.array([]); fNum=np.array([]); curv=np.array([]) #Counter, curvature at specific time (appended every loop)
r=np.array([]); #z=np.array([]) #radial distance, height (refreshed every loop)
frameCalcTotal=100

def calc_R(rc, zc): #Used for circle fitting later
	#distance of 2d points from center rc, zc. Returns a vector
	return np.sqrt((rReady-rc)**2 + (zReady-zc)**2) #c means center.
def f(c):
	#calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc)
	Ri = calc_R(*c) #what is the meaning of this. *c can be a set of arguments
	return Ri - Ri.mean()
	
#start loop over frames
for ts in u.trajectory: #Gives access to the coordinates over time
	if ts.frame == 1:
		print heads.positions
		print "Just printed (heads.positions) at frame #1. Check its shape"
	if ts.frame > (len(u.trajectory)-frameCalcTotal): 
	#can change 100 to other values (#frames to analyze)
		print "Frame Number:", ts.frame,"Progress: ", ts.frame-(len(u.trajectory)-frameCalcTotal), "/", frameCalcTotal
		count=np.append(count, ts.frame-frameCalcTotal) #column for calculation number
		fNum=np.append(fNum, ts.frame) #column for frame number
		#^Selecting inner region only, other edges may have unrealistic results^
		heads_pos = heads.positions #gets position of all head atoms
		#cutoff region. Do you need to do this now? later will cut off again
		heads_pos = heads_pos[(heads_pos[:,0]>85) & (heads_pos[:,0]<335)]
		#Code above switched from < to =<, to double count data at r=0. But may not do anything, since r is continuous, 0 is discrete.
		#Separating the columns into individual arrays.
		z=heads_pos[:,2]
		#Calculating radial distance
		r=np.sqrt((heads_pos[:,0] - cryst_com[0])**2 + (heads_pos[:,1] - cryst_com[1])**2) #[0] is x coordinate of COM. [2] is y.
		#Split data. 
		maxRadius=150 #won't fit to values with radius > maxRadius
		splitSelect = np.where(r<=maxRadius)#gives element numbers of positions smaller than 150
		fit=np.polyfit(r[splitSelect],z[splitSelect],2) #polynomial fit. Returns an array of coefficients.
		#Defining the polynomial function for plotting/ calculations.
		def split(x):
			y=fit[0]*x**2+fit[1]*x+fit[2]
			return y
		#Work on the upper leaf only.
		upperLeaf = np.where(z > split(r)) #binary array. List element# if z is larger than the split function.
		rUpperLeaf=r[upperLeaf]
		zUpperLeaf=z[upperLeaf]
		#Filter start/end points
		crystalWidth=50-25; minSelect=60+25; #v1=50, v2=25, v3=15, v4=10 angstroms
		maxSelect=minSelect+crystalWidth #min from folder name
		curvedSelection=np.where((rUpperLeaf>=minSelect) & (rUpperLeaf<=maxSelect)) #restricting radial distance region
		rReady=rUpperLeaf[curvedSelection]; zReady=zUpperLeaf[curvedSelection]
		#fit data to circle
		center_estimate = rReady.mean(), zReady.mean() #estimate initial
		center, ier = optimize.leastsq(f, center_estimate) #send estimate, and function in. Sum of squares of equations minimised.
		rC_2, zC_2 = center #solutions
		Ri = calc_R(rC_2, zC_2)
		print("The curvature is "+str(Ri.mean()))
		print(rC_2, zC_2, center)
		curv=np.append(curv, Ri.mean()) #curvature is average distance of points from center.
		#end of loop through each frame
np.savetxt('Fit_LeastSQ.txt', np.c_[count, fNum, curv]) #np.c_ translates slice objects to concatenation along 2nd axis.
print "Mean= ",curv.mean(), "Std=" , np.std(curv)