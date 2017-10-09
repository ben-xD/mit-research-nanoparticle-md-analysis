# -*- coding: utf-8 -*- #first comment specifies the encoding. 8-bit Unicode Transformation
"""This script fits a circle using both ODR and Sigmoid curve"""
import time,sys,math,random,os, numpy as np
from scipy import optimize, odr
from scipy.optimize import curve_fit, minimize_scalar
import matplotlib as mpl
mpl.use("agg") #make sure the backend is right to save a file to check after running script.
mpl.rcParams["agg.path.chunksize"]=20000 
mpl.rcParams["lines.linewidth"]=1 #rcParams: rc means configuration file or run commands parameters.  
#rcParams is a dictionary like variable
#backend is not case sensitive. Agg can save PNGs. Agg stands for antigrain geometry
import MDAnalysis
from MDAnalysis.coordinates.XTC import XTCReader
import matplotlib.pyplot as plt #the backend is specifically for pyplot.

if not os.path.exists("graphs/frameSnapshotODRSigm"):
	os.makedirs("graphs/frameSnapshotODRSigm")

#*****Importing data****************************************
GRO="../solution_eq.gro"#topology: which atoms, connections, constant properties (partial charge,..)
XTC="../solution_md.xtc"#trajectory: changing coordinate info
u=MDAnalysis.Universe(GRO,XTC)
print u.atoms #Universe.atoms is a AtomGroup and can be thought of as list consisting of Atom objects. 
print "Number of Frames = ", len(u.trajectory) ##frames or time step
#The Atom is the elementary and fundamental object in MDAnalysis.
#*****Selecting atom groups. *******************************
#provided by ASE docs = Atomic Simulation Environment
PO4 = u.select_atoms("name PO4")
heads = u.select_atoms("name PO4 NC3")
u_heads = heads.atoms[0:len(heads)/2] 
l_heads = heads.atoms[len(heads)/2:len(heads)] 
tails = u.select_atoms("name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B") 
crystal = u.select_atoms("resname SI")
cryst_com = crystal.center_of_mass()
print cryst_com
#Creating empty arrays**************************************
count=np.array([]); fNum=np.array([]) 
curv=np.array([]); #Counter, curvature at specific time (appended every loop)
curvSigm_FirstDevMax=np.array([]);
curvSigm_SecondDevMax=np.array([]);
curvSigm_CurvatureDepth=np.array([]);#change data type? dtype='int64'
r=np.array([]); #z=np.array([]) #radial distance, height (refreshed every loop)
frameCalcTotal=100
#Function definitions*********************************************
#Define implicitily, circle function for fitting
def circle(beta, x): #beta are the parameters
	return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2
def calc_R(rc, zc): #Used for circle fitting later
	#distance of 2d points from center rc, zc. Returns a vector
	return np.sqrt((rReady-rc)**2 + (zReady-zc)**2) #c means center.
#Sigmoid, and its derivatives to fit upper leaflet curvature
def sigmoid(x, x0, k, a, c): #defining the f() which goes into curve_fit
     y = a / (1 + np.exp(-k*(x-x0))) + c
     return y
def d_sigmoid(x, x0, k, a, c): 
    y = (a*(k*np.exp(-k*(x-x0)))) / (1 + np.exp(-k*(x-x0)))**2
    return y
def d2_sigmoid(x, x0, k, a, c): 
    y2 = ((1 - np.exp(math.ceil(-k*(x-x0)*100000/100000))) * (np.exp((math.ceil(-k*(x-x0)*100000/100000)))) * a * k**2) / (1 + np.exp((math.ceil(-k*(x-x0)*100000/100000))))**3
    return y2
d_fm = lambda x: -d_sigmoid(x, *pOpt)
d2_fm = lambda x: -d2_sigmoid(x, *pOpt)
#Defining the polynomial function for splitting upper leaf and lower leaf
def split(x):
	y=fit[0]*x**2+fit[1]*x+fit[2]
	return y
#**********************start loop over frames*********************
for ts in u.trajectory: #Gives access to the coordinates over time
	# if ts.frame == 1:
		# print heads.positions
		# print "Just printed heads.positions at frame #1. Check its shape"
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
		#Work on the upper leaf only.
		upperLeaf = np.where(z > split(r)) #binary array. List element# if z is larger than the split function.
		rUpperLeaf=r[upperLeaf]
		zUpperLeaf=z[upperLeaf]
		#Filter start/end points
		crystalWidth=50+25; minSelect=60; #v1=50, v2=25, v3=15, v4=10 angstroms
		maxSelect=minSelect+crystalWidth #min from folder name
		curvedSelection=np.where((rUpperLeaf>=minSelect) & (rUpperLeaf<=maxSelect)) #restricting radial distance region
		rReady=rUpperLeaf[curvedSelection]; zReady=zUpperLeaf[curvedSelection]
		#fit data to circle
		#initial guess for circle fit
		R_m = calc_R(rReady.mean(), zReady.mean()).mean()
		beta0 = [ rReady.mean(), zReady.mean(), R_m]#initial parameters
		#ODR********************
		lsc_data  = odr.Data(np.ma.row_stack([rReady, zReady]), y=1)#what does this line do? The data to fit.
		#row_stack stacks arrays in sequence vertically. 
		#scalar input into odr.Data for y= value indicates implicit model
		lsc_model = odr.Model(circle, implicit=True) #The Model class stores information about the function you wish to fit.
		lsc_odr   = odr.ODR(lsc_data, lsc_model, beta0) #The ODR class gathers all information and coordinates the running of the main fitting routine.
		lsc_out   = lsc_odr.run()
		rc2, zc2, R = lsc_out.beta #estimated parameter values
		Ri		  = calc_R(rc2,zc2)
		#End of ODR*************
		print("The Circle fit (ODR) calculated curvature is "+str(Ri.mean()))
		print(rc2, zc2)
		curv=np.append(curv, Ri.mean()) #curvature is average distance of points from center.
		#***Fit to Sigmoid***********************************
		initial_guess = [60,1e-2,100,130]#can improve initial guesss
		pOpt, pCov = curve_fit(sigmoid, rUpperLeaf, zUpperLeaf, p0=initial_guess)
		
		r_1 = minimize_scalar(d_fm, bounds=(-50, 200))#minimisation of scalar function of one variable
		r_2 = minimize_scalar(d2_fm, bounds=(-50, 200))
		#Saving data to columns
		curvSigm_FirstDevMax=np.append(curvSigm_FirstDevMax, d_sigmoid(r_1["x"], *pOpt))
		curvSigm_SecondDevMax=np.append(curvSigm_SecondDevMax, d2_sigmoid(r_2["x"], *pOpt))
		curvSigm_CurvatureDepth=np.append(curvSigm_CurvatureDepth, pOpt[2])
		print "Sigmoid First Derivative Max = ", d_sigmoid(r_1["x"], *pOpt)
		print "Sigmoid Second Derivative Max= ", d2_sigmoid(r_2["x"], *pOpt)
		print "Sigmoid Curvature Depth =", pOpt[2]
		#plot rReady and zReady to check data and circle center every once in a while
		if ts.frame%10==0:
			plt.figure(); plt.xlabel("Radial Distance / A"); plt.ylabel("Height / unit")
			plt.scatter(rReady,zReady)
			plt.scatter(rc2,zc2,c='r')
			plt.xlim(0, maxSelect+10)
			nameString="graphs/frameSnapshotODRSigm/No%i" % ts.frame
			plt.show(); plt.savefig(nameString)
		#end of loop through each frame
#removing broken values (NaN),	NaN=not a number	
#np.isnan returns true/false array when Nan/Not Nan
removingBroken1=np.logical_not(np.isnan(curvSigm_FirstDevMax))#logical_not: if true, returns false
curvSigm_FirstDevMaxClean = curvSigm_FirstDevMax[removingBroken1]
curvSigm_FirstDevMaxClean=curvSigm_FirstDevMax[removingBroken1]
curvSigm_SecondDevMaxClean=curvSigm_SecondDevMax[removingBroken1]
curvSigm_CurvatureDepthClean=curvSigm_CurvatureDepth[removingBroken1]
print("Removed ",len(removingBroken1)-np.sum(removingBroken1), "NaN data points.")
#removing >1 and 0. Can also use 0.2 instead of 0.0001. 
removingBroken2=np.where((curvSigm_FirstDevMaxClean<1) & (curvSigm_FirstDevMaxClean>0.001))
curvSigm_FirstDevMaxClean=curvSigm_FirstDevMaxClean[removingBroken2]
curvSigm_SecondDevMaxClean=curvSigm_SecondDevMaxClean[removingBroken2]
curvSigm_CurvatureDepthClean=curvSigm_CurvatureDepthClean[removingBroken2]
print("Then, removed ",(np.sum(np.where(removingBroken2==0))), " 'zero' or 'above 1' data points.")
#Printing information
#np.c_ translates slice objects to concatenation along 2nd axis.
print "Circle (ODR) Fit Mean = ",curv.mean(), "Circle Fit Std =" , np.std(curv)
print "Mean Sigmoid First Derivative Max= ", curvSigm_FirstDevMaxClean.mean()
print "Mean Sigmoid Second Derivative Max = ", curvSigm_SecondDevMaxClean.mean()
print "Mean Sigmoid Curvature Depth = ", curvSigm_CurvatureDepthClean.mean()
#*******************mean values as the last row******************
#count=np.append(count, 00)
#fNum=np.append(fNum, fNum.mean())
#curv=np.append(curv, curv.mean())
#curvSigm_FirstDevMax=np.append(curvSigm_FirstDevMax, curvSigm_FirstDevMax.mean())
#curvSigm_SecondDevMax=np.append(curvSigm_SecondDevMax, curvSigm_SecondDevMax.mean())
#curvSigm_CurvatureDepth=np.append(curvSigm_CurvatureDepth, curvSigm_CurvatureDepth.mean())
np.savetxt("Fit_ODR.txt", np.c_[count, fNum, curv,
curvSigm_FirstDevMax, curvSigm_SecondDevMax, curvSigm_CurvatureDepth])