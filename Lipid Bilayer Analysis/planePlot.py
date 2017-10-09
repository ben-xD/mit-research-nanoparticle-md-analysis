#!/usr/bin/env python

"""Learnt and created this code from scratch, with Alexis code as example."""

import time,sys,math,random,os, numpy as np, matplotlib.pyplot as plt
import MDAnalysis
from scipy.optimize import curve_fit, minimize_scalar, fmin
import sympy as sympy
from MDAnalysis.coordinates.XTC import XTCReader

GRO = '../solution_eq.gro' #GRO file is the topology: 
XTC = '../solution_md.xtc' #XTC file is the trajectory: 
u = MDAnalysis.Universe(GRO,XTC) #Reading data into python.

#Select atom groups 
PO4 = u.select_atoms("name PO4")
heads = u.select_atoms("name PO4 NC3")
u_heads = heads.atoms[0:len(heads)/2] 
l_heads = heads.atoms[len(heads)/2:len(heads)] 
tails = u.select_atoms("name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B") 
crystal = u.select_atoms("resname SI")
cryst_com = crystal.center_of_mass()
print cryst_com

print "Number of Frames is" + str(len(u.trajectory))

#pos_phos_traj = np.zeros((len(u.trajectory),total_phos, 3))
#print len(pos_phos_traj[0])

#(r,z) Coords for bilayer heads
r1 = np.array([])
r1F = np.array([]) #array for FULL PLANE plot data.
z1 = np.array([])
#(r,z) Coords for crystal 
r2 = np.array([])
z2 = np.array([])
#(r,z) Coords for lipid tails 
r3 = np.array([])
z3 = np.array([])

#----------------------------------------
#RADIAL COORDS
#Convert x,y,z to r,z (Cartesian to Radial distance in x,y only) to average curvature --> 
#Reference point is COM of crystal, Want to plot z vs r
for ts in u.trajectory:
	print ts.frame
	if ts.frame>200:
		print ts.frame, "again"
		heads_pos = heads.positions
		heads_pos_nor = heads_pos[(heads_pos[:,0]>85) & (heads_pos[:,0]<335)]
#	cryst_pos = crystal.positions
#	tails_pos = tails.positions
#	tails_com = tails.center_of_mass()
	#Lipid Heads: Initial 1/2 plane (condensed) representation.
		# r1 = np.append(r1, np.sqrt((heads_pos_nor[:,0] - cryst_com[0])**2 + (heads_pos_nor[:,1] - cryst_com[1])**2))
		# z1 = np.append(z1, heads_pos_nor[:,2]) #this is the interesting value
	
		#negative x values are left, positive x values are right. 
		#If positive, then the element position is returned. All returned positions are in a list.
		heads_larger=np.where((heads_pos_nor[:,0] - cryst_com[0]) >= 0)
		heads_smaller=np.where((heads_pos_nor[:,0] - cryst_com[0]) <= 0) #determines elements for left side of plane plot.
		#Code above switched from < to =<, to double count data at r=0. But may not do anything, since r is continuous, 0 is discrete.
		xM=heads_pos_nor[:,0] #separating column of x values into separate array. 
		yM=heads_pos_nor[:,1] #All rows, first (0) column.
		zM=heads_pos_nor[:,2] #0 = x, 1 = y, 2 = z. z is height.
		xL=xM[(heads_larger)] #getting all contents that are true/ larger
		yL=yM[(heads_larger)]#selecting the corresponding y values, to pair with x.
		zL=zM[(heads_larger)]
		xS=xM[(heads_smaller)] #separates x into its own array
		yS=yM[(heads_smaller)] #same for y, and z
		zS=zM[(heads_smaller)]
		#r is radial distance, being calculated. Initially r1 is empty.
		#Two lines below, are adding the radial distances of positive x positions to the array.
		r1 = np.append(r1, np.sqrt((xL - cryst_com[0])**2 + (yL - cryst_com[1])**2))
		r1F = np.append(r1F, np.sqrt((xL - cryst_com[0])**2 + (yL - cryst_com[1])**2))#FULL PLANE
		#Now appending radial distances of negative positioned atoms.
		r1 = np.append(r1, np.sqrt((xS - cryst_com[0])**2 + (yS - cryst_com[1])**2)) #add on the negative values to r1 and z1.
		r1F = np.append(r1F, -np.sqrt((xS - cryst_com[0])**2 + (yS - cryst_com[1])**2)) #keep them positive, for fitting!
		
		z1= np.append(z1, zL)#saving the z of all atoms in positive region first.
		z1= np.append(z1, zS)#in the same list, saving z of atoms in negative region.

	np.savetxt('heads_r_pos_FullPlane.txt', np.c_[r1F, z1])
	np.savetxt('heads_r_pos.txt', np.c_[r1, z1])