#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert gro file to txt for data analysis. Or just use proteinXYZ data..."""
#@author: Adrian Butterworth
import time,sys,math,random,os, numpy as np, matplotlib.pyplot as plt
import MDAnalysis
from scipy.optimize import curve_fit, minimize_scalar, fmin
import sympy as sympy
from MDAnalysis.coordinates.XTC import XTCReader

GRO = '../solution_eq_323.gro'#need to rename to protein GRO & XTC files
XTC = '../solution_md_323.xtc'
u = MDAnalysis.Universe(GRO,XTC) #Load data

#Select atom groups #need to reselect atom groups: need proteins, and separate from lipids
PO4 = u.select_atoms("name PO4")
heads = u.select_atoms("name PO4 NC3")
u_heads = heads.atoms[0:len(heads)/2] 
l_heads = heads.atoms[len(heads)/2:len(heads)] 
tails = u.select_atoms("name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B") 
crystal = u.select_atoms("resname SI")
cryst_com = crystal.center_of_mass()
# print cryst_com

#creating arrays
x=np.array([])
y=np.array([])
z=np.array([])
checkedFramesCount=100;
for ts in u.trajectory:
	print ts.frame
	if ts.frame>(len(u.trajectory)-checkedFramesCount):#Run data saving only if data is in last *100 frames
		heads_pos = heads.positions
		heads_pos_nor = heads_pos[(heads_pos[:,0]>85) & (heads_pos[:,0]<335)]#cutoff
		x=np.append(x, heads_pos_nor[:,0])
		y=np.append(y, heads_pos_nor[:,1])
		z=np.append(z, heads_pos_nor[:,2]) 
		
# xyz=np.stack((x,y,z),axis=1) #3 columns (x, y and z.) (every row is a 3D data point)
np.savetxt('xyzData.txt', np.c_[x,y,z]) 
		
