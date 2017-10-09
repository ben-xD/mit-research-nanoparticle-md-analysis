#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Contacts.py plots the frequency of contacts between lipid tail and water/antifreeze solvent."""
"""Created on Mon Jul 31 10:18: 2017 @author: Adrian Butterworth"""
""" """
#PROBLEMS? make sure matplotlibrc is present in dir

import math, os
import numpy as np, scipy.spatial.distance as scpd
import MDAnalysis as MDA
from MDAnalysis.analysis import contacts
import matplotlib as mpl
import matplotlib.pyplot as plt

#Ensure directory present
if not os.path.exists("graphs_contacts/"):
	os.makedirs("graphs_contacts/")

GRO="../solution_eq_323_1frame.gro"#lists atoms: Starting positions/ velocities
XTC="../solution_md_323_100frames.xtc"#Coordinates (one or more frames)
u = MDA.Universe(GRO,XTC)
print u.atoms #u.atoms is the fundamental object, a AtomGroup: list of Atom objects.

#Check ITP file for tails and solvent names. Keywords are case sensitive. 
#How do you ensure name selected correspond to tails! Check structure
#http://cgmartini.nl/index.php/force-field-parameters/lipids2/351-lipid.html?dir=PC&lipid=DOPC
tails = u.select_atoms("resname DOPC and name C1A D2A C3A C4A C1B D2B C3B C4B")
#DOPC stands for 1,2-Dioleoyl-sn-glycero-3-phosphocholine
#PC stands for phosphatidylcholine (PC): the first phospholipid to be discovered
#W is the molname/ resname for water atom/bead (Solvent)
sol = u.select_atoms("resname W WF") #WF=Antifreeze molname/ resname

#Creating arrays for calculations
xBinCount=10;yBinCount=10;
#USE BASH SCRIPT to select frames for analysis.
# frameCalcTotal=3; #100 frames correspond to 30ns. (dt=0.03ps, 1000frames/step)

#Creating arrays for final output data
#Keep adding bin values to the num_contacts_2D
num_contacts_2D_allframesSum = np.empty((xBinCount,yBinCount)) #Total number of contacts in bins.
print "Total Number of Frames = ", len(u.trajectory)
for ts in u.trajectory:
	print "Calculating ", ts.frame, " out of ", len(u.trajectory)

	#need to bin data each frame again to reselect beads
	tails_xyzpos=tails.positions
	tails_xpos=tails_xyzpos[0];tails_ypos=tails_xyzpos[1];tails_zpos=tails_xyzpos[2];
	sol_xyzpos=sol.positions
	sol_xpos=sol.positions[0];sol_ypos=sol.positions[1];sol_zpos=sol.positions[2];
	
	#Generating arrays for x and y positions every frame.
	x=np.linspace( np.min(tails_xpos), np.max(tails_xpos), xBinCount)#try fixed bin positions to remove noise.
	y=np.linspace( np.min(tails_ypos), np.max(tails_ypos), yBinCount)#This bin method will move with data boundaries.
	
	#loop over each bin in **X-DIRECTION**
	for counter, content in enumerate(x):#first row, loop over each column value
		#Selects atoms in the bin x region
		binSelectionX=np.where(tails_xyzpos[:,0]<x[counter+1] & tails_xyzpos[:,0]>=x[counter])
		#****************Be careful with np.where in 2D arrays.
		tails_xyz_Select=tails_xyzpos(binSelectionX)#***?can you do this for 2D array?
		#loop over each bin in **Y-DIRECTION**
		for counter2, content2 in enumerate(y):#Loop over each row
			#Selects atoms in the y region
			#****************Be careful with np.where in 2D arrays.
			binSelectionY=np.where(tails_xyz_Select[:,1]<y[counter2+1] & 
					tails_xyz_Select[:,1]>=y[counter2]) #this identifies the rows with valid values
			tails_xyz_Select2=tails_xyz_Select[binSelectionY,:]#or switch to :,binS if problems arise
			#need to select water molecules close by only to save time
			solClose=select.atoms("around 6 tails and resname W WF")#around 6 angstroms away from tails AtomGroup
			sol_xyzpos=solClose.positions

			#cdist calculates distances between EVERY 2 pairs (2 2D arrays given).
			distance=scpd.cdist(XA = tails_xyz_Select2, XB = sol_xyzpos)
			cutoff=6 #Angstroms. cutoff to count as IN CONTACT
			closeCount=sum(sum(distance<cutoff))#count number of positions smaller than cutoff
			num_contacts_2D_allframesSum[counter,counter2]=closeCount+num_contacts_2D_allframesSum[counter,counter2]
			#MDAnalysis count number of contacts. Need reference conformation or use radius cut!
				# contactFraction=MDA.analysis.contacts.Contacts(u,(tails,sol),refgroup=(,), method='radius_cut')
				# contactFraction=contacts.q1q2(u
				# contactFraction=contacts.Contacts(u, selection=(),method=is_any_closer,refgroup=())
					
#Prep data for output. All Frames means all frames are included in data.
#Average means divided by number of frames. The data still contains bins for x (1D) or xANDy (2D)
num_contacts_1D_allframesAverage = sum(num_contacts_2D_allframesSum)/frameCalcTotal #Sum 2D -> 1D So it can be plotted later
num_contacts_2D_allframesAverage = num_contacts_2D_allframesSum/frameCalcTotal #divide (average over all frames, but keep bins separate)

#Output
plt.figure(); plt.title("Number of Contacts between Solvent and Lipid (1D)")
plt.plot(x, num_contacts_1D_allframesAverage)
plt.xlabel("x Position (summed y values) / A");plt.ylabel("Frequency of Contacts")
plt.savefig("/graph_contacts/G1_1D_Solvent_Lipid_Contact_Count");

#Could do a heatmap for the 2D array.
plt.figure(); plt.xlabel("x Position / A");plt.ylabel("y Position / A")
plt.imshow(num_contacts_2D_allframesAverage);
plt.title("2D Contact Count (Heat Map)");
plt.savefig("/graph_contacts/G1_2D_Contact_Count_Heat_Map");