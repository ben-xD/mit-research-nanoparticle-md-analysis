#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""PLOTTING the data.:"""
#Load each calculated frameCount. Numpy: Add them, divide by total number of frames
#Plot pyplot heatmap of number of contacts in 2D (Check Tail_Solvent_Contact.py code)
#Bar chart for compressed into 1D.

"""Takes x,y,CountFrame1,CountFrame2... and plots it."""
"""Created on Mon Jul 31 10:18: 2017 @author: Adrian Butterworth"""
"""First take the average for each bin. Then can pcolormesh for heatmap."""
#PROBLEMS? make sure matplotlibrc is present in dir

import math, os, glob #glob is a unix style pathname pattern expansion
import numpy as np, scipy.spatial.distance as scpd
import MDAnalysis as MDA
from MDAnalysis.analysis import contacts
import matplotlib as mpl
import matplotlib.pyplot as plt

#Ensure directory present
if not os.path.exists("graphs/contacts/"):
	os.makedirs("graphs/contacts/")
	
#Load text and average them over all frames.
files_loaded = 0
total_list = np.array([0]) #let the broadcasting happen
all_files = glob.glob('framesCount/countOnlyFrame*.txt')
for frame_file in all_files:
    print(frame_file)
    total_list = total_list + np.genfromtxt(frame_file) #max_rows = 98
    files_loaded += 1

#remove 0s from total_list
# remove = np.where(total_list == 0)
# total_list = np.delete(total_list, remove)
# print(total_list, len(total_list), len(total_list)**0.5)
        
average_list = total_list / files_loaded
# average_list = np.append(average_list, [0,0])
print(len(average_list))
dim = int(len(average_list)**0.5)
average_2d_list = np.reshape(average_list, (dim, dim))
sum_x = sum(average_2d_list)

#Output
x =range(1,40,4)
plt.figure(1); plt.title("Number of Contacts between Solvent and Lipid (1D)")
plt.plot(x, sum_x)
plt.xlabel("x Position (summed y values) / A");plt.ylabel("Frequency of Contacts")
plt.savefig("graphs/contacts/G1.png");

#Could do a heatmap for the 2D array.
plt.figure(2); plt.xlabel("x Position / A");plt.ylabel("y Position / A")
plt.imshow(average_2d_list);
# pcolormesh();
plt.title("2D Contact Count (Heat Map)");
plt.savefig("graphs/contacts/G2.png");