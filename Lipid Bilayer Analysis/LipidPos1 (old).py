# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as mpl #first comment specifies the encoding. 8-bit Unicode Transformation Format
mpl.use("agg") #make sure the backend is right to save a file to check after running script.
mpl.rcParams["agg.path.chunksize"]=20000
mpl.rcParams["lines.linewidth"]=1
#backend is not case sensitive. Agg can save PNGs. Agg stands for antigrain geometry
import matplotlib.pyplot as plt #the backend is specifically for pyplot.
"""
Created on Tue Jul 18 09:08:20 2017

@author: Adrian
"""

r, z = np.genfromtxt("heads_r_pos_FullPlane.txt", unpack=True,dtype="float_")
#Need to remove lines, just have scatter points.


plt.plot(r,z,",b",label="Positions")
plt.xlabel("Radius from Center")
plt.ylabel("Height")
plt.title("Curvature of Lipid Bilayer")
plt.legend()
plt.show()
plt.savefig("Full Plane Split")

