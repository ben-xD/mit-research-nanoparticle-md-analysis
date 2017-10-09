Lipid Documents:
curvature  (Alexis) & curvature_v2: Extracts GRO data into upper leaflet, fits sigmoid, calculates curvature
analysis.sh (Alexis): 
distance.sh (Alexis): Creates gromacs index, runs gmx mindist, saves minimum distance

curvaturePlane (Alexis): Plots full plane lipid head heights vs positions. (negative positions possible)
curvaturePlane3: Least square fitting of curvature on every frames
curvaturePlane3odr_sigm: Fits the curve every frame.

curvatureBin_Tri: Binning lipid bilayer to allow more consistent triangle-fitting:curvature calculations.
curvatureBin_odr: Binning lipid bilayer to allow more consistent odr-fitting:curvature calculations.
curvatureBin_Tri_and_odr: Both of the above in one
curvatureBin_MultiTri: Fitting multiple triangles per data point
curvatureBin_LeastSQ: Not as good as other script.

curvatureInternet: Calculating curvature at every point by Splines

My progress:
1. Reverse engineering curvaturePlane to create workable data files
'heads_r_pos.txt' & 'heads_r_pos_FullPlane.txt'

2. LipidPos(123) to plot various graphs.

3. Curvature calculations without binning. (LEASTSQ) & (ODR and Sigmoid fitting)
