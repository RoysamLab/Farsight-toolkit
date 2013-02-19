MDL Skeletonization Code Instructions:

Xiaosong Yuan
11/12/2008

1. Prerequisites
Environment: MS Visual Studio
Libraries: Boost Graph Library (BGL) (Graph data structure / operations)
        LEDA ver 3.2.3 or above (A Library of Efficient Data Types and Algorithms)

Boost Graph Library can be downloaded:
   http://www.boost.org/doc/libs/1_37_0/libs/graph/doc/index.html
LEDA library can be downloaded: (with new version of BGL, LEDA library is not necessary)
   http://www.algorithmic-solutions.com/leda/index.htm



2. Code Components
(1) AnisoDiffuse: Image smoothing algorithm ¨C Anisotropic Diffusion
(2) ConnComntwFldfill: Connected component analysis algorithm with flood filling method (to prevent from memory outflow problem)
(3) ConnectComponents: Connected component analysis algorithm with scanline filling method
(4) Floodfill: Image floodfill algorithm
(5) Gen_VectorField: Vector field generation algorithm
(6) GradientVecField: Computing gradient vector field from intensity images
(7) MinSpanTree: MDL algorithm for graph (tree) generation
(8) Skel_Extrvalley3D: Computing high curvature seed points
(9) Skel_streamline: Producing 3D skeleton points from vector field and seed points
(10) volumeProcess: Image preprocessing algorithms (including image operations, crop, scale, etc.)


3. Code Description

ConnCompntwFldfill
Compute Connected Component with Flood filling method DepthFirstSearch may cause stack overflow for large datasets

ConnectComponents
Label the connected components of the input volume with zero background and remove the connected components with small number of voxels.
Input: original volume
Output: removed small objects

Floodfill
Flood filling accept a sequence of volumes. Input data contains either 0 or values greater or equal to 3.

Gen_VectorField
Compute the vector field of any volume objects
Input : Binary 3D volume with sizes.
Output: ASCII file with vector 3 components for all object voxels

GradientVecField
Compute the gradient vector field of any 3D density map
Input : 3D volume density map with any sizes.
Output: ASCII file with vector 3 components for all object voxels

MinSpanTree
Generate graph structure from a list of 3D points
Input format: .skel  (3D points)
Output format: .vtk  (3D graph format)

Skel_Extrvalley3D
Extract the ridges and valleys feature in 3D vector field
Input: 3D vector field
Output: 3D image-scalar field

Skel_streamline
Form streamlines from saddle points and seed points
Input: 1. 3D vector field
         2. Seed points
Output: skelelton file

VolumeProcess
Volume dataset processing accepts a sequence of volumes
Input parameters
      1. sizeExpand   
      2. preproess          


4. Code Usage - Command Scripts
(Use MBFsp5 as example)
# 1. set up dataset name
# 2. set up dim of dataset

cd volumeProcess\Release
volumeProcess  C:\Xiaosong\Data\MBFsp5\MBFsp5.308x512x49.raw  308 512 49  1  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk.308x512x49.raw  2
cd ..\..


cd ConnCompntwFldfill\Release
ConnCompntwFldfill  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk.308x512x49.raw  308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.raw  100
cd ..\..


cd AnisoDiffuse\Release
AnisoDiffuse  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.raw  308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.Aniso_k800d02t2.raw  0
cd ..\..


cd GradientVecField\Release
GradientVecField  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.Aniso_k800d02t2.raw  308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec

 
#GradientVecField  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.raw  308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.grdt.vec
cd ..\..


cd skel_Extrvalley3D\Release
skel_Extrvalley3D  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec  308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec_curv0.2.seed
cd ..\..

# Prepares
# 1. Need to set the number of seeds.
# 2. Need to run matlab to generate the vesselness map.

cd skel_streamline\Debug
skel_streamline  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec   308 512 49  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec_curv0.seed  C:\Xiaosong\Data\MBFsp5\MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt_div-1crt.05.step4k_sd.cv0.D4.skel
cd ..\..


cd MinSpanTree\Debug
#MinSpanTree C:\Xiaosong\Data\MBFsp5\  MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt_div-1crt.05.step4k_sd.cv0.D4.skel  MBFsp5_thrs2_dilatk_rmCC100.308x512x49.raw  308 512 49   C:\Xiaosong\Data\MBFsp5\out.vtk  out.txt  C:\Xiaosong\Data\MBFsp5\MBFsp5_thrs2_dilatk_rmCC100.308x512x49.Aniso_k800d02t2.raw
cd ..\..

# 3. Record number of critical points.
# 4. Set up the number of nodes in vtk files.



5. Intermediate Result Formats:
- filename.raw: 3D volume dataset with dimension indicated by the filename (x*y*z), with no header.
- filename.seed: text file with High curvature seed points, which contains a list of 3D points.
- filename.skel: text file with skeleton points, which contains a list of 3D points.
- filename.vtk: common vtk format, which contains MDL skeleton output.