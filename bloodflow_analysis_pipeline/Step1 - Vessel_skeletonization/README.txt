Code for skeletonization of 2-photon Gaussian Z-stacks, and subsequent visualization and calculation of segments of skeleton. 
Written by: Gokul Upadhyayula, Jiang Lan Fan, UC Berkeley

Developed and tested in MATLAB R2019a on Windows 10.
All dependencies are included. 

Skeletonization code based off of:
Kollmannsberger, Kerschnitzki et al., "The small world of osteocytes: connectomics of the lacuno-canalicular network in bone." New Journal of Physics 19:073019, 2017.
T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "Enhancement of Vascular Structures in 3D and 2D Angiographic Images", IEEE Transactions on Medical Imaging, 35(9), p. 2107-2118 (2016)
T. Jerman, F. Pernus, B. Likar, Z. Spiclin, "Blob Enhancement and Visualization for Improved Intracranial Aneurysm Detection", IEEE Transactions on Visualization and Computer Graphics, 22(6), p. 1705-1717 (2016)
with read/write libraries from:
F. Aguet, S. Upadhyayula, R. Gaudin, et al., "Membrane dynamics of dividing cells imaged by lattice light-sheet microscopy", Molecular Biology of the Cell 2016 Nov 7;27(22):3418-3435.

Instructions:

First run 
GU_extractVesselness_Skel('sampledata\kf31_0z_roi2_gauss_1umz_1p3238x.tif')
in Matlab, which will take a sample raw image stack of vessels and transform it into a skeleton. 
Newly created file "sampledata\kf31_0z_roi2_gauss_1umz_1p3238x_skelInterp.tif" is a 1x1x1-micron interpolated skeleton of the raw image stack.

Then, run nonlinear_transform_of_3d_graph_segment.m code.
This transforms skeleton tif stack into a 3D Matlab graph object, which can then be visualized in 3D and used to calculate lengths of graph segments. 
A kymograph of a vessel segment inputted into this code is non-linearly stretched into a new kymograph output based on the 3D profile of the segment.

Replace pointers to sample data to run on other data. Details are commented in the Matlab scripts.


