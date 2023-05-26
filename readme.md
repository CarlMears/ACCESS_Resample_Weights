# Read Me file for the Fortran Part of the resampling code

The code in the this reposistory calculates the resampling weights for resampling swath-geomtry 
radiometer measurements onto circular footprints at specified locations in the swath geometry.

The inputs are 
* source footprints for a range of scans
* target footprints
These are generated in python

The resampling calculation is broken into two parts.

### First, the "g" matrix is calculated.  
g is the overlap between pairs of source footprints.  Because
many footprints are too far apart to have any overlap, the matrix is fairly sparse.

The calculation is performed by calc_g_everywhere.f90 which builds to calc_g_everywhere.  The run time depends on the size of the 
source footprints.  For the larger footprints, it takes many hours.

Note that the g matrix is independent of the target size.

The output of this part is file containing the g matrix.  Because it is relatively sparse, it is stored as i(int32),j(int32),value(real64) triplets.  The matrix is symmetric, so the i,j orver is not important, but all values are stored (not just the upper or lower triangle).  

### Second, the resampling weights are computed.  
This is done by find_resampling_weights_precomute_g.f90 which build to find_resampling_weights_precompute_g.  The most costly part of this process is calculating 
overlap integrals between the target footprint and the various source footprints included.  Only 
source footprints within some distance of the target location are included.  This is currently 
set to 3X the target footprint size.