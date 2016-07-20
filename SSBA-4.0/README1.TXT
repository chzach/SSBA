Description
-----------

This is an implementation of a sparse Levenberg-Marquardt optimization
procedure and several bundle adjustment modules based on it. There are several
versions of bundle adjustment: probably the most interesting one are the
bundle_large*.cpp files in the Apps folder, which are used in my paper,
"Robust bundle adjustment revisited", ECCV 2014.

* bundle_large.cpp: contains the implementation for the IRLS, root psi, and
  naive lifting method (without the Schur complement trick). Which method to
  use is adjustable via boolean variables near the beginning of the file.

* bundle_large_triggs.cpp: uses the Triggs correction.

* bundle_large_lifted_schur.cpp: employs the lifted formulation for robust
  kernels made efficient via the Schur complement.

All programs expect a dataset file name as command line argument. The format
is the one used for the ``bundle adjustemnt in the large'' publication
(Agarwal et al., ECCV 2010). For convenience, one problem instance is found as
problem-257-65132-pre.txt in the Dataset folder.

Build and test instructions on a Unixoid box with Suitesparse/colamd installed
in a standard directory:

> mkdir build && cd build && cmake .. && make
> cd Apps
> ./bundle_large ../../Dataset/problem-257-65132-pre.txt
> ./bundle_large_triggs ../../Dataset/problem-257-65132-pre.txt
> ./bundle_large_lifted_schur ../../Dataset/problem-257-65132-pre.txt

You can adjust some parameters (metric or full BA, robust kernel, inlier
threshold) in Apps/bundle_large_common.h.

If you find this code useful, please refer to

@incollection{zach2014robust,
  title={Robust Bundle Adjustment Revisited},
  author={Zach, Christopher},
  booktitle={Computer Vision--ECCV 2014},
  pages={772--787},
  year={2014},
  publisher={Springer}
}

The preprint (with minor corrections) can be found in the Paper folder.

-Christopher Zach (christopher.m.zach@gmail.com)


Remainder of the description for the other BA applications
----------------------------------------------------------

1) Pure metric adjustment. Camera poses have 6 dof and 3D points have 3 dof.
2) Common, but adjustable intrinsic and distortion parameters. This is useful,
   if the set of images are taken with the same camera under constant zoom
   settings.
3) Variable intrinsics and distortion parameters for each view. This addresses
   the "community photo collection" setting, where each image is captured with
   a different camera and/or with varying zoom setting.

There are two demo applications in the Apps directory, bundle_common and
bundle_varying, which correspond to item 2) and 3) above.

The input data file for both applications is a text file with the following
numerical values:

First, the number of 3D points, views and 2D measurements:
<M> <N> <K>
Then, the values of the intrinsic matrix
    [ fx  skew  cx ]
K = [  0   fy   cy ]
    [  0    0    1 ],
and the distortion parameters according to the convention of the Bouget
toolbox:

  <fx> <skew> <cx> <fy> <cy> <k1> <k2> <p1> <p2>

For the bundle_varying application this is given <N> times, one for each
camera/view.
Then the <M> 3D point positions are given:

  <point-id> <X> <Y> <Z>

Note: the point-ids need not to be exactly from 0 to M-1, any (unique) ids
will do.
The camera poses are given subsequently:

  <view-id> <12 entries of the RT matrix>

There is a lot of confusion how to specify the orientation of cameras. We use
projection matrix notation, i.e. P = K [R|T], and a 3D point X in world
coordinates is transformed into the camera coordinate system by XX=R*X+T.

Finally, the <K> 2d image measurements (given in pixels) are provided:

  <view-id> <point-id> <x> <y> 1

See the example in the Dataset folder.


Performance

This software is able to perform successful loop closing for a video sequence
containing 1745 views, 37920 3D points and 627228 image measurements in about
16min on a 2.2 GHz Core 2. The footprint in memory was <700MB.


Requirements

Solving the augmented normal equation in the LM optimizer is done with LDL, a
Cholsky like decomposition method for sparse matrices (see
http://www.cise.ufl.edu/research/sparse/ldl). The appropriate column
reordering is done with COLAMD (see
http://www.cise.ufl.edu/research/sparse/colamd). Both packages are licensed
under the GNU LGPL.

This software was developed under Linux, but should compile equally well on
other operating systems.


News for SSBA 4.0

* Added support for the lifted approach.

News for SSBA 2.0

* Added a sparse LM implementation (struct ExtSparseLevenbergOptimizer)
  handling several least-squares terms in the cost function. This is useful
  when several types of measurements (e.g. image feature locations and GPS
  positions) are available. See Apps/bundle_ext_LM.cpp for a simple demo.

* Changed the default update rule for the damping parameter lambda to a
  simpler one (multiply and divide lambda by 10, depending on the cost
  function improvement). This seems to work better that the more complicated
  rule used before. 

* Fixed a trivial, but important bug in cost evaluation after the parameter
  update.

/*
Copyright (c) 2011-2014 Christopher Zach

This file is part of SSBA-4.0 (Simple Sparse Bundle Adjustment).

SSBA is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

SSBA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along
with SSBA. If not, see <http://www.gnu.org/licenses/>.
*/
