# Simple Spare Bundle Adjustment

## Description

This is an implementation of a sparse Levenberg-Marquardt optimization
procedure and several bundle adjustment modules based on it. There are several
versions of bundle adjustment: "Robust bundle adjustment revisited", ECCV 2014.

* bundle_large.cpp: contains the implementation for the IRLS, root psi, and naive lifting method (without the Schur complement trick). Which metod to use is adjustable via boolean variables near the beginning of the file.

* bundle_large_triggs.cpp: uses the Triggs correction.

* bundle_large_lifted_schur.cpp: employs the lifted formulation for roust kernels made efficient via the Schur complement.

* bundle_common.cpp: standard bundle implementation.

All programs expect a dataset file name as command line argument. The format is the one used for the ``bundle adjustemnt in the large'' publication (Agarwal et al., ECCV 2010). For convenience:
* One problem instance is found as problem-257-65132-pre.txt in the Dataset folder. This is for the bundle_large.cpp, bundle_large_ligted_triggs.cpp, and bundle_large_lifted_schur.cpp methods.
* One problem instance is found as common_dataset.txt in the Dataset directory, for the bundle_common.cpp method. This is also for the other implementations. 
* A result of the bundle_common.cpp method is in refined.txt file under the Dataset directory. 


## Building 
This program is build using Cmake. The [SSBA-4.0 implementation by chzach](https://github.com/chzach/SSBA), did not include the necessary [Suitesparse, COLAMD](http://faculty.cse.tamu.edu/davis/research.html) libraries. This repo has included the directives for this. To build the Suitesparse, COLAMD libries and the program files, follow: 
```
> mkdir build && cd build && cmake .. && make 
> cd Apps
```

### Troubleshooting
If an error occurs, it is mostlikely because of the CMakeLists.txt file. In order to correct for this, look and make sure the appropriate compiler and directories are being referenced. 

### Citation
I have made it so that the entire git project compiles within the directory. If you find this useful, please cite: 
```
@incollection{zach2014robust,
  title={Robust Bundle Adjustment Revisited},
  author={Zach, Christopher},
  booktitle={Computer Vision--ECCV 2014},
  pages={772--787},
  year={2014},
  publisher={Springer};
  note = {phatty added libraries July 2016}
}
```

## Header Format: 
This pertains to using the `bundle_common.cpp` like programs. The following is how the data should be formated. You can find an anotated example in the Dataset directory. For implementation, please make sure that there no blank lines or commments in your data file.

1. `<numPoints> <numCams> <num2Dmeasurements>`
2. `<focalXdirect> <skew> <cx> <focalYdirect> <cy> <k1> <k2> <p1> <p2>`
3. `<uniquePtId> <xCord> <yCord> <zCord>`
4. `<uniqueViewId> <12RTparameters(Row Major, left to right)>`
5. `<uniqueViewId> <uniquePtId> <camera2DxCord> <camera2DyCord> <1>`

## Performance:
This software is able to perform common_bundle.cpp adjustment on 10269 points, with 295 camera perspectives and 97798 2D projections, under the radial distortion pattern in 10.452 seconds; 2.7 GHz Intel Core i5, 8 GB 1867 MHz DDR3. 


### Updates:

* News for July 2016 update: Included the COLAMD, and Suitsparse library for compiling ease. 

* News for SSBA 4.0 :Added support for the lifted approach.

* News for SSBA 2.0: Added a sparse LM implementation (struct ExtSparseLevenbergOptimizer) handling several least-squares terms in the cost function. This is useful when several types of measurements (e.g. image feature locations and GPS positions) are available. See Apps/bundle_ext_LM.cpp for a simple demo. Changed the default update rule for the damping parameter lambda to a simpler one (multiply and divide lambda by 10, depending on the cost function improvement). This seems to work better that the more complicated rule used before. Fixed a trivial, but important bug in cost evaluation after the parameter update.

```
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
```
