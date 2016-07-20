#Description

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


