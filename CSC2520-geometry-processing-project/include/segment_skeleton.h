#ifndef SEGMENT_SKELETON_H
#define SEGMENT_SKELETON_H

#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "../include/Object.h"

// Compute the discrete angle defect at each vertex of a triangle mesh
// (`V`,`F`), that is, the _locally integrated_ discrete Gaussian
// curvature.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   D  #V list of angle defects in units radians
void segment_skeleton(
    Object & obj);
#endif
