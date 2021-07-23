#ifndef SAMPLE_OBJ_VERT_H
#define SAMPLE_OBJ_VERT_H

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
void sample_obj_vertices(
    int N_samples,
    Object & obj,
    Eigen::MatrixXd & NV);
#endif
