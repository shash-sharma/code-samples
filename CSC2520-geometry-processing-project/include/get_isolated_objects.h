#ifndef GET_ISOLATED_OBJECTS_H
#define GET_ISOLATED_OBJECTS_H

#include <iostream>
#include <vector>

#include <igl/squared_edge_lengths.h>
#include <igl/components.h>
#include <igl/adjacency_matrix.h>

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
void get_isolated_objects(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  std::vector<Object> & Objects,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & counts);
#endif
