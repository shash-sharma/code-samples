#include <iostream>
#include <vector>

#include <igl/squared_edge_lengths.h>
#include <igl/components.h>
#include <igl/adjacency_matrix.h>

#include <Eigen/Core>

#include "../include/get_N_obj.h"
#include "../include/Object.h"

void get_N_obj(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  int & N_obj,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & counts)
{

  	// First, get adjacency matrix
  	Eigen::SparseMatrix<double> A;
  	igl::adjacency_matrix(F, A);

  	// Now get connected components for vertices
  	igl::components(A, C, counts);

  	N_obj = counts.rows();

	// ------ Debugging ------
  	// for (int ii = 0; ii < C.rows(); ii++)
  	// 	std::cout << C(ii) << std::endl;

  	// std::cout << std::endl;

  	// for (int ii = 0; ii < counts.rows(); ii++)
  	// 	std::cout << counts(ii) << std::endl;
	// -----------------------


  	return;

}
