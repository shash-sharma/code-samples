#include <iostream>
#include <vector>

#include <igl/squared_edge_lengths.h>
#include <igl/components.h>
#include <igl/adjacency_matrix.h>

#include <Eigen/Core>

#include "../include/get_isolated_objects.h"
#include "../include/Object.h"

void get_isolated_objects(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  std::vector<Object> & Objects,
  Eigen::MatrixXd & C,
  Eigen::MatrixXd & counts)
{

  	int N_obj = Objects.size();
  	int N_vert = V.rows();
  	int N_face = F.rows();

  	// Temporary lists
  	std::vector<std::vector<int>> tempV (N_obj), tempF (N_obj);

  	// Populate the temporary list of vertices per Object
  	for (int ii = 0; ii < N_vert; ii++)
  		tempV[C(ii)].push_back(ii);

	// ------ Debugging ------
  	// for (int ii = 0; ii < N_obj; ii++)
  	// 	std::cout << tempV[ii].size() << std::endl;
	// -----------------------

  	// Populate the temporary list of faces per Object. Go through each face and figure out which Object it must belong to
  	for (int ii = 0; ii < N_face; ii++)
  	{
  		for (int jj = 0; jj < N_obj; jj++)
  		{
  			if (std::find(tempV[jj].begin(), tempV[jj].end(), F(ii, 0)) != tempV[jj].end() ||
  				std::find(tempV[jj].begin(), tempV[jj].end(), F(ii, 1)) != tempV[jj].end() ||
  				std::find(tempV[jj].begin(), tempV[jj].end(), F(ii, 2)) != tempV[jj].end())
			{
				tempF[jj].push_back(ii);
				break;
			}
  		}
  	}

	// ------ Debugging ------
  	// for (int ii = 0; ii < N_obj; ii++)
  	// 	std::cout << tempF[ii].size() << std::endl;
	// -----------------------

  	// Map the vertices into the vector of Objects in Eigen-friendly format
  	for (int ii = 0; ii < N_obj; ii++)
  	{
  		Eigen::Map<Eigen::VectorXi> V_temp(&(tempV[ii][0]), tempV[ii].size());
  		Objects[ii].V_ind = V_temp;

  		Objects[ii].V.resize(tempV[ii].size(), 3);
  		for (int jj = 0; jj < tempV[ii].size(); jj++)
  			Objects[ii].V.row(jj) = V.row(tempV[ii][jj]);
  	}

	// ------ Debugging ------
  	// for (int ii = 0; ii < N_obj; ii++)
  	// 	std::cout << Objects[ii].V.rows() << ", " << Objects[ii].V.cols() << std::endl;
  	// for (int ii = 0; ii < Objects[0].F.rows(); ii++)
  	// 	std::cout << Objects[0].F(ii, 0) << ", " << Objects[0].F(ii, 1) << ", " << Objects[0].F(ii, 2) << std::endl;
	// -----------------------

  	// Map the faces into the vector of Objects in Eigen-friendly format
  	for (int ii = 0; ii < N_obj; ii++)
  	{
  		Eigen::Map<Eigen::VectorXi> F_temp(&(tempF[ii][0]), tempF[ii].size());
  		Objects[ii].F_ind = F_temp;

  		Objects[ii].F.resize(tempF[ii].size(), 3);
  		for (int jj = 0; jj < tempF[ii].size(); jj++)
  			Objects[ii].F.row(jj) = F.row(tempF[ii][jj]);
  	}

	// ------ Debugging ------
  	// for (int ii = 0; ii < N_obj; ii++)
  	// 	std::cout << Objects[ii].F.rows() << ", " << Objects[ii].F.cols() << std::endl;
  	// for (int ii = 0; ii < Objects[0].F.rows(); ii++)
  	// 	std::cout << Objects[0].F(ii, 0) << ", " << Objects[0].F(ii, 1) << ", " << Objects[0].F(ii, 2) << std::endl;

	// -----------------------


	// ------ Debugging ------
	// int Nu = 0;
	// for (int ii = 0; ii < N_obj; ii++)
	// {
	// 	Nu += Objects[ii].V.rows() + Objects[ii].F.rows();
	// 	std::cout << Objects[ii].V.rows() << ", " << Objects[ii].F.rows() << std::endl;
	// }
	// std::cout << Nu << std::endl;
	// -----------------------



  	return;

}
