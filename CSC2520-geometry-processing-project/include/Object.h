#ifndef OBJECT_H
#define OBJECT_H

#include <vector>

#include <Eigen/Core>


class Object
{
  public:

	Eigen::VectorXi V_ind, F_ind, V_ind_samples, V_seg_ind, V_tagged;
	Eigen::MatrixXd V, V_samples, NV_lateral, V_seg_points;
	Eigen::MatrixXi F;
	std::vector<int> V_samples_sort_ind;
    

    // Constructors
    Object() {};
	
};

#endif




