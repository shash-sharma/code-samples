#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <Eigen/Core>

#include "../include/sample_obj_vertices.h"
#include "../include/Object.h"

void sample_obj_vertices(
    int N_samples,
    Object & obj,
    Eigen::MatrixXd & NV)
{

    // Get a random sampling of Object vertices
    srand(1);
    int min = 0, max = obj.V.rows();
    std::vector<int> samples_vec;
    for (int ii = 0; ii < N_samples; ii++)
        samples_vec.push_back(min + (rand() % static_cast<int>(max - min + 1)));

    Eigen::Map<Eigen::VectorXi> V_temp(&(samples_vec[0]), samples_vec.size());
    obj.V_ind_samples = V_temp;

    // Compute and store corresponding lateral normals
    obj.NV_lateral.resize(N_samples, 3);
    for (int ii = 0; ii < N_samples; ii++)
    {
        obj.NV_lateral.row(ii) = NV.row(obj.V_ind(obj.V_ind_samples(ii)));
        obj.NV_lateral(ii, 2) = 0.0; // Remove any z-component - we only care about lateral movement
    }

    // ------ Debugging ------
    // std::cout << min_V(0) << ", " << min_V(1) << ", " << min_V(2) << std::endl;
    // std::cout << max_V(0) << ", " << max_V(1) << ", " << max_V(2) << std::endl;
    // -----------------------


  	return;

}
