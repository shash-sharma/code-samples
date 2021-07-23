#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <igl/signed_distance.h>

#include <Eigen/Core>

#include "../include/skeletonize.h"
#include "../include/Object.h"

void skeletonize(
    Object & obj,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::SparseMatrix<double> & A)
{

    // Copy the global vertices over so we can modify at will
    Eigen::MatrixXd V_mod = V;

    // So that our signed distances are not affected by the z-coordinate, modify the z-coordinates to expand them outwards
    // To do this evenly upwards and downwards, find the max and min possible z-points.

    Eigen::VectorXd min_V = obj.V.colwise().minCoeff();
    Eigen::VectorXd max_V = obj.V.colwise().maxCoeff();

    double min_z = min_V(2);
    double max_z = max_V(2);

    // Expand outwards by a factor proportional to the height of this Object
    double factor = 3.0;
    double height = std::abs(max_z - min_z);
    double mid = 0.5*std::abs(max_z + min_z);

    // Adjust the corresponding vertices in the global vertex list copy
    double tol = 1.0e-15;
    for (int ii = 0; ii < obj.V_ind.rows(); ii++)
    {
        if (std::abs(V(obj.V_ind(ii), 2) - max_z) < tol)
            V_mod(obj.V_ind(ii), 2) += factor*height;
        else if (std::abs(V(obj.V_ind(ii), 2) - min_z) < tol)
            V_mod(obj.V_ind(ii), 2) -= factor*height;
    }

    // Set the sampled points to lie in the centre (z-wise) of the Object.
    // Then translate the sampled points in the direction of decreasing signed distance function.
    // Use a line search to do this. This is basically a gradient descent, except that we know the
    // descent direction already - it is given by the lateral per-vertex normals.
    // Choose a step size that corresponds to the desired level of accuracy.

    obj.V_samples.resize(obj.V_ind_samples.rows(), 3);
    double step = 0.5e-6;
    igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;

    for (int ii = 0; ii < obj.V_ind_samples.rows(); ii++)
    {
        obj.V_samples.row(ii) = V.row(obj.V_ind(obj.V_ind_samples(ii)));
        obj.V_samples(ii, 2) = mid;

        // Compute initial signed distance
        Eigen::VectorXd S0, S;
        Eigen::VectorXi I;
        Eigen::MatrixXd C, N;

        igl::signed_distance(obj.V_samples.row(ii), V_mod, F, sign_type, S0, I, C, N);

        // Now start the line search
        int flag = 0;
        while (flag == 0)
        {
            // Translate the current point
            obj.V_samples(ii, 0) += obj.NV_lateral(ii,0)*step;
            obj.V_samples(ii, 1) += obj.NV_lateral(ii,1)*step;

            // Compute new signed distance
            igl::signed_distance(obj.V_samples.row(ii), V_mod, F, sign_type, S, I, C, N);

            if (S(0) <= S0(0)) // We've gone too far. Go back and stop.
            {
                // Translate back to previous point
                obj.V_samples(ii, 0) -= obj.NV_lateral(ii,0)*step;
                obj.V_samples(ii, 1) -= obj.NV_lateral(ii,1)*step;

                flag = 1;
                break;
            }

            S0 = S;
        }
    } 

    // We want the skeletal points to be in order from one side of the structure to the other.

    // First, find the point closest to one corner of the structure. This must be one end of the skeleton.
    double dist = 1.0e16;
    obj.V_samples_sort_ind.resize(obj.V_ind_samples.rows());
    std::vector<int> temp_unsorted(obj.V_ind_samples.rows());
    int first_ind;
    for (int ii = 0; ii < obj.V_ind_samples.rows(); ii++)
    {
        temp_unsorted[ii] = ii;

        double dist_test = std::sqrt( (min_V(0) - obj.V_samples(ii, 0))*(min_V(0) - obj.V_samples(ii, 0)) + 
            (min_V(1) - obj.V_samples(ii, 1))*(min_V(1) - obj.V_samples(ii, 1)) + 
            (min_V(2) - obj.V_samples(ii, 2))*(min_V(2) - obj.V_samples(ii, 2)) );
        if (dist_test < dist)
        {
            dist = dist_test;
            obj.V_samples_sort_ind[0] = ii;
            first_ind = ii;
        }
    }
    temp_unsorted.erase(temp_unsorted.begin() + first_ind);

    // Now we have the starting point. So sort through the rest of the points.
    double norm_tol = 1.0e-2;
    for (int ii = 0; ii < obj.V_ind_samples.rows()-1; ii++)
    {
        dist = 1.0e16;
        int temp_ind;

        // Note - we have to make sure that the next closest point does not happen to be across a boundary.
        // To do this, make sure that the vector between the next nearest vertex has no component along the lateral normal
        for (int jj = 0; jj < temp_unsorted.size(); jj++)
        {
            Eigen::VectorXd dist_vec (3);
            dist_vec(0) = (obj.V_samples(obj.V_samples_sort_ind[ii], 0) - obj.V_samples(temp_unsorted[jj], 0));
            dist_vec(1) = (obj.V_samples(obj.V_samples_sort_ind[ii], 1) - obj.V_samples(temp_unsorted[jj], 1));
            dist_vec(2) = (obj.V_samples(obj.V_samples_sort_ind[ii], 2) - obj.V_samples(temp_unsorted[jj], 2));

            double dist_test = std::sqrt( dist_vec(0)*dist_vec(0) + 
                dist_vec(1)*dist_vec(1) + dist_vec(2)*dist_vec(2) );

            Eigen::VectorXd dist_vec_dot_N_lateral (3);
            dist_vec_dot_N_lateral (0) = dist_vec(0)*obj.NV_lateral(obj.V_samples_sort_ind[ii], 0);
            dist_vec_dot_N_lateral (1) = dist_vec(1)*obj.NV_lateral(obj.V_samples_sort_ind[ii], 1);
            dist_vec_dot_N_lateral (2) = dist_vec(2)*obj.NV_lateral(obj.V_samples_sort_ind[ii], 2);

            double dist_norm = std::abs(dist_vec_dot_N_lateral(0) + dist_vec_dot_N_lateral(1) + dist_vec_dot_N_lateral(2));

            if (dist_test < dist && dist_norm < norm_tol)
            {
                dist = dist_test;
                obj.V_samples_sort_ind[ii+1] = temp_unsorted[jj];
                temp_ind = jj;
            }
        }

        temp_unsorted.erase(temp_unsorted.begin() + temp_ind);
    }

    // Use the vertex adjacency matrix to find next nearest vertex, rather than physical distances - TODO
    // A(obj.V_ind(obj.V_ind_samples(obj.V_samples_sort_ind[0]))

    // Test sorting
    // for (int ii = 0; ii < obj.V_ind_samples.rows(); ii++)
    //     std::cout << obj.V_samples(obj.V_samples_sort_ind[ii], 1) << std::endl;


  	return;

}
