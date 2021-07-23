#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <igl/signed_distance.h>

#include <Eigen/Core>

#include "../include/get_min_signed_dist.h"
#include "../include/Object.h"

void get_min_signed_dist(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Object obj,
    Eigen::VectorXd & S,
    Eigen::MatrixXd & P)
{


    // First, create a grid over which to compute signed distances
    Eigen::VectorXd min_V = obj.V.colwise().minCoeff();
    Eigen::VectorXd max_V = obj.V.colwise().maxCoeff();

    // Decide number of grid points and grid spacing along each direction
    // double d = 1.0e-6;
    // int Nx = std::min((int)((max_V(0) - min_V(0))/d), 50);
    // int Ny = std::min((int)((max_V(1) - min_V(1))/d), 50);
    // int Nz = std::min((int)((max_V(2) - min_V(2))/d), 50);

    // int Nx = 50;
    // int Ny = 50;
    // int Nz = 5;
    // double dx = (max_V(0) - min_V(0))/((double)Nx-1.0);
    // double dy = (max_V(1) - min_V(1))/((double)Ny-1.0);
    // double dz = (max_V(2) - min_V(2))/((double)Nz-1.0);

    // Create list of query points (grid)
    // P.resize(Nx*Ny*Nz, 3);

    // // Create grid points
    // for (int ii = 0; ii < Nx; ii++)
    // {
    //     for (int jj = 0; jj < Ny; jj++)
    //     {
    //         for (int kk = 0; kk < Nz; kk++)
    //         {
    //             P((kk + jj*Nz + ii*(Nz*Ny)), 0) = min_V(0) + ((double)ii)*dx;
    //             P((kk + jj*Nz + ii*(Nz*Ny)), 1) = min_V(1) + ((double)jj)*dy;
    //             P((kk + jj*Nz + ii*(Nz*Ny)), 2) = min_V(2) + ((double)kk)*dz;
    //         }
    //     }
    // }


    P.resize(1, 3);
    P(0, 0) = -2.0e-6 + (min_V(0))/1.0;
    P(0, 1) = 8.0e-6 + (min_V(1))/1.0;
    P(0, 2) = (max_V(2) + min_V(2))/2.0;

// std::cout << "Here2" << std::endl;

    // Compute signed distances
    Eigen::VectorXi I;
    Eigen::MatrixXd C, N;

    // std::cout << "Here" << std::endl;

    igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL;
    igl::signed_distance(P, V, F, sign_type, S, I, C, N);

    // igl::SignedDistanceType sign_type = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
    // igl::signed_distance(P, obj.V, obj.F, sign_type, 
    //     std::numeric_limits<double>::min(), std::numeric_limits<double>::max(), S, I, C, N);


    if (S(0) != S(0))
        S(0) = 1.0e-15;
    std::cout << S(0) << std::endl;

    // double min_dist = S.minCoeff();
    // double thresh_dist = 0.05*min_dist;
    // Eigen::MatrixXd C_min;

    // for (int ii = 0; ii < S.rows(); ii++)
    // {
    //     if (S(ii) <= thresh_dist)

    // }

    // ------ Debugging ------
    // std::cout << min_V(0) << ", " << min_V(1) << ", " << min_V(2) << std::endl;
    // std::cout << max_V(0) << ", " << max_V(1) << ", " << max_V(2) << std::endl;
    // -----------------------


  	return;

}
