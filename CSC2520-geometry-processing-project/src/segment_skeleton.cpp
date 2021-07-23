#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <igl/signed_distance.h>

#include <Eigen/Core>

#include "../include/segment_skeleton.h"
#include "../include/Object.h"

void segment_skeleton(
    Object & obj)
{

    // Strategy:  start at the second point of the skeleton, so that there is at least one other point on either side.
    // Then, see if the points on either side of this point lie on the same line (within some angular tolerance).
    // Do this by looking at the lengths of the lines. Suppose B is the current point, A is the point on the left, C is the point on the right.
    // Then if the length of AC = 2AB = 2BC within some tolerance, they must be (almost) collinear.
    // If they're not collinear, then B constitutes a segment joint. Mark it as such and move on to the next point.
    // If you've travelled a certain distance without encountering a bend, then also mark that point as a segment.
    // Also, make sure that two joints are not too close to each other - this likely indicates false positives.

    std::vector<int> seg_ind_vec;

    // Settings - \todo make these function arguments
    double max_length = 200.0e-6;
    double current_length = 0.0;
    double false_positive_length = 50.0e-6;

    // max_length = 1.0e16;
    // false_positive_length = 0.0;

    for (int ii = 1; ii < obj.V_samples_sort_ind.size()-1; ii++)
    {
        double AC_xdiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 0) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 0);
        double AC_ydiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 1) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 1);
        double AC_zdiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 2) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 2);

        double AB_xdiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 0) - obj.V_samples(obj.V_samples_sort_ind[ii], 0);
        double AB_ydiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 1) - obj.V_samples(obj.V_samples_sort_ind[ii], 1);
        double AB_zdiff = obj.V_samples(obj.V_samples_sort_ind[ii-1], 2) - obj.V_samples(obj.V_samples_sort_ind[ii], 2);

        double BC_xdiff = obj.V_samples(obj.V_samples_sort_ind[ii], 0) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 0);
        double BC_ydiff = obj.V_samples(obj.V_samples_sort_ind[ii], 1) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 1);
        double BC_zdiff = obj.V_samples(obj.V_samples_sort_ind[ii], 2) - obj.V_samples(obj.V_samples_sort_ind[ii+1], 2);

        double AC = std::sqrt( (AC_xdiff)*(AC_xdiff) + (AC_ydiff)*(AC_ydiff) + (AC_zdiff)*(AC_zdiff) );
        double AB = std::sqrt( (AB_xdiff)*(AB_xdiff) + (AB_ydiff)*(AB_ydiff) + (AB_zdiff)*(AB_zdiff) );
        double BC = std::sqrt( (BC_xdiff)*(BC_xdiff) + (BC_ydiff)*(BC_ydiff) + (BC_zdiff)*(BC_zdiff) );

        double tol = 5.0e-1; // Percent
        current_length += AB;

        if (100.0*std::abs((AC - (AB + BC))/AC) > tol || current_length >= max_length) // Segment joint found
        {

            double dist_seg = 1.0e16;
            if (seg_ind_vec.size() > 0)
            {
                Eigen::MatrixXd dist_vec (1, 3);
                dist_vec = obj.V_samples.row(obj.V_samples_sort_ind[ii]) - obj.V_samples.row(seg_ind_vec[seg_ind_vec.size()-1]);
                dist_seg = std::sqrt( dist_vec(0, 0)*dist_vec(0, 0) + dist_vec(0, 1)*dist_vec(0, 1) + dist_vec(0, 2)*dist_vec(0, 2) );
            }

            if (dist_seg >= false_positive_length)
            {
                seg_ind_vec.push_back(obj.V_samples_sort_ind[ii]);
                current_length = 0.0;
            }
            
        }
    }

    Eigen::Map<Eigen::VectorXi> V_temp(&(seg_ind_vec[0]), seg_ind_vec.size());
    obj.V_seg_ind = V_temp;

    obj.V_seg_points.resize(seg_ind_vec.size(), 3);
    for (int ii = 0; ii < seg_ind_vec.size(); ii++)
        obj.V_seg_points.row(ii) = obj.V_samples.row(seg_ind_vec[ii]);

    // std::cout << obj.V_seg_points << std::endl;

  	return;

}
