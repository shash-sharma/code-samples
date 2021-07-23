#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "../include/segment_mesh.h"
#include "../include/Object.h"

void segment_mesh(
    Object & obj,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    std::vector<std::vector<int>> & VF)
{

    // Strategy: for each skeletal joint, we know the lateral normal vector that was driving that point inwards.
    // We also know that a unit vector along z would be perpendicular to this lateral normal.
    // With these two vectors known, we can easily define a plane for that joint.
    // Then, for each such plane, we can figure out which mesh vertices are on either side of that plane.
    // Each vertex can then be tagged as belonging to a segment. Any pre-tagged vertices can be ignored.

    obj.V_tagged.resize(obj.V.rows());
    obj.V_tagged.setZero(obj.V.rows());
    int tag1 = 1;
    int tag2 = 2;
    for (int ii = 0; ii < obj.V_seg_points.rows(); ii++)
    {
        // Compute the skeletal plane
        Eigen::Vector3d N_lat = obj.NV_lateral.row(obj.V_seg_ind(ii));
        Eigen::Vector3d N_vert (3);
        N_vert << 0.0, 0.0, 1.0;

        Eigen::Vector3d N_plane = N_lat.cross(N_vert);

        // Now we have a vector normal to the skeletal plane, and a point on that plane.
        // Thus we have the equation of the plane, which looks like ax + by + cz + d = 0.
        double a = N_plane(0);
        double b = N_plane(1);
        double c = N_plane(2);
        double d = -a*obj.V_seg_points(ii, 0) - b*obj.V_seg_points(ii, 1) - c*obj.V_seg_points(ii, 2);

        // Now go through all previously untagged mesh vertices and tag them as part of a segment depending on what side of the plane they're on
        for (int jj = 0; jj < obj.V.rows(); jj++)
        {
            // Plug this vertex into the equation of the plane
            double test_point = a*obj.V(jj, 0) + b*obj.V(jj, 1) + c*obj.V(jj, 2) + d;

            // Depending on the sign of this, this vertex is on one side of the plane or the other
            if (test_point <= 0.0)
                obj.V_tagged(jj) += tag1;
            else
                obj.V_tagged(jj) += tag2;
        }
    }

    obj.V_tagged = obj.V_tagged - (obj.V_tagged.minCoeff() - 1)*Eigen::VectorXi::Ones(obj.V.rows());

    // Now, tag the corresponding faces that belong to each segment - TODO
    // for (int ii = 0; ii < obj.V_tagged.size(); ii++)
    // {

    // }

  	return;

}
