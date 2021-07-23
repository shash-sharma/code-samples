#include "include/get_isolated_objects.h"
#include "include/get_N_obj.h"
#include "include/Object.h"
#include "include/get_min_signed_dist.h"
#include "include/sample_obj_vertices.h"
#include "include/skeletonize.h"
#include "include/segment_skeleton.h"
#include "include/segment_mesh.h"

#include <igl/avg_edge_length.h>
#include <igl/per_vertex_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/doublearea.h>
#include <igl/viewer/Viewer.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_matrix.h>

#include <Eigen/Core>

#include <iostream>
#include <vector>

// #include "cgal_wrappers.h"


int main (int argc, char *argv[])
{

  // Set up the viewer
  igl::viewer::Viewer viewer;
  std::cout<<R"(
S,s      Stretch, squish color axis range
0        Cycle through Objects to select one for visualizations
1        Step 1: Random sampling of points
2        Step 2: Skeletonize the sampled points
3        Step 3: Segment the skeleton
4        Step 4: Use skeletal "joints" to segment the mesh
5        Show full segmented mesh
)";

  int plot_obj = 0;

  // Scale for the color axis
  double scale = 100.0;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  // Load input mesh
  igl::read_triangle_mesh((argc>1?argv[1]:"../shared/data/test_case_0.obj"), V, F);

  // Find Object connectivity info and number of Objects
  int N_obj;
  Eigen::MatrixXd C, counts;
  get_N_obj(V, F, N_obj, C, counts);

  std::vector<Object> Objects (N_obj);

  // Get adjacency matrix for future convenience
  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);

  // Split geometry into distinct Objects
  get_isolated_objects(V, F, Objects, C, counts);

  // ------ Debugging ------
    // for (int ii = 0; ii < N_obj; ii++)
    //   std::cout << Objects[ii].V.rows() << ", " << Objects[ii].V.cols() << std::endl;
    // for (int ii = 0; ii < N_obj; ii++)
    //   std::cout << Objects[ii].F.rows() << ", " << Objects[ii].F.cols() << std::endl;
  // -----------------------

  // Get per-vertex normals
  igl::PerVertexNormalsWeightingType N_type = igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_UNIFORM;
  Eigen::MatrixXd NV;
  per_vertex_normals(V, F, N_type, NV);

  // Get random sampling of each Object's vertices and their corresponding lateral normals
  int N_samples = 30;
  for (int ii = 0; ii < N_obj; ii++)
    sample_obj_vertices(N_samples, Objects[ii], NV);

  // Associate vertices with triangles for later convenience
  std::vector<std::vector<int>> VF, VFi;
  igl::vertex_triangle_adjacency(V, F, VF, VFi);

  // Test random sampling of points
  Eigen::MatrixXd P;
  P.resize(Objects[plot_obj].V_ind_samples.rows(), 3);
  for (int ii = 0; ii < Objects[plot_obj].V_ind_samples.rows(); ii++)
    P.row(ii) = V.row(Objects[plot_obj].V_ind(Objects[plot_obj].V_ind_samples(ii)));


  const Eigen::RowVector3d orange (1.0, 0.5, 0.2);
  const Eigen::RowVector3d green (0.0, 1.0, 0.0);
  const Eigen::RowVector3d yellow (1.0, 0.9, 0.2);
  const Eigen::RowVector3d blue (0.0, 0.0, 1.0);
  const Eigen::RowVector3d grey (0.6, 0.6, 0.6);
  const Eigen::RowVector3d red (1.0, 0.0, 0.0);
  const Eigen::RowVector3d white (1.0, 1.0, 1.0);
  const Eigen::RowVector3d purple (1.0, 0.0, 1.0);

  std::vector<Eigen::RowVector3d> colours (8);
  colours[0] = orange;
  colours[1] = green;
  colours[2] = yellow;
  colours[3] = blue;
  colours[4] = grey;
  colours[5] = red;
  colours[6] = white;
  colours[7] = purple;

  for (int ii = 0; ii < N_obj; ii++)
  {
    // int ii = plot_obj;

    std::cout << "Processing object " << ii + 1 << " of " << N_obj << std::endl;

    // Move sampled points to Object centre to form a skeleton for each Object
    skeletonize(Objects[ii], V, F, A);

    // Test skeletonized points
    for (int ii = 0; ii < Objects[plot_obj].V_samples_sort_ind.size()-1; ii++)
      viewer.data.add_edges(
        Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii]), 
        Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii+1]), 
        Eigen::RowVector3d(1, 0, 0));

    // Segment the skeleton
    segment_skeleton(Objects[ii]);

    // Test skeleton joints
    // viewer.data.add_points(Objects[plot_obj].V_seg_points, Eigen::RowVector3d(0, 1, 0));

    // Segment the mesh by translating the skeleton joints on to the mesh
    segment_mesh(Objects[ii], V, F, VF);

    // Plot segmented vertices. Use a different colour for each segment.
    // int num_segments = Objects[plot_obj].V_tagged.maxCoeff();
    // for (int jj = 1; jj <= num_segments; jj++)
    //   for (int kk = 0; kk < Objects[plot_obj].V_tagged.size(); kk++)
    //     if (Objects[plot_obj].V_tagged(kk) == jj)
    //       viewer.data.add_points(Objects[plot_obj].V.row(kk), colours[jj-1]);

  }


  // Default view: different colour for each Object
  Eigen::VectorXd Z;
  Z.setZero(V.rows());

  for (int ii = 0; ii < Objects[plot_obj].V.rows(); ii++)
    Z(Objects[plot_obj].V_ind(ii)) = (double)100;

  viewer.data.set_mesh(V, F);
  viewer.core.invert_normals = true;
  viewer.core.show_lines = true;
  viewer.core.show_faces = true;

  const auto update = [&]()
  {
    Eigen::MatrixXd C;
    igl::parula(Z,-scale,scale,C);
    viewer.data.set_colors(C);
  };
  viewer.callback_key_pressed = 
    [&](igl::viewer::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
      case 'D':
      case 'd':
        viewer.core.show_overlay ^= 1;
        break;
      case '0':
        // Cycle through Objects to select one for visualizations
        plot_obj += 1;
        if (plot_obj == N_obj)
          plot_obj = 0;

        Z.setZero(V.rows());
        for (int ii = 0; ii < Objects[plot_obj].V.rows(); ii++)
          Z(Objects[plot_obj].V_ind(ii)) = (double)100;

        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = true;

        break;
      case '1':
        // Step 1: Random sampling of points
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = false;
        
        P.resize(Objects[plot_obj].V_ind_samples.rows(), 3);
        for (int ii = 0; ii < Objects[plot_obj].V_ind_samples.rows(); ii++)
          P.row(ii) = V.row(Objects[plot_obj].V_ind(Objects[plot_obj].V_ind_samples(ii)));

        viewer.data.add_points(P, Eigen::RowVector3d(1, 0, 0));

        break;
      case '2':
        // Step 2: Skeletonize the sampled points
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = false;

        viewer.data.add_points(Objects[plot_obj].V_samples, Eigen::RowVector3d(1, 0, 0));

        for (int ii = 0; ii < Objects[plot_obj].V_samples_sort_ind.size()-1; ii++)
          viewer.data.add_edges(
            Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii]), 
            Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii+1]), 
            Eigen::RowVector3d(1, 0, 0));

        break;
      case '3':
        // Step 3: Segment the skeleton
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = false;

        for (int ii = 0; ii < Objects[plot_obj].V_samples_sort_ind.size()-1; ii++)
          viewer.data.add_edges(
            Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii]), 
            Objects[plot_obj].V_samples.row(Objects[plot_obj].V_samples_sort_ind[ii+1]), 
            Eigen::RowVector3d(1, 0, 0));

        viewer.data.add_points(Objects[plot_obj].V_seg_points, Eigen::RowVector3d(0, 1, 0));

        break;
      case '4':
        // Step 4: Use skeletal "joints" to segment the mesh
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = false;

        // Plot segmented vertices. Use a different colour for each segment.
        for (int jj = 1; jj <= Objects[plot_obj].V_tagged.maxCoeff(); jj++)
          for (int kk = 0; kk < Objects[plot_obj].V_tagged.size(); kk++)
            if (Objects[plot_obj].V_tagged(kk) == jj)
              viewer.data.add_points(Objects[plot_obj].V.row(kk), colours[jj-1]);

        // Z.setZero(V.rows());
        // for (int ii = 0; ii < Objects[plot_obj].V.rows(); ii++)
        //   Z(Objects[plot_obj].V_ind(ii)) = Objects[plot_obj].V_tagged(ii)*100.0;

        break;
      case '5':
        // Show full segmented mesh
        viewer.data.clear();
        viewer.data.set_mesh(V, F);
        viewer.core.show_lines = true;
        viewer.core.show_faces = false;

        // Plot segmented vertices. Use a different colour for each segment.
        for (int ii = 0; ii < N_obj; ii++)
          for (int jj = 1; jj <= Objects[ii].V_tagged.maxCoeff(); jj++)
            for (int kk = 0; kk < Objects[ii].V_tagged.size(); kk++)
              if (Objects[ii].V_tagged(kk) == jj)
                viewer.data.add_points(Objects[ii].V.row(kk), colours[jj-1]);

        break;
      case 'S':
      case 's':
        scale *= key=='S' ? 2.0 : 0.5;
        std::cout<<"Color axis range: ["<<-scale<<","<<scale<<"]"<<std::endl;
        break;
      default:
        return false;
    }
    update();
    return true;
  };



  Eigen::MatrixXd lP(V.rows()*4,3);
  const double h = igl::avg_edge_length(V,F);
  Eigen::MatrixXd D1,D2;
  lP << V-0.5*h*D1, V+0.5*h*D1, V-0.5*h*D2, V+0.5*h*D2;
  Eigen::MatrixXi lE(2*V.rows(),2);
  Eigen::MatrixXd lC(2*V.rows(),3);

  for(int e = 0;e<V.rows();e++)
  {
    lE(e,0)          = e+0*V.rows();
    lE(e,1)          = e+1*V.rows();
    lE(V.rows()+e,0) = e+2*V.rows();
    lE(V.rows()+e,1) = e+3*V.rows();
    lC.row(         e) = orange;
    lC.row(V.rows()+e) = blue;
  }
  viewer.data.set_edges(lP,lE,lC);

  update();
  // viewer.core.show_lines = true;
  viewer.core.show_overlay = true;
  // viewer.data.face_based = false;
  viewer.launch();
  return EXIT_SUCCESS;
}
