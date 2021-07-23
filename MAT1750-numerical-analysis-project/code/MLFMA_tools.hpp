/**************************** MLFMM_tools.hpp ****************************

 * Library of routines for setup and management of MLFMA cubes.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 20, 2018

 *************************************************************************/


#ifndef MLFMATOOLS_H
#define MLFMATOOLS_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <complex>
#include <vector>
#include <map>
#include <stdlib.h>
#include "math.h"

// #include "mkl.h"
// #include <petscksp.h>
// #include <omp.h>

// #include "petsc_extensions.h"
#include "coordinate_system.hpp"
#include "mesh.hpp"


class Cube
{
 public:

 	// Cube side length
 	double w;

 	// Vertices
 	std::vector<vect3D> v = std::vector<vect3D> (8);

 	// Centre
 	vect3D centre;

 	// Is this the finest or coarsest level of the octree?
 	bool isLeaf = true, isMain = true;

 	// Octree level of this cube, and its unique global index, and its local index within its level
 	int level = -1, global_index = -1, local_index = -1;

 	// Octree parent of this cube
 	Cube *parent;

 	// Octree children of this cube
 	std::vector<Cube *> children;

 	// Connected, near, and far cubes
 	std::vector<Cube *> connected_cubes, near_cubes, far_cubes;

 	// Number of associated cubes
 	int N_connected, N_near, N_far;

 	// Indices to mesh elements that lie within this cube
 	std::vector<int> edges, triangles;

 	// Number of mesh elements
 	int N_edges, N_triangles;

 	void SetVertex(vect3D &_v, double _w = -1.0);
 	void SetVertex(std::vector<double> &_v, double _w = -1.0);
 	void SetVertex(double *_v, double _w = -1.0);

 	void SetVertices(std::vector<vect3D> &_v);
 	void SetVertices(vect3D *_v);

 	void PrintVertices();

 	void FindCentre();
 	void SetCentre(vect3D &_centre);
 	void SetCentre(std::vector<double> &_centre);
 	void SetCentre(double *_centre);

	void SetVerticesFromCentre(vect3D &_centre, double _w = -1.0);
	void SetVerticesFromCentre(std::vector<double> &_centre, double _w = -1.0);
	void SetVerticesFromCentre(double *_centre, double _w = -1.0);

};


class MLFMA_tools
{
 public:

 	Mesh *mesh;

	std::vector<vect3D> bounding_box = std::vector<vect3D> (8);

	std::vector<Cube *> all_cubes;
	std::vector<std::vector<Cube *>> cube_tree;

	int N_cubes, N_levels;

	double cube_threshold;

	void SetBoundingBox();
	void SetInitialCube();

	void RecursiveCubeDivision(Cube *prev_cube, std::vector<Cube *> &_all_cubes, double &threshold);
	void FinalizeCubeDivision();

	void InitializeSubCube(Cube *cube, vect3D &_v, double _w, int _level, Cube *_parent, int _local_index, int _global_index);
	void InitializeMainCube(vect3D &_v, double _w);

	bool isTriangleInCube(int triangle, Cube *cube);
	bool isEdgeInCube(int edge, Cube *cube);
	bool isCubeInBoundingBox(Cube *cube);

	void AssignTrianglesToCube(std::vector<int> &triangles, Cube *cube);
	void AssignEdgesToCube(std::vector<int> &edges, Cube *cube);

};


#endif
