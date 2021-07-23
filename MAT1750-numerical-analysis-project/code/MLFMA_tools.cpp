/**************************** MLFMM_tools.cpp ****************************

 * Library of routines for setup and management of MLFMA cubes.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 20, 2018

 *************************************************************************/


#include "MLFMA_tools.hpp"


/*! \brief Overloaded function to set the cube's vertex with smallest coordinates via an input vect3D object. The rest of the vertices are set accordingly, based on w.*/
void Cube::SetVertex(vect3D &_v, double _w)
{
	if (_w > 0.0)
		w = _w;

	for (int ii = 0; ii < v.size(); ii++)
		v[ii] = _v;

	// Cube base
	v[1].x += w;
	v[2].y += w;
	v[3].x += w;
	v[3].y += w;

	// Cube top
	v[4].z += w;
	v[5].x += w;
	v[5].z += w;
	v[6].y += w;
	v[6].z += w;
	v[7].x += w;
	v[7].y += w;
	v[7].z += w;

	FindCentre();

	return;
}


/*! \brief Overloaded function to set the cube's vertex with smallest coordinates via an input std::vector of size 3 (x,y,z). The rest of the vertices are set accordingly, based on w.*/
void Cube::SetVertex(std::vector<double> &_v, double _w)
{
	if (_w > 0.0)
		w = _w;

	vect3D _v_set(_v);
	SetVertex(_v_set);

	return;
}


/*! \brief Overloaded function to set the cube's vertex with smallest coordinates via an input array of size 3 (x,y,z). The rest of the vertices are set accordingly, based on w.*/
void Cube::SetVertex(double *_v, double _w)
{
	if (_w > 0.0)
		w = _w;

	vect3D _v_set(_v);
	SetVertex(_v_set);

	return;
}


/*! \brief Overloaded function to set cube vertices via an input std::vector of 8 vect3D objects.*/
void Cube::SetVertices(std::vector<vect3D> &_v)
{

	for (int ii = 0; ii < _v.size(); ii++)
		v[ii] = _v[ii];

	w = v[1].x - v[0].x;

	FindCentre();

	return;

}


/*! \brief Overloaded function to set cube vertices via an input array of 8 vect3D objects.*/
void Cube::SetVertices(vect3D *_v)
{

	for (int ii = 0; ii < 8; ii++)
		v[ii] = _v[ii];

	w = v[1].x - v[0].x;

	FindCentre();

	return;

}

/*! \brief Function to print cube vertices to terminal.*/
void Cube::PrintVertices()
{

	std::cout << "\n========================" << std::endl;
	std::cout << "Cube vertices:" << std::endl;
	for (int ii = 0; ii < v.size(); ii++)
		std::cout << v[ii] << std::endl;
	std::cout << "========================" << std::endl;

	return;

}


/*! \brief Function to find and set the centre coordinate of the cube.*/
void Cube::FindCentre()
{
	for (int ii = 0; ii < v.size(); ii++)
		centre += v[ii];

	centre /= 8.0;

	return;
}


/*! \brief Overloaded function to set cube centre via an input vect3D object.*/
void Cube::SetCentre(vect3D &_centre)
{
	centre = _centre;
	return;
}


/*! \brief Overloaded function to set cube centre via an input std::vector of size 3 (x,y,z).*/
void Cube::SetCentre(std::vector<double> &_centre)
{
	centre.x = _centre[0];
	centre.y = _centre[1];
	centre.z = _centre[2];

	return;
}


/*! \brief Overloaded function to set cube centre via an input array of size 3 (x,y,z).*/
void Cube::SetCentre(double *_centre)
{
	centre.x = _centre[0];
	centre.y = _centre[1];
	centre.z = _centre[2];

	return;
}


/*! \brief Overloaded function to set cube vertices upon being provided its centre and optionally its side length.*/
void Cube::SetVerticesFromCentre(vect3D &_centre, double _w)
{
	if (_w > 0.0)
		w = _w;

	SetCentre(_centre);

	// Set vertices
	_centre -= w/2.0;
	SetVertex(_centre);

	return;
}


/*! \brief Overloaded function to set cube vertices upon being provided its centre and optionally its side length.*/
void Cube::SetVerticesFromCentre(std::vector<double> &_centre, double _w)
{
	if (_w > 0.0)
		w = _w;

	vect3D _centre_set(_centre);
	SetVerticesFromCentre(_centre_set);

	return;
}


/*! \brief Overloaded function to set cube vertices upon being provided its centre and optionally its side length.*/
void Cube::SetVerticesFromCentre(double *_centre, double _w)
{
	if (_w > 0.0)
		w = _w;

	vect3D _centre_set(_centre);
	SetVerticesFromCentre(_centre_set);

	return;
}


void MLFMA_tools::InitializeSubCube(Cube *cube, vect3D &_v, double _w, int _level, Cube *_parent, int _local_index, int _global_index)
{
	cube->SetVertex(_v, _w);
	cube->level = _level;
	cube->parent = _parent;
	cube->local_index = _local_index;
	cube->global_index = _global_index;
	cube->isMain = false;

	AssignTrianglesToCube(_parent->triangles, cube);
	AssignEdgesToCube(_parent->edges, cube);

	all_cubes.push_back(cube);

	return;
}


void MLFMA_tools::InitializeMainCube(vect3D &_v, double _w)
{
	Cube *main_cube = new Cube;
	main_cube->SetVertex(_v, _w);
	main_cube->level = 0;
	main_cube->parent = NULL;
	main_cube->local_index = 0;
	main_cube->global_index = 0;
	main_cube->isMain = true;

	all_cubes.clear();
	all_cubes.push_back(main_cube);

	return;
}


/*! \brief Function to run recursive sub-division of a given outermost cube until the size of cubes hits the provided threshold.*/
void MLFMA_tools::RecursiveCubeDivision(Cube *prev_cube, std::vector<Cube *> &_all_cubes, double &threshold)
{

	double w = prev_cube->w;
	double w_new = w/2.0;

	// Test break condition
	if (w_new < threshold)
	{
		prev_cube->isLeaf = true;
		return;
	}

	// Strategy: starting from the vertex with smallest coordinates, we proceed to create 8 smaller cubes
	int level_new = prev_cube->level + 1;
	int global_index_new = (prev_cube->global_index)*8 + 1;
	int local_index_new = global_index_new - std::pow(8, prev_cube->level) - prev_cube->level;

	std::vector<Cube *> cubes (8);

	// First sub-cube
	cubes[0] = new Cube;
	InitializeSubCube(cubes[0], prev_cube->v[0], w_new, level_new, prev_cube, local_index_new, global_index_new);
	local_index_new++;
	global_index_new++;

	// Subsequent sub-cubes
	for (int ii = 1; ii < 8; ii++)
	{
		cubes[ii] = new Cube;
		InitializeSubCube(cubes[ii], cubes[0]->v[ii], w_new, level_new, prev_cube, local_index_new, global_index_new);
		local_index_new++;
		global_index_new++;
	}

	// Update children of previous cube
	prev_cube->children.resize(8);

	for (int ii = 1; ii < 8; ii++)
		prev_cube->children[ii] = _all_cubes[_all_cubes.size() - 8 + ii];


	// DEBUGGING
	for (int ii = 0; ii < 8; ii++)
	{
		std::cout << "\nCube " << cubes[ii]->global_index << ", local index " << cubes[ii]->local_index << ", level " << level_new << ":" << std::endl;
		// cubes[ii]->PrintVertices();
	}


	// Proceed to next level of recursion for each new sub-cube
	for (int ii = 0; ii < 8; ii++)
		RecursiveCubeDivision(cubes[ii], _all_cubes, threshold);


	return;
}


/*! \brief Function to be called after recursive cube division to update final object properties.*/
void MLFMA_tools::FinalizeCubeDivision()
{
	N_cubes = all_cubes.size();
	N_levels = all_cubes.back()->level + 1;

	cube_tree.clear();
	cube_tree.resize(N_levels);

	int start_index = 0;
	int end_index = 1;
	
	for (int ii = 0; ii < N_levels; ii++)
	{
		cube_tree[ii].insert(cube_tree[ii].begin(), all_cubes.begin() + start_index, all_cubes.begin() + end_index);

		start_index += std::pow(8, ii);
		end_index = start_index + std::pow(8, ii+1);
	}


	// DEBUGGING
	// for (int ii = 0; ii < N_levels; ii++)
	// {
	// 	std::cout << "Level " << ii << " has cubes " << std::flush;
	// 	for (int jj = 0; jj < cube_tree[ii].size(); jj++)
	// 		std::cout << cube_tree[ii][jj]->global_index << " " << std::flush;
	// 	std::cout << std::endl;
	// }

	return;
}



/*! \brief Function to check if a given mesh triangle is in a given cube.*/
bool MLFMA_tools::isTriangleInCube(int triangle, Cube *cube)
{
	// A triangle is in a cube if its centroid is within the cube. If the centroid is on a cube face or edge, then accept it as part of the cube. Duplicates can be removed later (TODO).
	
	vect3D centroid = mesh->triangles[triangle]->centroid;

	if (centroid >= cube->v[0] && centroid <= cube->v[7])
		return true;
	else
		return false;
}


/*! \brief Function to check if a given mesh edge is in a given cube.*/
bool MLFMA_tools::isEdgeInCube(int edge, Cube *cube)
{
	// An edge is in a cube if its midpoint is within the cube. If the midpoint is on a cube face or edge, then accept it as part of the cube. Duplicates can be removed later (TODO).
	
	vect3D v1 = mesh->edges[edge]->nodes[0]->p;
	vect3D v2 = mesh->edges[edge]->nodes[1]->p;
	vect3D centre = (v1 + v2)/2.0;

	if (centre >= cube->v[0] && centre <= cube->v[7])
		return true;
	else
		return false;
}


/*! \brief Function to check if a given cube is within the structure's bounding box.*/
bool MLFMA_tools::isCubeInBoundingBox(Cube *cube)
{
	if (cube->v[0] >= bounding_box[0] && cube->v[7] <= bounding_box[7])

	return false;
}


/*! \brief Find all triangles from a given list that are in the given cube, and assign them to the cube. Note: this clears the triangle indices currently stored in the cube.*/
void MLFMA_tools::AssignTrianglesToCube(std::vector<int> &triangles, Cube *cube)
{
	cube->triangles.resize(0);

	for (int ii = 0; ii < triangles.size(); ii++)
		if (isTriangleInCube(triangles[ii], cube))
			cube->triangles.push_back(triangles[ii]);

	cube->N_triangles = cube->triangles.size();

	return;
}


/*! \brief Find all edges from a given list that are in the given cube, and assign them to the cube. Note: this clears the edge indices currently stored in the cube.*/
void MLFMA_tools::AssignEdgesToCube(std::vector<int> &edges, Cube *cube)
{
	cube->edges.resize(0);

	for (int ii = 0; ii < edges.size(); ii++)
		if (isEdgeInCube(edges[ii], cube))
			cube->edges.push_back(edges[ii]);

	cube->N_edges = cube->edges.size();

	return;
}




