

#include <cmath>
#include <complex>
#include <string>
#include <fstream>		//file I/O
#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>      //std::setprecision
#include <chrono>
#include <ctime>

#include "MLFMA.hpp"
#include "MLFMA_tools.hpp"
#include "coordinate_system.hpp"



// Usage: 


void TestCubeFunctions()
{

	Cube test_cube;

	// Width of test cube
	double w = 5.0e-3;

	// Origin
	vect3D origin = {0.0, 0.0, 0.0};

	// Points to insert into test cube
	std::vector<vect3D> pts (8);

	for (int ii = 0; ii < pts.size(); ii++)
		pts[ii] = origin;

	// Cube base
	pts[0] = origin;
	pts[1].x += w;
	pts[2].y += w;
	pts[3].x += w;
	pts[3].y += w;

	// Cube top
	pts[4].z += w;
	pts[5].x += w;
	pts[5].z += w;
	pts[6].y += w;
	pts[6].z += w;
	pts[7].x += w;
	pts[7].y += w;
	pts[7].z += w;

	// Set cube vertices
	test_cube.SetVertices(pts);

	// Print cube vertices
	test_cube.PrintVertices();

	// Compute and print cube centre
	test_cube.FindCentre();
	std::cout << "Centre of cube: " << test_cube.centre << std::endl;

	// Reset cube vertices by setting a new centre and cube width
	w = 2.0;
	vect3D new_centre = {-1.0, -1.0, -1.0};
	test_cube.SetVerticesFromCentre(new_centre, w);

	// Print cube vertices
	test_cube.PrintVertices();
	std::cout << "Centre of cube: " << test_cube.centre << std::endl;


	return;

}


void TestCubeDivision()
{

	// Create initial cube
	double w = 1.0;
	vect3D origin = {0.0, 0.0, 0.0};

	MLFMA_tools test_obj;

	test_obj.InitializeMainCube(origin, w);

	double threshold = 0.25;
	test_obj.RecursiveCubeDivision(test_obj.all_cubes[0], test_obj.all_cubes, threshold);
	test_obj.FinalizeCubeDivision();

	return;

}


int main (int num, char** file)
{

	TestCubeFunctions();

	TestCubeDivision();

	vect3D a = {1.0, 0.2, 0.5};
	vect3D b = {1.0, 0.5, 0.0};

	std::cout << dot(a, b) << std::endl;
	std::cout << cross(a, b) << std::endl;
	std::cout << dot(a, cross(a, b)) << std::endl;
	std::cout << dot(b, cross(a, b)) << std::endl;

	return 0;

}






