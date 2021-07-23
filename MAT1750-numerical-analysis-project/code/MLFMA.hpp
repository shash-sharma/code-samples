/****************************** MLFMA.hpp *********************************

 * Library of routines for the multi-level fast multipole algorithm in 3D.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 20, 2018

 *************************************************************************/


#ifndef MLFMA_H
#define MLFMA_H


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <complex>
#include <vector>
#include <map>
#include <stdlib.h>
#include <ctime>
#include "math.h"

#include <petscksp.h>
// #include "mkl.h"
// #include <omp.h>

#include "petsc_extensions.h"
#include "constants.hpp"
#include "progress_bar.hpp"
#include "coordinate_system.hpp"
#include "mesh.hpp"
#include "MLFMA_tools.hpp"
#include "structure.hpp"
#include "post_process.hpp"

// ====== Solver ======

#include "MoM_integrals.hpp"


// Helper functions for the matrix-free iterative solver
PetscErrorCode usermult_external(Mat mat, Vec x, Vec y);


class MLFMA
{
 public:

	// Main structure
 	Structure *str;

	// Object meshes
 	std::vector<Mesh *> *objects;

	// Number of MLFMA instances / regions
	int N_mlfma = 1;

	// Pointer to all MLFMA regions
	MLFMA_tools *mtls;


	// Frequency data
	double f = 1.0e9;
	double omega = 2.0*M_PI*f;
	double k0 = omega*std::sqrt(eps0*mu0);
	

 	// ====== Near matrices ======

	// Single layer tangential operator for near interactions
	Mat L_near;


	void SetFrequency(double _f);
	
	void GenerateNearMatrices();
	void GenerateFarFunctions();

	void MatrixVectorProduct(Vec &x, Vec &y);
	
	void IterativeSolver_DeltaGap();

	void DirectSolver_DeltaGap();
	void DirectSolver(Mat &A, Mat &X, Mat &B);
	
	void ResetContainers();
	void ReleaseMemory();

	PetscErrorCode usermult(Mat mat, Vec x, Vec y);

};


#endif
