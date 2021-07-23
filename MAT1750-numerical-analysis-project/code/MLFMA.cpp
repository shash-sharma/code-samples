/****************************** MLFMA.cpp *********************************

 * Library of routines for the multi-level fast multipole algorithm in 3D.
 *
 * Author: Shashwat Sharma
 * Created on: Nov 26, 2018

 *************************************************************************/


#include "MLFMA.hpp"



/*! \brief Function to set the simulation frequency and run frequency-dependent initializations.*/
void MLFMA::SetFrequency(double _f)
{
	f = _f;
	omega = 2.0*M_PI*f;
	k0 = omega*std::sqrt(eps0*mu0);

	return;
}


/*! \brief Function to generate near region matrices for each MLFMA region.*/
void MLFMA::GenerateNearMatrices()
{

	MoM_integrator_setup setup;
	setup.k = k0;
	setup.omega = omega;
	
	mtls->GenerateNearMatrices(setup);
	mtls->GenerateFarFunctions(setup);

	mtls->InitializeMatVec(mtls->all_sphere_data);

	return;

}


/*! \brief Helper function for the matrix-free iterative solver.*/
PetscErrorCode MLFMA::usermult(Mat mat, Vec x, Vec y)
{
	MatrixVectorProduct(x, y);
	return 0;
}


/*! \brief Helper function for the matrix-free iterative solver.*/
PetscErrorCode usermult_external(Mat mat, Vec x, Vec y)
{
	MLFMA *solver;
	MatShellGetContext(mat, &solver);	
	solver->usermult(mat, x, y);

	return 0;
}


/*! \brief Matrix-vector product for the iterative solver, for a delta-gap excitation with RWG basis functions.*/
void MLFMA::MatrixVectorProduct(Vec &x, Vec &y)
{

	// ====== Near interactions ======

	MatMult(mtls->ST_near, x, y);


	// ====== Far interactions ======

	mtls->Aggregation_RWG(mtls->cube_tree.back(), mtls->all_sphere_data, x);
	mtls->Disaggregation_RWG(mtls->cube_tree.back(), mtls->all_sphere_data, y);


	return;

}


/*! \brief Function to solve an Ax = b system using GMRES in PETSc, with a delta-gap excitation.*/
void MLFMA::IterativeSolver_DeltaGap()
{

	int Nu = mtls->mesh->N_edges;
	
	// Solve the matrix with a delta gap excitation to test
	Vec b, x;
	VecCreateAndAllocate(b, Nu);
	VecCreateAndAllocate(x, Nu);
	
	VecSetValue(b, 0, 1.0e-9, INSERT_VALUES);
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);


  	// Initialize the iterative solver
	KSP ksp;
	PC pc;
	Mat A;
	MPI_Comm comm;
	int n_iter = 0, n_iter_avg = 0;
	comm = MPI_COMM_SELF;
	void *ctx = (void *)this;

	// Create the system matrix shell
	MatCreateShell(comm, Nu, Nu, Nu, Nu, ctx, &A);
	MatShellSetOperation(A, MATOP_MULT, (void(*)(void))usermult_external);

	// Set up the solver object
	KSPCreate(comm, &ksp);
	// KSPSetOperators(ksp, A, A);
	KSPSetOperators(ksp, A, mtls->ST_near);

	KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
	KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
	KSPGMRESSetRestart(ksp, 200);
	
	PetscReal rtol = 1.0e-6;
	PetscReal atol = 1.0e-16;
	int ksp_maxiter = 100;
	KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, ksp_maxiter);

	// KSPSetPCSide(ksp, PC_RIGHT);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCJACOBI);
	// PCSetType(pc, PCLU);
	
	KSPSetFromOptions(ksp);	
	KSPSetUp(ksp);

	
	std::cout << "Solving the system..." << std::endl;

	// Solve the system
	std::clock_t t1_KSPsolve = clock();
	KSPSolve(ksp, b, x);
	std::clock_t t2_KSPsolve = clock();
	double dt_KSPsolve = (t2_KSPsolve - t1_KSPsolve)/(double)CLOCKS_PER_SEC;

	std::cout << "KSP solve took " << dt_KSPsolve << " s." << std::endl;
	
	// Number of iterations run
	KSPGetIterationNumber(ksp, &n_iter);

	// Get convergence reason
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);

	

	PrintVecToFile(x, "./src/MLFMA/output/Vec_x_iter");


	std::complex<double> *Js, *E_inc;
	VecGetArray(x, &Js);
	VecGetArray(b, &E_inc);
	
	std::complex<double> Z_in = E_inc[0]/Js[0];
	std::cout << Z_in << std::endl;

	
	post_process post(mtls->mesh);
	// post.vis.VisualizeMesh("Mesh", false, true, VTK_DEFAULT, 1.0);

	post.VisualizeFaceScalarField("./src/MLFMA/output/Js_iter", Js, true, false);



	VecRestoreArray(x, &Js);
	VecRestoreArray(b, &E_inc);
	VecDestroy(&x);
	VecDestroy(&b);



	return;
	
}


/*! \brief Overloaded function to solve an Ax = b system using PETSc, with one right-hand side.*/
void MLFMA::DirectSolver_DeltaGap()
{
	
	// Solve the matrix with a delta gap excitation to test
	Vec b, x;
	VecCreateAndAllocate(b, mtls->mesh->N_edges);

	VecSetValue(b, 0, 1.0e-9, INSERT_VALUES);
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);


	std::cout << "Solving system..."  << std::flush;
	std::clock_t t1_sol = clock();

	MatSolve_SuperLU_dist(mtls->ST_near, x, b);

	std::clock_t t2_sol = clock();
	double dt_sol = (t2_sol - t1_sol)/(double)CLOCKS_PER_SEC;
	std::cout << "done in " << dt_sol << " s."  << std::endl;


	PrintVecToFile(x, "./src/MLFMA/output/Vec_x_direct");


	std::complex<double> *Js, *E_inc;
	VecGetArray(x, &Js);
	VecGetArray(b, &E_inc);
	
	std::complex<double> Z_in = E_inc[0]/Js[0];
	std::cout << Z_in << std::endl;



	
	post_process post(mtls->mesh);
	// post.vis.VisualizeMesh("Mesh", false, true, VTK_DEFAULT, 1.0);

	post.VisualizeFaceScalarField("./src/MLFMA/output/Js_direct", Js, true, false);


	VecRestoreArray(x, &Js);
	VecRestoreArray(b, &E_inc);
	VecDestroy(&x);
	VecDestroy(&b);


	return;
	
}


/*! \brief Overloaded function to solve an AX = B system using PETSc, for multiple right-hand sides.*/
void MLFMA::DirectSolver(Mat &A, Mat &X, Mat &B)
{

	MatMatSolve_SuperLU_dist(A, X, B);

	return;
	
}


/*! \brief This function is to be called in between frequency sweep points to reset selected data containers to prevent contamination, and undo frequency-dependent matrix scalings.*/
void MLFMA::ResetContainers()
{

	MatZeroEntries(L_near);

	return;

}


/*! \brief Function to release all memory at the end of a simulation.*/
void MLFMA::ReleaseMemory()
{

	MatDestroy(&L_near);

	return;

}
