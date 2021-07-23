#ifndef PETSC_EXTENSIONS_H
#define PETSC_EXTENSIONS_H

#include <cstring>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <set>

#include <petscksp.h>
#include <petscdraw.h>

const int PETSC_ROW_MAJOR = 1;
const int PETSC_COL_MAJOR = 2;


/*! \brief Data structure to store input settings and output statistics for iterative solvers.*/
struct KSP_data
{
	
	// Solver settings
	PetscReal rtol = 1.0e-4;
	PetscReal atol = 1.0e-16;
	int max_it = 1000;
	int restart = 1000;
	bool allocate_x = true;
	bool use_fischer_guess = false;
	int fischer_size = 200;
	
	bool compute_singular_values = false;
	bool plot_residual_norm = false;
	bool monitor_residual = false;
	bool report_stats = true;
	bool report_general = true;

	std::string plot_title = "";
	std::string residual_filename = "rrnorm.txt";

	// Preconditioner settings
	void *pc_ctx = NULL;
	PetscErrorCode (*pc_apply)(PC, Vec, Vec) = NULL;
	PCType pc_type = PCNONE;
	PCSide pc_side = PC_RIGHT;
	Mat *M = NULL;

	// Initial guess
	/* Vec *x0 = NULL; */
	Vec x0 = PETSC_NULL;

	// ILU settings
	bool use_ilu = false;
	int ilu_level = 2;
	int ilu_fillfactor = 2;

	// Output data
	int N_iters = -1;
	double s_max = -1.0, s_min = -1.0;
	double solve_time = -1.0;
	KSPConvergedReason reason;
	PetscDrawLG lgctx;
	std::vector<double> rrnorm;
	std::string print_reason;
	int na = -1;

	~KSP_data() { VecDestroy(&x0); };
	
};


/*! \brief Data structure for tridiagonal matrices and their LU factorization via Lapacke.*/
struct MatTridiag
{
	
	// Inputs
	int n;
	Vec dl = PETSC_NULL;
	Vec d = PETSC_NULL;
	Vec du = PETSC_NULL;

	// Factorization output
	Vec du2 = PETSC_NULL;
	std::vector<int> ipiv;
	bool factored = false;

	// Destructor
	~MatTridiag()
	{
		VecDestroy(&dl);
		VecDestroy(&d);
		VecDestroy(&du);
		VecDestroy(&du2);
	}

};


void PrintMatToFile(Mat Matrix, std::complex<double> scale, std::string file_name, bool verbose = true);
void PrintVecToFile(Vec Vector, std::string file_name, bool verbose = true);
void PrintVecToFile_Binary(Vec Vector, std::string file_name, bool verbose = true);
void LoadVecFromFile_Binary(Vec &Vector, std::string file_name, bool verbose = true);

void VecCreateAndAllocate(Vec &vec, int N, std::complex<double> val = 0.0);
void VecCreateAndAllocate(Vec *vec, int N, std::complex<double> val = 0.0);

void VecCreateFromCPP(Vec &out, std::vector<std::complex<double>> &in, int rows);
void VecCreateFromCPP(Vec &out, std::complex<double> *in, int rows);

void VecConcatenate(std::vector<Vec> &vecs, Vec &out);
	
void MatCreateAndAllocate_MATAIJ(Mat &mat, int rows, int cols, int nnz_const, int *nnz_per_row = NULL);
void MatCreateAndAllocate_MATDENSE(Mat &mat, int rows, int cols);
void MatCreateAndAllocate_MATTRIDIAG(MatTridiag &mat, int N, Vec &dl, Vec &d, Vec &du);

void MatCreate_MATDENSE(Mat &mat);
void MatAllocate_MATDENSE(Mat &mat, int rows, int cols);

void MatIdentity_MATDENSE(Mat &mat, int size);
void MatIdentity_MATAIJ(Mat &mat, int size);

void MatCreateFromCPP_MATDENSE(Mat &out, std::vector<std::complex<double>> &in, int rows, int cols, int ordering = PETSC_ROW_MAJOR);
void MatCreateFromCPP_MATDENSE(Mat &out, std::complex<double> *in, int rows, int cols, int ordering = PETSC_ROW_MAJOR);

void MatGetNNZ(Mat &mat, std::vector<int> &nnz);
void MatGetTridiagonal(Mat &mat, int n, Vec &dl, Vec &d, Vec &du);

void MatApplySparsityWindow(Mat A, Mat refmat, int nrows, int ncols, Mat *outmat);
void MatApplyWindowGivenNNZ(Mat A, int *nnz, std::vector<std::vector<int>> nnz_cols, int nrows, int ncols, Mat *outmat);
void MatApplyDistanceWindow_Edges(Mat &in, Mat &out, double threshold);
void MatApplyDistanceWindow_Triangles(Mat &in, Mat &out, double threshold);

void MatChopSparse(Mat A, double tol, int nrows, int ncols, Mat *outmat);

void MatScaleIfExists(Mat &mat, std::complex<double> val);
void MatAXPYIfExists(Mat &Y, std::complex<double> a, Mat &X);

void MatExtractDiagonal(Mat &in, Mat &out);
void MatExtractDiagonal(Mat &mat);
void MatExtractBlockDiagonal(Mat &in, Mat &out, int block_size);

void MatInvertDiagonal(Mat &in, Mat &out);
void MatInvertDiagonal(Mat &in);
void MatInvertDiagonal(Mat &in, Vec &out);

void MatInvert2x2BlockDiagonalMat(Mat &block00, Mat &block01, Mat &block10, Mat &block11);

void MatConcatenate(std::vector<Mat *> submatrices, int rows, int cols, Mat &out, std::vector<std::complex<double>> *a = NULL);
void MatAppend(Mat &mat, std::vector<Mat *> submatrices, int rows, int cols);

void MatFactorizeTridiag_Lapacke(MatTridiag &mat);
void MatSolveTridiag_Lapacke(MatTridiag &mat, Vec &x, Vec &b);

void MatInvert_SuperLU(Mat &mat, Mat &out);
void MatInvert_SuperLU_dist(Mat &mat, Mat &out);
void MatInvert_Petsc(Mat &mat, Mat &out);
void MatInvert_Lapacke(Mat &mat, Mat &out);

void ColVecsToMat_MATDENSE(Mat &mat, std::vector<Vec> &vecs);
void ColVecsToMat_MATAIJ(Mat &mat, std::vector<Vec> &vecs);
void MatToColVecs(Mat &mat, std::vector<Vec> &vecs);

void MatMatMult_IndexPairs(Mat &A, Mat &B, Mat &X, std::vector<std::set<int>> &idx, double *dt = NULL);

void MatFactorize_SuperLU(Mat &mat, Mat &mat_factored, double *dt = NULL, bool incomplete = false, int fill_factor = 1);
void MatFactorize_SuperLU_dist(Mat &mat, Mat &mat_factored, double *dt = NULL);
void MatFactorize_Petsc(Mat &mat, Mat &mat_factored, double *dt = NULL, bool incomplete = false, int fill_level = 1);

void MatMatSolve_SuperLU(Mat &A, Mat &X, Mat &B, bool factored = false, double *dt = NULL);
void MatMatSolve_SuperLU(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored = false, double *dt = NULL);

void MatMatSolve_SuperLU_dist(Mat &A, Mat &X, Mat &B, bool factored = false, double *dt = NULL);
void MatMatSolve_SuperLU_dist(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored = false, double *dt = NULL);

void MatMatSolve_Petsc(Mat &A, Mat &X, Mat &B, bool factored = false, double *dt = NULL);
void MatMatSolve_Petsc(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored = false, double *dt = NULL);

void MatMatSolve_Lapacke(Mat &A, Mat &X, Mat &B, double *dt = NULL);
void MatMatSolve_Lapacke(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, double *dt = NULL);

void MatMatSolve_IndexPairs(Mat &A, Mat &X, Mat &B, std::vector<std::set<int>> &idx, bool factored = false, double *dt = NULL);

void MatSolve_SuperLU(Mat &A, Vec &x, Vec &b, bool factored = false, bool allocate_x = true, double *dt = NULL);
void MatSolve_SuperLU_dist(Mat &A, Vec &x, Vec &b, bool factored = false, bool allocate_x = true, double *dt = NULL);
void MatSolve_Petsc(Mat &A, Vec &x, Vec &b, bool factored = false, bool allocate_x = true, double *dt = NULL);
void MatSolve_Lapacke(Mat &A, Vec &x, Vec &b, double *dt = NULL);

PetscErrorCode KSPStoreResidual(KSP ksp, int n, PetscReal rnorm, void *mctx);
void KSPInitialize_GMRES(KSP &ksp, Mat &A, KSP_data *setup = NULL);
void MatSolve_GMRES(Mat &A, Vec &x, Vec &b, KSP_data *setup = NULL, KSP *ksp_in = NULL);

void MatCreateFromShell(Mat &in, Mat &out);
void PCCreateFromShell(PC &pc, Mat &mat);

void MatSVD(Mat &mat, std::vector<double> &s);
void MatCondRank(Mat &mat, double *cond = NULL, int *rank = NULL, double *rank_rel = NULL, bool verbose = true);
void MatEig(Mat &mat, std::vector<std::complex<double>> &w);

void MatMatMask(Mat &A, Mat &B);
void MatMatPointwiseMult(Mat &A, Mat &B, Mat &C);
void MatAbs(Mat &A);

void MatEliminateRows(Mat &mat, std::vector<int> &rows);
void MatEliminateCols(Mat &mat, std::vector<int> &cols);
void MatInsertRows(Mat &mat, std::vector<int> &rows);
void MatInsertCols(Mat &mat, std::vector<int> &cols);

void MatCreateRowCopyMatrix(Mat &mat, int size, std::vector<int> &source_rows, std::vector<int> &target_rows);
void MatCreateRowAdditionMatrix(Mat &mat, int size, std::vector<int> &source_rows, std::vector<int> &target_rows, std::complex<double> val);

void MatMatMult_OverwriteA(Mat &A, Mat &B);
void MatMatTransposeMult_OverwriteA(Mat &A, Mat &B);
void MatMatMult_OverwriteB(Mat &A, Mat &B);
void MatTransposeMatMult_OverwriteB(Mat &A, Mat &B);
void MatMatTransposeMult_OverwriteB(Mat &A, Mat &B);
void MatPtAP_OverwriteA(Mat &A, Mat &P);
void MatBtAC_OverwriteA(Mat &A, Mat &B, Mat &C);
void MatRtAP(Mat &A, Mat &R, Mat &P, Mat &C);


std::string MatInfoToString(MatInfo &info);
Mat DuplicateTypeAndSize(Mat &mat);

void SaveDistMatToASCII(Mat Matrix, const std::string &file_name);
void SaveDistMatToBinary(Mat Matrix, const std::string &file_name);
void LoadDistMatFromBinary(Mat Matrix, const std::string &file_name);
void SaveDistVecToBinary(Vec vector, const std::string &file_name);
void LoadDistVecFromBinary(Vec vector, const std::string &file_name);

void Vec_Init_Dist(Vec &vec, int N);


// ====== Unit tests ======

void Test_CreateA(Mat &A);
void Test_MatInsertRows(std::string debug_folder = "./src/unit_tests/debug/");
void Test_MatInsertCols(std::string debug_folder = "./src/unit_tests/debug/");


#endif
