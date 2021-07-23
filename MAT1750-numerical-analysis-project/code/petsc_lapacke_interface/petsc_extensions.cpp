

#include <algorithm>
#include <stdexcept>

#include <lapacke.h>

#include "petsc_extensions.h"


void PrintMatToFile(Mat Matrix, std::complex<double> scale, std::string file_name, bool verbose)
{

	if (verbose)
		std::cout << "[EXPORT] Printing matrix to file..." << std::flush;

	int rows, cols;
	MatGetSize(Matrix, &rows, &cols);
	
	Mat temp = PETSC_NULL;
	MatDuplicate(Matrix, MAT_COPY_VALUES, &temp);
	
	MatScale(temp, scale);
	PetscViewer viewer;
	PetscViewerASCIIOpen(MPI_COMM_SELF, file_name.c_str(), &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	MatView(temp, viewer);

	MatDestroy(&temp);
	PetscViewerPopFormat(viewer);
	PetscViewerDestroy(&viewer);

	if (verbose)
		std::cout << "done." << std::endl;
	
	return;
}


void PrintVecToFile(Vec Vector, std::string file_name, bool verbose)
{

	if (verbose)
		std::cout << "[EXPORT] Printing vector to file..." << std::flush;
	
	PetscViewer viewer;
	PetscViewerASCIIOpen(MPI_COMM_SELF, file_name.c_str(), &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	VecView(Vector, viewer);
	PetscViewerPopFormat(viewer);
	PetscViewerDestroy(&viewer);

	if (verbose)
		std::cout << "done." << std::endl;

	return;
}


void PrintVecToFile_Binary(Vec Vector, std::string file_name, bool verbose)
{

	if (verbose)
		std::cout << "[EXPORT] Printing vector to file..." << std::flush;
	
	PetscViewer viewer;
	PetscViewerBinaryOpen(MPI_COMM_SELF, file_name.c_str(), FILE_MODE_WRITE, &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
	VecView(Vector, viewer);
	PetscViewerPopFormat(viewer);
	PetscViewerDestroy(&viewer);

	if (verbose)
		std::cout << "done." << std::endl;

	return;
}


void LoadVecFromFile_Binary(Vec &Vector, std::string file_name, bool verbose)
{

	if (verbose)
		std::cout << "[IMPORT] Importing vector from file..." << std::flush;

	PetscViewer viewer;
	PetscViewerBinaryOpen(MPI_COMM_SELF, file_name.c_str(), FILE_MODE_READ, &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_BINARY_MATLAB);
	VecLoad(Vector, viewer);
	PetscViewerPopFormat(viewer);
	PetscViewerDestroy(&viewer);

	if (verbose)
		std::cout << "done." << std::endl;

	return;
}


/*! \brief Allocate memory for a sequential Petsc vector. Vector entries are initialized to zeros.*/
void VecCreateAndAllocate(Vec &vec, int N, std::complex<double> val)
{

	if ((bool)vec)
		VecDestroy(&vec);
	VecCreate(MPI_COMM_SELF, &vec);
	VecSetType(vec, VECSEQ);
	VecSetSizes(vec, PETSC_DECIDE, N);
	VecSet(vec, val);
	VecAssemblyBegin(vec);
	VecAssemblyEnd(vec);

	return;

}


/*! \brief Allocate memory for a sequential Petsc vector. Vector entries are initialized to zeros.*/
void VecCreateAndAllocate(Vec *vec, int N, std::complex<double> val)
{

	if ((bool)(*vec))
		VecDestroy(vec);
	VecCreate(MPI_COMM_SELF, vec);
	VecSetType(*vec, VECSEQ);
	VecSetSizes(*vec, PETSC_DECIDE, N);
	VecSet(*vec, val);
	VecAssemblyBegin(*vec);
	VecAssemblyEnd(*vec);

	return;

}


/*! \brief Create a Petsc vector from a given std::vector object.*/
void VecCreateFromCPP(Vec &out, std::vector<std::complex<double>> &in, int rows)
{
	
	VecCreateAndAllocate(out, rows);

	std::vector<int> ix (rows);
	for (int ii = 0; ii < rows; ii++)
		ix[ii] = ii;
		
	VecSetValues(out, rows, &(ix[0]), &(in[0]), INSERT_VALUES);

	VecAssemblyBegin(out);
	VecAssemblyEnd(out);

	return;

}


/*! \brief Create a Petsc vector from a given raw C array.*/
void VecCreateFromCPP(Vec &out, std::complex<double> *in, int rows)
{
	
	VecCreateAndAllocate(out, rows);

	std::vector<int> ix (rows);
	for (int ii = 0; ii < rows; ii++)
		ix[ii] = ii;
		
	VecSetValues(out, rows, &(ix[0]), in, INSERT_VALUES);

	VecAssemblyBegin(out);
	VecAssemblyEnd(out);

	return;

}


/*! \brief Concatenate a series of given Petsc vectors.*/
void VecConcatenate(std::vector<Vec> &vecs, Vec &out)
{

	int size = 0;
	std::vector<std::vector<int>> idx (vecs.size()), idxl (vecs.size());
	std::vector<int> sizes (vecs.size());

	int kk = 0;
	for (std::size_t ii = 0; ii < vecs.size(); ii++)
	{
		if (vecs[ii] != PETSC_NULL)
		{		
			VecGetSize(vecs[ii], &sizes[ii]);
			size += sizes[ii];

			for (int jj = 0; jj < sizes[ii]; jj++)
			{
				idx[ii].push_back(kk);
				idxl[ii].push_back(jj);
				kk++;
			}
		}
	}

	VecCreateAndAllocate(out, size);

	for (std::size_t ii = 0; ii < vecs.size(); ii++)
	{
		if (vecs[ii] != PETSC_NULL)
		{
			std::complex<double> *vals;
			VecGetArray(vecs[ii], &vals);
			VecSetValues(out, sizes[ii], &idx[ii][0], vals, INSERT_VALUES);
			VecRestoreArray(vecs[ii], &vals);
		}
	}

	VecAssemblyBegin(out);
	VecAssemblyEnd(out);

	return;

}


/*! \brief Allocate memory for a sequential AIJ Petsc matrix. Matrix entries are initialized to zeros.*/
void MatCreateAndAllocate_MATAIJ(Mat &mat, int rows, int cols, int nnz_const, int *nnz_per_row)
{

	if ((bool)mat)
		MatDestroy(&mat);
	MatCreate(MPI_COMM_SELF, &mat);
	MatSetType(mat, MATAIJ);
	MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
	MatSeqAIJSetPreallocation(mat, nnz_const, nnz_per_row);
	MatZeroEntries(mat);

	return;

}


/*! \brief Allocate memory for a dense Petsc matrix. Matrix entries are initialized to zeros.*/
void MatCreateAndAllocate_MATDENSE(Mat &mat, int rows, int cols)
{

	if ((bool)mat)
		MatDestroy(&mat);
	MatCreate(MPI_COMM_SELF, &mat);
	MatSetType(mat, MATDENSE);
	MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
	MatSetUp(mat);
	MatZeroEntries(mat);

	return;

}


/*! \brief Create a tridiagonal matrix with given vectors of diagonals.*/
void MatCreateAndAllocate_MATTRIDIAG(MatTridiag &mat, int n, Vec &dl, Vec &d, Vec &du)
{

	mat.n = n;

	VecDuplicate(dl, &mat.dl);
	VecCopy(dl, mat.dl);

	VecDuplicate(d, &mat.d);
	VecCopy(d, mat.d);

	VecDuplicate(du, &mat.du);
	VecCopy(du, mat.du);

	mat.factored = false;

	return;

}


/*! \brief Create, but don't allocate, a dense Petsc matrix.*/
void MatCreate_MATDENSE(Mat &mat)
{

	if ((bool)mat)
		MatDestroy(&mat);
	MatCreate(MPI_COMM_SELF, &mat);
	MatSetType(mat, MATDENSE);

	return;

}


/*! \brief Allocate memory for a pre-created dense Petsc matrix. Matrix entries are initialized to zeros.*/
void MatAllocate_MATDENSE(Mat &mat, int rows, int cols)
{
	
	MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
	MatSetUp(mat);
	MatZeroEntries(mat);

	return;

}


/*! \brief Create a dense square identity matrix.*/
void MatIdentity_MATDENSE(Mat &mat, int size)
{

	MatCreateAndAllocate_MATDENSE(mat, size, size);

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	MatShift(mat, 1.0);

	return;

}


/*! \brief Create a sparse square identity matrix.*/
void MatIdentity_MATAIJ(Mat &mat, int size)
{

	MatCreateAndAllocate_MATAIJ(mat, size, size, 1, NULL);

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	MatShift(mat, 1.0);
	
	return;

}


/*! \brief Create a dense Petsc matrix from a given std::vector object, with respect to given matrix ordering.*/
void MatCreateFromCPP_MATDENSE(Mat &out, std::vector<std::complex<double>> &in, int rows, int cols, int ordering)
{
	
	MatCreateAndAllocate_MATDENSE(out, rows, cols);

	if (ordering == PETSC_ROW_MAJOR)
	{
		for (int ii = 0; ii < rows; ii++)
			for (int jj = 0; jj < cols; jj++)
				MatSetValues(out, 1, &ii, 1, &jj, &(in[ii*cols + jj]), INSERT_VALUES);
	}
	else if (ordering == PETSC_COL_MAJOR)
	{
		for (int jj = 0; jj < cols; jj++)
			for (int ii = 0; ii < rows; ii++)
				MatSetValues(out, 1, &ii, 1, &jj, &(in[ii + rows*jj]), INSERT_VALUES);
	}
	else
	{
		std::cout << "[ERROR] MatCreateFromCPP_MATDENSE(): Invalid matrix ordering requested." << std::endl;
		return;
	}

	MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY);


	return;

}


/*! \brief Create a dense Petsc matrix from a given raw C array, with respect to given matrix ordering.*/
void MatCreateFromCPP_MATDENSE(Mat &out, std::complex<double> *in, int rows, int cols, int ordering)
{
	
	MatCreateAndAllocate_MATDENSE(out, rows, cols);

	if (ordering == PETSC_ROW_MAJOR)
	{
		for (int ii = 0; ii < rows; ii++)
			for (int jj = 0; jj < cols; jj++)
				MatSetValues(out, 1, &ii, 1, &jj, &(in[ii*cols + jj]), INSERT_VALUES);
	}
	else if (ordering == PETSC_COL_MAJOR)
	{
		for (int jj = 0; jj < cols; jj++)
			for (int ii = 0; ii < rows; ii++)
				MatSetValues(out, 1, &ii, 1, &jj, &(in[ii + rows*jj]), INSERT_VALUES);
	}
	else
	{
		std::cout << "[ERROR] MatCreateFromCPP_MATDENSE(): Invalid matrix ordering requested." << std::endl;
		return;
	}

	MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY);


	return;

}


/*! \brief Get a vector containing the number of non-zero entries per row of a given matrix.*/
void MatGetNNZ(Mat &mat, std::vector<int> &nnz)
{

	int N_rows, N_cols;
	MatGetSize(mat, &N_rows, &N_cols);

	nnz.resize(N_rows);

	for (int ii = 0; ii < N_rows; ii++)
	{
		int ncols;
		MatGetRow(mat, ii, &ncols, NULL, NULL);
		nnz[ii] = ncols;
		MatRestoreRow(mat, ii, &ncols, NULL, NULL);
	}

	return;

}


/*! \brief Get vectors containing tridiagonal entries of a given square matrix.*/
void MatGetTridiagonal(Mat &mat, int N, Vec &dl, Vec &d, Vec &du)
{

	VecCreateAndAllocate(dl, N - 1, 0.0);
	VecCreateAndAllocate(d, N, 0.0);
	VecCreateAndAllocate(du, N - 1, 0.0);

	for (int ii = 0; ii < N; ii++)
	{
		
		int idx_dl = ii - 1, idx_d = ii, idx_du = ii + 1;
		std::complex<double> val_dl, val_d, val_du;
		
		if (ii > 0)
		{
			MatGetValues(mat, 1, &ii, 1, &idx_dl, &val_dl);
			VecSetValue(dl, idx_dl, val_dl, INSERT_VALUES);
		}

		MatGetValues(mat, 1, &ii, 1, &idx_d, &val_d);
		VecSetValue(d, idx_d, val_d, INSERT_VALUES);

		if (ii < N - 1)
		{
			MatGetValues(mat, 1, &ii, 1, &idx_du, &val_du);
			VecSetValue(du, ii, val_du, INSERT_VALUES);
		}

	}

	VecAssemblyBegin(dl);
	VecAssemblyEnd(dl);

	VecAssemblyBegin(d);
	VecAssemblyEnd(d);

	VecAssemblyBegin(du);
	VecAssemblyEnd(du);

	return;
	
}


/*! \brief Concatenate a series of (optionally scaled) matrices and assemble them into a large matrix. Matrices must be ordered in a row-major indexing scheme. The matrices to be assembled must have the proper sizes. A NULL entry is treated as a block of zeros. There must be at least one non-NULL (but possibly zero) matrix on each row and column of the final matrix.*/
void MatConcatenate(std::vector<Mat *> submatrices, int rows, int cols, Mat &out, std::vector<std::complex<double>> *a)
{

	if ((int)submatrices.size() != rows*cols)
		throw std::logic_error("[ERROR] MatConcatenate(): Number of submatrices does not match rows*cols.");

	if (a != NULL && (int)a->size() != rows*cols)
		throw std::logic_error("[ERROR] MatConcatenate(): Size of scaling vector does not match rows*cols.");

	Mat nested = PETSC_NULL;
	Mat *submatrices_copied = new Mat [rows*cols];

	// Loop through all indices and process each submatrix
	int nz = 0;
	for (int ii = 0; ii < rows*cols; ii++)
	{
		if (submatrices[ii] == NULL || *submatrices[ii] == PETSC_NULL)
		{
			submatrices_copied[ii] = NULL;
			nz++;
		}
		else
		{
			MatDuplicate(*(submatrices[ii]), MAT_COPY_VALUES, &(submatrices_copied[ii]));
			if (a != NULL)
				MatScale(submatrices_copied[ii], (*a)[ii]);
		}
	}

	// Create the full matrix in nested format
	if (nz != rows*cols)
	{
		MatCreateNest(MPI_COMM_SELF, rows, NULL, cols, NULL, submatrices_copied, &nested);

		// Clean-up
		for (int ii = 0; ii < rows*cols; ii++)
			MatDestroy(&submatrices_copied[ii]);

		// Convert the full matrix to AIJ, as output
		MatConvert(nested, MATAIJ, MAT_INITIAL_MATRIX, &out);
		MatDestroy(&nested);
	}
	else
		out = PETSC_NULL;

	delete [] submatrices_copied;
	
	return;

}


/*! \brief Append a series of matrices to the input matrix. The input matrix will be expanded in size as a result of this, and will be located at the (0,0) block. Matrices must be ordered in a row-major indexing scheme. The matrices to be assembled must have the proper sizes. A NULL entry is treated as a block of zeros. There must be at least one non-NULL (but possibly zero) matrix on each row and column of the final matrix.*/
void MatAppend(Mat &mat, std::vector<Mat *> submatrices, int rows, int cols)
{

	if ((int)submatrices.size() + 1 != rows*cols)
		throw std::logic_error("[ERROR] MatAppend(): Number of submatrices does not match rows*cols.");

	Mat copy = PETSC_NULL;
	MatDuplicate(mat, MAT_COPY_VALUES, &copy);

	submatrices.insert(submatrices.begin(), &copy);

	MatDestroy(&mat);
	MatConcatenate(submatrices, rows, cols, mat);
	MatDestroy(&copy);
	
	return;

}


/*! \brief Extract the diagonal of a given matrix and use it to generate a new matrix.*/
void MatExtractDiagonal(Mat &in, Mat &out)
{

	int nrows, ncols;
	MatGetSize(in, &nrows, &ncols);

	Vec temp = PETSC_NULL;
	VecCreateAndAllocate(temp, nrows, 0.0);

	MatGetDiagonal(in, temp);
	MatCreateAndAllocate_MATAIJ(out, nrows, ncols, 1.0, NULL);

	MatDiagonalSet(out, temp, INSERT_VALUES);
	VecDestroy(&temp);

	return;
	
}


/*! \brief Scale a matrix if it exists, otherwise do nothing.*/
void MatScaleIfExists(Mat &mat, std::complex<double> val)
{
	if (mat != PETSC_NULL)
		MatScale(mat, val);
	return;
}


/*! \brief Apply MatAXPY if both matrices exist, otherwise do nothing.*/
void MatAXPYIfExists(Mat &Y, std::complex<double> a, Mat &X)
{
	if (Y != PETSC_NULL && X != PETSC_NULL)
		MatAXPY(Y, a, X, DIFFERENT_NONZERO_PATTERN);
	return;
}


/*! \brief Extract the diagonal of a given matrix and use it to generate a new matrix in-place.*/
void MatExtractDiagonal(Mat &mat)
{

	int nrows, ncols;
	MatGetSize(mat, &nrows, &ncols);

	Vec temp = PETSC_NULL;
	VecCreateAndAllocate(temp, nrows, 0.0);

	MatGetDiagonal(mat, temp);
	MatDestroy(&mat);
	MatCreateAndAllocate_MATAIJ(mat, nrows, ncols, 1.0, NULL);

	MatDiagonalSet(mat, temp, INSERT_VALUES);
	VecDestroy(&temp);

	return;
	
}


/*! \brief Extract the diagonal elements of a given matrix, and assemble the output matrix as the inverse of that diagonal matrix.*/
void MatInvertDiagonal(Mat &in, Mat &out)
{

	int nrows, ncols;
	MatGetSize(in, &nrows, &ncols);

	Vec temp = PETSC_NULL;
	VecCreateAndAllocate(temp, nrows, 0.0);

	MatGetDiagonal(in, temp);
	VecReciprocal(temp);

	MatCreateAndAllocate_MATAIJ(out, nrows, ncols, 1.0, NULL);
	
	MatDiagonalSet(out, temp, INSERT_VALUES);
	VecDestroy(&temp);	

	return;
}


/*! \brief Extract the diagonal elements of a given matrix, and assemble the output matrix as the inverse of that diagonal matrix ***in place***.*/
void MatInvertDiagonal(Mat &in)
{

	int nrows, ncols;
	MatGetSize(in, &nrows, &ncols);

	Vec temp = PETSC_NULL;
	VecCreateAndAllocate(temp, nrows, 0.0);

	MatGetDiagonal(in, temp);
	VecReciprocal(temp);

	MatScale(in, 0.0);
	MatDiagonalSet(in, temp, INSERT_VALUES);
	VecDestroy(&temp);	

	return;
}


/*! \brief Extract the diagonal elements of a given matrix, and assemble the output vector as the inverse of that diagonal matrix.*/
void MatInvertDiagonal(Mat &in, Vec &out)
{

	int nrows, ncols;
	MatGetSize(in, &nrows, &ncols);

	VecCreateAndAllocate(out, nrows, 0.0);

	MatGetDiagonal(in, out);
	VecReciprocal(out);

	return;
}


/*! \brief Builds a Petsc matrix from the block diagonal entries of a given Petsc matrix. It is particularly useful for generating preconditioner blocks. Both matrices must be square and have the same size, and the size of the block must be provided as an input. If the size of the matrix is not exactly divisible by the size of the block, the last block of the matrix will have a different size.*/
void MatExtractBlockDiagonal(Mat &in, Mat &out, int block_size)
{
	int nrows, ncols;
	MatGetSize(in, &nrows, &ncols);

	// Validate block size
	int nblocks = nrows/block_size;
	int block_size_last = nrows%block_size;

	
	// Extract and copy over each block
	std::complex<double> val [block_size*block_size];
	for (int ii = 0; ii < nblocks; ii++)
	{
		// Create indices
		int idx [block_size];
		for (int jj = 0; jj < block_size; jj++)
			idx[jj] = jj + ii*block_size;

		MatGetValues(in, block_size, idx, block_size, idx, val);
		MatSetValues(out, block_size, idx, block_size, idx, val, INSERT_VALUES);
	}

	// Extract and copy over the leftover block
	if (block_size_last != 0)
	{
		// Create indices
		int idx [block_size_last];
		for (int jj = 0; jj < block_size_last; jj++)
			idx[jj] = jj + nblocks*block_size;

		MatGetValues(in, block_size_last, idx, block_size_last, idx, val);
		MatSetValues(out, block_size_last, idx, block_size_last, idx, val, INSERT_VALUES);
	}
	
	MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY);

	
	return;
}


/*! \brief Invert a 2x2 block matrix, where each block is a diagonal matrix. The blocks have to be provided as separate inputs. The output is provided ***in-place*** as separate blocks. NOTE: this function will obviously change the values in the input matrices! This function has been validated.*/
void MatInvert2x2BlockDiagonalMat(Mat &block00, Mat &block01, Mat &block10, Mat &block11)
{

	// Make sure the blocks are diagonal matrices
	MatExtractDiagonal(block00);
	MatExtractDiagonal(block01);
	MatExtractDiagonal(block10);
	MatExtractDiagonal(block11);

	// ------ Invert the (0, 0) block in place ------
	
	int nrows00, ncols00;
	MatGetSize(block00, &nrows00, &ncols00);
	
	// Extract diagonal entries of the 00 block
	Vec temp = PETSC_NULL;
	VecCreateAndAllocate(temp, nrows00, 0.0);
	MatGetDiagonal(block00, temp);

	// Invert diagonal entries
	VecReciprocal(temp);

	// Convert back to a matrix
	MatDiagonalSet(block00, temp, INSERT_VALUES);


	// ------ Create the Schur complement of the (0, 0) block ------

	PetscReal fill_ratio = 1.0;  // estimated fill ratio after matrix-matrix multiplication
	Mat temp_mat = PETSC_NULL;
	MatMatMatMult(block10, block00, block01, MAT_INITIAL_MATRIX, fill_ratio, &temp_mat);
	MatAYPX(temp_mat, -1.0, block11, SAME_NONZERO_PATTERN);


	// ------ Invert the Schur complement and place it in the (1, 1) block ------

	MatGetDiagonal(temp_mat, temp);
	VecReciprocal(temp);
	MatDiagonalSet(block11, temp, INSERT_VALUES);
	

	// ------ Assemble the (0, 1) block ------

	MatDestroy(&temp_mat);
	MatMatMatMult(block00, block01, block11, MAT_INITIAL_MATRIX, fill_ratio, &temp_mat);
	MatScale(temp_mat, -1.0);
	MatCopy(temp_mat, block01, SAME_NONZERO_PATTERN);


	// ------ Assemble an intermediate term of the (0, 0) block ------

	Mat block00_temp = PETSC_NULL;
	MatMatMatMult(block01, block10, block00, MAT_INITIAL_MATRIX, fill_ratio, &block00_temp);
	// Note: this is off by a minus sign, which will be fixed below.
	

	// ------ Assemble the (1, 0) block ------

	MatDestroy(&temp_mat);
	MatMatMatMult(block11, block10, block00, MAT_INITIAL_MATRIX, fill_ratio, &temp_mat);
	MatScale(temp_mat, -1.0);
	MatCopy(temp_mat, block10, SAME_NONZERO_PATTERN);


	// ------ Assemble the (0, 0) block ------

	MatAXPY(block00, -1.0, block00_temp, SAME_NONZERO_PATTERN);

	
	// ------ Cleanup ------

	MatDestroy(&temp_mat);
	MatDestroy(&block00_temp);
	VecDestroy(&temp);
	

	return;

}


/*! \brief Apply a window to a matrix, based the sparsity pattern of a given matrix. NOT TESTED - TODO \todo.*/
void MatApplySparsityWindow(Mat A, Mat refmat, int nrows, int ncols, Mat *outmat)
{

	// Extract nonzero structure of the reference matrix, and corresponding values of the input matrix
	int nnzcols;
	const int *cols_nz;
	int NumEntriesPerRow [nrows];

	for (int ii = 0; ii < nrows; ii++)
	{
		MatGetRow(refmat, ii, &nnzcols, NULL, NULL);
		NumEntriesPerRow[ii] = nnzcols;
	}

	// Preallocate memory for the ouput matrix
	MatCreate(MPI_COMM_SELF, outmat);
	MatSetType(*outmat, MATAIJ);
	MatSetSizes(*outmat, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols);
	MatSeqAIJSetPreallocation(*outmat, 0, NumEntriesPerRow);
	MatZeroEntries(*outmat);

	// Populate the windowed parts of the output matrix
	for (int ii = 0; ii < nrows; ii++)
	{
		MatGetRow(refmat, ii, NULL, &cols_nz, NULL);

		int idxm_t[1];
		idxm_t[0] = ii;
		const int *idxm = idxm_t;
		std::complex<double> A_vals[NumEntriesPerRow[ii]];
		MatGetValues(A, 1, idxm, NumEntriesPerRow[ii], cols_nz, A_vals);
		MatSetValues(*outmat, 1, idxm, NumEntriesPerRow[ii], cols_nz, A_vals, ADD_VALUES);
	}

	MatAssemblyBegin(*outmat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*outmat, MAT_FINAL_ASSEMBLY);

	return;
}


/*! \brief Apply a window to a matrix, based on a given sparsity pattern, such that all entries outside the window are zeroed out and compressed away. NOT TESTED - TODO \todo.*/
void MatApplyWindowGivenNNZ(Mat A, int *nnz, std::vector<std::vector<int>> nnz_cols, int nrows, int ncols, Mat *outmat)
{

	// Preallocate memory for the ouput matrix
	MatCreate(MPI_COMM_SELF, outmat);
	MatSetType(*outmat, MATAIJ);
	MatSetSizes(*outmat, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols);
	MatSeqAIJSetPreallocation(*outmat, 0, nnz);
	MatZeroEntries(*outmat);

	// Populate the windowed parts of the output matrix
	for (int ii = 0; ii < nrows; ii++)
	{

		int nnzcols = nnz[ii];
		const int *cols_nz = &nnz_cols[ii][0];

		int idxm_t[1];
		idxm_t[0] = ii;
		const int *idxm = idxm_t;
		std::complex<double> A_vals [nnz[ii]];
		MatGetValues(A, 1, idxm, nnz[ii], cols_nz, A_vals);
		MatSetValues(*outmat, 1, idxm, nnz[ii], cols_nz, A_vals, ADD_VALUES);
	}

	MatAssemblyBegin(*outmat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*outmat, MAT_FINAL_ASSEMBLY);

	return;
}


/*! \brief This function accepts a matrix of mesh edge interactions, and windows it to keep only interactions where the source and test edges fall within a given distance threshold. TODO \todo.*/
void MatApplyDistanceWindow_Edges(Mat &in, Mat &out, double threshold)
{
	// TODO
	
	return;

}


/*! \brief This function accepts a matrix of mesh triangle interactions, and windows it to keep only interactions where the source and test triangles fall within a given distance threshold. TODO \todo.*/
void MatApplyDistanceWindow_Triangles(Mat &in, Mat &out, double threshold)
{
	// TODO

	return;

}


/*! \brief Based on MatChop; removes values in a matrix whose absolute values are below some tolerance. The new matrix has a new sparsity pattern, corresponding to the entries that no longer exist. TODO \todo.*/
void MatChopSparse(Mat A, double tol, int nrows, int ncols, Mat *outmat)
{

}


/*! \brief Factorize a tridiagonal matrix using Lapacke's specialized routines.*/
void MatFactorizeTridiag_Lapacke(MatTridiag &mat)
{

	mat.ipiv.resize(mat.n);
	std::complex<double> *dl, *d, *du, *du2;
	
	VecCreateAndAllocate(mat.du2, mat.n-2, 0.0);

	VecGetArray(mat.dl, &dl);
	VecGetArray(mat.d, &d);
	VecGetArray(mat.du, &du);
	VecGetArray(mat.du2, &du2);

	int info = LAPACKE_zgttrf(mat.n, (lapack_complex_double *) dl, (lapack_complex_double *) d, (lapack_complex_double *) du, (lapack_complex_double *) du2, &mat.ipiv[0]);

	VecRestoreArray(mat.dl, &dl);
	VecRestoreArray(mat.d, &d);
	VecRestoreArray(mat.du, &du);
	VecRestoreArray(mat.du2, &du2);
	
	mat.factored = true;
	
	return;

}


/*! \brief Factorize (if not already factored) and solve a tridiagonal system, Ax = b, using Lapacke's specialized routines.*/
void MatSolveTridiag_Lapacke(MatTridiag &mat, Vec &x, Vec &b)
{

	if (!mat.factored)
		MatFactorizeTridiag_Lapacke(mat);

	std::complex<double> *dl, *d, *du, *du2, *x_ptr;

	// VecDestroy(&x);
	// VecDuplicate(b, &x);
	VecCopy(b, x);
	
	VecGetArray(mat.dl, &dl);
	VecGetArray(mat.d, &d);
	VecGetArray(mat.du, &du);
	VecGetArray(mat.du2, &du2);
	VecGetArray(x, &x_ptr);

	int into = LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', mat.n, 1, (lapack_complex_double *) dl, (lapack_complex_double *) d, (lapack_complex_double *) du, (lapack_complex_double *) du2, &mat.ipiv[0], (lapack_complex_double *) x_ptr, mat.n);

	VecRestoreArray(mat.dl, &dl);
	VecRestoreArray(mat.d, &d);
	VecRestoreArray(mat.du, &du);
	VecRestoreArray(mat.du2, &du2);
	VecRestoreArray(x, &x_ptr);
	
	return;

}


/*! \brief [Testing and debugging only] Invert a given sparse AIJ matrix after factorizing using SuperLU. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format.*/
void MatInvert_SuperLU(Mat &mat, Mat &out)
{
	
	Mat mat_factored = PETSC_NULL;
	MatFactorize_SuperLU(mat, mat_factored);

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	Mat RHS = PETSC_NULL;
	MatCreateAndAllocate_MATDENSE(RHS, rows, cols);
	MatShift(RHS, 1.0);

	MatCreateAndAllocate_MATDENSE(out, rows, cols);
	MatMatSolve(mat_factored, RHS, out);
	// MatConvert(out, MATAIJ, MAT_INPLACE_MATRIX, &out);

	MatDestroy(&RHS);
	MatDestroy(&mat_factored);

	return;

}


/*! \brief [Testing and debugging only] Invert a given sparse AIJ matrix after factorizing using SuperLU_dist. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format.*/
void MatInvert_SuperLU_dist(Mat &mat, Mat &out)
{
	
	Mat mat_factored = PETSC_NULL;
	MatFactorize_SuperLU_dist(mat, mat_factored);

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	Mat RHS = PETSC_NULL;
	MatCreateAndAllocate_MATDENSE(RHS, rows, cols);
	MatShift(RHS, 1.0);

	MatCreateAndAllocate_MATDENSE(out, rows, cols);
	MatMatSolve(mat_factored, RHS, out);
	// MatConvert(out, MATAIJ, MAT_INPLACE_MATRIX, &out);

	MatDestroy(&RHS);
	MatDestroy(&mat_factored);

	return;

}


/*! \brief [Testing and debugging only] Inverts a given sparse AIJ matrix after factorizing using Petsc's inbuilt factorization routines. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format.*/
void MatInvert_Petsc(Mat &mat, Mat &out)
{
	
	Mat mat_factored = PETSC_NULL;
	MatFactorize_Petsc(mat, mat_factored);

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	Mat RHS = PETSC_NULL;
	MatCreateAndAllocate_MATDENSE(RHS, rows, cols);
	MatShift(RHS, 1.0);

	MatCreateAndAllocate_MATDENSE(out, rows, cols);
	MatMatSolve(mat_factored, RHS, out);
	// MatConvert(out, MATAIJ, MAT_INPLACE_MATRIX, &out);

	MatDestroy(&RHS);
	MatDestroy(&mat_factored);

	return;

}


/*! \brief [Testing and debugging only] Invert a given dense matrix after factorizing using lapacke. The output matrix is returned in dense format because the inverse in general is expected to be dense.*/
void MatInvert_Lapacke(Mat &mat, Mat &out)
{
	
	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	Mat RHS = PETSC_NULL;
	MatCreateAndAllocate_MATDENSE(RHS, rows, cols);
	MatShift(RHS, 1.0);

	MatMatSolve_Lapacke(mat, out, RHS);

	MatDestroy(&RHS);

	return;

}


/*! \brief Concatenate a C++ vector of Petsc Vecs into a Petsc Mat, column-wise.*/
void ColVecsToMat_MATDENSE(Mat &mat, std::vector<Vec> &vecs)
{

	int rows;
	VecGetSize(vecs[0], &rows);
	int cols = vecs.size();
	MatCreateAndAllocate_MATDENSE(mat, rows, cols);

	std::vector<int> idxn (rows);
	for (int ii = 0; ii < rows; ii++)
		idxn[ii] = ii;

	for (int ii = 0; ii < cols; ii++)
	{
		std::complex<double> *vals;
		VecGetArray(vecs[ii], &vals);
		MatSetValues(mat, rows, &idxn[0], 1, &ii, vals, INSERT_VALUES);
		VecRestoreArray(vecs[ii], &vals);
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	return;

}


/*! \brief Concatenate a C++ vector of Petsc Vecs into a Petsc Mat, column-wise.*/
void ColVecsToMat_MATAIJ(Mat &mat, std::vector<Vec> &vecs)
{

	int rows;
	VecGetSize(vecs[0], &rows);
	int cols = vecs.size();
	MatCreateAndAllocate_MATAIJ(mat, rows, cols, cols, NULL);

	std::vector<int> idxn (rows);
	for (int ii = 0; ii < rows; ii++)
		idxn[ii] = ii;

	for (int ii = 0; ii < cols; ii++)
	{
		std::complex<double> *vals;
		VecGetArray(vecs[ii], &vals);
		MatSetValues(mat, rows, &idxn[0], 1, &ii, vals, INSERT_VALUES);
		VecRestoreArray(vecs[ii], &vals);
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	return;

}


/*! \brief Extract columns from a Petsc Mat and create a C++ vector of Petsc Vecs.*/
void MatToColVecs(Mat &mat, std::vector<Vec> &vecs)
{

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	std::vector<int> idxn (rows);
	for (int ii = 0; ii < rows; ii++)
		idxn[ii] = ii;

	vecs.clear();
	vecs.resize(cols);
	for (int ii = 0; ii < cols; ii++)
	{
		std::vector<std::complex<double>> vals (rows);
		MatGetValues(mat, rows, &idxn[0], 1, &ii, &vals[0]);

		vecs[ii] = PETSC_NULL;
		VecCreateAndAllocate(vecs[ii], rows);
		VecSetValues(vecs[ii], rows, &idxn[0], &vals[0], INSERT_VALUES);

		VecAssemblyBegin(vecs[ii]);
		VecAssemblyEnd(vecs[ii]);
	}

	return;

}


/*! \brief Computes specific elements of A*B without actually constructing the matrix. Note that idx corresponds to the desired indices of matrix X.*/
void MatMatMult_IndexPairs(Mat &A, Mat &B, Mat &X, std::vector<std::set<int>> &idx, double *dt)
{

	std::clock_t t1 = clock();

	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	int rows_B, cols_B;
	MatGetSize(B, &rows_B, &cols_B);

	// Allocate the output matrix
	std::vector<int> nnz (rows_B, 0);
	for (std::size_t ii = 0; ii < idx.size(); ii++)
		nnz[ii] = idx[ii].size();
	MatCreateAndAllocate_MATAIJ(X, rows_A, cols_B, 0, &nnz[0]);

	Mat BT;
	MatTranspose(B, MAT_INITIAL_MATRIX, &BT);

	// Iterate over target rows
	std::set<int>::iterator it;
	for (int ii = 0; ii < rows_A; ii++)
	{
		
		// Compute the row solution vector
		const std::complex<double> *vals_ii;
		const int *cols_ii;
		int ncols_ii;
		MatGetRow(A, ii, &ncols_ii, &cols_ii, &vals_ii);
		
		// Iterate over source rows
		for (it = idx[ii].begin(); it != idx[ii].end(); ++it)
		{
			int jj = *it;

			// Compute the column solution vector
			const std::complex<double> *vals_jj;
			const int *cols_jj;
			int ncols_jj;
			MatGetRow(BT, jj, &ncols_jj, &cols_jj, &vals_jj);

			std::complex<double> val = 0.0;
			for (int kk = 0; kk < ncols_jj; kk++)
			{
				int idx = std::distance(cols_ii, std::find(cols_ii, cols_ii + ncols_ii, cols_jj[kk]));
				if (idx < ncols_ii)
					val += vals_jj[kk]*vals_ii[idx];
			}

			MatRestoreRow(BT, jj, &ncols_jj, &cols_jj, &vals_jj);
			MatSetValue(X, ii, jj, val, INSERT_VALUES);
		}

		MatRestoreRow(A, ii, &ncols_ii, &cols_ii, &vals_ii);
		
	}
	
	MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY);

	
	// Clean-up
	MatDestroy(&BT);

	
	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


// ==================================================================================
// MatFactorize
// ==================================================================================

/*! \brief Factorize a given sparse AIJ matrix using SuperLU. Incomplete factorization is performed if "incomplete" is set to "true".*/
void MatFactorize_SuperLU(Mat &mat, Mat &mat_factored, double *dt, bool incomplete, int fill_factor)
{

	std::clock_t t1 = clock();

	MatDestroy(&mat_factored);

	MatFactorInfo info;
	MatFactorInfoInitialize(&info);

	IS row, col;
	MatGetOrdering(mat, MATORDERINGNATURAL, &row, &col);
	// Note: SuperLU does not care about the type of ordering set here

	MatFactorType ftype;
	
	if (incomplete)
	{
		ftype = MAT_FACTOR_ILU;
		PetscOptionsSetValue(NULL, "-mat_superlu_ilu_fillfactor", std::to_string(fill_factor).c_str());
	}
	else
		ftype = MAT_FACTOR_LU;
	
	MatGetFactor(mat, MATSOLVERSUPERLU, ftype, &mat_factored);
	
	MatLUFactorSymbolic(mat_factored, mat, row, col, &info);
	MatLUFactorNumeric(mat_factored, mat, &info);

	ISDestroy(&row);
	ISDestroy(&col);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;

	return;

}


/*! \brief Factorize a given sparse AIJ matrix using SuperLU_dist.*/
void MatFactorize_SuperLU_dist(Mat &mat, Mat &mat_factored, double *dt)
{

	std::clock_t t1 = clock();

	PC pc;
	
	PCCreate(MPI_COMM_SELF, &pc);
	PCSetOperators(pc, mat, mat);
	PCSetType(pc, PCLU);

	// For Petsc 3.7.6
	// PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU_DIST);
	// PCFactorSetUpMatSolverPackage(pc);

	// For Petsc 3.9.1
	PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);
	PCFactorSetUpMatSolverType(pc);

	PCFactorGetMatrix(pc, &mat_factored);
	PCSetUp(pc);

	// PCDestroy(&pc);
	
	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


// /*! \brief This function factorizes a given sparse AIJ matrix using Petsc's inbuilt factorization.*/
// void MatFactorize_Petsc(Mat &mat, Mat &mat_factored, double *dt, bool incomplete, int fill_level)
// {

// 	std::clock_t t1 = clock();
	
// 	PC pc;
	
// 	PCCreate(MPI_COMM_SELF, &pc);
// 	PCSetOperators(pc, mat, mat);

// 	PCType type;
// 	if (incomplete)
// 	{
// 		type = PCILU;
// 		PCFactorSetLevels(pc, fill_level);
// 	}
// 	else
// 		type = PCLU;

// 	PCSetType(pc, type);

// 	// For Petsc 3.7.6
// 	// PCFactorSetMatSolverPackage(pc, MATSOLVERPETSC);
// 	// PCFactorSetUpMatSolverPackage(pc);

// 	// For Petsc 3.9.1
// 	PCFactorSetMatSolverType(pc, MATSOLVERPETSC);
// 	PCFactorSetUpMatSolverType(pc);

// 	PCFactorGetMatrix(pc, &mat_factored);	
// 	PCSetUp(pc);

// 	// PCDestroy(&pc);

// 	std::clock_t t2 = clock();
// 	if (dt != NULL)
// 		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
// 	return;

// }


/*! \brief Factorizes a given sparse AIJ matrix using Petsc's inbuilt factorization.*/
void MatFactorize_Petsc(Mat &mat, Mat &mat_factored, double *dt, bool incomplete, int fill_level)
{

	std::clock_t t1 = clock();

	MatDestroy(&mat_factored);

	MatFactorInfo info;
	MatFactorInfoInitialize(&info);

	info.dtcol = 1;
	info.diagonal_fill = 1;

	IS row, col;
	MatGetOrdering(mat, MATORDERINGND, &row, &col); // This is what Petsc uses in PCLU
	// MatGetOrdering(mat, MATORDERINGNATURAL, &row, &col);

	MatFactorType ftype = MAT_FACTOR_LU;
	
	if (incomplete)
	{
		info.levels = fill_level;
		ftype = MAT_FACTOR_ILU;
	}
	
	MatGetFactor(mat, MATSOLVERPETSC, ftype, &mat_factored);
	
	MatLUFactorSymbolic(mat_factored, mat, row, col, &info);
	MatLUFactorNumeric(mat_factored, mat, &info);

	ISDestroy(&row);
	ISDestroy(&col);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;

	return;

}


// ==================================================================================
// MatMatSolve
// ==================================================================================

/*! \brief Solve AX = B for X, after factorizing A using SuperLU. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format. X should *not* have already been allocated.*/
void MatMatSolve_SuperLU(Mat &A, Mat &X, Mat &B, bool factored, double *dt)
{
	
	std::clock_t t1 = clock();

	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	int rows_B, cols_B;
	MatGetSize(B, &rows_B, &cols_B);
	
	MatCreateAndAllocate_MATDENSE(X, rows_A, cols_B);

	MatConvert(B, MATDENSE, MAT_INPLACE_MATRIX, &B);

	if (factored)
		MatMatSolve(A, B, X);
	else
	{
		Mat A_factored;
		MatFactorize_SuperLU(A, A_factored);
		MatMatSolve(A_factored, B, X);
		MatDestroy(&A_factored);
	}

	MatConvert(B, MATAIJ, MAT_INPLACE_MATRIX, &B);
	// MatConvert(X, MATAIJ, MAT_INPLACE_MATRIX, &X);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;

	return;

}


/*! \brief Solve AX = B for X, after factorizing A using SuperLU. Multiple right-hand sides are provided as an std vector of Petsc Vecs, and will get concatenated to solve simultaneously. The output will be extracted into columns of an std vector as well.*/
void MatMatSolve_SuperLU(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored, double *dt)
{
	
	std::clock_t t1 = clock();
 
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	Mat X_Mat = PETSC_NULL, B_Mat = PETSC_NULL;

	ColVecsToMat_MATDENSE(B_Mat, B);
	MatCreateAndAllocate_MATDENSE(X_Mat, rows, B.size());

	MatAssemblyBegin(X_Mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(X_Mat, MAT_FINAL_ASSEMBLY);

	if (factored)
		MatMatSolve(A, B_Mat, X_Mat);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_SuperLU(A, A_factored);
		MatMatSolve(A_factored, B_Mat, X_Mat);
		MatDestroy(&A_factored);
	}

	MatToColVecs(X_Mat, X);

	MatDestroy(&B_Mat);
	MatDestroy(&X_Mat);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Solve AX = B for X, after factorizing A using SuperLU_dist. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format. X should *not* have already been allocated.*/
void MatMatSolve_SuperLU_dist(Mat &A, Mat &X, Mat &B, bool factored, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);
	MatCreateAndAllocate_MATDENSE(X, rows, cols);

	MatConvert(B, MATDENSE, MAT_INPLACE_MATRIX, &B);

	if (factored)
		MatMatSolve(A, B, X);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_SuperLU_dist(A, A_factored);
		MatMatSolve(A_factored, B, X);
		MatDestroy(&A_factored);
	}

	MatConvert(B, MATAIJ, MAT_INPLACE_MATRIX, &B);
	// MatConvert(X, MATAIJ, MAT_INPLACE_MATRIX, &X);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Solve AX = B for X, after factorizing A using SuperLU_dist. Multiple right-hand sides are provided as an std vector of Petsc Vecs, and will get concatenated to solve simultaneously. The output will be extracted into columns of an std vector as well.*/
void MatMatSolve_SuperLU_dist(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	Mat X_Mat = PETSC_NULL, B_Mat = PETSC_NULL;

	ColVecsToMat_MATDENSE(B_Mat, B);
	MatCreateAndAllocate_MATDENSE(X_Mat, rows, B.size());

	MatAssemblyBegin(X_Mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(X_Mat, MAT_FINAL_ASSEMBLY);

	if (factored)
		MatMatSolve(A, B_Mat, X_Mat);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_SuperLU_dist(A, A_factored);
		MatMatSolve(A_factored, B_Mat, X_Mat);
		MatDestroy(&A_factored);
	}

	MatToColVecs(X_Mat, X);

	MatDestroy(&B_Mat);
	MatDestroy(&X_Mat);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Solves AX = B for X, after factorizing A using PETSc's default routines. The output matrix is returned in dense format because the inverse in general is expected to be dense. Uncomment the line near the end of this function to make it return in sparse AIJ format. X should *not* have already been allocated.*/
void MatMatSolve_Petsc(Mat &A, Mat &X, Mat &B, bool factored, double *dt)
{

	std::clock_t t1 = clock();

	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	int rows_B, cols_B;
	MatGetSize(B, &rows_B, &cols_B);
	
	MatCreateAndAllocate_MATDENSE(X, rows_A, cols_B);

	MatConvert(B, MATDENSE, MAT_INPLACE_MATRIX, &B);

	if (factored)
		MatMatSolve(A, B, X);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_Petsc(A, A_factored);
		MatMatSolve(A_factored, B, X);
		MatDestroy(&A_factored);
	}
	
	MatConvert(B, MATAIJ, MAT_INPLACE_MATRIX, &B);
	// MatConvert(X, MATAIJ, MAT_INPLACE_MATRIX, &X);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Solves AX = B for X, after factorizing A using PETSc's default routines. Multiple right-hand sides are provided as an std vector of Petsc Vecs, and will get concatenated to solve simultaneously. The output will be extracted into columns of an std vector as well.*/
void MatMatSolve_Petsc(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, bool factored, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	Mat X_Mat = PETSC_NULL, B_Mat = PETSC_NULL;

	ColVecsToMat_MATDENSE(B_Mat, B);
	MatCreateAndAllocate_MATDENSE(X_Mat, rows, B.size());

	MatAssemblyBegin(X_Mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(X_Mat, MAT_FINAL_ASSEMBLY);

	if (factored)
		MatMatSolve(A, B_Mat, X_Mat);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_Petsc(A, A_factored);
		MatMatSolve(A_factored, B_Mat, X_Mat);
		MatDestroy(&A_factored);
	}
	
	MatToColVecs(X_Mat, X);

	MatDestroy(&B_Mat);
	MatDestroy(&X_Mat);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}



/*! \brief Solves AX = B for X, using Lapacke's gesv directly. It is assumed that the input matrices are dense, and stored in row-major layout. If the input matrix is not dense, it will be converted to a dense matrix anyway.*/
void MatMatSolve_Lapacke(Mat &A, Mat &X, Mat &B, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	int rows_B, cols_B;
	MatGetSize(B, &rows_B, &cols_B);

	MatCreateAndAllocate_MATDENSE(X, rows_A, cols_B);
	
	// Ensure all matrices are dense
	MatConvert(A, MATDENSE, MAT_INPLACE_MATRIX, &A);
	MatConvert(B, MATDENSE, MAT_INPLACE_MATRIX, &B);

	// Lapacke overwrites the RHS matrix, so copy it over to the intended solution matrix
	MatCopy(B, X, DIFFERENT_NONZERO_PATTERN);
	
	// Extract raw pointers to matrix data
	std::complex<double> *A_data, *X_data;
	MatDenseGetArray(A, &A_data);
	MatDenseGetArray(X, &X_data);

	std::vector<int> ipiv (rows_A);
	// int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, rows_A, cols_B, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) X_data, cols_B);
	int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, rows_A, cols_B, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) X_data, rows_B);

	MatDenseRestoreArray(A, &A_data);
	MatDenseRestoreArray(X, &X_data);

	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Solves AX = B for X, using Lapacke's gesv directly. It is assumed that the input matrix is dense, and stored in row-major layout. If the input matrix is not dense, it will be converted to a dense matrix anyway.*/
void MatMatSolve_Lapacke(Mat &A, std::vector<Vec> &X, std::vector<Vec> &B, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	Mat B_Mat = PETSC_NULL;
	ColVecsToMat_MATDENSE(B_Mat, B);

	int rows_B, cols_B;
	MatGetSize(B_Mat, &rows_B, &cols_B);

	// Ensure all matrices are dense
	MatConvert(A, MATDENSE, MAT_INPLACE_MATRIX, &A);

	// Extract raw pointers to matrix data
	std::complex<double> *A_data, *B_data;
	MatDenseGetArray(A, &A_data);
	MatDenseGetArray(B_Mat, &B_data);

	std::vector<int> ipiv (rows_A);
	// int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, rows_A, cols_B, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) B_data, cols_B);
	int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, rows_A, cols_B, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) B_data, rows_B);

	MatDenseRestoreArray(A, &A_data);
	MatDenseRestoreArray(B_Mat, &B_data);

	MatToColVecs(B_Mat, X);
	MatDestroy(&B_Mat);

	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Function that computes specific elements of A^{-1} * B without actually constructing the matrix.*/
void MatMatSolve_IndexPairs(Mat &A, Mat &X, Mat &B, std::vector<std::set<int>> &idx, bool factored, double *dt)
{

	std::clock_t t1 = clock();

	Mat *A_factored = NULL;

	if (factored)
		A_factored = &A;
	else
		MatFactorize_SuperLU(A, *A_factored);

	int rows, cols;
	MatGetSize(A, &rows, &cols);

	// Allocate the output matrix
	std::vector<int> nnz (rows, 0);
	for (std::size_t ii = 0; ii < idx.size(); ii++)
		nnz[ii] = idx[ii].size();
	MatCreateAndAllocate_MATAIJ(X, rows, cols, 0, &nnz[0]);

	// Generate the index selector vector
	Vec ei;
	VecCreateAndAllocate(ei, cols, 0.0);

	Vec xi;
	VecCreateAndAllocate(xi, rows, 0.0);

	Mat BT;
	MatTranspose(B, MAT_INITIAL_MATRIX, &BT);

	// Iterate over target rows
	std::set<int>::iterator it;
	for (std::size_t ii = 0; ii < idx.size(); ii++)
	{
		
		// Set the current target row
		VecSet(ei, 0.0);
		VecSetValue(ei, ii, 1.0, INSERT_VALUES);
		VecAssemblyBegin(ei);
		VecAssemblyEnd(ei);
		
		// Compute the row solution vector
		MatSolveTranspose(*A_factored, ei, xi);

		std::complex<double> *xi_vec;
		VecGetArray(xi, &xi_vec);
		
		// Iterate over source rows
		for (it = idx[ii].begin(); it != idx[ii].end(); ++it)
		{
			int jj = *it;

			// Compute the column solution vector
			const std::complex<double> *vals;
			const int *cols;
			int ncols;
			MatGetRow(BT, jj, &ncols, &cols, &vals);

			std::complex<double> val = 0.0;
			for (int kk = 0; kk < ncols; kk++)
				val += vals[kk]*xi_vec[cols[kk]];

			MatRestoreRow(BT, jj, &ncols, &cols, &vals);
			MatSetValue(X, ii, jj, val, INSERT_VALUES);
		}

		VecRestoreArray(xi, &xi_vec);
		
	}
	
	MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(X, MAT_FINAL_ASSEMBLY);

	
	// Clean-up
	MatDestroy(&BT);
	VecDestroy(&ei);
	VecDestroy(&xi);

	
	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


// ==================================================================================
// MatSolve
// ==================================================================================

/*! \brief Function to solve Ax = b for x, after factorizing A using SuperLU. x should *not* have already been allocated.*/
void MatSolve_SuperLU(Mat &A, Vec &x, Vec &b, bool factored, bool allocate_x, double *dt)
{

	std::clock_t t1 = clock();
	
	if (allocate_x)
	{
		int rows, cols;
		MatGetSize(A, &rows, &cols);
		VecCreateAndAllocate(x, cols);
	}
	
	if (factored)
		MatSolve(A, b, x);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_SuperLU(A, A_factored);
		MatSolve(A_factored, b, x);
		MatDestroy(&A_factored);
	}

	
	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief Function to solve Ax = b for x, after factorizing A using SuperLU_dist. x should *not* have already been allocated.*/
void MatSolve_SuperLU_dist(Mat &A, Vec &x, Vec &b, bool factored, bool allocate_x, double *dt)
{

	std::clock_t t1 = clock();
	
	if (allocate_x)
	{
		int rows, cols;
		MatGetSize(A, &rows, &cols);
		VecCreateAndAllocate(x, cols);
	}
	
	if (factored)
		MatSolve(A, b, x);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_SuperLU_dist(A, A_factored);
		MatSolve(A_factored, b, x);
		MatDestroy(&A_factored);
	}


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief This function solves Ax = b for x, after factorizing A using Petsc's default routines. x should *not* have already been allocated.*/
void MatSolve_Petsc(Mat &A, Vec &x, Vec &b, bool factored, bool allocate_x, double *dt)
{

	std::clock_t t1 = clock();
	
	if (allocate_x)
	{
		int rows, cols;
		MatGetSize(A, &rows, &cols);
		VecCreateAndAllocate(x, cols);
	}

	if (factored)
		MatSolve(A, b, x);
	else
	{
		Mat A_factored = PETSC_NULL;
		MatFactorize_Petsc(A, A_factored);
		MatSolve(A_factored, b, x);
		MatDestroy(&A_factored);
	}


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


/*! \brief This function solves Ax = b for X, using Lapacke's gesv directly. It is assumed that the input matrix is dense, and stored in row-major layout. If the input matrix is not dense, it will be converted to a dense matrix anyway. Note: the input matrix is overwritten by the L and U factors.*/
void MatSolve_Lapacke(Mat &A, Vec &x, Vec &b, double *dt)
{

	std::clock_t t1 = clock();
	
	int rows_A, cols_A;
	MatGetSize(A, &rows_A, &cols_A);

	// VecDuplicate(b, &x);
	VecCopy(b, x);

	// Ensure all matrices are dense
	MatConvert(A, MATDENSE, MAT_INPLACE_MATRIX, &A);

	// Extract raw pointers to matrix data
	std::complex<double> *A_data, *b_data;
	MatDenseGetArray(A, &A_data);
	VecGetArray(x, &b_data);

	std::vector<int> ipiv (rows_A);
	// int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, rows_A, 1, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) b_data, 1);
	int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, rows_A, 1, (lapack_complex_double *) A_data, rows_A, &ipiv[0], (lapack_complex_double *) b_data, cols_A);

	MatDenseRestoreArray(A, &A_data);
	VecRestoreArray(x, &b_data);


	std::clock_t t2 = clock();
	if (dt != NULL)
		*dt = (t2 - t1)/(double)CLOCKS_PER_SEC;
	
	return;

}


// ==================================================================================
// KSP
// ==================================================================================

/*! \brief Function to store the estimated relative true residual norm at each iteration of an iterative solve. \todo This throws a bad_alloc error right now - fix it.*/
PetscErrorCode KSPStoreResidual(KSP ksp, int n, PetscReal rnorm, void *mctx)
{
	
	// The following snippet was adapted from https://www.mcs.anl.gov/petsc/petsc-current/src/ksp/ksp/interface/iterativ.c.html#KSPMonitorDefault

	// PetscReal truenorm, bnorm;
	// KSP_data *setup = static_cast<KSP_data *> (mctx);

	// KSPBuildResidual(ksp, PETSC_NULL, PETSC_NULL, &setup->resid);
	// VecNorm(resid, NORM_2, &truenorm);
	// VecDestroy(&setup->resid);
	// KSPGetRhs(ksp, &setup->rhs);
	// VecNorm(setup->rhs, NORM_2, &bnorm);
	// VecDestroy(&setup->rhs);

	// double rrnorm = (double)(truenorm/bnorm);
	// setup->rrnorm.push_back(rrnorm);
	
	return 0;
	
}


/*! \brief Function to initialize and set up the Petsc KSP solver. Separate initialization is useful for repeated solves.*/
void KSPInitialize_GMRES(KSP &ksp, Mat &A, KSP_data *setup)
{

	// FYI: Options database keys for monitoring the residual:
	// ksp_monitor_true_residual
	// ksp_monitor_singular_value
	// ksp_monitor_max
	// ksp_monitor_lg_true_residualnorm

	// KSPDestroy(&ksp);

	KSP_data default_setup;
	if (setup == NULL)
		setup = &default_setup;
	
	PC pc;
	MPI_Comm comm = MPI_COMM_SELF;

	// Set up the solver object
	KSPCreate(comm, &ksp);
	if (setup->M == NULL || std::string(setup->pc_type) == std::string(PCNONE))
		KSPSetOperators(ksp, A, A);
	else if (setup->M != NULL || std::string(setup->pc_type) != std::string(PCNONE))
		KSPSetOperators(ksp, A, *setup->M);
	// KSPSetReusePreconditioner(ksp, PETSC_TRUE);


	KSPSetType(ksp, KSPGMRES);

	
	// Apply solver settings
	if (setup->pc_side == PC_RIGHT)
		KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
	
	KSPGMRESSetRestart(ksp, setup->restart);
	KSPSetTolerances(ksp, setup->rtol, setup->atol, PETSC_DEFAULT, setup->max_it);

	// Apply preconditioner settings
	KSPSetPCSide(ksp, setup->pc_side);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, setup->pc_type);
	PCShellSetApply(pc, setup->pc_apply);
	PCShellSetContext(pc, setup->pc_ctx);
	
	// Residual monitoring
	if (setup->monitor_residual)
	{
		// Cancel any previous / default monitoring (besides run-time user settings)
		KSPMonitorCancel(ksp);

		PetscViewerAndFormat *vf;
		PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);

		// petsc 3.9.1
		KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorTrueResidualNorm, vf, NULL);
		
		// petsc 3.15
		// KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorTrueResidual, vf, NULL);
	}

	
	// Plot the residual norm, if desired
	if (setup->plot_residual_norm)
	{
		// petsc 3.9.1
		KSPMonitorLGTrueResidualNormCreate(comm, NULL, setup->plot_title.c_str(), 100, 100, 800, 600, &setup->lgctx);
		KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorLGTrueResidualNorm, setup->lgctx, NULL);
		
		// petsc 3.15
		// KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorTrueResidualDrawLG, NULL, NULL);
	}
		
	
	// Compute singular value estimates, if desired
	if (setup->compute_singular_values)
	{
		KSPSetComputeSingularValues(ksp, PETSC_TRUE);

		PetscViewerAndFormat *vf;
		PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);

		KSPMonitorSet(ksp, (PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*))KSPMonitorSingularValue, vf, NULL);
	}


	// Residual storage
	KSPSetResidualHistory(ksp, NULL, setup->max_it, PETSC_TRUE);

	// Allow for run-time user override of settings
	KSPSetFromOptions(ksp);

	// Set up the solver
	KSPSetUp(ksp);

	// Initial guess settings
	if (setup->use_fischer_guess)
		KSPSetUseFischerGuess(ksp, 1, setup->fischer_size);
	

	return;

}


/*! \brief This function iteratively solves Ax = b for x using PETSc's GMRES algorithm. The system matrix A may be a sparse or dense matrix, or a shell matrix for a matrix-free solve. In that case, the shell matrix must have already been fully defined through the PETSc functions usermult(), MatCreateShell() and MatShellSetOperation(). x should *not* have already been allocated. A custom initial guess vector x0 can optionally be provided through the setup data structure. Any preconditioning should be accounted for either by pre-multiplying the system matrix, or during the matrix-vector product, in case of a shell matrix.*/
void MatSolve_GMRES(Mat &A, Vec &x, Vec &b, KSP_data *setup, KSP *ksp_in)
{

	KSP_data default_setup;
	if (setup == NULL)
		setup = &default_setup;

	
	// Initialize the iterative solver
	KSP ksp = PETSC_NULL, *ksp_ptr = NULL;
	if (ksp_in == NULL)
	{
		KSPInitialize_GMRES(ksp, A, setup);
		ksp_ptr = &ksp;
	}
	else
		ksp_ptr = ksp_in;


	// Allocation of unknowns
	if (setup->allocate_x)
	{
		int rows;
		VecGetSize(b, &rows);
		VecCreateAndAllocate(x, rows);
	}

	
	// Initial guess settings
	if (setup->x0 == PETSC_NULL && !setup->use_fischer_guess)
		KSPSetInitialGuessNonzero(*ksp_ptr, PETSC_FALSE);
	else if (!setup->use_fischer_guess)
	{
		KSPSetInitialGuessNonzero(*ksp_ptr, PETSC_TRUE);
		VecCopy(setup->x0, x);
	}
	
	// KSPSetInitialGuessKnoll(*ksp_ptr, PETSC_TRUE);
	// KSPSetUseFischerGuess(ksp, 1, 300);


	// Solve the system
	std::clock_t t1_KSPsolve = clock();

	KSPSolve(*ksp_ptr, b, x);

	std::clock_t t2_KSPsolve = clock();
	setup->solve_time = (t2_KSPsolve - t1_KSPsolve)/(double)CLOCKS_PER_SEC;
	
	// Number of iterations that were run
	KSPGetIterationNumber(*ksp_ptr, &setup->N_iters);

	// Get convergence reason
	KSPGetConvergedReason(*ksp_ptr, &setup->reason);

	switch (setup->reason)
	{
	case 1: setup->print_reason = "CONVERGED: rtol_normal."; break;
	case 9: setup->print_reason = "CONVERGED: atol_normal."; break;
	case 2: setup->print_reason = "CONVERGED: rtol."; break;
	case 3: setup->print_reason = "CONVERGED: atol."; break;
	case -3: setup->print_reason = "DID NOT CONVERGE: Maximum iterations reached."; break;
	case -4: setup->print_reason = "DID NOT CONVERGE: Residual norm diverged."; break;
	case -9: setup->print_reason = "DID NOT CONVERGE: NaN or inf encountered."; break;
	case -11: setup->print_reason = "DID NOT CONVERGE: Preconditioner setup failed."; break;
	default: setup->print_reason = "DID NOT CONVERGE: Code " + std::to_string(setup->reason); break;
	}

	// Report iterative solve stats, if required
	if (setup->report_stats)
	{

		// Obtain the residual norm history

		// petsc 3.9.1
		double *rrnorm;

		// petsc 3.15
		// const double *rrnorm;

		KSPGetResidualHistory(*ksp_ptr, &rrnorm, &setup->na);

		// Normalize the residual norms
		double bnorm;
		VecNorm(b, NORM_2, &bnorm);
		setup->rrnorm.resize(setup->na);
		for (int ii = 0; ii < setup->na; ii++)
			setup->rrnorm[ii] = rrnorm[ii]/bnorm;

		std::cout << "\n================================================" << std::endl;
		std::cout << "Iterative solve completed in " << setup->solve_time << " s." << std::endl;
		std::cout << setup->print_reason << std::endl;
		std::cout << "Iterations: " << setup->N_iters << std::endl;
		std::cout << "Final relative residual norm: " << setup->rrnorm.back() << std::endl;
		std::cout << "================================================\n" << std::endl;
		
	}
	else if (setup->report_general)
	{
		std::cout << "\n================================================" << std::endl;
		std::cout << "Iterative solve completed in " << setup->solve_time << " s." << std::endl;
		std::cout << setup->print_reason << std::endl;
		std::cout << "Iterations: " << setup->N_iters << std::endl;
		std::cout << "================================================\n" << std::endl;
	}
	
	// Get max and min singular values, if required
	if (setup->compute_singular_values)
		KSPComputeExtremeSingularValues(*ksp_ptr, &setup->s_min, &setup->s_min);

	// Plot the residual norm, if desired
	if (setup->plot_residual_norm)
	{
		PetscDraw draw;
		PetscDrawLGGetDraw(setup->lgctx, &draw);
		PetscDrawSetPause(draw, -1);

		PetscDrawPause(draw);
		PetscDrawLGDestroy(&setup->lgctx);
	}

	if (ksp_in == NULL)
		KSPDestroy(ksp_ptr);

	
	return;

}


/*! \brief Given a shell matrix with a matrix-vector product defined via MATOP_MULT, this function explicitly builds the entire matrix using the mat-vec applied to a pair of index selection vectors. The returned matrix is of course dense.*/
void MatCreateFromShell(Mat &in, Mat &out)
{

	int N_rows, N_cols;
	MatGetSize(in, &N_rows, &N_cols);

	MatCreateAndAllocate_MATDENSE(out, N_rows, N_cols);

	// Generate the index selector vectors
	Vec ei = PETSC_NULL, ej = PETSC_NULL;
	VecCreateAndAllocate(ei, N_rows, 0.0);
	VecCreateAndAllocate(ej, N_cols, 0.0);

	Vec xj = PETSC_NULL;
	VecCreateAndAllocate(xj, N_rows, 0.0);

	// Traverse mat-vecs with the right-side selector vector
	for (int jj = 0; jj < N_cols; jj++)
	{

		// Set the current target column
		VecSet(ej, 0.0);
		VecSetValue(ej, jj, 1.0, INSERT_VALUES);
		VecAssemblyBegin(ej);
		VecAssemblyEnd(ej);
		
		// Compute the right-side mat-vec
		MatMult(in, ej, xj);

		// Traverse dots products with the left-side selector vector
		for (int ii = 0; ii < N_rows; ii++)
		{
				
			// The left-side dot product just requires selecting the correct element
			std::complex<double> val;
			VecGetValues(xj, 1, &ii, &val);
			MatSetValue(out, ii, jj, val, INSERT_VALUES);	

		}

	}

	MatAssemblyBegin(out, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(out, MAT_FINAL_ASSEMBLY);

	return;

}


/*! \brief Given a shell preconditioner object with its "apply" defined via PCShellSetApply and its constituent operators defined via PCSetOperators, this function explicitly builds the inverse preconditioner matrix by applying the preconditioner operator to a pair of index selection vectors. The returned matrix is of course dense.*/
void PCCreateFromShell(PC &pc, Mat &mat)
{

	Mat A = PETSC_NULL, P = PETSC_NULL;
	PCGetOperators(pc, &A, &P);

	int N_rows, N_cols;
	MatGetSize(A, &N_rows, &N_cols);

	MatCreateAndAllocate_MATDENSE(mat, N_rows, N_cols);

	// Generate the index selector vectors
	Vec ei = PETSC_NULL, ej = PETSC_NULL;
	VecCreateAndAllocate(ei, N_rows, 0.0);
	VecCreateAndAllocate(ej, N_cols, 0.0);

	Vec xj = PETSC_NULL;
	VecCreateAndAllocate(xj, N_rows, 0.0);

	// Traverse mat-vecs with the right-side selector vector
	for (int jj = 0; jj < N_cols; jj++)
	{

		// Set the current target column
		VecSet(ej, 0.0);
		VecSetValue(ej, jj, 1.0, INSERT_VALUES);
		VecAssemblyBegin(ej);
		VecAssemblyEnd(ej);
		
		// Compute the right-side mat-vec
		PCApply(pc, ej, xj);

		// Traverse dots products with the left-side selector vector
		for (int ii = 0; ii < N_rows; ii++)
		{
				
			// The left-side dot product just requires selecting the correct element
			std::complex<double> val;
			VecGetValues(xj, 1, &ii, &val);
			MatSetValue(mat, ii, jj, val, INSERT_VALUES);	

		}

	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	return;

}	


/*! \brief Function to compute the singular value decomposition of a Petsc-based matrix. Note that this will temporarily generate a dense copy of the input matrix.*/
void MatSVD(Mat &mat, std::vector<double> &s)
{

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	// Make a copy of the matrix so the original stays intact
	Mat A = PETSC_NULL;
	MatDuplicate(mat, MAT_COPY_VALUES, &A);

	// Ensure the matrix is dense
	MatConvert(A, MATDENSE, MAT_INPLACE_MATRIX, &A);

	// Extract raw pointers to matrix data
	std::complex<double> *A_data = NULL;
	MatDenseGetArray(A, &A_data);

	// Compute SVD
	std::complex<double> *u = NULL, *vt = NULL;
	int ldu = std::min(rows, cols), ldvt = cols;
	std::vector<double> superb (rows);
	s.resize(std::min(rows, cols));
	
	int info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR, 'N', 'N', rows, cols, (lapack_complex_double *) A_data, cols, &s[0], (lapack_complex_double *) u, ldu, (lapack_complex_double *) vt, ldvt, &superb[0]);

	if (info < 0)
		std::cout << "[WARNING] MatSVD(): LAPACKE_zgesvd returned " << info << "; parameter " << -info << " had an illegal value." << std::endl;
	else if (info > 0)
		std::cout << "[WARNING] MatSVD(): LAPACKE_zgesvd returned " << info << "; did not converge." << std::endl;
	
	MatDenseRestoreArray(A, &A_data);
	MatDestroy(&A);
	
	return;

}


/*! \brief Function to compute the rank and condition number (based on SVD) of a Petsc-based matrix. Note that this will temporarily generate a dense copy of the input matrix for SVD computation.*/
void MatCondRank(Mat &mat, double *cond, int *rank, double *rank_rel, bool verbose)
{

	if (verbose)
		std::cout << "============ MatCondRank() called ============" << std::endl;
	
	double t_cond, t_rank_rel;
	int t_rank;

	if (cond == NULL)
		cond = &t_cond;

	if (rank == NULL)
		rank = &t_rank;

	if (rank_rel == NULL)
		rank_rel = &t_rank_rel;

	
	double tol = 1.0e-16;
	std::vector<double> s;
	MatSVD(mat, s);

	const auto [min, max] = std::minmax_element(std::begin(s), std::end(s));
	
	*cond = (*max)/(*min);

	*rank = 0;
	for (std::size_t ii = 0; ii < s.size(); ii++)
		if (std::abs(s[ii])/std::abs(s[0]) > tol)
			*rank = *rank + 1;

	*rank_rel = ((double)(*rank))/((double)s.size());

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	if (verbose)
	{
		std::cout << "Size: " << rows << " x " << cols << std::endl;
		std::cout << "Condition number: " << *cond << std::endl;
		std::cout << "Rank: " << *rank << std::endl;
		std::cout << "Relative rank: " << *rank_rel << std::endl;
		std::cout << "==============================================" << std::endl;
	}
	
	return;

}


/*! \brief Function to compute the eigenvalue decomposition of a Petsc-based square matrix. Note that this will temporarily generate a dense copy of the input matrix.*/
void MatEig(Mat &mat, std::vector<std::complex<double>> &w)
{

	int rows, cols;
	MatGetSize(mat, &rows, &cols);

	if (rows != cols)
		throw std::logic_error("[ERROR] MatEig(): Matrix must be square.");

	int n = rows;

	// Make a copy of the matrix so the original stays intact
	Mat A = PETSC_NULL;
	MatDuplicate(mat, MAT_COPY_VALUES, &A);

	// Ensure the matrix is dense
	MatConvert(A, MATDENSE, MAT_INPLACE_MATRIX, &A);

	// Extract raw pointers to matrix data
	std::complex<double> *A_data = NULL;
	MatDenseGetArray(A, &A_data);

	// Compute eigenvalues
	int matrix_layout = LAPACK_ROW_MAJOR;	
	int lda = n, ldvl = n, ldvr = n;
	std::complex<double> *vl = NULL, *vr = NULL;	
	w.resize(n);

	int info = LAPACKE_zgeev(matrix_layout, 'N', 'N', n, (lapack_complex_double *) A_data, (lapack_int) lda, (lapack_complex_double *) &w[0], (lapack_complex_double *) vl, ldvl, (lapack_complex_double *) vr, ldvr);

	if (info < 0)
		std::cout << "[WARNING] MatEig(): LAPACKE_zgeev returned " << info << "; parameter " << -info << " had an illegal value." << std::endl;
	else if (info > 0)
		std::cout << "[WARNING] MatEig(): LAPACKE_zgeev returned " << info << "; did not converge." << std::endl;
	
	MatDenseRestoreArray(A, &A_data);
	MatDestroy(&A);
	
	return;

}


/*! \brief Function to mask matrix A based on nonzeros in matrix B. Could be more efficient. \todo*/
void MatMatMask(Mat &A, Mat &B)
{

	double tol = 1.0e-15;
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	for (int ii = 0; ii < rows; ii++)
	{
		for (int jj = 0; jj < cols; jj++)
		{
			std::complex<double> val_A = 0.0, val_B;
			MatGetValues(B, 1, &ii, 1, &jj, &val_B);

			if (std::abs(val_B) < tol)
				MatSetValues(A, 1, &ii, 1, &jj, &val_A, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	return;

}


/*! \brief Function to pointwise-multiply matrix A with matrix B. Could be more efficient. \todo*/
void MatMatPointwiseMult(Mat &A, Mat &B, Mat &C)
{

	double tol = 1.0e-15;
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
	MatSetOption(C, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	
	for (int ii = 0; ii < rows; ii++)
	{
		for (int jj = 0; jj < cols; jj++)
		{
			std::complex<double> val_A, val_B, val_C;
			MatGetValues(A, 1, &ii, 1, &jj, &val_A);
			MatGetValues(B, 1, &ii, 1, &jj, &val_B);
			val_C = val_A*val_B;

			// if (std::abs(val_C) > tol)
				MatSetValues(C, 1, &ii, 1, &jj, &val_C, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);

	return;

}


/*! \brief Function to set values of the input matrix to their element-wise absolute values. \todo Could be more efficient.*/
void MatAbs(Mat &A)
{

	double tol = 1.0e-15;
	
	int rows, cols;
	MatGetSize(A, &rows, &cols);

	Mat B = PETSC_NULL;
	MatDuplicate(A, MAT_COPY_VALUES, &B);
	MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

	for (int ii = 0; ii < rows; ii++)
	{
		for (int jj = 0; jj < cols; jj++)
		{
			std::complex<double> val;
			MatGetValues(B, 1, &ii, 1, &jj, &val);
			val = std::abs(val);

			// if (std::abs(val) > tol)
				MatSetValues(A, 1, &ii, 1, &jj, &val, INSERT_VALUES);
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	MatDestroy(&B);

	return;

}


/*! \brief Function to eliminate a given set of rows of a matrix, thus changing its dimensions.*/
void MatEliminateRows(Mat &mat, std::vector<int> &rows)
{

	if (rows.size() < 1)
		return;

	Mat orig;
	MatDuplicate(mat, MAT_COPY_VALUES, &orig);

	int rows_orig, cols_orig;
	MatGetSize(orig, &rows_orig, &cols_orig);

	int nnz_per_row [rows_orig];
	for (int ii = 0; ii < rows_orig; ii++)
		MatGetRow(orig, ii, &(nnz_per_row[ii]), NULL, NULL);

	int rows_new = rows_orig - rows.size();
	int nnz_per_row_new [rows_new];

	int jj = 0;
	for (int ii = 0; ii < rows_orig; ii++)
	{
		if (std::find(rows.begin(), rows.end(), ii) == rows.end())
		{
			nnz_per_row_new[jj] = nnz_per_row[ii];
			jj++;
		}
	}
	
	MatDestroy(&mat);
	MatCreateAndAllocate_MATAIJ(mat, rows_new, cols_orig, 0, nnz_per_row_new);

	for (int ii = 0; ii < rows_orig; ii++)
		MatRestoreRow(orig, ii, &(nnz_per_row[ii]), NULL, NULL);

	jj = 0;
	for (int ii = 0; ii < rows_orig; ii++)
	{
		if (std::find(rows.begin(), rows.end(), ii) == rows.end())
		{
			int nnz;
			const int *cols;
			const std::complex<double> *vals;
			MatGetRow(orig, ii, &nnz, &cols, &vals);
			MatSetValues(mat, 1, &jj, nnz, cols, vals, INSERT_VALUES);
			MatRestoreRow(orig, ii, &nnz, &cols, &vals);
			jj++;
		}
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	MatDestroy(&orig);

	return;

}


/*! \brief Function to eliminate a given set of columns of a matrix, thus changing its dimensions. Note that this function actually eliminates rows of the transpose of the input matrix.*/
void MatEliminateCols(Mat &mat, std::vector<int> &cols)
{

	if (cols.size() < 1)
		return;

	Mat mat_new;

	MatTranspose(mat, MAT_INITIAL_MATRIX, &mat_new);
	MatEliminateRows(mat_new, cols);

	MatDestroy(&mat);
	MatTranspose(mat_new, MAT_INITIAL_MATRIX, &mat);

	MatDestroy(&mat_new);
	
	return;

}


/*! \brief Function to insert rows of zeros into a matrix, thus changing its dimensions. The input vector should contain the indices that the new rows would have, after their insertion.*/
void MatInsertRows(Mat &mat, std::vector<int> &rows)
{

	if (rows.size() < 1)
		return;

	Mat orig;
	MatDuplicate(mat, MAT_COPY_VALUES, &orig);

	int rows_orig, cols_orig;
	MatGetSize(orig, &rows_orig, &cols_orig);

	int rows_new = rows_orig + rows.size();

	int nnz_per_row [rows_orig];
	for (int ii = 0; ii < rows_orig; ii++)
		MatGetRow(orig, ii, &(nnz_per_row[ii]), NULL, NULL);

	int nnz_per_row_new [rows_new];
	std::vector<int> idx_rows_orig (rows_orig);

	int jj = 0;
	for (int ii = 0; ii < rows_new; ii++)
	{
		if (std::find(rows.begin(), rows.end(), ii) == rows.end())
		{
			nnz_per_row_new[ii] = nnz_per_row[jj];
			idx_rows_orig[jj] = ii;
			jj++;
		}
		else
			nnz_per_row_new[ii] = 1;
	}
	
	MatDestroy(&mat);
	MatCreateAndAllocate_MATAIJ(mat, rows_new, cols_orig, 0, nnz_per_row_new);

	for (int ii = 0; ii < rows_orig; ii++)
		MatRestoreRow(orig, ii, &(nnz_per_row[ii]), NULL, NULL);

	std::complex<double> zero = 0.0;
	for (std::size_t ii = 0; ii < rows.size(); ii++)
		MatSetValues(mat, 1, &(rows[ii]), 0, &(rows[ii]), &zero, INSERT_VALUES);
	
	for (int ii = 0; ii < rows_orig; ii++)
	{
		int nnz;
		const int *cols;
		const std::complex<double> *vals;
		MatGetRow(orig, ii, &nnz, &cols, &vals);
		MatSetValues(mat, 1, &(idx_rows_orig[ii]), nnz, cols, vals, INSERT_VALUES);
		MatRestoreRow(orig, ii, &nnz, &cols, &vals);
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

	MatDestroy(&orig);

	return;

}


/*! \brief Function to insert columns of zeros into a matrix, thus changing its dimensions. The input vector should contain the indices that the new columns would have, after their insertion. Note that this function actually inserts rows into the transpose of the input matrix.*/
void MatInsertCols(Mat &mat, std::vector<int> &cols)
{

	if (cols.size() < 1)
		return;

	Mat mat_new;

	MatTranspose(mat, MAT_INITIAL_MATRIX, &mat_new);
	MatInsertRows(mat_new, cols);

	MatDestroy(&mat);
	MatTranspose(mat_new, MAT_INITIAL_MATRIX, &mat);

	MatDestroy(&mat_new);
	
	return;

}


/*! \brief Function to create a matrix which, when left-multiplied into an arbitrary matrix of compatible size, copies source rows of that matrix on to target rows. When its transpose is right-multiplied, it copies columns instead.*/
void MatCreateRowCopyMatrix(Mat &mat, int size, std::vector<int> &source_rows, std::vector<int> &target_rows)
{

	if (source_rows.size() != target_rows.size())
	{
		std::cout << "[Error] MatCreateRowCopyMatrix(): Number of source rows must equal number of target rows. Returning an identity matrix instead." << std::endl;
		MatIdentity_MATAIJ(mat, size);
		return;
	}

	MatCreateAndAllocate_MATAIJ(mat, size, size, 1 + source_rows.size(), NULL);

	std::complex<double> one = 1.0, zero = 0.0;
	
	for (int ii = 0; ii < size; ii++)
		MatSetValues(mat, 1, &ii, 1, &ii, &one, INSERT_VALUES);

	for (std::size_t ii = 0; ii < source_rows.size(); ii++)
	{
		MatSetValues(mat, 1, &target_rows[ii], 1, &target_rows[ii], &zero, INSERT_VALUES);
		MatSetValues(mat, 1, &target_rows[ii], 1, &source_rows[ii], &one, INSERT_VALUES);
	}

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);


	return;

}


/*! \brief Function to create a matrix which, when left-multiplied into an arbitrary matrix of compatible size, adds a scaled version of source rows of that matrix to target rows. When its transpose is right-multiplied, it adds columns instead.*/
void MatCreateRowAdditionMatrix(Mat &mat, int size, std::vector<int> &source_rows, std::vector<int> &target_rows, std::complex<double> val)
{

	if (source_rows.size() != target_rows.size())
	{
		std::cout << "[Error] MatCreateRowAdditionMatrix(): Number of source rows must equal number of target rows. Returning an identity matrix instead." << std::endl;
		MatIdentity_MATAIJ(mat, size);
		return;
	}

	MatCreateAndAllocate_MATAIJ(mat, size, size, 1 + source_rows.size(), NULL);

	std::complex<double> one = 1.0;
	
	for (int ii = 0; ii < size; ii++)
		MatSetValues(mat, 1, &ii, 1, &ii, &one, INSERT_VALUES);

	for (std::size_t ii = 0; ii < source_rows.size(); ii++)
		MatSetValues(mat, 1, &target_rows[ii], 1, &source_rows[ii], &val, INSERT_VALUES);

	MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);


	return;

}


/*! \brief Function to compute A = A*B using MatMatMult.*/
void MatMatMult_OverwriteA(Mat &A, Mat &B)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatMult_OverwriteA(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatMult_OverwriteA(): Matrix B is null; no action taken." << std::endl;
		return;
	}

	Mat C = PETSC_NULL;

	MatMatMult(A, B, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&A);
	MatDuplicate(C, MAT_COPY_VALUES, &A);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute A = A*B using MatMatMult.*/
void MatMatTransposeMult_OverwriteA(Mat &A, Mat &B)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatTransposeMult_OverwriteA(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatTransposeMult_OverwriteA(): Matrix B is null; no action taken." << std::endl;
		return;
	}
	
	Mat C = PETSC_NULL;

	MatMatTransposeMult(A, B, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&A);
	MatDuplicate(C, MAT_COPY_VALUES, &A);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute B = A*B using MatMatMult.*/
void MatMatMult_OverwriteB(Mat &A, Mat &B)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatMult_OverwriteB(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatMult_OverwriteB(): Matrix B is null; no action taken." << std::endl;
		return;
	}

	Mat C = PETSC_NULL;

	MatMatMult(A, B, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&B);
	MatDuplicate(C, MAT_COPY_VALUES, &B);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute B = A^T*B using MatTransposeMatMult.*/
void MatTransposeMatMult_OverwriteB(Mat &A, Mat &B)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatTransposeMatMult_OverwriteB(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatTransposeMatMult_OverwriteB(): Matrix B is null; no action taken." << std::endl;
		return;
	}
	
	Mat C = PETSC_NULL;

	MatTransposeMatMult(A, B, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&B);
	MatDuplicate(C, MAT_COPY_VALUES, &B);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute B = A*B^T using MatMatTransposeMult.*/
void MatMatTransposeMult_OverwriteB(Mat &A, Mat &B)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatTransposeMult_OverwriteB(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatMatTransposeMult_OverwriteB(): Matrix B is null; no action taken." << std::endl;
		return;
	}
	
	Mat C = PETSC_NULL;

	MatMatTransposeMult(A, B, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&B);
	MatDuplicate(C, MAT_COPY_VALUES, &B);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute A = P^T A P using MatMatMult (because Petsc's MatPtAP doesn't allow for mixed dense and sparse matrices yet).*/
void MatPtAP_OverwriteA(Mat &A, Mat &P)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatPtAP_OverwriteA(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (P == PETSC_NULL)
	{
		std::cout << "[WARNING] MatPtAP_OverwriteA(): Matrix P is null; no action taken." << std::endl;
		return;
	}
	
	Mat C = PETSC_NULL;
	MatMatMult(A, P, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&A);
	MatTransposeMatMult(P, C, MAT_INITIAL_MATRIX, 2.0, &A);

	MatDestroy(&C);

	return;

}


/*! \brief Function to compute A = B^T A C using MatMatMult.*/
void MatBtAC_OverwriteA(Mat &A, Mat &B, Mat &C)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatBtAC_OverwriteA(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (B == PETSC_NULL)
	{
		std::cout << "[WARNING] MatBtAC_OverwriteA(): Matrix B is null; no action taken." << std::endl;
		return;
	}
	if (C == PETSC_NULL)
	{
		std::cout << "[WARNING] MatBtAC_OverwriteA(): Matrix C is null; no action taken." << std::endl;
		return;
	}
	
	Mat D = PETSC_NULL;
	MatMatMult(A, C, MAT_INITIAL_MATRIX, 2.0, &D);
	MatDestroy(&A);
	MatTransposeMatMult(B, D, MAT_INITIAL_MATRIX, 2.0, &A);

	MatDestroy(&D);

	return;

}


/*! \brief Function to compute C = R^T A P using MatMatMult.*/
void MatRtAP(Mat &A, Mat &R, Mat &P, Mat &C)
{

	if (A == PETSC_NULL)
	{
		std::cout << "[WARNING] MatRtAP(): Matrix A is null; no action taken." << std::endl;
		return;
	}
	if (R == PETSC_NULL)
	{
		std::cout << "[WARNING] MatRtAP(): Matrix R is null; no action taken." << std::endl;
		return;
	}
	if (P == PETSC_NULL)
	{
		std::cout << "[WARNING] MatRtAP(): Matrix P is null; no action taken." << std::endl;
		return;
	}
	
	Mat temp = PETSC_NULL;
	MatMatMult(A, P, MAT_INITIAL_MATRIX, 2.0, &temp);
	MatTransposeMatMult(R, temp, MAT_INITIAL_MATRIX, 2.0, &C);
	MatDestroy(&temp);
	
	return;
	
}


std::string MatInfoToString(MatInfo &info)
{
    std::stringstream str_str;
    str_str << std::setw(10);
    str_str << std::left << "Assemblies:\t" << std::right <<  info.assemblies << std::endl;
    str_str << std::left <<  "Block Size:\t" << std::right <<  info.block_size << std::endl;
    str_str << std::left <<  "Factor Mallocs:\t" << std::right <<  info.factor_mallocs<< std::endl;
    str_str << std::left <<  "Fill Ratio Given:\t" << std::right <<  info.fill_ratio_given << "\t\t\t";
    str_str << std::left <<  "Fill Ratio Needed:\t" << std::right <<  info.fill_ratio_needed << std::endl;
    str_str << std::left <<  "Mallocs:\t" << std::right <<  info.mallocs << "\t\t\t\t\t";
    str_str << std::left <<  "Memory:\t" << std::right <<  info.memory<< std::endl;
    str_str << std::left <<  "NZ Allocated:\t" << std::right <<  info.nz_allocated << "\t\t";
    str_str << std::left <<  "NZ Unneeded:\t" << std::right <<  info.nz_unneeded << "\t\t";
    str_str << std::left <<  "NZ Used:\t" << std::right <<  info.nz_used << std::endl;
    return str_str.str();
}


// =================================================================================
// Distributed matrices
// =================================================================================

void SaveDistMatToASCII(Mat Matrix, const std::string &file_name)
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(MPI_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(Matrix, viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
}


void SaveDistMatToBinary(Mat Matrix, const std::string &file_name)
{
    PetscViewer viewer;
    PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &viewer);
    MatView(Matrix, viewer);
    PetscViewerDestroy(&viewer);
}


void LoadDistMatFromBinary(Mat Matrix, const std::string &file_name)
{
    PetscViewer viewer;
    PetscViewerBinaryOpen(MPI_COMM_WORLD, file_name.c_str(), FILE_MODE_READ, &viewer);
    MatLoad(Matrix, viewer);
    PetscViewerDestroy(&viewer);
}


void SaveDistVecToBinary(Vec vector, const std::string &file_name)
{
	PetscViewer viewer;
	MPI_Comm comm;
	PetscObjectGetComm((PetscObject)vector, &comm);
	PetscViewerBinaryOpen(comm, file_name.c_str(), FILE_MODE_WRITE, &viewer);
	VecView(vector, viewer);
	PetscViewerDestroy(&viewer);
}

void LoadDistVecFromBinary(Vec vector, const std::string &file_name)
{
	PetscViewer viewer;
	MPI_Comm comm;
	PetscObjectGetComm((PetscObject)vector, &comm);
	PetscViewerBinaryOpen(comm, file_name.c_str(), FILE_MODE_READ, &viewer);
	VecLoad(vector, viewer);
	PetscViewerDestroy(&viewer);
}


/*! \brief Overloaded Petsc-enhancing wrapper function to allocate memory for a distributed Petsc
 * vector. Vector entries are initialized to zeros.*/
void Vec_Init_Dist(Vec &vec, int N)
{
	VecCreate(MPI_COMM_WORLD, &vec);
	VecSetType(vec, VECMPI);
	VecSetSizes(vec, PETSC_DECIDE, N);
	VecSet(vec, 0.0);
	VecAssemblyBegin(vec);
	VecAssemblyEnd(vec);
}


Mat DuplicateTypeAndSize(Mat &mat)
{

    Mat duplicant;

    PetscInt g_rows, g_cols, l_rows, l_cols;
    MatGetSize(mat, &g_rows, &g_cols);
    MatGetLocalSize(mat, &l_rows, &l_cols);

    MatCreate(MPI_COMM_WORLD, &duplicant);

    MatSetSizes(duplicant, l_rows, l_cols, g_rows, g_cols);
    return duplicant;
}


// =================================================================================
// Unit tests
// =================================================================================

/*! \brief Function to create a test matrix A for various Petsc unit tests.*/
void Test_CreateA(Mat &A)
{

	// Test matrix A:
	// 
	// A = 1 0 4 4 7
	//     0 0 3 2 0
	//     6 7 0 2 0
	//     6 5 5 0 1
	//     8 2 1 0 0
	// size(A) = 5 rows x 5 cols

	int A_rows = 5, A_cols = 5;
	int row, col;
	std::complex<double> val;
	MatCreateAndAllocate_MATAIJ(A, A_rows, A_cols, 5, NULL);

	// Row 0
	
	row = 0; col = 0; val = 1.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 0; col = 2; val = 4.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 0; col = 3; val = 4.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 0; col = 4; val = 7.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	// Row 1
	
	row = 1; col = 2; val = 3.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 1; col = 3; val = 2.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	// Row 2
	
	row = 2; col = 0; val = 6.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 2; col = 1; val = 7.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 2; col = 3; val = 2.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	// Row 3
	
	row = 3; col = 0; val = 6.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);   

	row = 3; col = 1; val = 5.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);   

	row = 3; col = 2; val = 5.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);   

	row = 3; col = 4; val = 1.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	// Row 4
	
	row = 4; col = 0; val = 8.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 4; col = 1; val = 2.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	row = 4; col = 2; val = 1.0;
	MatSetValues(A, 1, &row, 1, &col, &val, INSERT_VALUES);

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	return;
	
}


void Test_MatInsertRows(std::string debug_folder)
{

	Mat A;
	Test_CreateA(A);

	std::vector<int> idx = {0, 3, 4, 6, 9};

	// Expected output:
	// 
	// A = 0 0 0 0 0 <
	//     1 0 4 4 7
	//     0 0 3 2 0
	//     0 0 0 0 0 <
	//     0 0 0 0 0 <
	//     6 7 0 2 0
	//     0 0 0 0 0 <
	//     6 5 5 0 1
	//     8 2 1 0 0
	//     0 0 0 0 0 <
	// size(A) = 10 rows x 5 cols

	PrintMatToFile(A, 1.0, debug_folder + "Test_MatInsertRows_A_before.debug");
	
	MatInsertRows(A, idx);

	PrintMatToFile(A, 1.0, debug_folder + "Test_MatInsertRows_A_after.debug");

	return;

}


void Test_MatInsertCols(std::string debug_folder)
{

	Mat A;
	Test_CreateA(A);

	std::vector<int> idx = {0, 3, 4, 6, 9};

	// Expected output:
	// 
	// A = 0 1 0 0 0 4 0 4 7 0
	//     0 0 0 0 0 3 0 2 0 0
	//     0 6 7 0 0 0 0 2 0 0
	//     0 6 5 0 0 5 0 0 1 0
	//     0 8 2 0 0 1 0 0 0 0
	//     ^     ^ ^   ^     ^
	// size(A) = 5 rows x 10 cols

	PrintMatToFile(A, 1.0, debug_folder + "Test_MatInsertCols_A_before.debug");
	
	MatInsertCols(A, idx);

	PrintMatToFile(A, 1.0, debug_folder + "Test_MatInsertCols_A_after.debug");

	return;

}



