/***************************VECTORMATH_H******************************
 * VectorMath.h contains function declaration for elementary math
 * operations used in this program.
 *
 *
 * Author: Utkarsh R. Patel
 ***************************************************************/



#ifndef VECTORMATH_H
#define VECTORMATH_H

#include <complex>


		/// Calculates L2-norm of a vector. 
		template <typename T>
		T VectorNorm(T * v, int N = 3){
			T result = 0.0;
			for(int ii = 0; ii < N; ii++)
				result += v[ii]*v[ii];

			return sqrt(result);
			}  		

		
		/// Calculates sum of two vectors.
		template<typename T>
		T * VectorAdd(T *v1, T *v2, T *result, int N = 3) 
		{
			for(int ii = 0; ii < N; ii++)
			{
				result[ii] = v1[ii] + v2[ii];
			}
			return result;
		}
		

		/// Calculates subtraction of two vectors.
		template<typename T>
		T * VectorMinus(T *v1, T *v2, T *result, int N = 3) 
		{
			for(int ii = 0; ii < N; ii++)
			{
				result[ii] = v1[ii] - v2[ii];
			}
			return result;
		}
		
		/// Calculates result of a vector multiplied by a scalar (for double arguments)
		template<typename T>
		T * VectorMultiplyScalar(T v1, T *v2, T * result, int N = 3)
		{
			for(int ii = 0; ii < N; ii++)
			{
				result[ii] = v1*v2[ii];
			}
			return result;
		}

		/// Calculates dot product of two vectors.
		template<typename T>
		T VectorDot (T * v1, T *v2, int N = 3 )
		{
			T result = 0.0;
			for(int ii = 0; ii < N; ii++)
			{
				result += v1[ii]*v2[ii];
			}
			return result;
		}

		/// Calculates cross product of two vectors
		template<typename T>
		T * VectorCross(T * v1, T *v2, T * result)
		{
			result[0] = v1[1]*v2[2] - v1[2]*v2[1];
			result[1] = v1[2]*v2[0] - v1[0]*v2[2];
			result[2] = v1[0]*v2[1] - v1[1]*v2[0];

			return result;
		}



		/// Prints vectors onto terminal
		template<typename T>
		void VectorPrint(T * v1)
		{
			std::cout << "(" << v1[0] << "," << v1[1] << "," << v1[2] << ")" << std::endl;
			return;
		}


#endif
