/*! \file    DenseMat.hpp
 *  \brief   Operations about small dense mat
 *  \author  Shizhe Li
 *  \date    Oct/24/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DENSEMAT_HEADER__
#define __DENSEMAT_HEADER__

// Standard header files
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

extern "C" {

////// BLAS functions

/// Scales a vector by a constant.
void dscal_(const int* n, const double* alpha, double* x, const int* incx);

/// Forms the dot product of two vectors.
double ddot_(const int* n, double* a, const int* inca, double* b, const int* incb);

/// Copies a vector, src, to a vector, dst.
int dcopy_(const int* n, const double* src, const int* incx, double* dst,
           const int* incy);

/// Constant times a vector plus a vector.
int daxpy_(const int* n, const double* alpha, const double* x, const int* incx,
           double* y, const int* incy);

/// Computes the Euclidean norm of a vector.
double dnrm2_(const int* n, double* x, const int* incx);

/// Computes the sum of the absolute values of a vector.
double dasum_(const int* n, double* x, const int* incx);

/// Finds the index of element having max absolute value.
int idamax_(const int* n, double* x, const int* incx);

/// Performs matrix-matrix operations C : = alpha * op(A) * op(B) + beta * C.
int dgemm_(const char* transa, const char* transb, const int* m, const int* n,
           const int* k, const double* alpha, const double* A, const int* lda,
           const double* B, const int* ldb, const double* beta, double* C,
           const int* ldc);

////// LAPACK functions

/// Computes the solution to system of linear equations A * X = B for general matrices.
int dgesv_(const int* n, const int* nrhs, double* A, const int* lda, int* ipiv,
           double* b, const int* ldb, int* info);

/// Computes the solution to system of linear equations A * X = B for symm matrices.
int dsysv_(const char* uplo, const int* n, const int* nrhs, double* A, const int* lda,
           int* ipiv, double* b, const int* ldb, double* work, const int* lwork,
           int* info);
}

void Dcopy(const int& N, double* dst, const double* src);

/// Computes the L1-norm of a vector.
double Dnorm1(const int& N, double* x);

/// Computes the L2-norm of a vector.
double Dnorm2(const int& N, double* x);

/// Scales a vector by a constant.
void Dscalar(const int& n, const double& alpha, double* x);

/// Constant times a vector plus a vector.
void Daxpy(const int& n, const double& alpha, const double* x, double* y);

/// Computes C' = alpha B'A' + beta C', all matrices are column-major.
void DaABpbC(const int& m, const int& n, const int& k, const double& alpha,
             const double* A, const double* B, const double& beta, double* C);

/// Computes y = a A x + b y.
void DaAxpby(const int& m, const int& n, const double& a, const double* A,
             const double* x, const double& b, double* y);

/// Calls dgesv to solve the linear system for general matrices.
void LUSolve(const int& nrhs, const int& N, double* A, double* b, int* pivot);

/// Calls dsysy to solve the linear system for symm matrices.
void SYSSolve(const int& nrhs, const char* uplo, const int& N, double* A, double* b, int* pivot);


/// Prints a double vector.
void PrintDX(const int& N, const double* x);

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/24/2021      Create file                          */
/*  Chensong Zhang      Jan/16/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/