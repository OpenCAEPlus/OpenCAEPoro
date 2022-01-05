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

#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace std;

extern "C" {
	//BLAS functions
	void dscal_(const int* n, const double* alpha, double* x, const int* incx);
	double ddot_(const int* n, double* a, const int* inca, double* b, const int* incb);
	int dcopy_(const int* n, const double* src, const int* incx, double* dst, const int* incy);
	int daxpy_(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);
	double dnrm2_(const int* n, double* x, const int* incx);
	double dasum_(const int* n, double* x, const int* incx);
	int idamax_(const int* n, double* x, const int* incx);
	int dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
		const double* alpha, const double* A, const int* lda, const double* B, const int* ldb,
		const double* beta, double* C, const int* ldc);

	//LAPACK functions
	int dgesv_(const int* n, const int* nrhs, double* A, const int* lda, int* ipiv, double* b, const int* ldb, int* info);
	int dsysv_(const char* uplo, const int* n, const int* nrhs, double* A, const int* lda, int* ipiv,
				double* b, const int* ldb, double* work, const int* lwork, int* info);
}

double Dnorm1(const int& N, double* x);

double Dnorm2(const int& N, double* x);

void Dscalar(const int& n, const double& alpha, double* x);

void Daxpy(const int& n, const double& alpha, const double* x, double* y);

void DaABpbC(const int& m, const int& n, const int& k, const double& alpha, const double* A, const double* B, const double& beta, double* C);

void DaAxpby(const int& m, const int& n, const double& a, const double* A, const double* x, const double& b, double* y);

void LUSolve(const int& N, double* A, double* b, int* pivot);

void SYSSolve(const char* uplo, const int& N, double* A, double* b, int* pivot);

void PrintDX(const int& N, const double* x);

#endif


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/24/2021      Create file                          */
/*----------------------------------------------------------------------------*/