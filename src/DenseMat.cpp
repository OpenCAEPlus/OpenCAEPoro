/*! \file    DenseMat.cpp
 *  \brief   Dense matrix-vector operations
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "DenseMat.hpp"


void Dcopy(const int& N, double* dst, const double* src)
{
    const int incx = 1, incy = 1;
    dcopy_(&N, src, &incx, dst, &incy);
}

// WARNING: absolute sum!
double Dnorm1(const int& N, double* x)
{
    const int incx = 1;
    return dasum_(&N, x, &incx);
}

double Dnorm2(const int& N, double* x)
{
    const int incx = 1;
    return dnrm2_(&N, x, &incx);
}

void Dscalar(const int& n, const double& alpha, double* x)
{
    // x = a x
    const int incx = 1;
    dscal_(&n, &alpha, x, &incx);
}

void Daxpy(const int& n, const double& alpha, const double* x, double* y)
{
    // y= ax +y
    const int incx = 1, incy = 1;
    daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void DaABpbC(const int& m, const int& n, const int& k, const double& alpha,
             const double* A, const double* B, const double& beta, double* C)
{
    /*  C' = alpha B'A' + beta C'
     *  A: m x k
     *  B: k x n
     *  C: m x n
     *  all column majored matrices, no tranpose
     */

    const char transa = 'N', transb = 'N';
    dgemm_(&transa, &transb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n);
}

void DaAxpby(const int& m, const int& n, const double& a, const double* A,
             const double* x, const double& b, double* y)
{
    /*  y= aAx+by
     */
    for (int i = 0; i < m; i++) {
        y[i] = b * y[i];
        for (int j = 0; j < n; j++) {
            y[i] += a * A[i * n + j] * x[j];
        }
    }
}

void LUSolve(const int& nrhs, const int& N, double* A, double* b, int* pivot)
{
    int info;

    dgesv_(&N, &nrhs, A, &N, pivot, b, &N, &info);

    if (info < 0) {
        cout << "Wrong Input !" << endl;
    } else if (info > 0) {
        cout << "Singular Matrix !" << endl;
    }
}

void SYSSolve(const int& nrhs, const char* uplo, const int& N, double* A, double* b, int* pivot)
{
    int    info;
    double work[1] = {0};
    int    lwork   = 1;

    dsysv_(uplo, &N, &nrhs, A, &N, pivot, b, &N, work, &lwork, &info);
    if (info < 0) {
        cout << "Wrong Input !" << endl;
    } else if (info > 0) {
        cout << "Singular Matrix !" << endl;
    }
}

void PrintDX(const int& N, const double* x)
{
    for (int i = 0; i < N; i++) {
        cout << setprecision(9) << x[i] << "    ";
    }
    cout << endl;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/
