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

void Dscalar(int n, double alpha, double* x)
{
    // x = a x
    const int incx = 1;
    dscal_(&n, &alpha, x, &incx);
}

void Daxpy(int n, double alpha, const double* x, double* y)
{
    // y= ax +y
    const int incx = 1, incy = 1;
    daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void DaABpbC(int m, int n, int k, double alpha, const double* A, const double* B,
             double beta, double* C)
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

void DaAxpby(int m, int n, double a, const double* A, const double* x, double b,
             double* y)
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

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/
