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

void CalEigenSY(const int& N, float* A, float* w, float* work, const int& lwork)
{
    int  info;
    int  iwork[1] = {0};
    int  liwork   = 1;
    char uplo{'U'};
    char Nonly{'N'};

    ssyevd_(&Nonly, &uplo, &N, A, &N, w, work, &lwork, iwork, &liwork, &info);
    if (info > 0) {
        cout << "failed to compute eigenvalues!" << endl;
    }
}

// void MinEigenS(const int& N, float* a, float* w)
//{
//     MKL_INT info = LAPACKE_ssyevd(LAPACK_COL_MAJOR, 'N', 'U', N, a, N, w);
//     if (info > 0) {
//         cout << "failed to compute eigenvalues!" << endl;
//     }
// }

void Dcopy(const int& N, double* dst, const double* src)
{
    const int incx = 1, incy = 1;
    dcopy_(&N, src, &incx, dst, &incy);
}

double Ddot(int n, double* a, double* b)
{
    const int inca = 1, incb = 1;
    return ddot_(&n, a, &inca, b, &incb);
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

void DaABpbC(const int&    m,
             const int&    n,
             const int&    k,
             const double& alpha,
             const double* A,
             const double* B,
             const double& beta,
             double*       C)
{
    /*  C' = alpha B'A' + beta C'
     *  A: m x k
     *  B: k x n
     *  C: m x n
     *  all column majored matrices, no tranpose
     *  A' in col-order in Fortran = A in row-order in C/Cpp
     */

    const char transa = 'N', transb = 'N';
    dgemm_(&transa, &transb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n);
}

void myDABpC(const int&    m,
             const int&    n,
             const int&    k,
             const double* A,
             const double* B,
             double*       C)
{
    // C = AB + C
    // A: m*n  B:n*k  C:m*k
    // all matrix are row majored matrices

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            for (int m = 0; m < n; m++) {
                C[i * k + j] += A[i * n + m] * B[m * k + j];
            }
        }
    }
}

void myDABpCp(const int&    m,
              const int&    n,
              const int&    k,
              const double* A,
              const double* B,
              double*       C,
              const int*    flag,
              const int     N)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) C[i * k + j] += A[i * n + p] * B[p * k + j];
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + 3 + p * (N + 1) + m] *
                                        B[(3 + p * (N + 1) + m) * k + j];
                    }
                }
            }
        }
    }
}

void myDABpCp1(const int&    m,
               const int&    n,
               const int&    k,
               const double* A,
               const double* B,
               double*       C,
               const int*    flag,
               const int     N)
{

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int s = 0;
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) {
                    C[i * k + j] += A[i * n + p] * B[s * k + j];
                    s++;
                }
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + 3 + p * (N + 1) + m] * B[s * k + j];
                        s++;
                    }
                }
            }
        }
    }
}

void myDABpCp2(const int&    m,
               const int&    n,
               const int&    k,
               const double* A,
               const double* B,
               double*       C,
               const int*    flag,
               const int     N)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            int s = 0;
            for (int p = 0; p < 3; p++) {
                if (flag[p] != 0) {
                    C[i * k + j] += A[i * n + s] * B[p * k + j];
                    s++;
                }
            }
            for (int p = 0; p < 2; p++) {
                if (flag[p] != 0) {
                    for (int m = 0; m < N; m++) {
                        C[i * k + j] += A[i * n + s] * B[(3 + p * (N + 1) + m) * k + j];
                        s++;
                    }
                }
            }
        }
    }
}

void DaAxpby(const int&    m,
             const int&    n,
             const double& a,
             const double* A,
             const double* x,
             const double& b,
             double*       y)
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

int SYSSolve(const int&  nrhs,
             const char* uplo,
             const int&  N,
             double*     A,
             double*     b,
             int*        pivot,
             double*     work,
             const int&  lwork)
{
    int info;

    dsysv_(uplo, &N, &nrhs, A, &N, pivot, b, &N, work, &lwork, &info);
    if (info < 0) {
        cout << "Wrong Input !" << endl;
    } else if (info > 0) {
        cout << "Singular Matrix !" << endl;
    }

    return info;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/21/2021      Create file                          */
/*----------------------------------------------------------------------------*/
