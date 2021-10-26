#include "DenseMat.hpp"

void Dscalar(int n, double alpha, double* x) {
    //x = a x
    const int incx = 1;
    dscal_(&n, &alpha, x, &incx);
}


void Daxpy(int n, double alpha, const double* x, double* y) {
    //y= ax +y
    const int incx = 1, incy = 1;
    daxpy_(&n, &alpha, x, &incx, y, &incy);
}


void DaABpbC(int m, int n, int k, double alpha, const double* A, const double* B, double beta, double* C) {
    /*  C' = alpha B'A' + beta C'
     *  A: m x k
     *  B: k x n
     *  C: m x n
     *  all column majored matrices, no tranpose
     */

    const char transa = 'N', transb = 'N';
    dgemm_(&transa, &transb, &n, &m, &k, &alpha, B, &n, A, &k, &beta, C, &n);
}