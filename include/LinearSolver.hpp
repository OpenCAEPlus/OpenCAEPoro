/*! \file    LinearSolver.hpp
 *  \brief   Lienar solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoro header files
#include "LinearSystem.hpp"
#include "OCPConst.hpp"

// Fasp header files
#include "fasp.h"
extern "C" {
#include "fasp_functs.h"
}

using namespace std;

/// A template class designed to stores and solves linear system derived from discrete
/// method. the maxtrix is stored in the form of row-segmented CSRX internaly, whose
/// sparsity pattern is almost the same as neighbor in Connection_BB. if martix is block
/// matrix, then the T will be a vector.
template <typename T> class Solver
{
    friend class OpenCAEPoro;
    friend class Connection_BB;
    friend class Well;

public:
    /// Allocate memory for linear system where possible maximum numbers of row are
    /// used.
    void AllocateMem(const OCP_USI& dimMax);
    /// Allocate memory for each row of matrix where the most terrible condition was
    /// considered.
    void AllocateColValMem();
    /// read solver param from input file which is supplied by user.
    void SetupParam(const string& dir, const string& file);
    /// output the solution to fileX.
    void PrintfSolution(const string& fileX) const;

    // FASP
    /// initialize the sover param for FASP.
    void InitParam_Fasp();
    /// read the sover param for FASP.
    void ReadParam_Fasp();
    /// convert the internal matrix structure into the the format required by FASP.
    void AssembleMat_Fasp();
    /// solve the linear system by FASP and return the status.
    int FaspSolve();
    /// free the matrix used for FASP.
    void Free_Fasp() { fasp_dcsr_free(&A_Fasp); };
    /// output the mat and rhs to fileA and fileb.
    void PrintfMatCSR(const string& fileA, const string& fileb) const;
    /// clear the internal matrix.
    void ClearData();
    /// return the solution.
    vector<T>& GetSolution() { return u; }

    /// check if NAN or INF occurs, used in debug mode.
    void CheckVal() const;

private:
    // internal mat structure.
    /// the maximum dimens matrix could have, it's fixed, always memory saving.
    /// it's used to allocate memory for the mat at the beginning of simulation.
    OCP_USI maxDim;
    OCP_USI dim; ///< the actual dimens of mat, it may changed all the time but always a
                 ///< little less than maxDim.
    /// the maximum possible capacity of each row of the mat.
    /// it's just a little greater than actual sizes, so it's very memory saving.
    /// it's used to allocate memory for the mat at the beginning of simulation.
    vector<USI> rowCapacity;
    vector<vector<OCP_USI>>
        colId; ///< column-index of nonzero entry in the format of row-segmented.
    vector<vector<T>> val; ///< values of nonzero entry in the format of row-segmented.
    vector<USI> diagPtr;   ///< the ith entry indicates the location of diagal entry of
                           ///< the ith row in val.
    /// an auxiliary variable used to help Setup entry in diagnal line.
    /// it will only be used when matrix is being assembling.
    vector<T> diagVal;
    vector<T> b; ///< right hands of linear system.
    vector<T> u; ///< solutiom of linear system.

    // for FASP solver
    string  solveDir;
    string  solveFile;
    dCSRmat A_Fasp;
    dvector b_Fasp;
    dvector x_Fasp;

    input_param inParam; // parameters from input files
    ITS_param   itParam; // parameters for itsolver

    AMG_param amgParam; // parameters for AMG
    ILU_param iluParam; // parameters for ILU
    SWZ_param swzParam; // parameters for Schwarz method

    // for other solver
};

template <typename T> void Solver<T>::AllocateMem(const OCP_USI& dimMax)
{
    maxDim = dimMax;
    rowCapacity.resize(maxDim, 0);
    colId.resize(maxDim);
    val.resize(maxDim);
    diagPtr.resize(maxDim, 0);
    diagVal.resize(maxDim, 0);
    b.resize(maxDim, 0);
    u.resize(maxDim, 0);
}

template <typename T> void Solver<T>::AllocateColValMem()
{
    for (OCP_USI n = 0; n < maxDim; n++) {
        colId[n].reserve(rowCapacity[n]);
        val[n].reserve(rowCapacity[n]);
    }
}

template <typename T> void Solver<T>::ClearData()
{
    for (OCP_USI i = 0; i < maxDim; i++) {
        colId[i].clear();
        val[i].clear();
    }
    diagPtr.assign(maxDim, 0);
    diagVal.assign(maxDim, 0);
    b.assign(maxDim, 0);
    u.assign(maxDim, 0);
}

template <typename T> void Solver<T>::AssembleMat_Fasp()
{
    // b & x
    b_Fasp.row = dim;
    b_Fasp.val = b.data();
    x_Fasp.row = dim;
    x_Fasp.val = u.data();
    // A
    OCP_USI nnz = 0;
    for (OCP_USI i = 0; i < dim; i++) {
        nnz += colId[i].size();
    }
    A_Fasp = fasp_dcsr_create(dim, dim, nnz);
    // IA

    OCP_USI count = 0;
    for (OCP_USI i = 1; i < dim + 1; i++) {
        USI nnz_Row  = colId[i - 1].size();
        A_Fasp.IA[i] = A_Fasp.IA[i - 1] + nnz_Row;

        for (USI j = 0; j < nnz_Row; j++) {
            A_Fasp.JA[count]  = colId[i - 1][j];
            A_Fasp.val[count] = val[i - 1][j];
            count++;
        }
    }
}

template <typename T> int Solver<T>::FaspSolve()
{
    int status = FASP_SUCCESS;

    const int print_level  = inParam.print_level;
    const int solver_type  = inParam.solver_type;
    const int precond_type = inParam.precond_type;
    const int output_type  = inParam.output_type;

    if (output_type) {
        const char* outputfile = "out/Solver.out";
        printf("Redirecting outputs to file: %s ...\n", outputfile);
        freopen(outputfile, "w", stdout); // open a file for stdout
    }

    // Preconditioned Krylov methods
    if (solver_type >= 1 && solver_type <= 20) {

        // Using no preconditioner for Krylov iterative methods
        if (precond_type == PREC_NULL) {
            status = fasp_solver_dcsr_krylov(&A_Fasp, &b_Fasp, &x_Fasp, &itParam);
        }

        // Using diag(A) as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dcsr_krylov_diag(&A_Fasp, &b_Fasp, &x_Fasp, &itParam);
        }

        // Using AMG as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
            status = fasp_solver_dcsr_krylov_amg(&A_Fasp, &b_Fasp, &x_Fasp, &itParam,
                                                 &amgParam);
        }

        // Using ILU as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_ILU) {
            if (print_level > PRINT_NONE) fasp_param_ilu_print(&iluParam);
            status = fasp_solver_dcsr_krylov_ilu(&A_Fasp, &b_Fasp, &x_Fasp, &itParam,
                                                 &iluParam);
        }

        // Undefined iterative methods
        else {
            printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);
            status = ERROR_SOLVER_PRECTYPE;
        }

    }

    // AMG as the iterative solver
    else if (solver_type == SOLVER_AMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
        fasp_solver_amg(&A_Fasp, &b_Fasp, &x_Fasp, &amgParam);
    }

    // Full AMG as the iterative solver
    else if (solver_type == SOLVER_FMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgParam);
        fasp_solver_famg(&A_Fasp, &b_Fasp, &x_Fasp, &amgParam);
    }

#if WITH_MUMPS // use MUMPS directly
    else if (solver_type == SOLVER_MUMPS) {
        status = fasp_solver_mumps(&A_Fasp, &b_Fasp, &x_Fasp, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_SuperLU // use SuperLU directly
    else if (solver_type == SOLVER_SUPERLU) {
        status = fasp_solver_superlu(&A_Fasp, &b_Fasp, &x_Fasp, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#if WITH_UMFPACK // use UMFPACK directly
    else if (solver_type == SOLVER_UMFPACK) {
        // Need to sort the matrix A for UMFPACK to work
        dCSRmat A_trans = fasp_dcsr_create(A_Fasp.row, A_Fasp.col, A_Fasp.nnz);
        fasp_dcsr_transz(&A_Fasp, NULL, &A_trans);
        fasp_dcsr_sort(&A_trans);
        status = fasp_solver_umfpack(&A_trans, &b_Fasp, &x_Fasp, print_level);
        fasp_dcsr_free(&A_trans);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

#ifdef WITH_PARDISO // use PARDISO directly
    else if (solver_type == SOLVER_PARDISO) {
        fasp_dcsr_sort(&A_Fasp);
        status = fasp_solver_pardiso(&A_Fasp, &b_Fasp, &x_Fasp, print_level);
        if (status >= 0) status = 1; // Direct solver returns 1
    }
#endif

    else {
        printf("### ERROR: Wrong solver type %d!!!\n", solver_type);
        status = ERROR_SOLVER_TYPE;
    }

    if (status < 0) {
        printf("### ERROR: Solver failed! Exit status = %d.\n\n", status);
    }

    if (output_type) fclose(stdout);
    return status;
}

template <typename T>
void Solver<T>::PrintfMatCSR(const string& fileA, const string& fileb) const
{
    string FileA = solveDir + fileA;
    string Fileb = solveDir + fileb;
    // csr format
    ofstream outA(FileA);
    if (!outA.is_open()) cout << "Can not open " << FileA << endl;
    outA << A_Fasp.row << endl;
    OCP_USI nnz0 = A_Fasp.nnz;
    for (OCP_USI i = 0; i < A_Fasp.row + 1; i++) outA << A_Fasp.IA[i] + 1 << endl;
    for (OCP_USI i = 0; i < nnz0; i++) outA << A_Fasp.JA[i] + 1 << endl;
    for (OCP_USI i = 0; i < nnz0; i++) outA << A_Fasp.val[i] << endl;
    outA.close();
    // out rhs
    ofstream outb(Fileb);
    if (!outb.is_open()) cout << "Can not open " << Fileb << endl;
    outb << b_Fasp.row << endl;
    for (OCP_USI i = 0; i < b_Fasp.row; i++) outb << b_Fasp.val[i] << endl;
    outb.close();
}

template <typename T> void Solver<T>::PrintfSolution(const string& fileU) const
{
    string   FileU = solveDir + fileU;
    ofstream outu(FileU);
    if (!outu.is_open()) cout << "Can not open " << FileU << endl;
    outu << dim << endl;
    for (OCP_USI i = 0; i < dim; i++) outu << u[i] << endl;
    outu.close();
}

template <typename T> void Solver<T>::SetupParam(const string& dir, const string& file)
{
    solveDir  = dir;
    solveFile = file;
#ifdef __SOLVER_FASP__
    ReadParam_Fasp();
#endif //  __SOLVER_FASP__
}

template <typename T> void Solver<T>::CheckVal() const
{
    for (OCP_USI n = 0; n < dim; n++) {
        for (auto v : val[n]) {
            if (!isfinite(v)) {
                cout << "###ERROR   ";
                ERRORcheck("NAN or INF in MAT");
                exit(0);
            }
        }
    }
}

template <typename T> void Solver<T>::ReadParam_Fasp()
{
    string file = solveDir + solveFile;
    InitParam_Fasp(); // Set default solver parameters
    std::ifstream ifs(file);
    if (!ifs.is_open()) {
        std::cout << "The input file " << file << " is missing!" << std::endl;
        file = solveDir + "../conf/csr.dat";
        ifs.open(file);
        if (!ifs.is_open()) {
            std::cout << "The input file " << file << " is missing!" << std::endl;
            std::cout << "Using the default parameters of FASP" << std::endl;
        } else {
            ifs.close();
            std::cout << "Using the input file " << file << std::endl;
            fasp_param_input(file.data(), &inParam);
        }
    } else {
        ifs.close(); // if file has been opened, close it first
        fasp_param_input(file.data(), &inParam);
    }
    fasp_param_init(&inParam, &itParam, &amgParam, &iluParam, &swzParam);
}

template <typename T> void Solver<T>::InitParam_Fasp()
{
    // Input/output
    inParam.print_level = PRINT_MIN;
    inParam.output_type = 0;

    // Problem information
    inParam.solver_type  = SOLVER_VFGMRES;
    inParam.precond_type = PREC_AMG;
    inParam.stop_type    = STOP_REL_RES;

    // Solver parameters
    inParam.itsolver_tol   = 1e-4;
    inParam.itsolver_maxit = 100;
    inParam.restart        = 30;

    // ILU method parameters
    inParam.ILU_type    = ILUk;
    inParam.ILU_lfil    = 0;
    inParam.ILU_droptol = 0.001;
    inParam.ILU_relax   = 0;
    inParam.ILU_permtol = 0.0;

    // Schwarz method parameters
    inParam.SWZ_mmsize    = 200;
    inParam.SWZ_maxlvl    = 2;
    inParam.SWZ_type      = 1;
    inParam.SWZ_blksolver = SOLVER_DEFAULT;

    // AMG method parameters
    inParam.AMG_type                = CLASSIC_AMG;
    inParam.AMG_levels              = 20;
    inParam.AMG_cycle_type          = V_CYCLE;
    inParam.AMG_smoother            = SMOOTHER_GS;
    inParam.AMG_smooth_order        = CF_ORDER;
    inParam.AMG_presmooth_iter      = 1;
    inParam.AMG_postsmooth_iter     = 1;
    inParam.AMG_relaxation          = 1.0;
    inParam.AMG_coarse_dof          = 500;
    inParam.AMG_coarse_solver       = 0;
    inParam.AMG_tol                 = 1e-6;
    inParam.AMG_maxit               = 1;
    inParam.AMG_ILU_levels          = 0;
    inParam.AMG_SWZ_levels          = 0;
    inParam.AMG_coarse_scaling      = OFF;
    inParam.AMG_amli_degree         = 1;
    inParam.AMG_nl_amli_krylov_type = 2;

    // Classical AMG specific
    inParam.AMG_coarsening_type      = 1;
    inParam.AMG_interpolation_type   = 1;
    inParam.AMG_max_row_sum          = 0.9;
    inParam.AMG_strong_threshold     = 0.3;
    inParam.AMG_truncation_threshold = 0.2;
    inParam.AMG_aggressive_level     = 0;
    inParam.AMG_aggressive_path      = 1;

    // Aggregation AMG specific
    inParam.AMG_aggregation_type   = PAIRWISE;
    inParam.AMG_quality_bound      = 8.0;
    inParam.AMG_pair_number        = 2;
    inParam.AMG_strong_coupled     = 0.25;
    inParam.AMG_max_aggregation    = 9;
    inParam.AMG_tentative_smooth   = 0.67;
    inParam.AMG_smooth_filter      = ON;
    inParam.AMG_smooth_restriction = ON;
}

#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/