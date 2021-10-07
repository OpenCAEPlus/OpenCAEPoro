/*! \file    Solver.hpp
 *  \brief   Solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/07/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

#include "MAT.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "fasp.h"

extern "C" {
#include "fasp_functs.h"
}

/// A template class designed to stores and solves linear system derived from discrete method.
/// the maxtrix is stored in the form of row-segmented CSRX internaly, whose sparsity pattern is almost the same as Neighbor in Connection_BB.
/// if martix is block matrix, then the T will be a vector.
template <typename T> class Solver
{
    friend class OpenCAEPoro;
    friend class Connection_BB;
    friend class Well;

public:
    /// allocate memory for linear system where possible maximum numbers of row are used.
    void allocate(const OCP_USI& dimMax);
    /// allocate memory for each row of matrix where the most terrible condition was considered.
    void allocateColVal();
    /// read solver param from input file which is supplied by user.
    void setupParam(const string& dir, const string& file);
    /// output the solution to fileX.
    void showSolution(const string& fileX) const;

    // FASP
    /// initialize the sover param for FASP.
    void init_param_Fasp();
    /// read the sover param for FASP.
    void read_param_Fasp();
    /// convert the internal matrix structure into the the format required by FASP.
    void assemble_Fasp();
    /// solve the linear system by FASP and return the status.
    int  faspsolve();
    /// free the matrix used for FASP.
    void free_Fasp() { fasp_dcsr_free(&A_Fasp); };
    /// output the mat and rhs to fileA and fileb.
    void showMat_CSR(const string& fileA, const string& fileb) const;
    /// clear the internal matrix.
    void clearData();
    /// return the solution.
    vector<T>& getSol() { return u; }

    /// check if NAN or INF occurs, used in debug mode.
    void checkVal() const;

private:
    // internal mat structure.
    /// the maximum dimens matrix could have, it's fixed, always memory saving.
    /// it's used to allocate memory for the mat at the beginning of simulation. 
    OCP_USI     MaxDim;     
    OCP_USI     Dim; ///< the actual dimens of mat, it may changed all the time but always a little less than MaxDim.
    /// the maximum possible capacity of each row of the mat.
    /// it's just a little greater than actual sizes, so it's very memory saving.
    /// it's used to allocate memory for the mat at the beginning of simulation. 
    vector<USI>            RowCapacity; 
    vector<vector<OCP_USI>> ColId; ///< column-index of nonzero entry in the format of row-segmented.
    vector<vector<T>>   Val;    ///< values of nonzero entry in the format of row-segmented.
    vector<USI>              DiagPtr;   ///< the ith entry indicates the location of diagal entry of the ith row in Val.
    /// an auxiliary variable used to help setup entry in diagnal line.
    /// it will only be used when matrix is being assembling.
    vector<T>                DiagVal; 
    vector<T> b;        ///< right hands of linear system.
    vector<T> u;        ///< solutiom of linear system.


    // for FASP solver
    string  SolveDir;
    string  SolveFile;
    dCSRmat A_Fasp;
    dvector b_Fasp;
    dvector x_Fasp;

    input_param inparam; // parameters from input files
    ITS_param   itparam; // parameters for itsolver

    AMG_param amgparam; // parameters for AMG
    ILU_param iluparam; // parameters for ILU
    SWZ_param swzparam; // parameters for Schwarz method

    // for other solver
};

template <typename T> void Solver<T>::allocate(const OCP_USI& dimMax)
{
    MaxDim = dimMax;
    RowCapacity.resize(MaxDim, 0);
    ColId.resize(MaxDim);
    Val.resize(MaxDim);
    DiagPtr.resize(MaxDim, 0);
    DiagVal.resize(MaxDim, 0);
    b.resize(MaxDim, 0);
    u.resize(MaxDim, 0);
}

template <typename T> void Solver<T>::allocateColVal()
{
    for (OCP_USI n = 0; n < MaxDim; n++) {
        ColId[n].reserve(RowCapacity[n]);
        Val[n].reserve(RowCapacity[n]);
    }
}

template <typename T> void Solver<T>::clearData()
{
    for (OCP_USI i = 0; i < MaxDim; i++) {
        ColId[i].clear();
        Val[i].clear();
    }
    DiagPtr.assign(MaxDim, 0);
    DiagVal.assign(MaxDim, 0);
    b.assign(MaxDim, 0);
    u.assign(MaxDim, 0);
}

template <typename T> void Solver<T>::assemble_Fasp()
{
    // b & x
    b_Fasp.row = Dim;
    b_Fasp.val = b.data();
    x_Fasp.row = Dim;
    x_Fasp.val = u.data();
    // A
    OCP_USI nnz = 0;
    for (OCP_USI i = 0; i < Dim; i++) {
        nnz += ColId[i].size();
    }
    A_Fasp = fasp_dcsr_create(Dim, Dim, nnz);
    // IA

    OCP_USI count = 0;
    for (OCP_USI i = 1; i < Dim + 1; i++) {
        USI nnz_Row  = ColId[i - 1].size();
        A_Fasp.IA[i] = A_Fasp.IA[i - 1] + nnz_Row;

        for (USI j = 0; j < nnz_Row; j++) {
            A_Fasp.JA[count]  = ColId[i - 1][j];
            A_Fasp.val[count] = Val[i - 1][j];
            count++;
        }
    }
}


template <typename T> int Solver<T>::faspsolve()
{
    int status = FASP_SUCCESS;

    const int print_level  = inparam.print_level;
    const int solver_type  = inparam.solver_type;
    const int precond_type = inparam.precond_type;
    const int output_type  = inparam.output_type;

    if (output_type) {
        const char* outputfile = "out/Solver.out";
        printf("Redirecting outputs to file: %s ...\n", outputfile);
        freopen(outputfile, "w", stdout); // open a file for stdout
    }

    // Preconditioned Krylov methods
    if (solver_type >= 1 && solver_type <= 20) {

        // Using no preconditioner for Krylov iterative methods
        if (precond_type == PREC_NULL) {
            status = fasp_solver_dcsr_krylov(&A_Fasp, &b_Fasp, &x_Fasp, &itparam);
        }

        // Using diag(A) as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dcsr_krylov_diag(&A_Fasp, &b_Fasp, &x_Fasp, &itparam);
        }

        // Using AMG as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
            status = fasp_solver_dcsr_krylov_amg(&A_Fasp, &b_Fasp, &x_Fasp, &itparam,
                                                 &amgparam);
        }

        // Using ILU as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_ILU) {
            if (print_level > PRINT_NONE) fasp_param_ilu_print(&iluparam);
            status = fasp_solver_dcsr_krylov_ilu(&A_Fasp, &b_Fasp, &x_Fasp, &itparam,
                                                 &iluparam);
        }

        // Undefined iterative methods
        else {
            printf("### ERROR: Wrong preconditioner type %d!!!\n", precond_type);
            status = ERROR_SOLVER_PRECTYPE;
        }

    }

    // AMG as the iterative solver
    else if (solver_type == SOLVER_AMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
        fasp_solver_amg(&A_Fasp, &b_Fasp, &x_Fasp, &amgparam);
    }

    // Full AMG as the iterative solver
    else if (solver_type == SOLVER_FMG) {
        if (print_level > PRINT_NONE) fasp_param_amg_print(&amgparam);
        fasp_solver_famg(&A_Fasp, &b_Fasp, &x_Fasp, &amgparam);
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

template <typename T> void Solver<T>::showMat_CSR(const string& fileA, const string& fileb) const
{
    string FileA = SolveDir + fileA;
    string Fileb = SolveDir + fileb;
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

template <typename T> void Solver<T>::showSolution(const string& fileU) const
{
    string   FileU = SolveDir + fileU;
    ofstream outu(FileU);
    if (!outu.is_open()) cout << "Can not open " << FileU << endl;
    outu << Dim << endl;
    for (OCP_USI i = 0; i < Dim; i++) outu << u[i] << endl;
    outu.close();
}

template <typename T> void Solver<T>::setupParam(const string& dir, const string& file)
{
    SolveDir  = dir;
    SolveFile = file;
#ifdef __SOLVER_FASP__
    read_param_Fasp();
#endif //  __SOLVER_FASP__
}

template <typename T> void Solver<T>::checkVal() const
{
    for (OCP_USI n = 0; n < Dim; n++) {
        for (auto v : Val[n]) {
            if (!isfinite(v)) {
                cout << "###ERROR   ";
                ERRORcheck("NAN or INF in MAT");
                exit(0);
            }
        }
    }
}

template <typename T> void Solver<T>::read_param_Fasp()
{
    string file = SolveDir + SolveFile;
    init_param_Fasp(); // Set default solver parameters
    std::ifstream ifs(file);
    if (!ifs.is_open()) {
        std::cout << "The input file " << file << " is missing!" << std::endl;
        file = SolveDir + "../conf/csr.dat";
        ifs.open(file);
        if (!ifs.is_open()) {
            std::cout << "The input file " << file << " is missing!" << std::endl;
            std::cout << "Using the default parameters of FASP" << std::endl;
        } else {
            ifs.close();
            std::cout << "Using the input file " << file << std::endl;
            fasp_param_input(file.data(), &inparam);
        }
    } else {
        ifs.close(); // if file has been opened, close it first
        fasp_param_input(file.data(), &inparam);
    }
    fasp_param_init(&inparam, &itparam, &amgparam, &iluparam, &swzparam);
}

template <typename T> void Solver<T>::init_param_Fasp()
{
    // Input/output
    inparam.print_level = PRINT_MIN;
    inparam.output_type = 0;

    // Problem information
    inparam.solver_type  = SOLVER_VFGMRES;
    inparam.precond_type = PREC_AMG;
    inparam.stop_type    = STOP_REL_RES;

    // Solver parameters
    inparam.itsolver_tol   = 1e-4;
    inparam.itsolver_maxit = 100;
    inparam.restart        = 30;

    // ILU method parameters
    inparam.ILU_type    = ILUk;
    inparam.ILU_lfil    = 0;
    inparam.ILU_droptol = 0.001;
    inparam.ILU_relax   = 0;
    inparam.ILU_permtol = 0.0;

    // Schwarz method parameters
    inparam.SWZ_mmsize    = 200;
    inparam.SWZ_maxlvl    = 2;
    inparam.SWZ_type      = 1;
    inparam.SWZ_blksolver = SOLVER_DEFAULT;

    // AMG method parameters
    inparam.AMG_type                = CLASSIC_AMG;
    inparam.AMG_levels              = 20;
    inparam.AMG_cycle_type          = V_CYCLE;
    inparam.AMG_smoother            = SMOOTHER_GS;
    inparam.AMG_smooth_order        = CF_ORDER;
    inparam.AMG_presmooth_iter      = 1;
    inparam.AMG_postsmooth_iter     = 1;
    inparam.AMG_relaxation          = 1.0;
    inparam.AMG_coarse_dof          = 500;
    inparam.AMG_coarse_solver       = 0;
    inparam.AMG_tol                 = 1e-6;
    inparam.AMG_maxit               = 1;
    inparam.AMG_ILU_levels          = 0;
    inparam.AMG_SWZ_levels          = 0;
    inparam.AMG_coarse_scaling      = OFF;
    inparam.AMG_amli_degree         = 1;
    inparam.AMG_nl_amli_krylov_type = 2;

    // Classical AMG specific
    inparam.AMG_coarsening_type      = 1;
    inparam.AMG_interpolation_type   = 1;
    inparam.AMG_max_row_sum          = 0.9;
    inparam.AMG_strong_threshold     = 0.3;
    inparam.AMG_truncation_threshold = 0.2;
    inparam.AMG_aggressive_level     = 0;
    inparam.AMG_aggressive_path      = 1;

    // Aggregation AMG specific
    inparam.AMG_aggregation_type   = PAIRWISE;
    inparam.AMG_quality_bound      = 8.0;
    inparam.AMG_pair_number        = 2;
    inparam.AMG_strong_coupled     = 0.25;
    inparam.AMG_max_aggregation    = 9;
    inparam.AMG_tentative_smooth   = 0.67;
    inparam.AMG_smooth_filter      = ON;
    inparam.AMG_smooth_restriction = ON;
}

#endif