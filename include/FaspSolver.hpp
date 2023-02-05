/*! \file    FaspSolver.hpp
 *  \brief   Declaration of classes interfacing to the FASP solvers
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FASPSOLVER_HEADER__
#define __FASPSOLVER_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// faspsolver header files
extern "C" {
#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
}

// faspcpr header files
#if WITH_FASPCPR
extern "C" {
#include "faspcpr.h"
#include "faspcpr_functs.h"
}
#endif

// fasp4blkoil header files
#if WITH_FASP4BLKOIL
extern "C" {
#include "fasp4blkoil.h"
#include "fasp4blkoil_functs.h"
}
#endif

// fasp4cuda header files
// Note: It should not inside extern "C" {} !
#if WITH_FASP4CUDA
#include "fasp4cuda.h"
#include "fasp4cuda_functs.h"
#endif

// OpenCAEPoro header files
#include "LinearSolver.hpp"

using namespace std;

// Standard preconditioner types
#define PC_NULL  60 ///< None:  no preconditioner
#define PC_FASP1 61 ///< FASP1: MSP, default for FIM from 2020
#define PC_FASP2 62 ///< FASP2: MSP, experimental only
#define PC_FASP3 63 ///< FASP3: MSP, monolithic preconditioner
#define PC_FASP4 64 ///< FASP4: MSP, default for FIM from 2015
#define PC_FASP5 65 ///< FASP5: MSP, experimental only
#define PC_DIAG  68 ///< DIAG:  diagonal preconditioner
#define PC_BILU  69 ///< BILU:  block ILU preconditioner

// Sharing-setup preconditioner types
#define PC_FASP1_SHARE 71 ///< Sharing setup stage for PC_FASP1, use with caution
#define PC_FASP4_SHARE 74 ///< Sharing setup stage for PC_FASP4, use with caution
#define RESET_CONST    35 ///< Sharing threshold for PC_FASP1_SHARE, PC_FASP4_SHARE

/// Basic FASP solver class.
class FaspSolver : public LinearSolver
{
public:
    /// Set FASP parameters.
    void SetupParam(const string& dir, const string& file) override;

    /// Get number of iterations used by iterative solver.
    USI GetNumIters() const override { return itParam.maxit; }

public:
    string      solveDir;  ///< Current work dir
    string      solveFile; ///< Relative path of fasp file
    input_param inParam;   ///< Parameters from input files
    ITS_param   itParam;   ///< Parameters for iterative method
    AMG_param   amgParam;  ///< Parameters for AMG method
    ILU_param   iluParam;  ///< Parameters for ILU method
    SWZ_param   swzParam;  ///< Parameters for Schwarz method
};

/// Scalar solvers in CSR format from FASP.
class ScalarFaspSolver : public FaspSolver
{
    friend class LinearSystem;

private:
    /// Allocate memory for the linear system.
    void Allocate(const vector<USI>& rowCapacity,
                  const OCP_USI&     maxDim,
                  const USI&         blockDim) override;

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>&     colId,
                     const vector<vector<OCP_DBL>>& val,
                     const OCP_USI&                 dim,
                     const USI&                     blockDim,
                     vector<OCP_DBL>&               rhs,
                     vector<OCP_DBL>&               u) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

private:
    dCSRmat A; ///< Matrix for scalar-value problems
    dvector b; ///< Right-hand side for scalar-value problems
    dvector x; ///< Solution for scalar-value problems
};

/// Vector solvers in BSR format from FASP.
class VectorFaspSolver : public FaspSolver
{
    friend class LinearSystem;

private:
    /// Allocate memory for the linear system.
    void Allocate(const vector<USI>& rowCapacity,
                  const OCP_USI&     maxDim,
                  const USI&         blockDim) override;

    /// Initialize the Params for linear solver.
    void InitParam() override;

    /// Assemble coefficient matrix.
    void AssembleMat(const vector<vector<USI>>&     colId,
                     const vector<vector<OCP_DBL>>& val,
                     const OCP_USI&                 dim,
                     const USI&                     blockDim,
                     vector<OCP_DBL>&               rhs,
                     vector<OCP_DBL>&               u) override;

    /// Solve the linear system.
    OCP_INT Solve() override;

    /// Apply decoupling to the linear system.
    void Decoupling(dBSRmat* Absr,
                    dvector* b,
                    dBSRmat* Asc,
                    dvector* fsc,
                    ivector* order,
                    double*  Dmatvec,
                    int      decouple_type);

private:
    dBSRmat A; ///< Matrix for vector-value problems
    dvector b; ///< Right-hand side for vector-value problems
    dvector x; ///< Solution for vector-value problems

    dBSRmat Asc;   ///< Scaled matrix for vector-value problems
    dvector fsc;   ///< Scaled right-hand side for vector-value problems
    ivector order; ///< User-defined ordering for smoothing process

    vector<OCP_DBL> Dmat; ///< Decoupling matrices
};

#endif // __FASPSOLVER_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/08/2022      Update Doxygen                       */
/*  Chensong Zhang      Jan/19/2022      Set FASP4BLKOIL as optional          */
/*  Li Zhao             Apr/04/2022      Set FASP4CUDA   as optional          */
/*----------------------------------------------------------------------------*/
