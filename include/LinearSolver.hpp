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

#ifndef __LINEARSOLVER_HEADER__
#define __LINEARSOLVER_HEADER__

// Standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

// faspsolver header files
extern "C" {
#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
}

// fasp4blkoil header files
extern "C" {
#include "fasp4blkoil.h"
#include "fasp4blkoil_functs.h"
}

// OpenCAEPoro header files
#include "DenseMat.hpp"
#include "LinearSystem.hpp"
#include "OCPConst.hpp"

// Preconditioner types
#define PC_NULL 60  ///< None: no preconditioner
#define PC_FASP1 61 ///< FASP1 preconditioner: default for FIM from 2020
#define PC_FASP2 62 ///< FASP2 preconditioner: experimental
#define PC_FASP3 63 ///< FASP3 preconditioner: monolithic
#define PC_FASP4 64 ///< FASP4 preconditioner: default for FIM from 2015
#define PC_FASP5 65 ///< FASP5 preconditioner: experimental
#define PC_DIAG 68  ///< DIAG preconditioner
#define PC_BILU 69  ///< BILU preconditioner

// Preconditioners types with shared setup
#define PC_FASP1_SHARE 71 ///< Sharing the setup stage for PC_FASP1
#define PC_FASP4_SHARE 74 ///< Sharing the setup stage for PC_FASP4
#define RESET_CONST 35    ///< Sharing threshold for PC_FASP1_SHARE and PC_FASP4_SHARE

using namespace std;

/// A template class designed to stores and solves linear system derived from discrete
/// method. the maxtrix is stored in the form of row-segmented CSRX internaly, whose
/// sparsity pattern is almost the same as neighbor in BulkConn.
class LinearSolver
{
    friend class OpenCAEPoro;
    friend class BulkConn;
    friend class Well;

public:
    /// Allocate memory for linear system with max possible number of rows.
    void AllocateRowMem(const OCP_USI& dimMax, const USI& nb);
    /// Allocate memory for each matrix row with max possible number of columns.
    void AllocateColMem();
    /// Allocate memory for scalar-value problems in FASP.
    void AllocateFasp();
    /// Allocate memory for vector-value problems in FASP.
    void AllocateBFasp();
    /// Enlarge row capacity. // TODO: maybe should be called EnlargeRowCap?
    void RowCapPlus(const OCP_USI& row, const USI& n) { rowCapacity[row] += n; }

    /// Setup solver param for scalar-value problems from input file.
    void SetupParam(const string& dir, const string& file);
    /// Setup solver param for vector-value problems from input file.
    void SetupParamB(const string& dir, const string& file);

    /// Return the solution.
    vector<OCP_DBL>& GetSolution() { return u; }
    /// Check whether NAN or INF occurs in solution, used in debug mode.
    void CheckVal() const;
    /// Output the mat and rhs to fileA and fileb. //TODO: output to some obj?
    void OutputLinearSystem(const string& fileA, const string& fileb) const;
    /// Output the solution to a disk file name.
    void OutputSolution(const string& filename) const;

    //---------------------------------//
    // FASP: for scalar-value problems //
    //---------------------------------//

    /// Initialize the solver param for FASP with default values.
    void InitParam_Fasp();
    /// Read the solver param for FASP from a file.
    void ReadParam_Fasp(); // TODO: use a file name for input param?
    /// Convert the internal matrix into FASP format.
    void AssembleMat_Fasp(); // TODO: maybe called AssembleFaspMat? Why no RHS (copy)?
    /// Solve the linear system using FASP and return the status.
    int FaspSolve();
    /// Clear the internal matrix data for scalar-value problems.
    void ClearData();

    //----------------------------------//
    // BFASP: for vector-value problems //
    //----------------------------------//

    /// Initialize the solver param for BFASP with default values.
    void InitParam_BFasp();
    /// Read the solver param for BFASP from a file.
    void ReadParam_BFasp();
    /// Convert the internal matrix into BFASP format.
    void AssembleMat_BFasp();
    /// Convert the internal right-hand side into BFASP format.
    void AssembleRhs_BFasp(const vector<OCP_DBL>& rhs);
    /// Solve the linear system using BFASP and return the status.
    int BFaspSolve();
    /// Clear the internal matrix data for vector-value problems.
    void ClearDataB();

private: // TODO: Maybe should not be here?
    void decoupling(dBSRmat* Absr, dvector* b, int scal_type, dBSRmat* Asc,
                    dvector* fsc, ivector* order, double* Dmatvec, int decouple_type);
    void decouple_abf(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_anl(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_truetrans_alg(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_truetrans(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_quasi(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_trueabf(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_rowsum(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_abftrue(dBSRmat* A, REAL* diaginv, dBSRmat* Asc);
    void decouple_true_scale(dBSRmat* A, REAL* diaginv, dBSRmat* B);
    void decouple_rotate(dBSRmat* A, REAL* diaginv, dBSRmat* B);

private:
    // Used for internal mat structure.
    USI blockDim;  ///< Dimens of small block matrix.
    USI blockSize; ///< Size of small block matrix. // TODO: Is it blockDim*blockDim?

    // maxDim: fixed and used to allocate memory at the beginning of simulation;
    // dim:    might change during simulation but always less than maxDim.
    OCP_USI maxDim; ///< Maximal possible dimension of matrix.
    OCP_USI dim;    ///< Actual dimension of matrix.

    // The following values are stored for each row. Among them, rowCapacity is the max
    // possible capacity of each row of the matrix. It is just a little bigger than the
    // actual size and is used to allocate memory at the beginning of simulation.
    // diagVal is an auxiliary variable used to help setup entries in diagnal line and
    // it will only be used when matrix is assembled.
    vector<USI>             rowCapacity; ///< Maximal capacity of each row.
    vector<vector<OCP_USI>> colId;       ///< Column indices of nonzero entry.
    vector<USI>             diagPtr;     ///< Indices of diagonal entries.
    vector<vector<OCP_DBL>> val;         ///< Nonzero values.
    vector<OCP_DBL>         diagVal;     ///< Diagonal values
    vector<OCP_DBL>         b;           ///< Right-hand side of linear system.
    vector<OCP_DBL>         u;           ///< Solution of linear system.

    // for FASP solver
    string  solveDir;
    string  solveFile;
    dCSRmat A_Fasp;
    dvector b_Fasp;
    dvector x_Fasp;

    // for Block FASP solver
    dBSRmat         A_BFasp;
    dvector         b_BFasp;
    dvector         x_BFasp;
    dBSRmat         Asc;
    dvector         fsc;
    ivector         order;
    vector<OCP_DBL> Dmat;

    input_param inParam;  ///< Parameters from input files.
    ITS_param   itParam;  ///< Parameters for iterative method.
    AMG_param   amgParam; ///< Parameters for AMG method.
    ILU_param   iluParam; ///< Parameters for ILU method.
    SWZ_param   swzParam; ///< Parameters for Schwarz method.
};

#endif /* end if __LINEARSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Nov/09/2021      Rewrite Doxygen                      */
/*----------------------------------------------------------------------------*/