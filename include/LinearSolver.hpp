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

// OpenCAEPoro header files
#include "LinearSystem.hpp"
#include "OCPConst.hpp"
#include "DenseMat.hpp"

// Fasp header files
#include "fasp.h"
extern "C" {
#include "fasp_functs.h"
#include "fasp4blkoil.h"
#include "fasp4blkoil_functs.h"
}

#define  PC_NULL         60 // no preconditioner
#define  PC_FASP1        61 // default preconditioner for FIM from 2020
#define  PC_FASP2        62
#define  PC_FASP3        63
#define  PC_FASP4        64 // default preconditioner for FIM from 2015
#define  PC_FASP5        65
#define  PC_DIAG         68 // diagonal preconditioner
#define  PC_BILU         69 // block incomplete factorization preconditioner
#define  PC_FASP1_SHARE  71 // Share the setup stage for PC_FASP1
#define  PC_FASP4_SHARE  74 // Share the setup stage for PC_FASP4
#define  RESET_CONST     35 // Sharing threshold for PC_FASP1_SHARE and PC_FASP4_SHARE

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
    /// Allocate memory for linear system where possible maximum numbers of row are
    /// used.
    void AllocateRowMem(const OCP_USI& dimMax, const USI& nb);
    /// Allocate memory for each row of matrix where the most terrible condition was
    /// considered.
    void AllocateColMem();
    void AllocateFasp();
    void AllocateBFasp();

    void RowCapPlus(const OCP_USI& row, const USI& n) { rowCapacity[row] += n; }
    /// Read solver param from input file which is supplied by user.
    void SetupParam(const string& dir, const string& file);
    /// Read solver param for Block matrix from input file which is supplied by user.
    void SetupParamB(const string& dir, const string& file);
    /// output the solution to fileX.
    void PrintfSolution(const string& fileX) const;

    // FASP
    /// initialize the sover param for FASP.
    void InitParam_Fasp();
    /// read the sover param for Block FASP.
    void ReadParam_Fasp();
    /// convert the internal matrix structure into the the format required by FASP.
    void AssembleMat_Fasp();
    /// solve the linear system by FASP and return the status.
    int FaspSolve();
    // Block Fasp
    void InitParam_BFasp();
    void ReadParam_BFasp();
    void AssembleMat_BFasp();
    void AssembleRhs_BFasp(const vector<OCP_DBL>& rhs);
    int BFaspSolve();
    void decoupling(dBSRmat* Absr, dvector* b, int scal_type,
        dBSRmat* Asc, dvector* fsc, ivector* order,
        double* Dmatvec, int decouple_type);
    void decouple_abf(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_anl(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_truetrans_alg(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_truetrans(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_quasi(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_trueabf(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_rowsum(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_abftrue(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* Asc);

    void decouple_true_scale(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* B);

    void decouple_rotate(dBSRmat* A,
        REAL* diaginv,
        dBSRmat* B);

    /// output the mat and rhs to fileA and fileb.
    void PrintfMatCSR(const string& fileA, const string& fileb) const;
    /// clear the internal matrix.
    void ClearData();
    /// return the solution.
    vector<OCP_DBL>& GetSolution() { return u; }

    /// check if NAN or INF occurs, used in debug mode.
    void CheckVal() const;

private:
    // internal mat structure.
    USI blockDim;  ///< Dimens of small block matrix.
    USI blockSize; ///< Size of small block matrix.
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
    vector<vector<OCP_DBL>> val; ///< values of nonzero entry in the format of row-segmented.
    vector<USI> diagPtr;   ///< the ith entry indicates the location of diagal entry of
                           ///< the ith row in val.
    /// an auxiliary variable used to help Setup entry in diagnal line.
    /// it will only be used when matrix is being assembling.
    vector<OCP_DBL> diagVal;
    vector<OCP_DBL> b; ///< right hands of linear system.
    vector<OCP_DBL> u; ///< solutiom of linear system.

    // for FASP solver
    string  solveDir;
    string  solveFile;
    dCSRmat A_Fasp;
    dvector b_Fasp;
    dvector x_Fasp;
    // for Block FASP solver
    dBSRmat A_BFasp;
    dvector b_BFasp;
    dvector x_BFasp;
    dBSRmat Asc;
    dvector fsc;
    ivector order;
    vector<OCP_DBL> Dmat;

    input_param inParam; // parameters from input files
    ITS_param   itParam; // parameters for itsolver

    AMG_param amgParam; // parameters for AMG
    ILU_param iluParam; // parameters for ILU
    SWZ_param swzParam; // parameters for Schwarz method

    // for other solver
};



#endif   /* end if __LINEARSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/