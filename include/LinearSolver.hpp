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

// Fasp header files
#include "fasp.h"
extern "C" {
#include "fasp_functs.h"
}

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
    vector<OCP_DBL>& GetSolution() { return u; }

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