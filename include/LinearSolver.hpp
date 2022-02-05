/*! \file    LinearSolver.hpp
 *  \brief   LinearSolver class declaration
 *  \author  Shizhe Li
 *  \date    Nov/22/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __LINEARSOLVER_HEADER__
#define __LINEARSOLVER_HEADER__

// Standard header files
#include <string>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"

using namespace std;

/// Virtual base class for linear solvers.
class LinearSolver
{
public:
    /// Read the params for linear solvers from an input file.
    virtual void SetupParam(const string& dir, const string& file) = 0;

    /// Initialize the params for linear solvers.
    virtual void InitParam() = 0;

    /// Allocate maximum memory for linear solvers.
    virtual void Allocate(const vector<USI>& rowCapacity, const OCP_USI& maxDim,
                          const USI& blockDim) = 0;

    /// Assemble matrix for linear solver from the internal matrix data.
    virtual void AssembleMat(const vector<vector<USI>>&     colId,
                             const vector<vector<OCP_DBL>>& val, const OCP_USI& dim,
                             const USI& blockDim, vector<OCP_DBL>& rhs,
                             vector<OCP_DBL>& u) = 0;

    /// Solve the linear system and return the number of iterations.
    virtual OCP_INT Solve(vector<OCP_DBL>& u) = 0;

    /// Get number of iterations.
    virtual USI GetNumIters() const = 0;
};

#endif // __LINEARSOLVER_HEADER__

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*  Chensong Zhang      Jan/18/2022      Update Doxygen                       */
/*----------------------------------------------------------------------------*/