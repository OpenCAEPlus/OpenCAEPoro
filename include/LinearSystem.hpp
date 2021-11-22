/*! \file    LinearSystem.hpp
 *  \brief   Linear solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __LINEARSYSTEM_HEADER__
#define __LINEARSYSTEM_HEADER__

// Standard header files
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>


// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "FaspSolver.hpp"


using namespace std;

/// Linear solvers for discrete systems.
//  Note: The matrix is stored in the form of row-segmented CSRx internaly, whose
//  sparsity pattern is almost the same as neighbor in BulkConn.
class LinearSystem
{
    friend class OpenCAEPoro;
    friend class BulkConn;
    friend class Well;

public:

    /// Allocate memory for linear system with max possible number of rows.
    void AllocateRowMem(const OCP_USI& dimMax, const USI& nb);
    /// Allocate memory for each matrix row with max possible number of columns.
    void AllocateColMem();
    /// Enlarge row capacity
    void EnlargeRowCap(const OCP_USI& row, const USI& n) { rowCapacity[row] += n; }
    /// Assign Rhs
    void AssembleRhs(const vector<OCP_DBL>& rhs);
    /// Clear the internal matrix data for scalar-value problems.
    void ClearData();
    

    /// Return the solution.
    vector<OCP_DBL>& GetSolution() { return u; }
    /// Check whether NAN or INF occurs in solution, used in debug mode.
    void CheckVal() const;
    /// Output the mat and rhs to fileA and fileb. //TODO: output to some obj?
    void OutputLinearSystem(const string& fileA, const string& fileb) const;
    /// Output the solution to a disk file name.
    void OutputSolution(const string& filename) const;
 

    // Linear Solver
    /// Setup LinearSolver
    void SetupLinearSolver(const USI& i, const string& dir, const string& file);
    /// Allocate memory for Linear Solver
    void AllocateLinearSolver() { LS->Allocate(rowCapacity, maxDim, blockDim); }
    /// Assemble Mat for Linear Solver
    void AssembleMatLinearSolver() { LS->AssembleMat(colId, val, dim, blockDim, b, u); }
    /// Solve the Linear System
    OCP_INT Solve() { return LS->Solve(u); }

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

    string                  solveDir;    ///< Current workdir.

    LinearSolver*           LS;

};

#endif /* end if __LINEARSOLVER_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Nov/09/2021      Remove decoupling methods            */
/*  Chensong Zhang      Nov/22/2021      renamed to LinearSystem              */
/*----------------------------------------------------------------------------*/