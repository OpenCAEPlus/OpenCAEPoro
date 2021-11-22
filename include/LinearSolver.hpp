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
#include <vector>
#include <string>

// OpenCAEPoro header files
#include "OCPConst.hpp"


using namespace std;

/// virtual Basic class for linear solver
class LinearSolver
{
public:
	/// Read Linear Solver Param form input file
	virtual void SetupParam(const string& dir, const string& file) = 0;
	/// Initialize the Params for Linear Solver
	virtual void InitParam() = 0;
	/// Allocate Maximum memory for Linear Solver according to Internal matrix
	virtual void Allocate(const vector<USI>& rowCapacity, const OCP_USI& maxDim, const USI& blockDim) = 0;
	/// Assemble Matrix for Linear Solver from Internal matrix
	virtual void AssembleMat(const vector<vector<USI>>& colId, const vector<vector<OCP_DBL>>& val,
		const OCP_USI& dim, const USI& blockDim, vector<OCP_DBL>& rhs, vector<OCP_DBL>& u) = 0;
	/// Solve the Linear System and return the Iterations and Solutions.
	virtual OCP_INT Solve(vector<OCP_DBL>& u) = 0;
};




#endif

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Nov/22/2021      Create file                          */
/*----------------------------------------------------------------------------*/