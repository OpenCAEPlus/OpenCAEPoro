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


#ifndef __FASPSOLVER_HEADER__
#define __FASPSOLVER_HEADER__


// Standard header files
#include <vector>
#include <string>
#include <fstream>
#include <iostream>


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
#include "LinearSolver.hpp"


// Preconditioner types
#define PC_NULL 60        ///< None: no preconditioner
#define PC_FASP1 61       ///< FASP1 preconditioner: default for FIM from 2020
#define PC_FASP2 62       ///< FASP2 preconditioner: experimental
#define PC_FASP3 63       ///< FASP3 preconditioner: monolithic
#define PC_FASP4 64       ///< FASP4 preconditioner: default for FIM from 2015
#define PC_FASP5 65       ///< FASP5 preconditioner: experimental
#define PC_DIAG 68        ///< DIAG preconditioner
#define PC_BILU 69        ///< BILU preconditioner
#define PC_FASP1_SHARE 71 ///< Sharing the setup stage for PC_FASP1
#define PC_FASP4_SHARE 74 ///< Sharing the setup stage for PC_FASP4
#define RESET_CONST 35    ///< Sharing threshold for PC_FASP1_SHARE and PC_FASP4_SHARE


using namespace std;

/// Basic virtual Fasp solver
class FaspSolver : public LinearSolver
{
public:
	void SetupParam(const string& dir, const string& file) override;

public:
	string      solveDir;  ///< Current work dir
	string      solveFile; ///< Relative path of fasp file.
	input_param inParam;   ///< Parameters from input files.
	ITS_param   itParam;   ///< Parameters for iterative method.
	AMG_param   amgParam;  ///< Parameters for AMG method.
	ILU_param   iluParam;  ///< Parameters for ILU method.
	SWZ_param   swzParam;  ///< Parameters for Schwarz method.
};


/// General Fasp Solver
class ScalarFaspSolver : public FaspSolver
{
	friend class LinearSystem;

private:

	/// Allocate memory for scalar-value problems in FASP.
	void Allocate(const vector<USI>& rowCapacity, const OCP_USI& maxDim, const USI& blockDim) override;
	/// Initialize the Params for Fasp
	void InitParam() override;
	/// Assemble Matrix for Fasp
	void AssembleMat(const vector<vector<USI>>& colId, const vector<vector<OCP_DBL>>& val,
					const OCP_USI& dim, const USI& blockDim, vector<OCP_DBL>& rhs, vector<OCP_DBL>& u) override;
	/// Solve the Linear system 
	OCP_INT Solve(vector<OCP_DBL>& u) override;
	

private:

	dCSRmat         A;  ///< Matrix for scalar-value problems.
	dvector         b;  ///< Right-hand side for scalar-value problems.
	dvector         x;  ///< Solution for scalar-value problems.

};


/// Block Fasp Solver
class VectorFaspSolver : public FaspSolver
{
	friend class LinearSystem;

private:
	/// Allocate memory for scalar-value problems in BFASP.
	void Allocate(const vector<USI>& rowCapacity, const OCP_USI& maxDim, const USI& blockDim) override;
	/// Initialize the Params for BFasp
	void InitParam() override;
	/// Assemble Matrix for BFasp
	void AssembleMat(const vector<vector<USI>>& colId, const vector<vector<OCP_DBL>>& val,
		const OCP_USI& dim, const USI& blockDim, vector<OCP_DBL>& rhs, vector<OCP_DBL>& u) override;
	/// Solve the Linear system
	OCP_INT Solve(vector<OCP_DBL>& u) override;
	/// Decouple the matrix
	void Decoupling(dBSRmat* Absr, dvector* b, dBSRmat* Asc, dvector* fsc,
		ivector* order, double* Dmatvec, int decouple_type);
private:

	dBSRmat         A; ///< Matrix for vector-value problems.
	dvector         b; ///< Right-hand side for vector-value problems.
	dvector         x; ///< Solution for vector-value problems.
	dBSRmat         Asc;     ///< Scaled matrix for vector-value problems.
	dvector         fsc;     ///< Scaled right-hand side for vector-value problems.
	ivector         order;   ///< User-defined ordering for smoothing process.
	vector<OCP_DBL> Dmat;    // TODO: What is this for???
};



#endif