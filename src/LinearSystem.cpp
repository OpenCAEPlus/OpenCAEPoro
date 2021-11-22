/*! \file    LinearSolver.cpp
 *  \brief   Linear solver for scalar-value problems
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "LinearSystem.hpp"

void LinearSystem::AllocateRowMem(const OCP_USI& dimMax, const USI& nb)
{
    blockSize = nb * nb;
    blockDim  = nb;
    maxDim    = dimMax;
    rowCapacity.resize(maxDim);
    colId.resize(maxDim);
    val.resize(maxDim);
    diagPtr.resize(maxDim);
    diagVal.resize(maxDim * blockSize);
    b.resize(maxDim * blockDim);
    u.resize(maxDim * blockDim);
}

void LinearSystem::AllocateColMem()
{
    for (OCP_USI n = 0; n < maxDim; n++) {
        colId[n].reserve(rowCapacity[n]);
        val[n].reserve(rowCapacity[n] * blockSize);
    }
}


void LinearSystem::ClearData()
{
    for (OCP_USI i = 0; i < maxDim; i++) {
        colId[i].clear();
        val[i].clear();
    }
    // diagPtr.assign(maxDim, 0);
    diagVal.assign(maxDim * blockSize, 0);
    b.assign(maxDim * blockDim, 0);
    // In fact, for linear system the current solution is a good initial solution for
    // next step, so u will not be set to zero. u.assign(maxDim, 0);
}


void LinearSystem::AssembleRhs(const vector<OCP_DBL>& rhs)
{
    OCP_USI nrow = dim * blockDim;
    for (OCP_USI i = 0; i < nrow; i++) {
        b[i] = rhs[i];
    }
}


void LinearSystem::OutputLinearSystem(const string& fileA, const string& fileb) const
{
    string FileA = solveDir + fileA;
    string Fileb = solveDir + fileb;

    
}

void LinearSystem::OutputSolution(const string& fileU) const
{
    string   FileU = solveDir + fileU;
    ofstream outu(FileU);
    if (!outu.is_open()) cout << "Can not open " << FileU << endl;
    outu << dim << endl;
    OCP_USI nrow = dim * blockDim;
    for (OCP_USI i = 0; i < nrow; i++) outu << u[i] << endl;
    outu.close();
}

void LinearSystem::CheckVal() const
{
    for (OCP_USI n = 0; n < dim; n++) {
        for (auto v : val[n]) {
            if (!isfinite(v)) {
                OCP_ABORT("NAN or INF in MAT");
            }
        }
    }
}


/// Setup LinearSolver
void LinearSystem::SetupLinearSolver(const USI& i, const string& dir, const string& file)
{
    solveDir = dir;
    switch (i)
    {
    case 1:
        // Fasp
        LS = new ScalarFaspSolver;
        break;
    case 2:
        // Blcok Fasp
        LS = new VectorFaspSolver;
        break;
    default:
        OCP_ABORT("Wrong Linear Solver type!");
        break;
    }
    LS->SetupParam(dir, file);
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/