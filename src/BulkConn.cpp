/*! \file    BulkConn.cpp
 *  \brief   BulkConn class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// Standard header files
#include <cassert>
#include <cmath>
#include <ctime>

// OpenCAEPoro header files
#include "BulkConn.hpp"

/////////////////////////////////////////////////////////////////////
// General
/////////////////////////////////////////////////////////////////////

/// It should be called after Grid and Bulk Setup.
void BulkConn::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;

    numConn = 0;
    numBulk = myGrid.activeGridNum;

    neighbor.resize(numBulk);
    selfPtr.resize(numBulk);
    neighborNum.resize(numBulk);

    vector<GPair> tmp1, tmp2;
    USI           len;
    OCP_USI       bIdb;

    for (OCP_USI n = 0; n < myGrid.numGrid; n++) {
        const GB_Pair& GBtmp = myGrid.activeMap_G2B[n];

        if (GBtmp.IsAct()) {
            bIdb = GBtmp.GetId();

            // Get rid of inactive neighbor
            tmp1 = myGrid.gNeighbor[n];
            len  = tmp1.size();
            tmp2.clear();
            for (USI i = 0; i < len; i++) {
                const GB_Pair& GBtmp2 = myGrid.activeMap_G2B[tmp1[i].id];

                if (GBtmp2.IsAct()) {
                    tmp1[i].id = GBtmp2.GetId();
                    tmp2.push_back(tmp1[i]);
                }
            }
            // Add Self
            tmp2.push_back(GPair(bIdb, 0.0));
            // Sort: Ascending
            sort(tmp2.begin(), tmp2.end(), GPair::lessG);
            // Find SelfPtr and Assign to neighbor and area
            len = tmp2.size();
            for (USI i = 0; i < len; i++) {
                neighbor[bIdb].push_back(tmp2[i].id);
                if (tmp2[i].id == bIdb) {
                    selfPtr[bIdb] = i;
                }
            }
            for (USI j = selfPtr[bIdb] + 1; j < len; j++) {
                area.push_back(tmp2[j].area);
                iteratorConn.push_back(BulkPair(bIdb, tmp2[j].id));
            }
            neighborNum[bIdb] = len;
        }
    }

    numConn = iteratorConn.size();

    // PrintConnectionInfoCoor(myGrid);
}

void BulkConn::CalIteratorConn()
{
    OCP_FUNCNAME;

    iteratorConn.reserve(numConn);
    // generate iterator for BB from iteratorConn
    for (OCP_USI bId = 0; bId < numBulk; bId++) {
        USI beginIt = selfPtr[bId] + 1;
        USI nbc     = neighborNum[bId];

        for (USI c = beginIt; c < nbc; c++) {
            OCP_USI eId = neighbor[bId][c];
            iteratorConn.push_back(BulkPair(bId, eId));
        }
    }

    assert(iteratorConn.size() == numConn);
}

/// This method should be called only once at the beginning.
void BulkConn::AllocateMat(LinearSystem& MySolver) const
{
    OCP_FUNCNAME;

    for (OCP_USI n = 0; n < numBulk; n++) {
        MySolver.rowCapacity[n] += neighborNum[n];
    }
}

void BulkConn::SetupMatSparsity(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    myLS.dim = numBulk;
    for (OCP_USI n = 0; n < numBulk; n++) {
        myLS.colId[n].assign(neighbor[n].begin(), neighbor[n].end());
        myLS.diagPtr[n] = selfPtr[n];
    }
}

void BulkConn::UpdateLastStep()
{
    OCP_FUNCNAME;

    lastUpblock          = upblock;
    lastUpblock_Rho      = upblock_Rho;
    lastUpblock_Trans    = upblock_Trans;
    lastUpblock_Velocity = upblock_Velocity;
}

void BulkConn::Reset()
{
    OCP_FUNCNAME;

    upblock          = lastUpblock;
    upblock_Rho      = lastUpblock_Rho;
    upblock_Trans    = lastUpblock_Trans;
    upblock_Velocity = lastUpblock_Velocity;
}

void BulkConn::CheckDiff() const
{
    OCP_FUNCNAME;

    // upblock
    OCP_DBL tmp;
    cout << setprecision(18);
    for (OCP_USI c = 0; c < numConn; c++) {
        tmp = fabs(upblock[c] - lastUpblock[c]);
        if (tmp < TINY) {
            cout << ">> Difference in upblock index at \t" << tmp << "\n";
        }
        tmp = fabs(upblock_Rho[c] - lastUpblock_Rho[c]);
        if (tmp < TINY) {
            cout << ">> Difference in upblock Rho at \t" << tmp << "\n";
        }
        tmp = fabs(upblock_Trans[c] - lastUpblock_Trans[c]);
        if (tmp < TINY) {
            cout << ">> Difference in upblock Trans at \t" << tmp << "\n";
        }
        tmp = fabs(upblock_Velocity[c] - lastUpblock_Velocity[c]);
        if (tmp < TINY) {
            cout << ">> Difference in upblock Velocity at \t" << tmp << "\n";
        }
    }
}

void BulkConn::PrintConnectionInfo(const Grid& myGrid) const
{
    for (OCP_USI i = 0; i < numBulk; i++) {
        cout << "(" << myGrid.activeMap_B2G[i] << ")"
             << "\t";

        for (auto& v : neighbor[i]) {
            cout << myGrid.activeMap_B2G[v] << "\t";
        }
        cout << "[" << selfPtr[i] << "]";
        cout << "\t" << neighborNum[i];
        cout << "\n";
    }

    for (OCP_USI i = 0; i < numConn; i++) {
        cout << myGrid.activeMap_B2G[iteratorConn[i].BId] 
            << "\t" << myGrid.activeMap_B2G[iteratorConn[i].EId] << "\n";
    }
}

void BulkConn::PrintConnectionInfoCoor(const Grid& myGrid) const
{
    OCP_USI bIdg, eIdg;
    USI I, J, K;
    cout << "BulkConn : " << numConn << endl;
    for (OCP_USI c = 0; c < numConn; c++) {
        bIdg = myGrid.activeMap_B2G[iteratorConn[c].BId];
        eIdg = myGrid.activeMap_B2G[iteratorConn[c].EId];
        myGrid.GetIJKGrid(I, J, K, bIdg);
        cout << "(" << setw(3) << I << "," << setw(3) << J << "," << setw(3) << K << ")    ";
        cout << setw(6) << bIdg;
        cout << "       ";
        myGrid.GetIJKGrid(I, J, K, eIdg);
        cout << "(" << setw(3) << I << "," << setw(3) << J << "," << setw(3) << K << ")    ";
        cout << setw(6) << eIdg;
        cout << setw(20) << setprecision(8) << fixed << area[c] * CONV2;

        cout << endl;
    }
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void BulkConn::AllocateAuxIMPEC(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    upblock_Trans.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);
    lastUpblock.resize(numConn * np);
    lastUpblock_Rho.resize(numConn * np);
    lastUpblock_Trans.resize(numConn * np);
    lastUpblock_Velocity.resize(numConn * np);
}

void BulkConn::AssembleMatIMPEC(LinearSystem& myLS, const Bulk& myBulk,
                                const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    // accumulate term
    OCP_DBL Vp0, Vp, vf, vfp, P;
    OCP_DBL cr = myBulk.rockC1;
    for (OCP_USI n = 0; n < numBulk; n++) {
        Vp0 = myBulk.rockVpInit[n];
        Vp  = myBulk.rockVp[n];
        vfp = myBulk.vfp[n];
        P   = myBulk.P[n];
        vf  = myBulk.vf[n];

        OCP_DBL temp    = cr * Vp0 - vfp;
        myLS.diagVal[n] = temp;
        myLS.b[n]       = temp * P + dt * (vf - Vp);
    }

    // flux term
    OCP_USI bId, eId, uId;
    USI     np = myBulk.numPhase;
    USI     nc = myBulk.numCom;
    OCP_DBL valupi, valdowni;
    OCP_DBL valup, rhsup, valdown, rhsdown;

    // Be careful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId     = iteratorConn[c].BId;
        eId     = iteratorConn[c].EId;
        valup   = 0;
        rhsup   = 0;
        valdown = 0;
        rhsdown = 0;

        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];
            if (!myBulk.phaseExist[uId * np + j]) continue;

            valupi   = 0;
            valdowni = 0;

            for (USI i = 0; i < nc; i++) {
                valupi +=
                    myBulk.vfi[bId * nc + i] * myBulk.xij[uId * np * nc + j * nc + i];
                valdowni +=
                    myBulk.vfi[eId * nc + i] * myBulk.xij[uId * np * nc + j * nc + i];
            }
            OCP_DBL dD   = myBulk.depth[bId] - myBulk.depth[eId];
            OCP_DBL dPc  = myBulk.Pc[bId * np + j] - myBulk.Pc[eId * np + j];
            OCP_DBL temp = myBulk.xi[uId * np + j] * upblock_Trans[c * np + j] * dt;
            valup += temp * valupi;

            // if (!isfinite(valup)) {
            //    cout << "###ERROR   ";
            //    ERRORcheck("NAN or INF in MAT");
            //}

            valdown += temp * valdowni;
            temp *= upblock_Rho[c * np + j] * GRAVITY_FACTOR * dD - dPc;
            rhsup += temp * valupi;
            rhsdown -= temp * valdowni;
        }

        USI diagptr = selfPtr[bId];
        if (bId != lastbId) {
            // new bulk
            assert(myLS.val[bId].size() == diagptr);
            myLS.val[bId].push_back(myLS.diagVal[bId]);
            lastbId = bId;
        }

        //if (!isfinite(valup) || !isfinite(valdown)) {
        //    cout << "NAN or INF in MAT" << endl;
        //}

        myLS.val[bId][diagptr] += valup;
        myLS.val[bId].push_back(-valup);
        myLS.val[eId].push_back(-valdown);
        myLS.diagVal[eId] += valdown;
        myLS.b[bId] += rhsup;
        myLS.b[eId] += rhsdown;
    }

    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == selfPtr[n]) myLS.val[n].push_back(myLS.diagVal[n]);
    }
}

OCP_DBL BulkConn::CalCFLIMPEC(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI     np   = myBulk.numPhase;
    OCP_DBL cfl  = 0;
    OCP_DBL temp = 0;
    for (OCP_USI c = 0; c < numConn; c++) {

        for (USI j = 0; j < np; j++) {
            OCP_USI uId = upblock[c * np + j];

            if (myBulk.phaseExist[uId * np + j]) { // uId -> uId * np + j  fix bugs.
                temp = fabs(upblock_Velocity[c * np + j]) * dt;
                temp /= myBulk.vj[uId * np + j];
                if (cfl < temp) cfl = temp;
            }
        }
    }
    return cfl;
}

void BulkConn::CalCFL01IMPEC(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI np = myBulk.numPhase;
    for (OCP_USI c = 0; c < numConn; c++) {

        for (USI j = 0; j < np; j++) {
            OCP_USI uId = upblock[c * np + j];

            if (myBulk.phaseExist[uId * np + j]) {
                myBulk.cfl[uId * np + j] += fabs(upblock_Velocity[c * np + j]) * dt;
            }
        }
    }
}

void BulkConn::CalFluxIMPEC(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    // calculate a step flux using iteratorConn
    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j;
    OCP_DBL Pbegin, Pend, rho;
    USI     np = myBulk.numPhase;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId         = iteratorConn[c].BId;
        eId         = iteratorConn[c].EId;
        OCP_DBL Akd = area[c];

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            bool exbegin = myBulk.phaseExist[bId_np_j];
            bool exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.P[eId];
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                Pbegin = myBulk.P[bId];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]          = bId;
                upblock_Trans[c * np + j]    = 0;
                upblock_Velocity[c * np + j] = 0;
                // if (!isfinite(upblock_Trans[c * np + j])) {
                //    cout << "###ERROR   ";
                //    ERRORcheck("NAN or INF in MAT");
                //}
                continue;
            }

            uId          = bId;
            bool    exup = exbegin;
            OCP_DBL dP   = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                         (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            }

            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;

            // (fix bugs) if exup is false, then trans and vec are zero.
            if (exup) {
                OCP_USI uId_np_j = uId * np + j;
                OCP_DBL trans =
                    CONV1 * CONV2 * Akd * myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j];
                upblock_Trans[c * np + j]    = trans;
                upblock_Velocity[c * np + j] = trans * dP;
            } else {
                upblock_Trans[c * np + j]    = 0;
                upblock_Velocity[c * np + j] = 0;
            }
        }
    }
}

void BulkConn::MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;

    for (OCP_USI c = 0; c < numConn; c++) {
        OCP_USI bId = iteratorConn[c].BId;
        OCP_USI eId = iteratorConn[c].EId;

        for (USI j = 0; j < np; j++) {
            OCP_USI uId      = upblock[c * np + j];
            OCP_USI uId_np_j = uId * np + j;

            if (!myBulk.phaseExist[uId_np_j]) continue;

            OCP_DBL phaseVelocity = upblock_Velocity[c * np + j];
            for (USI i = 0; i < nc; i++) {
                OCP_DBL dNi = dt * phaseVelocity * myBulk.xi[uId_np_j] *
                              myBulk.xij[uId_np_j * nc + i];
                myBulk.Ni[eId * nc + i] += dNi;
                myBulk.Ni[bId * nc + i] -= dNi;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void BulkConn::AllocateAuxFIM(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    lastUpblock.resize(numConn * np);
}

void BulkConn::AssembleMat_FIM(LinearSystem& myLS, const Bulk& myBulk,
                               const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI np     = myBulk.numPhase;
    const USI nc     = myBulk.numCom;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockC1 * myBulk.rockVpInit[n] - myBulk.vfp[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }
        for (USI i = 0; i < bsize; i++) {
            myLS.diagVal[n * bsize + i] = bmat[i];
        }
    }
    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);

    OCP_USI bId, eId, uId;
    OCP_USI uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        Akd    = CONV1 * CONV2 * area[c];
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        for (USI j = 0; j < np; j++) {
            uId      = upblock[c * np + j];
            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            dP = myBulk.Pj[bId * np + j] - myBulk.Pj[eId * np + j] -
                 myBulk.rho[bId * np + j] * dGamma;
            xi     = myBulk.xi[uId_np_j];
            kr     = myBulk.kr[uId_np_j];
            mu     = myBulk.mu[uId_np_j];
            muP    = myBulk.muP[uId_np_j];
            xiP    = myBulk.xiP[uId_np_j];
            rhoP   = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[uId_np_j * nc + i];

                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;
                tmp = transIJ * (-rhoP * dGamma);
                tmp += xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (bId == uId) {
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpE[(i + 1) * ncol] += tmp;
                }

                // Saturation -- Second var
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] +=
                        transIJ * myBulk.dPcj_dS[bId * np * np + j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] +=
                        transIJ * myBulk.dPcj_dS[eId * np * np + j * np + k];
                    tmp = Akd * xij * xi / mu *
                          myBulk.dKr_dS[uId * np * np + j * np + k] * dP;
                    if (bId == uId) {
                        dFdXsB[(i + 1) * ncol2 + k] += tmp;
                    } else {
                        dFdXsE[(i + 1) * ncol2 + k] += tmp;
                    }
                }
                // Cij -- Second var
                for (USI k = 0; k < nc; k++) {
                    rhox = myBulk.rhox[uId_np_j * nc + k];
                    xix  = myBulk.xix[uId_np_j * nc + k];
                    mux  = myBulk.mux[uId_np_j * nc + k];
                    tmp  = -transIJ * rhox * dGamma;
                    tmp += xij * transJ * xix * dP;
                    tmp += -transIJ * mux / mu * dP;
                    if (k == i) {
                        tmp += xi * transJ * dP;
                    }
                    if (bId == uId) {
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    } else {
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                }
            }
        }

        USI diagptr = selfPtr[bId];

        if (bId != lastbId) {
            // new bulk
            assert(myLS.val[bId].size() == diagptr * bsize);
            OCP_USI id = bId * bsize;
            myLS.val[bId].insert(myLS.val[bId].end(), myLS.diagVal.data() + id,
                                 myLS.diagVal.data() + id + bsize);

            lastbId = bId;
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin
        // Add
        for (USI i = 0; i < bsize; i++) {
            myLS.val[bId][diagptr * bsize + i] += bmat[i];
        }
        // End
        // Insert
        Dscalar(bsize, -1, bmat.data());
        myLS.val[eId].insert(myLS.val[eId].end(), bmat.begin(), bmat.end());

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin
        // Insert
        myLS.val[bId].insert(myLS.val[bId].end(), bmat.begin(), bmat.end());
        // Add
        Dscalar(bsize, -1, bmat.data());
        for (USI i = 0; i < bsize; i++) {
            myLS.diagVal[eId * bsize + i] += bmat[i];
        }
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == selfPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                               myLS.diagVal.data() + n * bsize + bsize);
    }
}

void BulkConn::CalFluxFIM(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    OCP_USI   bId, eId, uId;
    OCP_USI   bId_np_j, eId_np_j;
    OCP_DBL   Pbegin, Pend, rho;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId         = iteratorConn[c].BId;
        eId         = iteratorConn[c].EId;
        OCP_DBL Akd = area[c];

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            bool exbegin = myBulk.phaseExist[bId_np_j];
            bool exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.P[eId];
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                Pbegin = myBulk.P[bId];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j] = bId;
                continue;
            }

            uId        = bId;
            OCP_DBL dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                         (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock[c * np + j] = uId;
        }
    }
}

void BulkConn::CalResFIM(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    const USI np  = myBulk.numPhase;
    const USI nc  = myBulk.numCom;
    const USI len = nc + 1;
    OCP_USI   bId, eId, uId, bIdb;

    // Accumalation Term
    for (OCP_USI n = 0; n < numBulk; n++) {

        bId  = n * len;
        bIdb = n * nc;

        res[bId] = myBulk.rockVp[n] - myBulk.vf[n];
        for (USI i = 0; i < nc; i++) {
            res[bId + 1 + i] = myBulk.Ni[bIdb + i] - myBulk.lNi[bIdb + i];
        }
    }

    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL Pbegin, Pend, rho, dP;
    OCP_DBL tmp, dNi;
    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            bool exbegin = myBulk.phaseExist[bId_np_j];
            bool exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend   = myBulk.P[eId];
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                Pbegin = myBulk.P[bId];
                Pend   = myBulk.Pj[eId_np_j];
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j] = bId;
                // upblock_Rho[c * np + j] = rho;
                continue;
            }

            uId = bId;
            dP  = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                 (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            // upblock_Rho[c * np + j] = rho;
            upblock[c * np + j] = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            tmp = dt * CONV1 * CONV2 * area[c] * myBulk.xi[uId_np_j] *
                  myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
            for (USI i = 0; i < nc; i++) {
                dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                res[bId * len + 1 + i] += dNi;
                res[eId * len + 1 + i] -= dNi;
            }
        }
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/