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

#include <cassert>
#include <cmath>
#include <ctime>

#include "BulkConn.hpp"

// Active Conn & Active Bulk

void BulkConn::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    InitSize(myBulk);
    CalConn(myGrid, myBulk.numPhase);
    CalIteratorConn();
    CalArea(myGrid, myBulk);
}

void BulkConn::AllocateConnIMPEC(const USI& np)
{
    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    upblock_Trans.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);
    lastUpblock.resize(numConn * np);
    lastUpblock_Rho.resize(numConn * np);
    lastUpblock_Trans.resize(numConn * np);
    lastUpblock_Velocity.resize(numConn * np);
}


void BulkConn::AllocateConnFIM(const USI& np)
{
    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    // useless
    upblock_Trans.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);

}

void BulkConn::InitSize(const Bulk& myBulk)
{
    numConn = 0;
    numBulk = myBulk.numBulk;

    neighbor.resize(numBulk);
    selfPtr.resize(numBulk);
    neighborNum.resize(numBulk);
}

void BulkConn::CalConn(const Grid& myGrid, const USI& np)
{
    USI     nx   = myGrid.nx;
    USI     ny   = myGrid.ny;
    USI     nz   = myGrid.nz;
    OCP_USI nxny = nx * ny;

    // bIdb : begin id in bulk
    // bIdg : begin id in grid
    // eIdb : end id in bulk
    bool    activity;
    OCP_USI eIdb;
    for (OCP_USI bIdb = 0; bIdb < numBulk; bIdb++) {
        OCP_USI bIdg = myGrid.activeMap_B2G[bIdb];
        USI     k    = bIdg / nxny;
        USI     j    = (bIdg - k * nxny) / nx;
        USI     i    = bIdg - k * nxny - j * nx;

        USI count = 0;
        // up
        if (k > 0) {
            activity = myGrid.activeMap_G2B[bIdg - nxny].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg - nxny].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // back
        if (j > 0) {
            activity = myGrid.activeMap_G2B[bIdg - nx].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg - nx].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // left
        if (i > 0) {
            activity = myGrid.activeMap_G2B[bIdg - 1].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg - 1].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        numConn += count;
        // self
        neighbor[bIdb].push_back(bIdb);
        selfPtr[bIdb] = count;
        count++;

        // right
        if (i < nx - 1) {
            activity = myGrid.activeMap_G2B[bIdg + 1].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg + 1].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // front
        if (j < ny - 1) {
            activity = myGrid.activeMap_G2B[bIdg + nx].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg + nx].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        // down
        if (k < nz - 1) {
            activity = myGrid.activeMap_G2B[bIdg + nxny].GetAct();
            if (activity) {
                eIdb = myGrid.activeMap_G2B[bIdg + nxny].GetId();
                neighbor[bIdb].push_back(eIdb);
                count++;
            }
        }

        neighborNum[bIdb] = count;
    }
}

void BulkConn::CalIteratorConn()
{
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

void BulkConn::CalArea(const Grid& myGrid, const Bulk& myBulk)
{
    // calculate efficient area of interface of bulk to bulk
    // using iteratorConn
    area.reserve(numConn);

    for (OCP_USI n = 0; n < numConn; n++) {
        OCP_USI bIdb = iteratorConn[n].BId;
        OCP_USI eIdb = iteratorConn[n].EId;
        // area.push_back(1);
        area.push_back(CalAkd(myGrid, myBulk, bIdb, eIdb));
    }

    assert(area.size() == numConn);
}

OCP_DBL BulkConn::CalAkd(const Grid& myGrid, const Bulk& myBulk,
                              const OCP_USI& bIdb, const OCP_USI& eIdb) const
{
    OCP_USI bIdg = myGrid.activeMap_B2G[bIdb];
    OCP_USI eIdg = myGrid.activeMap_B2G[eIdb];
    OCP_INT diff = eIdg - bIdg;
    if (diff < 0) diff *= -1;

    if (diff == 1) {
        // x - direction
        OCP_DBL T1 = myBulk.rockKx[bIdb] * myBulk.ntg[bIdb] * myBulk.dy[bIdb] *
                     myBulk.dz[bIdb] / myBulk.dx[bIdb];
        OCP_DBL T2 = myBulk.rockKx[eIdb] * myBulk.ntg[eIdb] * myBulk.dy[eIdb] *
                     myBulk.dz[eIdb] / myBulk.dx[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else if (diff - myGrid.nx == 0) {
        // y - direction
        OCP_DBL T1 = myBulk.rockKy[bIdb] * myBulk.ntg[bIdb] * myBulk.dz[bIdb] *
                     myBulk.dx[bIdb] / myBulk.dy[bIdb];
        OCP_DBL T2 = myBulk.rockKy[eIdb] * myBulk.ntg[eIdb] * myBulk.dz[eIdb] *
                     myBulk.dx[eIdb] / myBulk.dy[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else if (diff - myGrid.nx * myGrid.ny == 0) {
        // z - direction  ----  no ntg
        OCP_DBL T1 =
            myBulk.rockKz[bIdb] * myBulk.dx[bIdb] * myBulk.dy[bIdb] / myBulk.dz[bIdb];
        OCP_DBL T2 =
            myBulk.rockKz[eIdb] * myBulk.dx[eIdb] * myBulk.dy[eIdb] / myBulk.dz[eIdb];
        return (2 / (1 / T1 + 1 / T2));
    } else {
        cout << "bIdg = " << bIdg << "eIdg = " << eIdg << endl;
        OCP_ABORT("Wrong bIdg and eIdg in function");
    }
}

// Connection function, no matter what grid is

OCP_DBL BulkConn::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    USI     np   = myBulk.numPhase;
    OCP_DBL cfl  = 0;
    OCP_DBL temp = 0;
    for (OCP_USI c = 0; c < numConn; c++) {

        for (USI j = 0; j < np; j++) {
            OCP_USI uId = upblock[c * np + j];

            if (myBulk.phaseExist[uId * np + j]) {  // uId -> uId * np + j  fix bugs.
                temp = fabs(upblock_Velocity[c * np + j]) * dt;
                temp /= myBulk.vj[uId * np + j];
                if (cfl < temp) cfl = temp;
            }
        }
    }
    return cfl;
}

void BulkConn::CalCFL01(const Bulk& myBulk, const OCP_DBL& dt) const
{
    USI     np = myBulk.numPhase;
    for (OCP_USI c = 0; c < numConn; c++) {

        for (USI j = 0; j < np; j++) {
            OCP_USI uId = upblock[c * np + j];

            if (myBulk.phaseExist[uId * np + j]) {
                myBulk.cfl[uId * np + j] += fabs(upblock_Velocity[c * np + j]) * dt;
            }
        }
    }
}

void BulkConn::CalFlux(const Bulk& myBulk)
{
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
                upblock[c * np + j] = bId;
                upblock_Trans[c * np + j] = 0;
                upblock_Velocity[c * np + j] = 0;
                //if (!isfinite(upblock_Trans[c * np + j])) {
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
                upblock_Trans[c * np + j] = trans;
                upblock_Velocity[c * np + j] = trans * dP;
            } else {
                upblock_Trans[c * np + j] = 0;
                upblock_Velocity[c * np + j] = 0;
            }
        }
    }
}

void BulkConn::MassConserve(Bulk& myBulk, const OCP_DBL& dt) const
{
    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;

    for (OCP_USI c = 0; c < numConn; c++) {
        OCP_USI bId = iteratorConn[c].BId;
        OCP_USI eId = iteratorConn[c].EId;

        for (USI j = 0; j < np; j++) {
            OCP_USI uId = upblock[c * np + j];
            OCP_USI uId_np_j = uId * np + j;

            if (!myBulk.phaseExist[uId_np_j]) continue;

            OCP_DBL phaseVelocity = upblock_Velocity[c * np + j];
            for (USI i = 0; i < nc; i++) {
                OCP_DBL dNi = dt * phaseVelocity * myBulk.xi[uId_np_j] *
                              myBulk.cij[uId_np_j * nc + i];
                myBulk.Ni[eId * nc + i] += dNi;
                myBulk.Ni[bId * nc + i] -= dNi;
            }
        }
    }
}

void BulkConn::AllocateMat(LinearSolver& MySolver) const
{
    for (OCP_USI n = 0; n < numBulk; n++) {
        MySolver.rowCapacity[n] += neighborNum[n];
    }
}

void BulkConn::InitAssembleMat(LinearSolver& mySolver) const
{
    mySolver.dim = numBulk;
    for (OCP_USI n = 0; n < numBulk; n++) {
        mySolver.colId[n].assign(neighbor[n].begin(), neighbor[n].end());
        mySolver.diagPtr[n] = selfPtr[n];
    }
}

void BulkConn::AssembleMat_IMPEC(LinearSolver& mySolver, const Bulk& myBulk,
                                      const OCP_DBL& dt) const
{
    // accumulate term
    OCP_DBL Vp0, Vp, vf, vfp, P;
    OCP_DBL cr = myBulk.rockC1;
    for (OCP_USI n = 0; n < numBulk; n++) {
        Vp0 = myBulk.rockVpInit[n];
        Vp  = myBulk.rockVp[n];
        vfp = myBulk.vfp[n];
        P   = myBulk.P[n];
        vf  = myBulk.vf[n];

        OCP_DBL temp        = cr * Vp0 - vfp;
        mySolver.diagVal[n] = temp;
        mySolver.b[n]       = temp * P + dt * (vf - Vp);
    }

    // check
    // ofstream outb("testb.dat");
    // if (!outb.is_open())
    //	cout << "Can not open " << "testb.dat" << endl;
    // outb << mySolver.dim << endl;
    // for (int i = 0; i < mySolver.dim; i++)
    //	outb << mySolver.b[i] << endl;
    // outb.close();

    // flux term
    OCP_USI bId, eId, uId;
    USI     np = myBulk.numPhase;
    USI     nc = myBulk.numCom;
    OCP_DBL valupi, valdowni;
    OCP_DBL valup, rhsup, valdown, rhsdown;
    // OCP_USI    lastbId = -1;
    // Becareful when first bulk has no neighbors!
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
                    myBulk.vfi[bId * nc + i] * myBulk.cij[uId * np * nc + j * nc + i];
                valdowni +=
                    myBulk.vfi[eId * nc + i] * myBulk.cij[uId * np * nc + j * nc + i];
            }
            OCP_DBL dD   = myBulk.depth[bId] - myBulk.depth[eId];
            OCP_DBL dPc  = myBulk.Pc[bId * np + j] - myBulk.Pc[eId * np + j];
            OCP_DBL temp = myBulk.xi[uId * np + j] * upblock_Trans[c * np + j] * dt;
            valup += temp * valupi;

            //if (!isfinite(valup)) {
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
            assert(mySolver.val[bId].size() == diagptr);
            mySolver.val[bId].push_back(mySolver.diagVal[bId]);
            lastbId = bId;
        }

        mySolver.val[bId][diagptr] += valup;
        mySolver.val[bId].push_back(-valup);
        mySolver.val[eId].push_back(-valdown);
        mySolver.diagVal[eId] += valdown;
        mySolver.b[bId] += rhsup;
        mySolver.b[eId] += rhsdown;
    }

    // Add the rest of diag value
    // important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (mySolver.val[n].size() == selfPtr[n])
            mySolver.val[n].push_back(mySolver.diagVal[n]);
    }
}

void BulkConn::AssembleMat_FIM(LinearSolver& mySolver, const Bulk& myBulk, const OCP_DBL& dt) const
{
    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;
    USI ncol = nc + 1;
    USI ncol2 = np*nc + np;
    USI bsize = ncol * ncol;
    USI bsize2 = ncol * ncol2;

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
            mySolver.diagVal[n * bsize + i] = bmat[i];
        }
    }
    // flux term
    OCP_DBL Akd;
    OCP_DBL transJ,transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);

    OCP_USI bId, eId, uId;
    OCP_USI uId_np_j;
    OCP_DBL kr, mu, xi, cij, rhoP, xiP, muP, rhoC, xiC, muC;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        dFdXpB.assign(bsize, 0);
        dFdXpE.assign(bsize, 0);
        dFdXsB.assign(bsize2, 0);
        dFdXsE.assign(bsize2, 0);
        Akd = CONV1 * CONV2 * area[c];
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];
            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j])    continue;
            dP = myBulk.Pj[bId * np + j] - myBulk.Pj[eId * np + j]
                - myBulk.rho[bId * np + j] * dGamma;
            xi = myBulk.xi[uId_np_j];
            kr = myBulk.kr[uId_np_j];
            mu = myBulk.mu[uId_np_j];   
            muP = myBulk.muP[uId_np_j];
            xiP = myBulk.xiP[uId_np_j];
            rhoP = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;
            
            for (USI i = 0; i < nc; i++) {
                cij = myBulk.cij[uId_np_j * nc + i];

                transIJ = cij * xi * transJ;
                
                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;
                tmp = transIJ * (-rhoP * dGamma);
                tmp += cij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (bId == uId) {
                    dFdXpB[(i + 1) * ncol] += tmp;
                }
                else {
                    dFdXpE[(i + 1) * ncol] += tmp;
                }

                // Saturation -- Second var
                for (USI k = 0; k < np; k++) {
                    dFdXsB[(i + 1) * ncol2 + k] += transIJ * myBulk.dPcj_dS[bId * np * np + j * np + k];
                    dFdXsE[(i + 1) * ncol2 + k] += transIJ * myBulk.dPcj_dS[eId * np * np + j * np + k];
                    tmp = Akd * cij * xi / mu * myBulk.dKr_dS[uId * np * np + j * np + k] * dP;
                    if (bId == uId) {
                        dFdXsB[(i + 1) * ncol2 + k] += tmp;
                    }
                    else {
                        dFdXsE[(i + 1) * ncol2 + k] += tmp;
                    }
                }
                // Cij -- Second var
                for (USI k = 0; k < nc; k++) {
                    rhoC = myBulk.rhoC[uId_np_j * nc + k];
                    xiC = myBulk.xiC[uId_np_j * nc + k];
                    muC = myBulk.muC[uId_np_j * nc + k];
                    tmp = -transIJ * rhoC * dGamma;
                    tmp += cij * transJ * xiC * dP;
                    tmp += -transIJ * muC / mu * dP;
                    if (k == i) {
                        tmp += xi * transJ * dP;
                    }
                    if (bId == uId) {
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                    else {
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                }
            }
        }

        USI diagptr = selfPtr[bId];
        
        if (bId != lastbId) {
            // new bulk
            assert(mySolver.val[bId].size() == diagptr * bsize);  
            OCP_USI id = bId * bsize;
            mySolver.val[bId].insert(mySolver.val[bId].end(),
                mySolver.diagVal.data() + id, mySolver.diagVal.data() + id + bsize);
            
            lastbId = bId;
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * bsize2], 1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin
        // Add
        for (USI i = 0; i < bsize; i++) {
            mySolver.val[bId][diagptr * bsize + i] += bmat[i];
        }
        // End
        // Insert
        Dscalar(bsize, -1, bmat.data());
        mySolver.val[eId].insert(mySolver.val[eId].end(), bmat.begin(), bmat.end());

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * bsize2], 1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin
        // Insert
        mySolver.val[bId].insert(mySolver.val[bId].end(), bmat.begin(), bmat.end());
        // Add
        Dscalar(bsize, -1, bmat.data());
        for (USI i = 0; i < bsize; i++) {
            mySolver.diagVal[eId * bsize + i] += bmat[i];
        }
    }
    // Add the rest of diag value
    // important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (mySolver.val[n].size() == selfPtr[n] * bsize)
            mySolver.val[n].insert(mySolver.val[n].end(),
                mySolver.diagVal.data() + n * bsize, mySolver.diagVal.data() + n * bsize + bsize);
    }
}

void BulkConn::CalResFIM(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;
    USI len = nc + 1;
    OCP_USI bId, eId, uId, bIdb;
    // Accumalation Term
    for (OCP_USI n = 0; n < numBulk; n++) {
        bId = n * len;
        bIdb = n * nc;
        res[bId] = myBulk.rockVp[n] - myBulk.vf[n];
        for (USI i = 0; i < nc; i++) {
            res[bId + 1 + i] = myBulk.Ni[bIdb + i] - myBulk.lNi[bIdb + i];
        }
    }

    //cout << "----------------------------------------" << endl;
    //for (OCP_USI i = 0; i < numBulk * (nc + 1); i++) {
    //    cout << res[i] << endl;
    //}

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
            bool exend = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend = myBulk.Pj[eId_np_j];
                rho = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            }
            else if (exbegin && (!exend)) {
                Pbegin = myBulk.Pj[bId_np_j];
                Pend = myBulk.P[eId];
                rho = myBulk.rho[bId_np_j];
            }
            else if ((!exbegin) && (exend)) {
                Pbegin = myBulk.P[bId];
                Pend = myBulk.Pj[eId_np_j];
                rho = myBulk.rho[eId_np_j];
            }
            else {
                upblock[c * np + j] = bId;
                upblock_Rho[c * np + j] = rho;
                continue;
            }

            uId = bId;
            bool    exup = exbegin;
            dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
                exup = exend;
            }
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j] = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;

            tmp = dt * CONV1 * CONV2 * area[c] * myBulk.xi[uId_np_j] * myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
            for (USI i = 0; i < nc; i++) {
                dNi = tmp * myBulk.cij[uId_np_j * nc + i];
                /*res[bId * len + 1 + i] -= dNi;
                res[eId * len + 1 + i] += dNi;*/
                res[bId * len + 1 + i] += dNi;
                res[eId * len + 1 + i] -= dNi;
            }
        }
    }
}

void BulkConn::UpdateLastStep()
{
    lastUpblock = upblock;
    lastUpblock_Rho = upblock_Rho;
    lastUpblock_Trans = upblock_Trans;
    lastUpblock_Velocity = upblock_Velocity;

}

void BulkConn::Reset()
{
    upblock = lastUpblock;
    upblock_Rho = lastUpblock_Rho;
    upblock_Trans = lastUpblock_Trans;
    upblock_Velocity = lastUpblock_Velocity;
}

void BulkConn::CheckDiff() const
{
    // upblock
    OCP_DBL tmp;
    cout << setprecision(18);
    for (OCP_USI c = 0; c < numConn; c++) {
        tmp = fabs(upblock[c] - lastUpblock[c]);
        if (tmp != 0.0) {
            cout << "Difference in upblock\t" << tmp << "\n";
        }
        tmp = fabs(upblock_Rho[c] - lastUpblock_Rho[c]);
        if ( tmp != 0.0) {
            cout << "Difference in upblock_Rho\t" << tmp << "\n";
        }
        tmp = fabs(upblock_Trans[c] - lastUpblock_Trans[c]);
        if ( tmp != 0.0) {
            cout << "Difference in upblock_Trans\t" << tmp << "\n";
        }
        tmp = fabs(upblock_Velocity[c] - lastUpblock_Velocity[c]);
        if ( tmp != 0.0) {
            cout << "Difference in upblock_Velocity\t" << tmp << "\n";
        }
    }
}


void BulkConn::GetConnectionInfo() const
{
    for (OCP_USI i = 0; i < numBulk; i++) {
        cout << "(" << i << ")"
             << "\t";

        for (auto v : neighbor[i]) {
            cout << v << "\t";
        }
        cout << "[" << selfPtr[i] << "]";
        cout << "\t" << neighborNum[i];
        cout << "\n";
    }

    for (OCP_USI i = 0; i < numConn; i++) {
        cout << iteratorConn[i].BId << "\t" << iteratorConn[i].EId << "\n";
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