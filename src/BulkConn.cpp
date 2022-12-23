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
void BulkConn::Setup(const Grid& myGrid)
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
        const GB_Pair& GBtmp = myGrid.map_All2Act[n];

        if (GBtmp.IsAct()) {
            bIdb = GBtmp.GetId();

            // Get rid of inactive neighbor
            tmp1 = myGrid.gNeighbor[n];
            len  = tmp1.size();
            tmp2.clear();
            for (USI i = 0; i < len; i++) {
                const GB_Pair& GBtmp2 = myGrid.map_All2Act[tmp1[i].id];

                if (GBtmp2.IsAct()) {
                    tmp1[i].id = GBtmp2.GetId();
                    tmp2.push_back(tmp1[i]);
                }
            }
            // Add Self
            tmp2.push_back(GPair(bIdb, 0, 0.0, 0.0));
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
                iteratorConn.push_back(BulkPair(bIdb, tmp2[j].id, tmp2[j].direction,
                                                tmp2[j].areaB, tmp2[j].areaE));
            }
            neighborNum[bIdb] = len;
        }
    }

    numConn = iteratorConn.size();

    // PrintConnectionInfoCoor(myGrid);
}


void BulkConn::CalAkd(const Bulk& myBulk)
{
    OCP_USI bId, eId;
    OCP_DBL areaB, areaE;
    OCP_DBL T1, T2;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        areaB = iteratorConn[c].areaB;
        areaE = iteratorConn[c].areaE;
        switch (iteratorConn[c].direction)
        {
        case 1:
            T1 = myBulk.ntg[bId] * myBulk.rockKx[bId] *  areaB;
            T2 = myBulk.ntg[eId] * myBulk.rockKx[eId] *  areaE;
            break;
        case 2:
            T1 = myBulk.ntg[bId] * myBulk.rockKy[bId] *  areaB;
            T2 = myBulk.ntg[eId] * myBulk.rockKy[eId] *  areaE;
            break;
        case 3:
            T1 = myBulk.rockKz[bId] * areaB;
            T2 = myBulk.rockKz[eId] * areaE;
            break;
        default:
            OCP_ABORT("Wrong Direction!");
        }
        iteratorConn[c].area = 1 / (1 / T1 + 1 / T2);
    }
}


void BulkConn::SetupWellBulk_K(Bulk& myBulk) const
{
    // For K = 1 now, defaulted

    USI len = myBulk.wellBulkId.size();
    for (USI n = 0; n < len; n++) {
        for (auto& v : neighbor[n]) {
            USI      clen = myBulk.wellBulkId.size();
            OCP_BOOL flag = OCP_FALSE;
            for (USI i = 0; i < clen; i++) {
                if (v == myBulk.wellBulkId[n]) {
                    flag = OCP_TRUE;
                    break;
                }
            }
            if (!flag) {
                myBulk.wellBulkId.push_back(v);
            }
        }
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


void BulkConn::PrintConnectionInfo(const Grid& myGrid) const
{
    for (OCP_USI i = 0; i < numBulk; i++) {
        cout << "(" << myGrid.map_Act2All[i] << ")" << "\t";

        for (auto& v : neighbor[i]) {
            cout << myGrid.map_Act2All[v] << "\t";
        }
        cout << "[" << selfPtr[i] << "]";
        cout << "\t" << neighborNum[i];
        cout << "\n";
    }

    for (OCP_USI i = 0; i < numConn; i++) {
        cout << myGrid.map_Act2All[iteratorConn[i].bId] << "\t"
             << myGrid.map_Act2All[iteratorConn[i].eId] << "\n";
    }
}

void BulkConn::PrintConnectionInfoCoor(const Grid& myGrid) const
{
    OCP_USI bIdg, eIdg;
    OCP_USI bIdb, eIdb;
    USI     I, J, K;
    USI     sp = myGrid.GetNumDigitIJK();
    cout << "BulkConn : " << numConn << endl;
    for (OCP_USI c = 0; c < numConn; c++) {
        bIdb = iteratorConn[c].bId;
        eIdb = iteratorConn[c].eId;
        bIdg = myGrid.map_Act2All[bIdb];
        eIdg = myGrid.map_Act2All[eIdb];
        myGrid.GetIJKGrid(I, J, K, bIdg);
        cout << GetIJKformat(I, J, K, sp) << "   ";
        cout << setw(6) << bIdg;
        cout << "    ";
        cout << setw(6) << bIdb;
        cout << "    ";
        myGrid.GetIJKGrid(I, J, K, eIdg);
        cout << GetIJKformat(I, J, K, sp) << "   ";
        cout << setw(6) << eIdg;
        cout << "    ";
        cout << setw(6) << eIdb;
        cout << setw(20) << setprecision(8) << fixed << iteratorConn[c].area * CONV2;

        cout << endl;
    }
}

void BulkConn::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;
    fill(myBulk.cfl.begin(), myBulk.cfl.end(), 0);

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

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void BulkConn::AllocateIMPEC_IsoT(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    upblock_Trans.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);

    lupblock.resize(numConn * np);
    lupblock_Rho.resize(numConn * np);
    lupblock_Trans.resize(numConn * np);
    lupblock_Velocity.resize(numConn * np);
}

void BulkConn::AssembleMatIMPEC(LinearSystem&  myLS,
                                const Bulk&    myBulk,
                                const OCP_DBL& dt) const
{
    OCP_FUNCNAME;
    myLS.AddDim(numBulk);

    // accumulate term
    OCP_DBL Vpp, Vp, vf, vfP, P;
    for (OCP_USI n = 0; n < numBulk; n++) {
        vf  = myBulk.vf[n];
        vfP = myBulk.vfP[n];
        P   = myBulk.lP[n];
        Vpp = myBulk.rockVntg[n] * myBulk.poroP[n];
        Vp  = myBulk.rockVp[n];

        myLS.NewDiag(n, Vpp - vfP);
        myLS.AddRhs(n, (Vpp - vfP) * P + dt * (vf - Vp));
    }

    // flux term
    OCP_USI bId, eId;
    OCP_USI uId_np_j;
    OCP_DBL valupi, valdowni;
    OCP_DBL valup, rhsup, valdown, rhsdown;
    OCP_DBL dD, tmp;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;

    // Be careful when first bulk has no neighbors!
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        dD  = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        valup   = 0;
        rhsup   = 0;
        valdown = 0;
        rhsdown = 0;

        for (USI j = 0; j < np; j++) {
            uId_np_j = upblock[c * np + j] * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;

            valupi   = 0;
            valdowni = 0;

            for (USI i = 0; i < nc; i++) {
                valupi   += myBulk.vfi[bId * nc + i] * myBulk.xij[uId_np_j * nc + i];
                valdowni += myBulk.vfi[eId * nc + i] * myBulk.xij[uId_np_j * nc + i];
            }
            
            tmp     =  myBulk.xi[uId_np_j] * upblock_Trans[c * np + j] * dt;
            valup   += tmp * valupi;
            valdown += tmp * valdowni;
            tmp     *= upblock_Rho[c * np + j] * dD 
                    - (myBulk.Pc[bId * np + j] - myBulk.Pc[eId * np + j]);
            rhsup   += tmp * valupi;
            rhsdown -= tmp * valdowni;
        }
        myLS.AddDiag(bId, valup);
        myLS.AddDiag(eId, valdown);
        myLS.NewOffDiag(bId, eId, -valup);
        myLS.NewOffDiag(eId, bId, -valdown);
        myLS.AddRhs(bId, rhsup);
        myLS.AddRhs(eId, rhsdown);
    }
}


void BulkConn::CalFluxIMPEC(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    // calculate a step flux using iteratorConn
    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j;
    OCP_BOOL exbegin, exend, exup;
    OCP_DBL  rho, dP, Akd;

    const USI np = myBulk.numPhase;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            exbegin = myBulk.phaseExist[bId_np_j];
            exend   = myBulk.phaseExist[eId_np_j];
            
            if ((exbegin) && (exend)) {
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]          = bId;
                upblock_Trans[c * np + j]    = 0;
                upblock_Velocity[c * np + j] = 0;
                continue;
            }

            dP   = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[bId]) 
                 - (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[eId]);           
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            }
            else {
                uId  = bId;
                exup = exbegin;
            }

            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;

            if (exup) {
                upblock_Trans[c * np + j]    = Akd * myBulk.kr[uId * np + j] / myBulk.mu[uId * np + j];
                upblock_Velocity[c * np + j] = upblock_Trans[c * np + j] * dP;
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

    OCP_USI bId, eId, uId;
    OCP_USI uId_np_j;
    OCP_DBL phaseVelocity, dNi;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;

        for (USI j = 0; j < np; j++) {
            uId      = upblock[c * np + j];
            uId_np_j = uId * np + j;

            if (!myBulk.phaseExist[uId_np_j]) continue;

            phaseVelocity = upblock_Velocity[c * np + j];
            for (USI i = 0; i < nc; i++) {
                dNi = dt * phaseVelocity * myBulk.xi[uId_np_j] *
                              myBulk.xij[uId_np_j * nc + i];
                myBulk.Ni[eId * nc + i] += dNi;
                myBulk.Ni[bId * nc + i] -= dNi;
            }
        }
    }
}

void BulkConn::ResetIMPEC()
{
    OCP_FUNCNAME;

    upblock             = lupblock;
    upblock_Rho         = lupblock_Rho;
    upblock_Trans       = lupblock_Trans;
    upblock_Velocity    = lupblock_Velocity;
}

void BulkConn::UpdateLastStepIMPEC()
{
    OCP_FUNCNAME;

    lupblock            = upblock;
    lupblock_Rho        = upblock_Rho;
    lupblock_Trans      = upblock_Trans;
    lupblock_Velocity   = upblock_Velocity;
}


/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void BulkConn::AllocateFIM_IsoT(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
}

void BulkConn::AssembleMat_FIM(LinearSystem&  myLS,
                               const Bulk&    myBulk,
                               const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    myLS.AddDim(numBulk);

    const USI np     = myBulk.numPhase;
    const USI nc     = myBulk.numCom;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockVntg[n] * myBulk.poroP[n] - myBulk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }
        myLS.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);

    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL phaseExistBj, phaseExistEj;
    OCP_DBL  kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL  dP, dGamma;
    OCP_DBL  tmp;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        for (USI j = 0; j < np; j++) {
            uId      = upblock[c * np + j];
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;

            if (!myBulk.phaseExist[uId_np_j]) continue;

            phaseExistBj = myBulk.phaseExist[bId_np_j];
            phaseExistEj = myBulk.phaseExist[eId_np_j];

            dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                 upblock_Rho[c * np + j] * dGamma;
            xi     = myBulk.xi[uId_np_j];
            kr     = myBulk.kr[uId_np_j];
            mu     = myBulk.mu[uId_np_j];
            muP    = myBulk.muP[uId_np_j];
            xiP    = myBulk.xiP[uId_np_j];
            rhoP   = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistEj) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistBj) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }
                // Second var
                if (bId == uId) {
                    // Saturation
                    for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[bId_np_j * np + k];
                        tmp =
                            Akd * xij * xi / mu * myBulk.dKr_dS[uId_np_j * np + k] * dP;
                        dFdXsB[(i + 1) * ncol2 + k] += tmp;
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * myBulk.dPcj_dS[eId_np_j * np + k];
                    }
                    // Cij
                    if (!phaseExistEj) {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] +=
                                -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }
                } else {
                    // Saturation
                    for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[bId_np_j * np + k];
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * myBulk.dPcj_dS[eId_np_j * np + k];
                        tmp =
                            Akd * xij * xi / mu * myBulk.dKr_dS[uId_np_j * np + k] * dP;
                        dFdXsE[(i + 1) * ncol2 + k] += tmp;
                    }
                    // Cij
                    if (!phaseExistBj) {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                                -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }
                }
            }
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        myLS.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        myLS.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        myLS.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        myLS.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}


/// rho = (rho1 + rho2)/2
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
    OCP_DBL rho, dP;
    OCP_DBL tmp, dNi;
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]     = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }

            uId = bId;
            dP  = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                 (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            // dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);
            // dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho;
            // cout << setprecision(6) << scientific << dP << "   "
            //     << myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho << "   "
            //     << dP - (myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho)
            //     << endl;
            if (dP < 0) {
                uId = eId;
            }
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            tmp = dt * Akd * myBulk.xi[uId_np_j] * myBulk.kr[uId_np_j] /
                  myBulk.mu[uId_np_j] * dP;

            for (USI i = 0; i < nc; i++) {
                dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                res[bId * len + 1 + i] += dNi;
                res[eId * len + 1 + i] -= dNi;
            }
        }
    }
}

/// rho = (S1*rho1 + S2*rho2)/(S1+S2)
void BulkConn::CalFluxFIMS(const Grid& myGrid, const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    OCP_USI   bId, eId, uId;
    OCP_USI   bId_np_j, eId_np_j;
    OCP_DBL   Pbegin, Pend, rho;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho = (myBulk.rho[bId_np_j] * myBulk.S[bId_np_j] +
                       myBulk.rho[eId_np_j] * myBulk.S[eId_np_j]) /
                      (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
            } else if (exbegin && (!exend)) {
                rho = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]     = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }
            Pbegin     = myBulk.Pj[bId_np_j];
            Pend       = myBulk.Pj[eId_np_j];
            uId        = bId;
            OCP_DBL dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                         (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock[c * np + j]     = uId;
            upblock_Rho[c * np + j] = rho;
        }
    }
}

void BulkConn::CalResFIMS(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
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
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {
                rho = (myBulk.rho[bId_np_j] * myBulk.S[bId_np_j] +
                       myBulk.rho[eId_np_j] * myBulk.S[eId_np_j]) /
                      (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
            } else if (exbegin && (!exend)) {
                rho = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]     = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }
            Pbegin = myBulk.Pj[bId_np_j];
            Pend   = myBulk.Pj[eId_np_j];
            uId    = bId;
            dP     = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                 (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            tmp = dt * Akd * myBulk.xi[uId_np_j] * myBulk.kr[uId_np_j] /
                  myBulk.mu[uId_np_j] * dP;

            for (USI i = 0; i < nc; i++) {
                dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                res[bId * len + 1 + i] += dNi;
                res[eId * len + 1 + i] -= dNi;
            }
        }
    }
}


/////////////////////////////////////////////////////////////////////
// FIM(new)
/////////////////////////////////////////////////////////////////////

void BulkConn::AssembleMat_FIM_new(LinearSystem&  myLS,
                                   const Bulk&    myBulk,
                                   const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    myLS.AddDim(numBulk);

    const USI np      = myBulk.numPhase;
    const USI nc      = myBulk.numCom;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockVntg[n] * myBulk.poroP[n] - myBulk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }
        myLS.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        USI jxB = 0;
        USI jxE = 0;
        ncolB   = 0;
        ncolE   = 0;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            phasedS_B[j]   = myBulk.pSderExist[bId * np + j];
            phasedS_E[j]   = myBulk.pSderExist[eId * np + j];
            if (phasedS_B[j]) jxB++;
            if (phasedS_E[j]) jxE++;
            pVnumComB[j] = myBulk.pVnumCom[bId * np + j];
            pVnumComE[j] = myBulk.pVnumCom[eId * np + j];
            ncolB += pVnumComB[j];
            ncolE += pVnumComE[j];
        }
        ncolB += jxB;
        ncolE += jxE;

        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                 upblock_Rho[c * np + j] * dGamma;
            xi     = myBulk.xi[uId_np_j];
            kr     = myBulk.kr[uId_np_j];
            mu     = myBulk.mu[uId_np_j];
            muP    = myBulk.muP[uId_np_j];
            xiP    = myBulk.xiP[uId_np_j];
            rhoP   = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    } else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                }
            }
            jxB += pVnumComB[j];
            jxE += pVnumComE[j];
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP],
                1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        myLS.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        myLS.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP],
                1, bmat.data());

        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        myLS.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        myLS.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}


void BulkConn::AssembleMat_FIM_newS(LinearSystem&  myLS,
                                    const Bulk&    myBulk,
                                    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    myLS.AddDim(numBulk);

    const USI np      = myBulk.numPhase;
    const USI nc      = myBulk.numCom;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockVntg[n] * myBulk.poroP[n] - myBulk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }

        myLS.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<USI>      pEnumComB(np, 0);
    vector<USI>      pEnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;
    OCP_DBL wghtb, wghte;


    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        const USI npB = myBulk.phaseNum[bId];
        ncolB         = npB;
        const USI npE = myBulk.phaseNum[eId];
        ncolE         = npE;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            pEnumComB[j]   = myBulk.pVnumCom[bId * np + j];
            pEnumComE[j]   = myBulk.pVnumCom[eId * np + j];
            ncolB += pEnumComB[j];
            ncolE += pEnumComE[j];
        }

        USI jxB = npB;
        USI jxE = npE;
        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pEnumComB[j];
                jxE += pEnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                 upblock_Rho[c * np + j] * dGamma;
            xi     = myBulk.xi[uId_np_j];
            kr     = myBulk.kr[uId_np_j];
            mu     = myBulk.mu[uId_np_j];
            muP    = myBulk.muP[uId_np_j];
            xiP    = myBulk.xiP[uId_np_j];
            rhoP   = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * dGamma *
                        (-myBulk.rhoP[bId_np_j] * myBulk.S[bId_np_j]) /
                        (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * dGamma *
                        (-myBulk.rhoP[eId_np_j] * myBulk.S[eId_np_j]) /
                        (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0;
                        wghte = 0;
                        if (j1 == j && phaseExistE[j]) {
                            tmp = -dGamma / ((myBulk.S[bId_np_j] + myBulk.S[eId_np_j]) *
                                             (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]));
                            wghtb = tmp * myBulk.rho[bId_np_j] * myBulk.S[eId_np_j];
                            wghte = tmp * myBulk.rho[eId_np_j] * myBulk.S[bId_np_j];
                        }

                        if (phaseExistB[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * (myBulk.dPcj_dS[bId_np_j * np + j1] + wghtb);
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phaseExistE[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * (myBulk.dPcj_dS[eId_np_j * np + j1] + wghte);
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pEnumComB[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    } else {
                        wghtb = myBulk.S[bId_np_j] /
                                (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        wghte = myBulk.S[eId_np_j] /
                                (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] * wghtb;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * dGamma * myBulk.rhox[eId_np_j * nc + k] *
                                wghte;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0;
                        wghte = 0;
                        if (j1 == j && phaseExistB[j]) {
                            tmp = -dGamma / ((myBulk.S[bId_np_j] + myBulk.S[eId_np_j]) *
                                             (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]));
                            wghtb = tmp * myBulk.rho[bId_np_j] * myBulk.S[eId_np_j];
                            wghte = tmp * myBulk.rho[eId_np_j] * myBulk.S[bId_np_j];
                        }

                        if (phaseExistB[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * (myBulk.dPcj_dS[bId_np_j * np + j1] + wghtb);
                            j1SB++;
                        }
                        if (phaseExistE[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * (myBulk.dPcj_dS[eId_np_j * np + j1] + wghte);
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pEnumComE[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    } else {
                        wghtb = myBulk.S[bId_np_j] /
                                (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        wghte = myBulk.S[eId_np_j] /
                                (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] * wghte;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * dGamma * myBulk.rhox[bId_np_j * nc + k] *
                                wghtb;
                        }
                        // WARNING !!!
                        if (i < pEnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                }
            }
            jxB += pEnumComB[j];
            jxE += pEnumComE[j];
        }

        // Assemble
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP],
                1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        myLS.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        myLS.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP],
                1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        myLS.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        myLS.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}

void BulkConn::AssembleMat_FIM_new_n(LinearSystem&  myLS,
                                     const Bulk&    myBulk,
                                     const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    myLS.AddDim(numBulk);

    const USI np      = myBulk.numPhase;
    const USI nc      = myBulk.numCom;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> tmpb(ncol, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockVntg[n] * myBulk.poroP[n] - myBulk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }

        myLS.NewDiag(n, bmat);
        myLS.AddRhs(n * ncol, -myBulk.resPc[n]);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    USI              ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        const USI npB = myBulk.phaseNum[bId];
        const USI npE = myBulk.phaseNum[eId];

        USI jxB = 0;
        USI jxE = 0;
        ncolB   = 0;
        ncolE   = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            phasedS_B[j]   = myBulk.pSderExist[bId * np + j];
            phasedS_E[j]   = myBulk.pSderExist[eId * np + j];
            if (phasedS_B[j]) jxB++;
            if (phasedS_E[j]) jxE++;
            pVnumComB[j] = myBulk.pVnumCom[bId * np + j];
            pVnumComE[j] = myBulk.pVnumCom[eId * np + j];
            ncolB += pVnumComB[j];
            ncolE += pVnumComE[j];
        }
        ncolB += jxB;
        ncolE += jxE;

        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
                continue;
            }

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP       = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                       upblock_Rho[c * np + j] * dGamma;
            xi     = myBulk.xi[uId_np_j];
            kr     = myBulk.kr[uId_np_j];
            mu     = myBulk.mu[uId_np_j];
            muP    = myBulk.muP[uId_np_j];
            xiP    = myBulk.xiP[uId_np_j];
            rhoP   = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;

            for (USI i = 0; i < nc; i++) {
                xij     = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                } else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                } else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0;
                USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistE[j]) {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] +=
                                xi * transJ * dP / myBulk.nj[uId_np_j];
                    } else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] +=
                                xi * transJ * dP / myBulk.nj[uId_np_j];
                    }
                } else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + j1SB] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                            j1SB++;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + j1SE] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                  myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                            j1SE++;
                        }
                    }
                    // Cij
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] +=
                                xi * transJ * dP / myBulk.nj[uId_np_j];
                    } else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix  = myBulk.xix[uId_np_j * nc + k];
                            mux  = myBulk.mux[uId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] +=
                                xi * transJ * dP / myBulk.nj[uId_np_j];
                    }
                }
            }
            jxB += pVnumComB[j];
            jxE += pVnumComE[j];
        }

        // Assemble rhs
        // Begin
        if (npB > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolB, -1.0, dFdXsB.data(),
                    &myBulk.res_n[myBulk.resIndex[bId]], 1.0, &tmpb[0]);
            myLS.AddRhs(bId, tmpb);
        }

        // Assemble mat
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP],
                1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        myLS.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        myLS.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // Assemble rhs
        // End
        if (npE > 2) {
            fill(tmpb.begin(), tmpb.end(), 0.0);
            DaAxpby(ncol, ncolE, -1.0, dFdXsE.data(),
                    &myBulk.res_n[myBulk.resIndex[eId]], 1.0, &tmpb[0]);
            myLS.AddRhs(eId, tmpb);
        }

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP],
                1, bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        myLS.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        myLS.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
}

/////////////////////////////////////////////////////////////////////
// AIMc
/////////////////////////////////////////////////////////////////////

void BulkConn::SetupFIMBulk(Bulk& myBulk, const OCP_BOOL& NRflag) const
{
    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;

    myBulk.bulkTypeAIM.Init();

    OCP_USI  bIdp, bIdc;
    OCP_BOOL flag;

    for (OCP_USI n = 0; n < numBulk; n++) {
        bIdp = n * np;
        bIdc = n * nc;
        flag = OCP_FALSE;
        // CFL
        for (USI j = 0; j < np; j++) {
            if (myBulk.cfl[bIdp + j] > 0.8) {
                flag = OCP_TRUE;
                break;
            }
        }
        // Volume error
        if (!flag) {
            if ((fabs(myBulk.vf[n] - myBulk.rockVp[n]) / myBulk.rockVp[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
        }

        // NR Step
        if (!flag && NRflag) {
            // dP
            if (fabs(myBulk.dPNR[n] / myBulk.P[n]) > 1E-3) {
                flag = OCP_TRUE;
            }
            // dNi
            if (!flag) {
                for (USI i = 0; i < myBulk.numCom; i++) {
                    if (fabs(myBulk.dNNR[bIdc + i] / myBulk.Ni[bIdc + i]) > 1E-3) {
                        flag = OCP_TRUE;
                        break;
                    }
                }
            }
        }

        if (flag) {
            // find it
            // myBulk.map_Bulk2FIM[n] = 1;
            for (auto& v : neighbor[n]) {
                // n is included also
                myBulk.bulkTypeAIM.SetBulkType(v, 1);
            }
        }
    }

    // add WellBulk's 2-neighbor as Implicit bulk
    for (auto& p : myBulk.wellBulkId) {
        for (auto& v : neighbor[p]) {
            for (auto& v1 : neighbor[v]) myBulk.bulkTypeAIM.SetBulkType(v1, 1);
        }
    }
}

void BulkConn::AllocateAIMc_IsoT(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);

    lupblock.resize(numConn * np);
    lupblock_Rho.resize(numConn * np);
    lupblock_Velocity.resize(numConn * np);
}

void BulkConn::AssembleMat_AIMc(LinearSystem&  myLS,
                                const Bulk&    myBulk,
                                const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    myLS.AddDim(numBulk);

    const USI np      = myBulk.numPhase;
    const USI nc      = myBulk.numCom;
    const USI ncol    = nc + 1;
    const USI ncol2   = np * nc + np;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI n = 0; n < numBulk; n++) {
        bmat[0] = myBulk.rockVntg[n] * myBulk.poroP[n] - myBulk.vfP[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }

        myLS.NewDiag(n, bmat);
    }

    // flux term
    OCP_DBL          Akd;
    OCP_DBL          transJ, transIJ;
    vector<OCP_DBL>  dFdXpB(bsize, 0);
    vector<OCP_DBL>  dFdXpE(bsize, 0);
    vector<OCP_DBL>  dFdXsB(bsize2, 0);
    vector<OCP_DBL>  dFdXsE(bsize2, 0);
    OCP_DBL*         dFdXpI;
    OCP_DBL*         dFdXsI;
    vector<OCP_BOOL> phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistE(np, OCP_FALSE);
    vector<OCP_BOOL> phaseExistI(np, OCP_FALSE);
    OCP_BOOL         phaseExistU;
    vector<OCP_BOOL> phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_E(np, OCP_FALSE);
    vector<OCP_BOOL> phasedS_I(np, OCP_FALSE);
    vector<USI>      pVnumComB(np, 0);
    vector<USI>      pVnumComE(np, 0);
    vector<USI>      pVnumComI(np, 0);
    OCP_USI          ncolB, ncolE, ncolI;

    OCP_USI  bId, eId, uId;
    OCP_USI  bId_np_j, eId_np_j, uId_np_j, ibId_np_j;
    OCP_DBL  kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL  dP, dGamma;
    OCP_DBL  tmp;
    OCP_BOOL bIdFIM, eIdFIM, uIdFIM;


    for (OCP_USI c = 0; c < numConn; c++) {
        bId    = iteratorConn[c].bId;
        eId    = iteratorConn[c].eId;
        Akd    = CONV1 * CONV2 * iteratorConn[c].area;
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);
        bIdFIM = eIdFIM = OCP_FALSE;
        if (myBulk.bulkTypeAIM.IfFIMbulk(bId)) bIdFIM = OCP_TRUE;
        if (myBulk.bulkTypeAIM.IfFIMbulk(eId)) eIdFIM = OCP_TRUE;

        if (!bIdFIM && !eIdFIM) {
            // both are explicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = myBulk.phaseExist[bId * np + j];
                phaseExistE[j] = myBulk.phaseExist[eId * np + j];
                uId            = upblock[c * np + j];
                phaseExistU    = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    continue;
                }
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                xi       = myBulk.xi[uId_np_j];
                kr       = myBulk.kr[uId_np_j];
                mu       = myBulk.mu[uId_np_j];
                transJ   = Akd * xi * kr / mu;
                for (USI i = 0; i < nc; i++) {
                    xij     = myBulk.xij[uId_np_j * nc + i];
                    transIJ = xij * transJ;
                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;
                }
            }
        } else if ((bIdFIM && !eIdFIM) || (!bIdFIM && eIdFIM)) {

            // one is explicit, one is implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            ncolI   = 0;
            USI jxI = 0;
            if (bIdFIM) {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = myBulk.phaseExist[bId * np + j];
                    phasedS_I[j]   = myBulk.pSderExist[bId * np + j];
                    if (phasedS_I[j]) jxI++;
                    pVnumComI[j] = myBulk.pVnumCom[bId * np + j];
                    ncolI += pVnumComI[j];
                }
                dFdXpI = &dFdXpB[0];
                dFdXsI = &dFdXsB[0];
            } else {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = myBulk.phaseExist[eId * np + j];
                    phasedS_I[j]   = myBulk.pSderExist[eId * np + j];
                    if (phasedS_I[j]) jxI++;
                    pVnumComI[j] = myBulk.pVnumCom[eId * np + j];
                    ncolI += pVnumComI[j];
                }
                dFdXpI = &dFdXpE[0];
                dFdXsI = &dFdXsE[0];
            }
            ncolI += jxI;

            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = myBulk.phaseExist[bId * np + j];
                phaseExistE[j] = myBulk.phaseExist[eId * np + j];
                uId            = upblock[c * np + j];
                phaseExistU    = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    continue;
                }

                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                dP       = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                     upblock_Rho[c * np + j] * dGamma;
                xi     = myBulk.xi[uId_np_j];
                kr     = myBulk.kr[uId_np_j];
                mu     = myBulk.mu[uId_np_j];
                transJ = Akd * kr / mu;
                uIdFIM = (uId == bId ? bIdFIM : eIdFIM);

                if (bIdFIM)
                    ibId_np_j = bId_np_j;
                else
                    ibId_np_j = eId_np_j;

                if (uIdFIM) {
                    // upblock is implicit
                    OCP_DBL rhoWgt = 1.0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;

                    muP  = myBulk.muP[ibId_np_j];
                    xiP  = myBulk.xiP[ibId_np_j];
                    rhoP = myBulk.rhoP[ibId_np_j];
                    for (USI i = 0; i < nc; i++) {
                        xij     = myBulk.xij[ibId_np_j * nc + i];
                        transIJ = xij * xi * transJ;

                        // Pressure -- Primary var
                        dFdXpB[(i + 1) * ncol] += transIJ;
                        dFdXpE[(i + 1) * ncol] -= transIJ;

                        tmp = transIJ * (-rhoP * dGamma) * rhoWgt;
                        tmp += xij * transJ * xiP * dP;
                        tmp += -transIJ * muP / mu * dP;
                        dFdXpI[(i + 1) * ncol] += tmp;

                        USI j1S = 0;
                        // Saturation -- Second var
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_I[j1]) {
                                dFdXsI[(i + 1) * ncolI + j1S] +=
                                    transIJ * myBulk.dPcj_dS[ibId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      myBulk.dKr_dS[ibId_np_j * np + j1] * dP;
                                dFdXsI[(i + 1) * ncolI + j1S] += tmp;
                                j1S++;
                            }
                        }
                        // Cij
                        for (USI k = 0; k < pVnumComI[j]; k++) {
                            rhox = myBulk.rhox[ibId_np_j * nc + k] * rhoWgt;
                            xix  = myBulk.xix[ibId_np_j * nc + k];
                            mux  = myBulk.mux[ibId_np_j * nc + k];
                            tmp  = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComI[j])
                            dFdXsI[(i + 1) * ncolI + jxI + i] += xi * transJ * dP;
                    }
                    jxI += pVnumComI[j];
                } else {
                    // downblock is implicit
                    OCP_DBL rhoWgt = 0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;
                    rhoP = myBulk.rhoP[ibId_np_j];

                    for (USI i = 0; i < nc; i++) {
                        xij     = myBulk.xij[uId_np_j * nc + i];
                        transIJ = xij * xi * transJ;

                        // Pressure -- Primary var
                        dFdXpB[(i + 1) * ncol] += transIJ;
                        dFdXpE[(i + 1) * ncol] -= transIJ;

                        tmp = transIJ * (-rhoP * dGamma) * rhoWgt;
                        dFdXpI[(i + 1) * ncol] += tmp;

                        USI j1S = 0;
                        // Saturation -- Second var
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_I[j1]) {
                                dFdXsI[(i + 1) * ncolI + j1S] +=
                                    transIJ * myBulk.dPcj_dS[ibId_np_j * np + j1];
                                j1S++;
                            }
                        }
                        // Cij
                        for (USI k = 0; k < pVnumComI[j]; k++) {
                            rhox = myBulk.rhox[ibId_np_j * nc + k] * rhoWgt;
                            tmp  = -transIJ * rhox * dGamma;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                    }
                    jxI += pVnumComI[j];
                }
            }
        } else {
            // both are implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            USI jxB = 0;
            USI jxE = 0;
            ncolB   = 0;
            ncolE   = 0;

            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = myBulk.phaseExist[bId * np + j];
                phaseExistE[j] = myBulk.phaseExist[eId * np + j];
                phasedS_B[j]   = myBulk.pSderExist[bId * np + j];
                phasedS_E[j]   = myBulk.pSderExist[eId * np + j];
                if (phasedS_B[j]) jxB++;
                if (phasedS_E[j]) jxE++;
                pVnumComB[j] = myBulk.pVnumCom[bId * np + j];
                pVnumComE[j] = myBulk.pVnumCom[eId * np + j];
                ncolB += pVnumComB[j];
                ncolE += pVnumComE[j];
            }
            ncolB += jxB;
            ncolE += jxE;

            for (USI j = 0; j < np; j++) {
                uId = upblock[c * np + j];

                phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    jxB += pVnumComB[j];
                    jxE += pVnumComE[j];
                    continue;
                }

                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                dP       = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                     upblock_Rho[c * np + j] * dGamma;
                xi     = myBulk.xi[uId_np_j];
                kr     = myBulk.kr[uId_np_j];
                mu     = myBulk.mu[uId_np_j];
                muP    = myBulk.muP[uId_np_j];
                xiP    = myBulk.xiP[uId_np_j];
                rhoP   = myBulk.rhoP[uId_np_j];
                transJ = Akd * kr / mu;

                for (USI i = 0; i < nc; i++) {
                    xij     = myBulk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    tmp = xij * transJ * xiP * dP;
                    tmp += -transIJ * muP / mu * dP;
                    if (!phaseExistE[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpB[(i + 1) * ncol] += tmp;
                    } else if (!phaseExistB[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpE[(i + 1) * ncol] += tmp;
                    } else {
                        dFdXpB[(i + 1) * ncol] +=
                            transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                        dFdXpE[(i + 1) * ncol] +=
                            transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                        if (bId == uId) {
                            dFdXpB[(i + 1) * ncol] += tmp;
                        } else {
                            dFdXpE[(i + 1) * ncol] += tmp;
                        }
                    }

                    // Second var
                    USI j1SB = 0;
                    USI j1SE = 0;
                    if (bId == uId) {
                        // Saturation
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_B[j1]) {
                                dFdXsB[(i + 1) * ncolB + j1SB] +=
                                    transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                                dFdXsB[(i + 1) * ncolB + j1SB] += tmp;
                                j1SB++;
                            }
                            if (phasedS_E[j1]) {
                                dFdXsE[(i + 1) * ncolE + j1SE] -=
                                    transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                                j1SE++;
                            }
                        }
                        // Cij
                        if (!phaseExistE[j]) {
                            for (USI k = 0; k < pVnumComB[j]; k++) {
                                rhox = myBulk.rhox[uId_np_j * nc + k];
                                xix  = myBulk.xix[uId_np_j * nc + k];
                                mux  = myBulk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        } else {
                            for (USI k = 0; k < pVnumComB[j]; k++) {
                                rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                                xix  = myBulk.xix[uId_np_j * nc + k];
                                mux  = myBulk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                                dFdXsE[(i + 1) * ncolE + jxE + k] +=
                                    -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 *
                                    dGamma;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        }
                    } else {
                        // Saturation
                        for (USI j1 = 0; j1 < np; j1++) {
                            if (phasedS_B[j1]) {
                                dFdXsB[(i + 1) * ncolB + j1SB] +=
                                    transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                                j1SB++;
                            }
                            if (phasedS_E[j1]) {
                                dFdXsE[(i + 1) * ncolE + j1SE] -=
                                    transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                                tmp = Akd * xij * xi / mu *
                                      myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                                dFdXsE[(i + 1) * ncolE + j1SE] += tmp;
                                j1SE++;
                            }
                        }
                        // Cij
                        if (!phaseExistB[j]) {
                            for (USI k = 0; k < pVnumComE[j]; k++) {
                                rhox = myBulk.rhox[uId_np_j * nc + k];
                                xix  = myBulk.xix[uId_np_j * nc + k];
                                mux  = myBulk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComE[j])
                                dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                        } else {
                            for (USI k = 0; k < pVnumComE[j]; k++) {
                                rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                                xix  = myBulk.xix[uId_np_j * nc + k];
                                mux  = myBulk.mux[uId_np_j * nc + k];
                                tmp  = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                                dFdXsB[(i + 1) * ncolB + jxB + k] +=
                                    -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 *
                                    dGamma;
                            }
                            // WARNING !!!
                            if (i < pVnumComE[j])
                                dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                        }
                    }
                }
                jxB += pVnumComB[j];
                jxE += pVnumComE[j];
            }
        }

        if (bIdFIM && !eIdFIM)
            ncolB = ncolI;
        else if (!bIdFIM && eIdFIM)
            ncolE = ncolI;

        // Assemble
        bmat = dFdXpB;
        if (bIdFIM) {
            DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(),
                    &myBulk.dSec_dPri[bId * lendSdP], 1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - Begin -- add
        myLS.AddDiag(bId, bmat);
        // End - Begin -- insert
        Dscalar(bsize, -1, bmat.data());
        myLS.NewOffDiag(eId, bId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        if (eIdFIM) {
            DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(),
                    &myBulk.dSec_dPri[eId * lendSdP], 1, bmat.data());
        }
        Dscalar(bsize, dt, bmat.data());
        // Begin - End -- insert
        myLS.NewOffDiag(bId, eId, bmat);
        // End - End -- add
        Dscalar(bsize, -1, bmat.data());
        myLS.AddDiag(eId, bmat);

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0])) {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif
    }
}

void BulkConn::CalResAIMc(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    // IMPORTANT!!!
    // in AIMc for IMPEC Bulk, P was updated in each Newton Step, but Pj didn't.
    // So here Pj must be replaced with P + Pcj, otherwise wrong results will arise

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
    OCP_DBL rho, dP;
    OCP_DBL tmp, dNi;
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].bId;
        eId = iteratorConn[c].eId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];
            OCP_BOOL exup = exbegin;

            if ((exbegin) && (exend)) {
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j]          = bId;
                upblock_Velocity[c * np + j] = 0;
                upblock_Rho[c * np + j]      = 0;
                continue;
            }

            uId           = bId;           
            dP            = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                 (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            }
            
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;
            uId_np_j = uId * np + j;

            if (exup) {
                upblock_Velocity[c * np + j] =
                    Akd * myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
            } else {
                upblock_Velocity[c * np + j] = 0;
                continue;
            }

            tmp = dt * myBulk.xi[uId_np_j] * upblock_Velocity[c * np + j];
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