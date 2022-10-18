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
                iteratorConn.push_back(BulkPair(bIdb, tmp2[j].id, tmp2[j].area));
            }
            neighborNum[bIdb] = len;
        }
    }

    numConn = iteratorConn.size();

    // PrintConnectionInfoCoor(myGrid);
}


void BulkConn::SetupWellBulk_K(Bulk& myBulk) const
{ 
    // For K = 1 now, defaulted

    USI len = myBulk.wellBulkId.size();
    for (USI n = 0; n < len; n++) {
        for (auto& v : neighbor[n]) {
            USI clen = myBulk.wellBulkId.size();
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
    OCP_USI bIdb, eIdb;
    USI I, J, K;
    USI sp = myGrid.GetNumDigitIJK();
    cout << "BulkConn : " << numConn << endl;
    for (OCP_USI c = 0; c < numConn; c++) {
        bIdb = iteratorConn[c].BId;
        eIdb = iteratorConn[c].EId;
        bIdg = myGrid.activeMap_B2G[bIdb];
        eIdg = myGrid.activeMap_B2G[eIdb];
        myGrid.GetIJKGrid(I, J, K, bIdg);
        cout << GetIJKformat(to_string(I), to_string(J), to_string(K), sp) << "   ";
        cout << setw(6) << bIdg;
        cout << "    ";
        cout << setw(6) << bIdb;
        cout << "    ";
        myGrid.GetIJKGrid(I, J, K, eIdg);
        cout << GetIJKformat(to_string(I), to_string(J), to_string(K), sp) << "   ";        
        cout << setw(6) << eIdg;
        cout << "    ";
        cout << setw(6) << eIdb;
        cout << setw(20) << setprecision(8) << fixed << iteratorConn[c].area * CONV2;

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
        vf = myBulk.vf[n];
        vfp = myBulk.vfp[n];
        P = myBulk.lP[n];
        Vp0 = myBulk.rockVpInit[n];
        Vp  = myBulk.rockVp[n];
                
        OCP_DBL temp    = cr * Vp0 - vfp;
        myLS.diagVal[n] = temp;
        myLS.b[n] = temp * P + dt * (vf - Vp);
        // myLS.b[n] = temp * P + (vf - Vp);
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

        USI diagptr = myLS.diagPtr[bId];
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
        if (myLS.val[n].size() == myLS.diagPtr[n]) myLS.val[n].push_back(myLS.diagVal[n]);
    }
}


void BulkConn::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
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

    //static USI myiter = 0;
    //myiter++;

    // calculate a step flux using iteratorConn
    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j;
    OCP_DBL Pbegin, Pend, rho;
    USI     np = myBulk.numPhase;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId         = iteratorConn[c].BId;
        eId         = iteratorConn[c].EId;
        OCP_DBL Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];

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
            OCP_BOOL    exup = exbegin;
            OCP_DBL dP   = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                         (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId  = eId;
                exup = exend;
            }

            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j]     = uId;

            if (exup) {
                OCP_USI uId_np_j = uId * np + j;
                OCP_DBL trans =  Akd * myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j];
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
    upblock_Rho.resize(numConn * np);
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
    for (USI i = 1; i < ncol; i++) {
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
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_BOOL    phaseExistBj, phaseExistEj;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;

    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
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
                xij = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;
                                           
                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;
               
                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistEj) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                }
                else if (!phaseExistBj) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                }
                else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }
                // Second var
				if (bId == uId) {
                    // Saturation 
					for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[bId_np_j * np + k];                       
						tmp = Akd * xij * xi / mu *
							myBulk.dKr_dS[uId_np_j * np + k] * dP;
						dFdXsB[(i + 1) * ncol2 + k] += tmp;
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * myBulk.dPcj_dS[eId_np_j * np + k];
					}
                    // Cij                    
                    if (!phaseExistEj) {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        dFdXsB[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }                   
				}
				else {
                    // Saturation
					for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[bId_np_j * np + k];
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * myBulk.dPcj_dS[eId_np_j * np + k];
						tmp = Akd * xij * xi / mu *
							myBulk.dKr_dS[uId_np_j * np + k] * dP;
						dFdXsE[(i + 1) * ncol2 + k] += tmp;
					}
                    // Cij                   
                    if (!phaseExistBj) {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        dFdXsE[(i + 1) * ncol2 + np + j * nc + i] += xi * transJ * dP;
                    }                    
				}
            }
        }

        USI diagptr = myLS.diagPtr[bId];

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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif

    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
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
    OCP_DBL   dGamma, rho, dP;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId         = iteratorConn[c].BId;
        eId         = iteratorConn[c].EId;
        // dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];           
            if ((exbegin) && (exend)) {                              
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j] = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }
            uId        = bId;
            
            dP = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock[c * np + j] = uId;
            upblock_Rho[c * np + j] = rho;
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
    OCP_DBL dGamma, rho, dP;
    OCP_DBL tmp, dNi;
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend   = myBulk.phaseExist[eId_np_j];
           
            if ((exbegin) && (exend)) {               
                rho    = (myBulk.rho[bId_np_j] + myBulk.rho[eId_np_j]) / 2;
            } else if (exbegin && (!exend)) {
                rho    = myBulk.rho[bId_np_j];
            } else if ((!exbegin) && (exend)) {
                rho    = myBulk.rho[eId_np_j];
            } else {
                upblock[c * np + j] = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }

            uId = bId;           
            dP  = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                 (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            //dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);
            //dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho;
            //cout << setprecision(6) << scientific << dP << "   "
            //    << myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho << "   "
            //    << dP - (myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] - dGamma * rho)
            //    << endl;
            if (dP < 0) {
                uId = eId;
            }
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j] = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            tmp = dt * Akd * myBulk.xi[uId_np_j] *
                  myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;

            for (USI i = 0; i < nc; i++) {
                dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                res[bId * len + 1 + i] += dNi;
                res[eId * len + 1 + i] -= dNi;
            }
        }
    }
}


/// rho = (S1*rho1 + S2*rho2)/(S1+S2)
void BulkConn::CalFluxFIMS(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    OCP_USI   bId, eId, uId;
    OCP_USI   bId_np_j, eId_np_j;
    OCP_DBL   Pbegin, Pend, rho;

    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend = myBulk.phaseExist[eId_np_j];
           
            if ((exbegin) && (exend)) {               
                rho = (myBulk.rho[bId_np_j]* myBulk.S[bId_np_j] + myBulk.rho[eId_np_j] * myBulk.S[eId_np_j])
                    / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
            }
            else if (exbegin && (!exend)) {
                rho = myBulk.rho[bId_np_j];
            }
            else if ((!exbegin) && (exend)) {
                rho = myBulk.rho[eId_np_j];
            }
            else {
                upblock[c * np + j] = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }
            Pbegin = myBulk.Pj[bId_np_j];
            Pend = myBulk.Pj[eId_np_j];
            uId = bId;
            OCP_DBL dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock[c * np + j] = uId;
            upblock_Rho[c * np + j] = rho;
        }
    }
}

void BulkConn::CalResFIMS(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI len = nc + 1;
    OCP_USI   bId, eId, uId, bIdb;

    // Accumalation Term
    for (OCP_USI n = 0; n < numBulk; n++) {

        bId = n * len;
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
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        for (USI j = 0; j < np; j++) {
            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;

            OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
            OCP_BOOL exend = myBulk.phaseExist[eId_np_j];

            if ((exbegin) && (exend)) {               
                rho = (myBulk.rho[bId_np_j] * myBulk.S[bId_np_j] + myBulk.rho[eId_np_j] * myBulk.S[eId_np_j])
                    / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
            }
            else if (exbegin && (!exend)) {
                rho = myBulk.rho[bId_np_j];
            }
            else if ((!exbegin) && (exend)) {
                rho = myBulk.rho[eId_np_j];
            }
            else {
                upblock[c * np + j] = bId;
                upblock_Rho[c * np + j] = 0;
                continue;
            }
            Pbegin = myBulk.Pj[bId_np_j];
            Pend = myBulk.Pj[eId_np_j];
            uId = bId;
            dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
            if (dP < 0) {
                uId = eId;
            }
            upblock_Rho[c * np + j] = rho;
            upblock[c * np + j] = uId;

            uId_np_j = uId * np + j;
            if (!myBulk.phaseExist[uId_np_j]) continue;
            tmp = dt * Akd * myBulk.xi[uId_np_j] *
                myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;

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


void BulkConn::AssembleMat_FIM_new(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
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
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistE(np, OCP_FALSE);
    OCP_BOOL            phaseExistU;
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_E(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<USI>     pVnumComE(np, 0);
    USI             ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;


    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        USI jxB = 0;  USI jxE = 0;
        ncolB = 0;      ncolE = 0;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            phasedS_B[j] = myBulk.pSderExist[bId * np + j];
            phasedS_E[j] = myBulk.pSderExist[eId * np + j];
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
            dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                upblock_Rho[c * np + j] * dGamma;
            //dP = (myBulk.Pj[bId_np_j] - GRAVITY_FACTOR * upblock_Rho[c * np + j] * myBulk.depth[bId]) -
            //    (myBulk.Pj[eId_np_j] - GRAVITY_FACTOR * upblock_Rho[c * np + j] * myBulk.depth[eId]);
            xi = myBulk.xi[uId_np_j];
            kr = myBulk.kr[uId_np_j];
            mu = myBulk.mu[uId_np_j];
            muP = myBulk.muP[uId_np_j];
            xiP = myBulk.xiP[uId_np_j];
            rhoP = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;


            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[uId_np_j * nc + i];
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
                USI j1SB = 0; USI j1SE = 0;
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }                    
                }
                else {
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
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

        USI diagptr = myLS.diagPtr[bId];

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
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + n * bsize + bsize);
    }
}


void BulkConn::AssembleMat_FIM_new1(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    // less than AssembleMat_FIM_new() but loops can be execute with multithreads

    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
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
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistE(np, OCP_FALSE);
    OCP_BOOL            phaseExistU;
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_E(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<USI>     pVnumComE(np, 0);
    vector<USI>     jsB(np, 0);
    vector<USI>     jsE(np, 0);
    vector<USI>     jxB(np, 0);
    vector<USI>     jxE(np, 0);
    USI             ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;


    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        jxB[0] = 0;     jxE[0] = 0;
        fill(jsB.begin(), jsB.end(), 0);
        fill(jsE.begin(), jsE.end(), 0);
        ncolB = 0;      ncolE = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            phasedS_B[j] = myBulk.pSderExist[bId * np + j];
            phasedS_E[j] = myBulk.pSderExist[eId * np + j];
            if (phasedS_B[j]) {
                jxB[0]++;
                for (USI j1 = j + 1; j1 < np; j1++) jsB[j1]++;
            }
            if (phasedS_E[j]) {
                jxE[0]++;
                for (USI j1 = j + 1; j1 < np; j1++) jsE[j1]++;
            }
            pVnumComB[j] = myBulk.pVnumCom[bId * np + j];
            pVnumComE[j] = myBulk.pVnumCom[eId * np + j];
            ncolB += pVnumComB[j];
            ncolE += pVnumComE[j];
        }
        ncolB += jxB[0];
        ncolE += jxE[0];
        for (USI j = 0; j < np - 1; j++) {
            jxB[j + 1] = jxB[j] + pVnumComB[j];
            jxE[j + 1] = jxE[j] + pVnumComE[j];
            
        }

        for (USI j = 0; j < np; j++) {
            uId = upblock[c * np + j];

            phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
            if (!phaseExistU) {
                continue;
            }              

            bId_np_j = bId * np + j;
            eId_np_j = eId * np + j;
            uId_np_j = uId * np + j;
            dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                upblock_Rho[c * np + j] * dGamma;
            xi = myBulk.xi[uId_np_j];
            kr = myBulk.kr[uId_np_j];
            mu = myBulk.mu[uId_np_j];
            muP = myBulk.muP[uId_np_j];
            xiP = myBulk.xiP[uId_np_j];
            rhoP = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;


            for (USI i = 0; i < nc; i++) {

                xij = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                }
                else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                }
                else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + jsB[j1]] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsB[(i + 1) * ncolB + jsB[j1]] += tmp;
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + jsE[j1]] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];                            
                        }
                    }
                    // Cij
                    const USI tmpId = pVnumComB[j];
                    if (!phaseExistE[j]) {                        
                        for (USI k = 0; k < tmpId; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB[j] + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < tmpId)
                            dFdXsB[(i + 1) * ncolB + jxB[j] + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < tmpId; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB[j] + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE[j] + k] += -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < tmpId)
                            dFdXsB[(i + 1) * ncolB + jxB[j] + i] += xi * transJ * dP;
                    }
                }
                else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {
                        if (phasedS_B[j1]) {
                            dFdXsB[(i + 1) * ncolB + jsB[j1]] +=
                                transIJ * myBulk.dPcj_dS[bId_np_j * np + j1];
                        }
                        if (phasedS_E[j1]) {
                            dFdXsE[(i + 1) * ncolE + jsE[j1]] -=
                                transIJ * myBulk.dPcj_dS[eId_np_j * np + j1];
                            tmp = Akd * xij * xi / mu *
                                myBulk.dKr_dS[uId_np_j * np + j1] * dP;
                            dFdXsE[(i + 1) * ncolE + jsE[j1]] += tmp;
                        }
                    }
                    // Cij
                    const USI tmpId = pVnumComE[j];
                    if (!phaseExistB[j]) {
                        for (USI k = 0; k < tmpId; k++) {
                            rhox = myBulk.rhox[uId_np_j * nc + k];
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE[j] + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < tmpId)
                            dFdXsE[(i + 1) * ncolE + jxE[j] + i] += xi * transJ * dP;
                    }
                    else {
                        for (USI k = 0; k < tmpId; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE[j] + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB[j] + k] += -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < tmpId)
                            dFdXsE[(i + 1) * ncolE + jxE[j] + i] += xi * transJ * dP;
                    }
                }
            }
        }

        USI diagptr = myLS.diagPtr[bId];

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
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + n * bsize + bsize);
    }
}


void BulkConn::AssembleMat_FIM_newS(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
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
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistE(np, OCP_FALSE);
    OCP_BOOL            phaseExistU;
    vector<USI>     pEnumComB(np, 0);
    vector<USI>     pEnumComE(np, 0);
    USI             ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;
    OCP_DBL wghtb, wghte;

    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        const USI npB = myBulk.phaseNum[bId] + 1;    ncolB = npB;
        const USI npE = myBulk.phaseNum[eId] + 1;    ncolE = npE;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            pEnumComB[j] = myBulk.pVnumCom[bId * np + j];
            pEnumComE[j] = myBulk.pVnumCom[eId * np + j];
            ncolB += pEnumComB[j];
            ncolE += pEnumComE[j];
        }


        USI jxB = npB;  USI jxE = npE;
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
            dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                upblock_Rho[c * np + j] * dGamma;
            xi = myBulk.xi[uId_np_j];
            kr = myBulk.kr[uId_np_j];
            mu = myBulk.mu[uId_np_j];
            muP = myBulk.muP[uId_np_j];
            xiP = myBulk.xiP[uId_np_j];
            rhoP = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;


            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;
                            
                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                }
                else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                }
                else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * dGamma * (-myBulk.rhoP[bId_np_j] * myBulk.S[bId_np_j])
                        / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * dGamma * (-myBulk.rhoP[eId_np_j] * myBulk.S[eId_np_j])
                        / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }

                // Second var
                USI j1SB = 0; USI j1SE = 0;
                if (bId == uId) {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0; wghte = 0;
                        if (j1 == j && phaseExistE[j]) {
                            tmp = -dGamma / ((myBulk.S[bId_np_j] + myBulk.S[eId_np_j]) * (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]));
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }
                    else {
                        wghtb = myBulk.S[bId_np_j] / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        wghte = myBulk.S[eId_np_j] / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] * wghtb;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += -transIJ * dGamma * myBulk.rhox[eId_np_j * nc + k] * wghte;
                        }
                        // WARNING !!!
                        if (i < pEnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                    }                   
                }
                else {
                    // Saturation
                    for (USI j1 = 0; j1 < np; j1++) {

                        wghtb = 0; wghte = 0;
                        if (j1 == j && phaseExistB[j]) {
                            tmp = -dGamma / ((myBulk.S[bId_np_j] + myBulk.S[eId_np_j]) * (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]));
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pEnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                    }
                    else {
                        wghtb = myBulk.S[bId_np_j] / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        wghte = myBulk.S[eId_np_j] / (myBulk.S[bId_np_j] + myBulk.S[eId_np_j]);
                        for (USI k = 0; k < pEnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] * wghte;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += -transIJ * dGamma * myBulk.rhox[bId_np_j * nc + k] * wghtb;
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

        USI diagptr = myLS.diagPtr[bId];

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
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + n * bsize + bsize);
    }
}


void BulkConn::AssembleMat_FIM_new_n(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < ncol; i++) {
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
        myLS.b[n * ncol] -= myBulk.resPc[n];
    }

    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistE(np, OCP_FALSE);
    OCP_BOOL            phaseExistU;
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_E(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<USI>     pVnumComE(np, 0);
    USI             ncolB, ncolE;

    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;


    // Becareful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
        fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
        fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
        fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

        const USI npB = myBulk.phaseNum[bId] + 1;
        const USI npE = myBulk.phaseNum[eId] + 1;

        USI jxB = 0;  USI jxE = 0;
        ncolB = 0;      ncolE = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[bId * np + j];
            phaseExistE[j] = myBulk.phaseExist[eId * np + j];
            phasedS_B[j] = myBulk.pSderExist[bId * np + j];
            phasedS_E[j] = myBulk.pSderExist[eId * np + j];
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
            dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                upblock_Rho[c * np + j] * dGamma;
            //dP = myBulk.Pj[bId * np + j] - myBulk.Pj[eId * np + j] -
            //    myBulk.rho[uId * np + j] * dGamma;
            xi = myBulk.xi[uId_np_j];
            kr = myBulk.kr[uId_np_j];
            mu = myBulk.mu[uId_np_j];
            muP = myBulk.muP[uId_np_j];
            xiP = myBulk.xiP[uId_np_j];
            rhoP = myBulk.rhoP[uId_np_j];
            transJ = Akd * kr / mu;


            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[uId_np_j * nc + i];
                transIJ = xij * xi * transJ;

                // Pressure -- Primary var
                dFdXpB[(i + 1) * ncol] += transIJ;
                dFdXpE[(i + 1) * ncol] -= transIJ;

                
                tmp = xij * transJ * xiP * dP;
                tmp += -transIJ * muP / mu * dP;
                if (!phaseExistE[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpB[(i + 1) * ncol] += tmp;
                }
                else if (!phaseExistB[j]) {
                    tmp += transIJ * (-rhoP * dGamma);
                    dFdXpE[(i + 1) * ncol] += tmp;
                }
                else {
                    dFdXpB[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                    dFdXpE[(i + 1) * ncol] +=
                        transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                    if (bId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                }
                
                // Second var
                USI j1SB = 0; USI j1SE = 0;
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP / myBulk.nj[uId_np_j];
                    }
                    else {
                        for (USI k = 0; k < pVnumComB[j]; k++) {
                            rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            dFdXsE[(i + 1) * ncolE + jxE + k] += -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComB[j])
                            dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP / myBulk.nj[uId_np_j];
                    }                   
                }
                else {
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
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP / myBulk.nj[uId_np_j];
                    }
                    else {
                        for (USI k = 0; k < pVnumComE[j]; k++) {
                            rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                            xix = myBulk.xix[uId_np_j * nc + k];
                            mux = myBulk.mux[uId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            tmp -= transIJ * dP / myBulk.nj[uId_np_j];
                            dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            dFdXsB[(i + 1) * ncolB + jxB + k] += -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
                        }
                        // WARNING !!!
                        if (i < pVnumComE[j])
                            dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP / myBulk.nj[uId_np_j];
                    }                   
                }
            }
            jxB += pVnumComB[j];
            jxE += pVnumComE[j];
        }

        USI diagptr = myLS.diagPtr[bId];

        if (bId != lastbId) {
            // new bulk
            assert(myLS.val[bId].size() == diagptr * bsize);
            OCP_USI id = bId * bsize;
            myLS.val[bId].insert(myLS.val[bId].end(), myLS.diagVal.data() + id,
                myLS.diagVal.data() + id + bsize);

            lastbId = bId;
        }

        // Assemble rhs
        // Begin
        if (npB > 2) {
            DaAxpby(ncol, ncolB, -1.0, dFdXsB.data(),
                &myBulk.res_n[myBulk.resIndex[bId]], 1.0, &myLS.b[bId * ncol]);

            //cout << "----------- " << bId << " ------------" << endl;
            //PrintDX(ncol, &myLS.b[bId * ncol]);
        }

        // Assemble mat
        bmat = dFdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // Assemble rhs
        // End
        if (npE > 2) {
            DaAxpby(ncol, ncolE, -1.0, dFdXsE.data(),
                &myBulk.res_n[myBulk.resIndex[eId]], 1.0, &myLS.b[eId * ncol]);

            //cout << "----------- " << eId << " ------------" << endl;
            //PrintDX(ncol, &myLS.b[eId * ncol]);
        }

        // End
        bmat = dFdXpE;
        DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP], 1,
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or INF in bmat !");
        }
#endif
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + n * bsize + bsize);
    }
}


/////////////////////////////////////////////////////////////////////
// AIMt
/////////////////////////////////////////////////////////////////////

void BulkConn::SetupFIMBulk(Bulk& myBulk, const OCP_BOOL& NRflag) const
{
    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;

    myBulk.FIMBulk.clear();
    fill(myBulk.map_Bulk2FIM.begin(), myBulk.map_Bulk2FIM.end(), -1);

    OCP_USI bIdp, bIdc;
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
            // cout << "[" << n << "]" << endl;
            // dP
            if (fabs(myBulk.dPNR[n] / myBulk.P[n]) > 1E-3) {        
                flag = OCP_TRUE;
                //cout << scientific << "P  " << myBulk.dPNR[n] / myBulk.P[n] << "   ";
                //cout << myBulk.dPNR[n] << "   " << myBulk.P[n] << "   ";
                //cout << endl;
            }
            // dNi
            if (!flag) {
                for (USI i = 0; i < myBulk.numCom; i++) {
                    if (fabs(myBulk.dNNR[bIdc + i] / myBulk.Ni[bIdc + i]) > 1E-3) {
                        flag = OCP_TRUE;
                        //cout << scientific << "Ni[" << i << "]  " << myBulk.dNNR[bIdc + i] / myBulk.Ni[bIdc + i]<< "   ";
                        //cout << myBulk.dNNR[bIdc + i] << "   " << myBulk.Ni[bIdc + i] << "   ";
                        //cout << endl;
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
                myBulk.map_Bulk2FIM[v] = 1;
            }
        }
    }

    // add WellBulk
    for (auto& p : myBulk.wellBulkId) {
        for (auto& v : neighbor[p]) {
            for (auto& v1 : neighbor[v])
                myBulk.map_Bulk2FIM[v1] = 1;
        }
    }
    USI iter = 0;
    for (OCP_USI n = 0; n < myBulk.numBulk; n++) {
        if (myBulk.map_Bulk2FIM[n] > 0) {
            myBulk.map_Bulk2FIM[n] = iter;
            iter++;
            myBulk.FIMBulk.push_back(n);
        }
    }
    myBulk.numFIMBulk = myBulk.FIMBulk.size();

    //myBulk.numFIMBulk = numBulk;
    //myBulk.FIMBulk.resize(numBulk);
    //for (OCP_USI n = 0; n < numBulk; n++) {
    //    myBulk.map_Bulk2FIM[n] = n;
    //    myBulk.FIMBulk[n] = n;
    //}      
}

void BulkConn::AddFIMBulk(Bulk& myBulk)
{
    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;

    fill(myBulk.map_Bulk2FIM.begin(), myBulk.map_Bulk2FIM.end(), -1);

    OCP_USI bIdp, bIdc;
    OCP_BOOL flag;

    for (OCP_USI n = 0; n < numBulk; n++) {
        bIdp = n * np;
        bIdc = n * nc;
        flag = OCP_FALSE;
        // cfl
        //for (USI j = 0; j < np; j++) {
        //    if (myBulk.cfl[bIdp + j] > 0.8) {
        //        flag = OCP_TRUE;
        //        break;
        //    }
        //}
        // Ni
        if (!flag) {
            for (USI i = 0; i < nc; i++) {
                if (myBulk.Ni[bIdc + i] < 0) {
                    flag = OCP_TRUE;
                    break;
                }
            }
        }
        // Volume error
        if (!flag) {
            if ((fabs(myBulk.vf[n] - myBulk.rockVp[n]) / myBulk.rockVp[n]) > 0.001) {
                flag = OCP_TRUE;
            }
        }

        if (flag) {
            // find it
            for (auto& v : neighbor[n]) {
                // n is included also
                myBulk.map_Bulk2FIM[v] = 1;
            }
        }
    }
    // add last FIMBulk
    for (USI i = 0; i < myBulk.numFIMBulk; i++) {
        myBulk.map_Bulk2FIM[myBulk.FIMBulk[i]] = 1;
    }

    // add WellBulk
    for (auto& p : myBulk.wellBulkId) {
        for (auto& v : neighbor[p]) {
            for (auto& v1 : neighbor[v])
                myBulk.map_Bulk2FIM[v1] = 1;
        }
    }

    myBulk.FIMBulk.clear();
    USI iter = 0;
    for (OCP_USI n = 0; n < myBulk.numBulk; n++) {
        if (myBulk.map_Bulk2FIM[n] > 0) {
            myBulk.map_Bulk2FIM[n] = iter;
            iter++;
            myBulk.FIMBulk.push_back(n);
        }
    }
    myBulk.numFIMBulk = myBulk.FIMBulk.size();

}

void BulkConn::SetupFIMBulkBoundAIMs(Bulk& myBulk)
{
    USI iter = 0;
    OCP_USI n;
    for (USI fn = 0; fn < myBulk.numFIMBulk; fn++) {
        n = myBulk.FIMBulk[fn];
        for (auto& v : neighbor[n]) {
            if (myBulk.map_Bulk2FIM[v] < 0) {
                myBulk.FIMBulk.push_back(v);
                myBulk.map_Bulk2FIM[v] = myBulk.numFIMBulk + iter;  
                iter++;
            }
        }
    }
}

void BulkConn::AllocateAuxAIMt()
{
    neighborNumGacc.resize(numBulk, 0);

    for (OCP_USI n = 1; n < numBulk - 1; n++) {
        neighborNumGacc[n] += neighborNum[n - 1] - selfPtr[n - 1] - 1;
        neighborNumGacc[n + 1] += neighborNumGacc[n];
    }
    neighborNumGacc[numBulk - 1] = numConn;
}

void BulkConn::SetupMatSparsityAIMt(LinearSystem& myLS, const Bulk& myBulk) const
{
    for (USI bIde = 0; bIde < myBulk.numFIMBulk; bIde++) {
        OCP_USI bId = myBulk.FIMBulk[bIde];
        USI iter = 0;
        for (auto& v : neighbor[bId]) {
            if (myBulk.map_Bulk2FIM[v] > -1) {
                myLS.colId[bIde].push_back(myBulk.map_Bulk2FIM[v]);               
                if (v == bId) {
                    myLS.diagPtr[bIde] = iter;
                }
                iter++;
            }           
        }
    }
    myLS.dim = myBulk.numFIMBulk;
}

void BulkConn::AssembleMat_AIMt(LinearSystem& myLS, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    vector<OCP_DBL> bmat(bsize, 0);

    // Accumulation term
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {
        OCP_USI n = myBulk.FIMBulk[fn];
        bmat[0] = myBulk.rockC1 * myBulk.rockVpInit[n] - myBulk.vfp[n];
        for (USI i = 0; i < nc; i++) {
            bmat[i + 1] = -myBulk.vfi[n * nc + i];
        }
        for (USI i = 0; i < bsize; i++) {
            myLS.diagVal[fn * bsize + i] = bmat[i];
        }
    }
    // flux term
    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);

    OCP_USI bId, eId, uId;  // index in bulks
    OCP_USI minId, maxId;
    OCP_INT bIde, eIde, uIde;   // index in equations
    OCP_USI c;              // index in flux
    OCP_USI uId_np_j, uIde_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;
    OCP_BOOL flagFIM;
    USI count = 0;

    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {

        count = 0;
        bIde = fn;
        bId = myBulk.FIMBulk[fn];
        for (auto& v : neighbor[bId]) {
            if (v == bId)
                continue;
            eId = v;
            // find the index of flux: c
            minId = (bId > eId ? eId : bId);
            maxId = (bId < eId ? eId : bId);
            c = neighborNumGacc[minId];

            for (USI i = selfPtr[minId] + 1; i < neighborNum[minId]; i++) {               
                if (neighbor[minId][i] == maxId) {
                    break;
                }
                c++;
            }
            Akd = CONV1 * CONV2 * iteratorConn[c].area;
            // cout << iteratorConn[c].BId << "   " << iteratorConn[c].EId << endl;

            eIde = myBulk.map_Bulk2FIM[eId];
            if (eIde < 0)  flagFIM = OCP_FALSE;
            else           flagFIM = OCP_TRUE; 

            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);
            dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);

            for (USI j = 0; j < np; j++) {
                uId = upblock[c * np + j];
                uId_np_j = uId * np + j;
                if (!myBulk.phaseExist[uId_np_j]) continue;
                dP = myBulk.Pj[bId * np + j] - myBulk.Pj[eId * np + j] -
                    myBulk.rho[bId * np + j] * dGamma;
                xi = myBulk.xi[uId_np_j];
                kr = myBulk.kr[uId_np_j];
                mu = myBulk.mu[uId_np_j];
                transJ = Akd * kr / mu;

                uIde = (bId == uId ? bIde : eIde);
                if (uIde > -1) {
                    uIde_np_j = uIde * np + j;
                    muP = myBulk.muP[uIde_np_j];
                    xiP = myBulk.xiP[uIde_np_j];
                    rhoP = myBulk.rhoP[uIde_np_j];
                }
                
                for (USI i = 0; i < nc; i++) {
                    xij = myBulk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    if (flagFIM) {
                        dFdXpE[(i + 1) * ncol] -= transIJ;
                    }
                    if (uIde > -1) {
                        tmp = transIJ * (-rhoP * dGamma);
                        tmp += xij * transJ * xiP * dP;
                        tmp += -transIJ * muP / mu * dP;
                        if (bId == uId) {
                            dFdXpB[(i + 1) * ncol] += tmp;
                        }
                        else if (flagFIM) {
                            dFdXpE[(i + 1) * ncol] += tmp;
                        }
                    }
                    
                    // Saturation -- Second var
                    for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[bIde * np * np + j * np + k];
                        if (flagFIM) {
                            dFdXsE[(i + 1) * ncol2 + k] -=
                                transIJ * myBulk.dPcj_dS[eIde * np * np + j * np + k];
                        }
                        if (uIde > -1) {
                            tmp = Akd * xij * xi / mu *
                                myBulk.dKr_dS[uIde * np * np + j * np + k] * dP;

                            if (bId == uId) {
                                dFdXsB[(i + 1) * ncol2 + k] += tmp;
                            }
                            else if (flagFIM) {
                                dFdXsE[(i + 1) * ncol2 + k] += tmp;
                            }
                        }                                              
                    }

                    // Cij -- Second var
                    if (uIde > -1) {
                        for (USI k = 0; k < nc; k++) {
                            rhox = myBulk.rhox[uIde_np_j * nc + k];
                            xix = myBulk.xix[uIde_np_j * nc + k];
                            mux = myBulk.mux[uIde_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            if (k == i) {
                                tmp += xi * transJ * dP;
                            }
                            if (bId == uId) {
                                dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            }
                            else if (flagFIM) {
                                dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                            }
                        }
                    }                    
                }
            }
            USI diagptr = myLS.diagPtr[bIde];

            // Assemble
            // Begin
            bmat = dFdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &myBulk.dSec_dPri[bIde * bsize2], 1,
                bmat.data());
            Dscalar(bsize, dt, bmat.data());

            if (count < diagptr) {
                // Add to diagVal
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[bIde * bsize + i] += bmat[i];
                }
            }
            else if (count == diagptr) {
                OCP_USI id = bIde * bsize;
                myLS.val[bIde].insert(myLS.val[bIde].end(), myLS.diagVal.data() + id,
                    myLS.diagVal.data() + id + bsize);
                // Add to val
                for (USI i = 0; i < bsize; i++) {
                    myLS.val[bIde][diagptr * bsize + i] += bmat[i];
                }
                count++;
            }
            else {
                // Add to val
                for (USI i = 0; i < bsize; i++) {
                    myLS.val[bIde][diagptr * bsize + i] += bmat[i];
                }
            }

            // End
            if (flagFIM) {
                bmat = dFdXpE;
                DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &myBulk.dSec_dPri[eIde * bsize2], 1,
                    bmat.data());
                Dscalar(bsize, dt, bmat.data());
                // Begin
                // Insert
                myLS.val[bIde].insert(myLS.val[bIde].end(), bmat.begin(), bmat.end());
                count++;
            }
        }
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < myBulk.numFIMBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + (n + 1) * bsize);
    }
}

void BulkConn::CalResAIMt(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI len = nc + 1;
    OCP_USI   bId, eId, uId, bIdc;
    OCP_INT   bIde;
    
    // Accumalation Term
    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {

        OCP_USI n = myBulk.FIMBulk[fn];
        bIde = fn * len;
        bId = n * len;
        bIdc = n * nc;

        res[bIde] = myBulk.rockVp[n] - myBulk.vf[n];
        for (USI i = 0; i < nc; i++) {
            res[bIde + 1 + i] = myBulk.Ni[bIdc + i] - myBulk.lNi[bIdc + i];
        }
    }

    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_USI maxId, minId, c;
    OCP_DBL Pbegin, Pend, rho, dP;
    OCP_DBL tmp, dNi;
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.

    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {
        bIde = fn;
        bId = myBulk.FIMBulk[fn];
        for (auto& v : neighbor[bId]) {
            if (v == bId)
                continue;
            eId = v;

            // find the index of flux: c
            minId = (bId > eId ? eId : bId);
            maxId = (bId < eId ? eId : bId);
            c = neighborNumGacc[minId];

            for (USI i = selfPtr[minId] + 1; i < neighborNum[minId]; i++) {
                if (neighbor[minId][i] == maxId) {
                    break;
                }
                c++;
            }
            Akd = CONV1 * CONV2 * iteratorConn[c].area;

            for (USI j = 0; j < np; j++) {
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;

                OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
                OCP_BOOL exend = myBulk.phaseExist[eId_np_j];

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
                    // upblock_Rho[c * np + j] = rho;
                    continue;
                }

                uId = bId;
                dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                    (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
                if (dP < 0) {
                    uId = eId;
                }
                // upblock_Rho[c * np + j] = rho;
                upblock[c * np + j] = uId;

                uId_np_j = uId * np + j;
                if (!myBulk.phaseExist[uId_np_j]) continue;
                tmp = dt * Akd * myBulk.xi[uId_np_j] *
                    myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
                for (USI i = 0; i < nc; i++) {
                    dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                    res[bIde * len + 1 + i] += dNi;
                }
            }
        }
    }
}

void BulkConn::CalResAIMs(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI len = nc + 1;
    OCP_USI   bId, eId, uId, bIdc;

    // Accumalation Term
    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {

        OCP_USI n = myBulk.FIMBulk[fn];
        bId = n * len;
        bIdc = n * nc;

        res[bId] = myBulk.rockVp[n] - myBulk.vf[n];
        for (USI i = 0; i < nc; i++) {
            res[bId + 1 + i] = myBulk.Ni[bIdc + i] - myBulk.lNi[bIdc + i];
        }
    }

    OCP_USI bId_np_j, eId_np_j, uId_np_j;
    OCP_USI maxId, minId, c;
    OCP_DBL Pbegin, Pend, rho, dP;
    OCP_DBL tmp, dNi;
    OCP_DBL Akd;
    // Flux Term
    // Calculate the upblock at the same time.

    for (OCP_USI fn = 0; fn < myBulk.numFIMBulk; fn++) {
        bId = myBulk.FIMBulk[fn];
        for (auto& v : neighbor[bId]) {
            if (v == bId)
                continue;
            eId = v;

            // find the index of flux: c
            minId = (bId > eId ? eId : bId);
            maxId = (bId < eId ? eId : bId);
            c = neighborNumGacc[minId];

            for (USI i = selfPtr[minId] + 1; i < neighborNum[minId]; i++) {
                if (neighbor[minId][i] == maxId) {
                    break;
                }
                c++;
            }
            Akd = CONV1 * CONV2 * iteratorConn[c].area;

            for (USI j = 0; j < np; j++) {
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;

                OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
                OCP_BOOL exend = myBulk.phaseExist[eId_np_j];

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
                    continue;
                }

                uId = bId;
                dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
                    (Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
                if (dP < 0) {
                    uId = eId;
                }

                uId_np_j = uId * np + j;
                if (!myBulk.phaseExist[uId_np_j]) continue;
                tmp = dt * Akd * myBulk.xi[uId_np_j] *
                    myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
                for (USI i = 0; i < nc; i++) {
                    dNi = tmp * myBulk.xij[uId_np_j * nc + i];
                    res[bId * len + 1 + i] += dNi;
                }
            }
        }
    }
}

void BulkConn::AssembleMat_AIMs(LinearSystem& myLS, vector<OCP_DBL>& res, const Bulk& myBulk,
    const OCP_DBL& dt) const
{
    // accumulate term
    OCP_DBL Vp0, Vp, vf, vfp;
    OCP_DBL cr = myBulk.rockC1;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    vector<OCP_DBL> bmat(bsize, 0);
    for (USI i = 1; i < nc + 1; i++) {
        bmat[i * ncol + i] = 1;
    }
    vector<OCP_DBL> bmatTmp(bmat);

    for (OCP_USI n = 0; n < numBulk; n++) {
        Vp0 = myBulk.rockVpInit[n];
        Vp = myBulk.rockVp[n];
        vfp = myBulk.vfp[n];
        
        if (myBulk.map_Bulk2FIM[n] < 0 || myBulk.map_Bulk2FIM[n] >= myBulk.numFIMBulk) {
            // IMPEC bulk
            bmat[0] = cr * Vp0 - vfp;
            vf = myBulk.vf[n];
            // Method 1
            res[n * ncol] = bmat[0] * (myBulk.lP[n] - myBulk.P[n]) + dt * (vf - Vp);
            // Method 2
            // res[n * ncol] = bmat[0] * myBulk.P[n] + dt * (vf - Vp);
            for (USI i = 0; i < bsize; i++) {
                myLS.diagVal[n * bsize + i] = bmat[i];
            }
        }
        else {
            // FIM Bulk
            bmatTmp[0] = cr * Vp0 - vfp;
            for (USI i = 0; i < nc; i++) {
                bmatTmp[i + 1] = -myBulk.vfi[n * nc + i];
            }
            for (USI i = 0; i < bsize; i++) {
                myLS.diagVal[n * bsize + i] = bmatTmp[i];
            }
        }
    }

    // flux term
    OCP_USI bId, eId, uId;
    OCP_INT bIde, eIde, uIde;
    OCP_DBL valupi, valdowni;
    OCP_DBL valup, rhsup, valdown, rhsdown;
    OCP_USI uId_np_j, uIde_np_j;
    OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
    OCP_DBL dP, dGamma;
    OCP_DBL tmp;
    OCP_BOOL bIdFIM, eIdFIM;
    OCP_BOOL otherFIM;
    OCP_USI FIMbId, FIMeId, FIMbIde, FIMeIde;
    OCP_USI IMPECbId, IMPECeId;
    USI diagptr;

    OCP_DBL         Akd;
    OCP_DBL         transJ, transIJ;
    vector<OCP_DBL> dFdXpB(bsize, 0);
    vector<OCP_DBL> dFdXpE(bsize, 0);
    vector<OCP_DBL> dFdXsB(bsize2, 0);
    vector<OCP_DBL> dFdXsE(bsize2, 0);
    vector<OCP_DBL> IMPECbmat(bsize, 0);

    // Be careful when first bulk has no neighbors!
    OCP_USI lastbId = iteratorConn[0].EId;
    for (OCP_USI c = 0; c < numConn; c++) {
        bId = iteratorConn[c].BId;
        eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;

        bIde = myBulk.map_Bulk2FIM[bId];
        eIde = myBulk.map_Bulk2FIM[eId];

        bIdFIM = eIdFIM = OCP_FALSE;
        if (bIde > -1 && bIde < myBulk.numFIMBulk)  bIdFIM = OCP_TRUE;
        if (eIde > -1 && eIde < myBulk.numFIMBulk)  eIdFIM = OCP_TRUE;

        if (bIdFIM || eIdFIM) {
            // There exist at least one FIM bulk
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            FIMbId = bId;
            FIMeId = eId;
            FIMbIde = bIde;
            FIMeIde = eIde;
            otherFIM = eIdFIM;
            if (!bIdFIM) {
                FIMbId = eId;
                FIMeId = bId;
                FIMbIde = eIde;
                FIMeIde = bIde;
                otherFIM = bIdFIM;
            }
            dGamma = GRAVITY_FACTOR * (myBulk.depth[FIMbId] - myBulk.depth[FIMeId]);
            for (USI j = 0; j < np; j++) {
                uId = upblock[c * np + j];
                uId_np_j = uId * np + j;
                uIde = myBulk.map_Bulk2FIM[uId];
                uIde_np_j = uIde * np + j;

                if (!myBulk.phaseExist[uId_np_j]) continue;
                dP = myBulk.Pj[FIMbId * np + j] - myBulk.Pj[FIMeId * np + j] -
                    myBulk.rho[FIMbId * np + j] * dGamma;
                xi = myBulk.xi[uId_np_j];
                kr = myBulk.kr[uId_np_j];
                mu = myBulk.mu[uId_np_j];
                muP = myBulk.muP[uIde_np_j];
                xiP = myBulk.xiP[uIde_np_j];
                rhoP = myBulk.rhoP[uIde_np_j];
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
                    if (FIMbId == uId) {
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }

                    // Saturation -- Second var
                    for (USI k = 0; k < np; k++) {
                        dFdXsB[(i + 1) * ncol2 + k] +=
                            transIJ * myBulk.dPcj_dS[FIMbIde * np * np + j * np + k];
                        dFdXsE[(i + 1) * ncol2 + k] -=
                            transIJ * myBulk.dPcj_dS[FIMeIde * np * np + j * np + k];
                        tmp = Akd * xij * xi / mu *
                            myBulk.dKr_dS[uIde * np * np + j * np + k] * dP;
                        if (FIMbId == uId) {
                            dFdXsB[(i + 1) * ncol2 + k] += tmp;
                        }
                        else {
                            dFdXsE[(i + 1) * ncol2 + k] += tmp;
                        }
                    }
                    // Cij -- Second var
                    for (USI k = 0; k < nc; k++) {
                        rhox = myBulk.rhox[uIde_np_j * nc + k];
                        xix = myBulk.xix[uIde_np_j * nc + k];
                        mux = myBulk.mux[uIde_np_j * nc + k];
                        tmp = -transIJ * rhox * dGamma;
                        tmp += xij * transJ * xix * dP;
                        tmp += -transIJ * mux / mu * dP;
                        if (k == i) {
                            tmp += xi * transJ * dP;
                        }
                        if (FIMbId == uId) {
                            dFdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                        else {
                            dFdXsE[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                        }
                    }
                }
            }
            USI diagptr = myLS.diagPtr[bId];

            // insert diag
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
            DaABpbC(ncol, ncol, ncol2, 1, dFdXsB.data(), &myBulk.dSec_dPri[FIMbIde * bsize2], 1,
                bmat.data());
            Dscalar(bsize, dt, bmat.data());
            // Begin
            // Add
            if (FIMbId < FIMeId) {
                diagptr = myLS.diagPtr[FIMbId];
                for (USI i = 0; i < bsize; i++) {
                    myLS.val[FIMbId][diagptr * bsize + i] += bmat[i];
                }
            }
            else {
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[FIMbId * bsize + i] += bmat[i];
                }
            }
            
            
            if (otherFIM) {
                // End
                // Insert
                Dscalar(bsize, -1, bmat.data());
                myLS.val[FIMeId].insert(myLS.val[FIMeId].end(), bmat.begin(), bmat.end());                
            }
            
            // End
            bmat = dFdXpE;
            DaABpbC(ncol, ncol, ncol2, 1, dFdXsE.data(), &myBulk.dSec_dPri[FIMeIde * bsize2], 1,
                bmat.data());
            Dscalar(bsize, dt, bmat.data());
            // Begin
            // Insert
            if (!otherFIM) {
                // delete der about Ni
                for (USI i = 0; i < ncol; i++) {
                    for (USI j = 1; j < ncol; j++) {
                        bmat[i * ncol + j] = 0;
                    }
                }
            }
            myLS.val[FIMbId].insert(myLS.val[FIMbId].end(), bmat.begin(), bmat.end());
            if (otherFIM) {
                // Add
                Dscalar(bsize, -1, bmat.data());
                if (FIMbId < FIMeId) {
                    for (USI i = 0; i < bsize; i++) {
                        myLS.diagVal[FIMeId * bsize + i] += bmat[i];
                    }
                }
                else {
                    diagptr = myLS.diagPtr[FIMeId];
                    for (USI i = 0; i < bsize; i++) {
                        myLS.val[FIMeId][diagptr * bsize + i] += bmat[i];
                    }
                }               
            }
        }

        if (!bIdFIM || !eIdFIM) {
            // There exist at least one IMPEC bulk
            valup = 0;
            rhsup = 0;
            valdown = 0;
            rhsdown = 0;

            IMPECbId = bId;
            IMPECeId = eId;
            otherFIM = eIdFIM;
            if (bIdFIM) {
                IMPECbId = eId;
                IMPECeId = bId;
                otherFIM = bIdFIM;
            }

            dGamma = GRAVITY_FACTOR* (myBulk.depth[IMPECbId] - myBulk.depth[IMPECeId]);
            for (USI j = 0; j < np; j++) {
                uId = upblock[c * np + j];
                if (!myBulk.phaseExist[uId * np + j]) continue;

                valupi = 0;
                valdowni = 0;

                for (USI i = 0; i < nc; i++) {
                    valupi +=
                        myBulk.vfi[IMPECbId * nc + i] * myBulk.xij[uId * np * nc + j * nc + i];
                    valdowni +=
                        myBulk.vfi[IMPECeId * nc + i] * myBulk.xij[uId * np * nc + j * nc + i];
                }
                OCP_DBL dPc = myBulk.Pc[IMPECbId * np + j] - myBulk.Pc[IMPECeId * np + j];
                dP = myBulk.Pj[IMPECbId * np + j] - myBulk.Pj[IMPECeId * np + j] -
                    upblock_Rho[c * np + j] * dGamma;
                OCP_DBL temp = myBulk.xi[uId * np + j] * upblock_Trans[c * np + j] * dt;
                
                valup += temp * valupi;
                valdown += temp * valdowni;
                // Method 1
                temp *= (-dP);
                // Method 2
                // temp *= upblock_Rho[c * np + j] * dGamma - dPc;
                rhsup += temp * valupi;
                rhsdown -= temp * valdowni;
            }

            // Add diag
            USI diagptr = myLS.diagPtr[bId];
            if (bId != lastbId) {
                // new bulk
                assert(myLS.val[bId].size() == diagptr * bsize);
                OCP_USI id = bId * bsize;
                myLS.val[bId].insert(myLS.val[bId].end(), myLS.diagVal.data() + id,
                    myLS.diagVal.data() + id + bsize);
                lastbId = bId;
            }

            // Begin
            if (IMPECbId < IMPECeId) {
                diagptr = myLS.diagPtr[IMPECbId];
                myLS.val[IMPECbId][diagptr * bsize] += valup;
            }
            else {
                myLS.diagVal[IMPECbId * bsize] += valup;
            }
            
            IMPECbmat[0] = (-valup);
            myLS.val[IMPECbId].insert(myLS.val[IMPECbId].end(), IMPECbmat.begin(), IMPECbmat.end());

            // End
            if (!otherFIM) {
                IMPECbmat[0] = (-valdown);
                myLS.val[IMPECeId].insert(myLS.val[IMPECeId].end(), IMPECbmat.begin(), IMPECbmat.end());
                
                if (IMPECbId < IMPECeId) {
                    myLS.diagVal[IMPECeId * bsize] += valdown;
                }
                else {
                    diagptr = myLS.diagPtr[IMPECeId];
                    myLS.val[IMPECeId][diagptr * bsize] += valup;
                }                
            }            


            // rhs
            res[IMPECbId * ncol] += rhsup;
            if (!otherFIM) {
                res[IMPECeId * ncol] += rhsdown;
            }
            
        }
    }
    // Add the rest of diag value. Important!
    for (OCP_USI n = 0; n < numBulk; n++) {
        if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
            myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
                myLS.diagVal.data() + n * bsize + bsize);
    }
}


void BulkConn::AllocateAuxAIMc(const USI& np)
{
    OCP_FUNCNAME;

    upblock.resize(numConn * np);
    upblock_Rho.resize(numConn * np);
    upblock_Velocity.resize(numConn * np);
}


void BulkConn::AssembleMat_AIMc(LinearSystem& myLS, const Bulk& myBulk, const OCP_DBL& dt) const
{
    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

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
    OCP_DBL* dFdXpI;
    OCP_DBL* dFdXsI;
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistE(np, OCP_FALSE);
    vector<OCP_BOOL>    phaseExistI(np, OCP_FALSE);
    OCP_BOOL            phaseExistU;
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_E(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_I(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<USI>     pVnumComE(np, 0);
    vector<USI>     pVnumComI(np, 0);
    OCP_USI     ncolB, ncolE, ncolI;


    OCP_USI bId, eId, uId;
    OCP_USI bId_np_j, eId_np_j, uId_np_j, ibId_np_j;
	OCP_DBL kr, mu, xi, xij, rhoP, xiP, muP, rhox, xix, mux;
	OCP_DBL dP, dGamma;
	OCP_DBL tmp;
    OCP_BOOL bIdFIM, eIdFIM, uIdFIM;
	USI diagptr;

	// Be careful when first bulk has no neighbors!
	OCP_USI lastbId = iteratorConn[0].EId;
	for (OCP_USI c = 0; c < numConn; c++) {
		bId = iteratorConn[c].BId;
		eId = iteratorConn[c].EId;
        Akd = CONV1 * CONV2 * iteratorConn[c].area;
        dGamma = GRAVITY_FACTOR * (myBulk.depth[bId] - myBulk.depth[eId]);
        bIdFIM = eIdFIM = OCP_FALSE;
        if (myBulk.map_Bulk2FIM[bId] > -1) bIdFIM = OCP_TRUE;
        if (myBulk.map_Bulk2FIM[eId] > -1) eIdFIM = OCP_TRUE;

        if (!bIdFIM && !eIdFIM) {
            // both are explicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            for (USI j = 0; j < np; j++) {               
                phaseExistB[j] = myBulk.phaseExist[bId * np + j];
                phaseExistE[j] = myBulk.phaseExist[eId * np + j];
                uId = upblock[c * np + j];
                phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);                 
                if (!phaseExistU) {
                    continue;
                }               
                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                xi = myBulk.xi[uId_np_j];
                kr = myBulk.kr[uId_np_j];
                mu = myBulk.mu[uId_np_j];
                transJ = Akd * xi * kr / mu;
                for (USI i = 0; i < nc; i++) {
                    xij = myBulk.xij[uId_np_j * nc + i];
                    transIJ = xij * transJ;
                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;
                }
            }
        }
        else if ((bIdFIM && !eIdFIM) || (!bIdFIM && eIdFIM)) {

            // one is explicit, one is implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            ncolI = 0;
            USI jxI = 0;
            if (bIdFIM) {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = myBulk.phaseExist[bId * np + j];
                    phasedS_I[j] = myBulk.pSderExist[bId * np + j];
                    if (phasedS_I[j]) jxI++;
                    pVnumComI[j] = myBulk.pVnumCom[bId * np + j];
                    ncolI += pVnumComI[j];
                }
                dFdXpI = &dFdXpB[0];
                dFdXsI = &dFdXsB[0];
            }
            else {
                for (USI j = 0; j < np; j++) {
                    phaseExistI[j] = myBulk.phaseExist[eId * np + j];
                    phasedS_I[j] = myBulk.pSderExist[eId * np + j];
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
                uId = upblock[c * np + j];
                phaseExistU = (uId == bId ? phaseExistB[j] : phaseExistE[j]);
                if (!phaseExistU) {
                    continue;
                }

                bId_np_j = bId * np + j;
                eId_np_j = eId * np + j;
                uId_np_j = uId * np + j;
                dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                    upblock_Rho[c * np + j] * dGamma;
                xi = myBulk.xi[uId_np_j];
                kr = myBulk.kr[uId_np_j];
                mu = myBulk.mu[uId_np_j];
                transJ = Akd * kr / mu;
                uIdFIM = (uId == bId ? bIdFIM : eIdFIM);

                if (bIdFIM) ibId_np_j = bId_np_j;
                else        ibId_np_j = eId_np_j;

                if (uIdFIM) {
                    // upblock is implicit
                    OCP_DBL rhoWgt = 1.0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;

                    muP = myBulk.muP[ibId_np_j];
                    xiP = myBulk.xiP[ibId_np_j];
                    rhoP = myBulk.rhoP[ibId_np_j];
                    for (USI i = 0; i < nc; i++) {
                        xij = myBulk.xij[ibId_np_j * nc + i];
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
                            xix = myBulk.xix[ibId_np_j * nc + k];
                            mux = myBulk.mux[ibId_np_j * nc + k];
                            tmp = -transIJ * rhox * dGamma;
                            tmp += xij * transJ * xix * dP;
                            tmp += -transIJ * mux / mu * dP;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                        // WARNING !!!
                        if (i < pVnumComI[j])
                            dFdXsI[(i + 1) * ncolI + jxI + i] += xi * transJ * dP;
                    }
                    jxI += pVnumComI[j];
                }
                else {
                    // downblock is implicit
                    OCP_DBL rhoWgt = 0;
                    if (phaseExistB[j] && phaseExistE[j]) rhoWgt = 0.5;
                    rhoP = myBulk.rhoP[ibId_np_j];

                    for (USI i = 0; i < nc; i++) {
                        xij = myBulk.xij[uId_np_j * nc + i];
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
                            tmp = -transIJ * rhox * dGamma;
                            dFdXsI[(i + 1) * ncolI + jxI + k] += tmp;
                        }
                    }
                    jxI += pVnumComI[j];
                }
            }
        }
        else {
            // both are implicit
            fill(dFdXpB.begin(), dFdXpB.end(), 0.0);
            fill(dFdXpE.begin(), dFdXpE.end(), 0.0);
            fill(dFdXsB.begin(), dFdXsB.end(), 0.0);
            fill(dFdXsE.begin(), dFdXsE.end(), 0.0);

            USI jxB = 0;        USI jxE = 0;
            ncolB = 0;          ncolE = 0;

            for (USI j = 0; j < np; j++) {
                phaseExistB[j] = myBulk.phaseExist[bId * np + j];
                phaseExistE[j] = myBulk.phaseExist[eId * np + j];
                phasedS_B[j] = myBulk.pSderExist[bId * np + j];
                phasedS_E[j] = myBulk.pSderExist[eId * np + j];
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
                dP = myBulk.Pj[bId_np_j] - myBulk.Pj[eId_np_j] -
                    upblock_Rho[c * np + j] * dGamma;
                xi = myBulk.xi[uId_np_j];
                kr = myBulk.kr[uId_np_j];
                mu = myBulk.mu[uId_np_j];
                muP = myBulk.muP[uId_np_j];
                xiP = myBulk.xiP[uId_np_j];
                rhoP = myBulk.rhoP[uId_np_j];
                transJ = Akd * kr / mu;


                for (USI i = 0; i < nc; i++) {
                    xij = myBulk.xij[uId_np_j * nc + i];
                    transIJ = xij * xi * transJ;

                    // Pressure -- Primary var
                    dFdXpB[(i + 1) * ncol] += transIJ;
                    dFdXpE[(i + 1) * ncol] -= transIJ;

                    tmp = xij * transJ * xiP * dP;
                    tmp += -transIJ * muP / mu * dP;
                    if (!phaseExistE[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpB[(i + 1) * ncol] += tmp;
                    }
                    else if (!phaseExistB[j]) {
                        tmp += transIJ * (-rhoP * dGamma);
                        dFdXpE[(i + 1) * ncol] += tmp;
                    }
                    else {
                        dFdXpB[(i + 1) * ncol] +=
                            transIJ * (-myBulk.rhoP[bId_np_j] * dGamma) / 2;
                        dFdXpE[(i + 1) * ncol] +=
                            transIJ * (-myBulk.rhoP[eId_np_j] * dGamma) / 2;
                        if (bId == uId) {
                            dFdXpB[(i + 1) * ncol] += tmp;
                        }
                        else {
                            dFdXpE[(i + 1) * ncol] += tmp;
                        }
                    }

                    // Second var
                    USI j1SB = 0; USI j1SE = 0;
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
                                xix = myBulk.xix[uId_np_j * nc + k];
                                mux = myBulk.mux[uId_np_j * nc + k];
                                tmp = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        }
                        else {
                            for (USI k = 0; k < pVnumComB[j]; k++) {
                                rhox = myBulk.rhox[bId_np_j * nc + k] / 2;
                                xix = myBulk.xix[uId_np_j * nc + k];
                                mux = myBulk.mux[uId_np_j * nc + k];
                                tmp = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += -transIJ * myBulk.rhox[eId_np_j * nc + k] / 2 * dGamma;
                            }
                            // WARNING !!!
                            if (i < pVnumComB[j])
                                dFdXsB[(i + 1) * ncolB + jxB + i] += xi * transJ * dP;
                        }
                    }
                    else {
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
                                xix = myBulk.xix[uId_np_j * nc + k];
                                mux = myBulk.mux[uId_np_j * nc + k];
                                tmp = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                            }
                            // WARNING !!!
                            if (i < pVnumComE[j])
                                dFdXsE[(i + 1) * ncolE + jxE + i] += xi * transJ * dP;
                        }
                        else {
                            for (USI k = 0; k < pVnumComE[j]; k++) {
                                rhox = myBulk.rhox[eId_np_j * nc + k] / 2;
                                xix = myBulk.xix[uId_np_j * nc + k];
                                mux = myBulk.mux[uId_np_j * nc + k];
                                tmp = -transIJ * rhox * dGamma;
                                tmp += xij * transJ * xix * dP;
                                tmp += -transIJ * mux / mu * dP;
                                dFdXsE[(i + 1) * ncolE + jxE + k] += tmp;
                                dFdXsB[(i + 1) * ncolB + jxB + k] += -transIJ * myBulk.rhox[bId_np_j * nc + k] / 2 * dGamma;
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

        USI diagptr = myLS.diagPtr[bId];
        if (bId != lastbId) {
            // new bulk
            assert(myLS.val[bId].size() == diagptr * bsize);
            OCP_USI id = bId * bsize;
            myLS.val[bId].insert(myLS.val[bId].end(), myLS.diagVal.data() + id,
                myLS.diagVal.data() + id + bsize);

            lastbId = bId;
        }

        if (bIdFIM && !eIdFIM) ncolB = ncolI;
        else if (!bIdFIM && eIdFIM) ncolE = ncolI;

        // Assemble
        bmat = dFdXpB;
        if (bIdFIM) {
            DaABpbC(ncol, ncol, ncolB, 1, dFdXsB.data(), &myBulk.dSec_dPri[bId * lendSdP], 1,
                bmat.data());
        }
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

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

        // End
        bmat = dFdXpE;
        if (eIdFIM) {
            DaABpbC(ncol, ncol, ncolE, 1, dFdXsE.data(), &myBulk.dSec_dPri[eId * lendSdP], 1,
                bmat.data());
        }       
        Dscalar(bsize, dt, bmat.data());
        // Begin
        // Insert
        myLS.val[bId].insert(myLS.val[bId].end(), bmat.begin(), bmat.end());
        // Add
        Dscalar(bsize, -1, bmat.data());
        for (USI i = 0; i < bsize; i++) {
            myLS.diagVal[eId * bsize + i] += bmat[i];
        }

#ifdef OCP_NANCHECK
        if (!CheckNan(bmat.size(), &bmat[0]))
        {
            OCP_ABORT("INF or NAN in bmat !");
        }
#endif

	}
	// Add the rest of diag value. Important!
	for (OCP_USI n = 0; n < numBulk; n++) {
		if (myLS.val[n].size() == myLS.diagPtr[n] * bsize)
			myLS.val[n].insert(myLS.val[n].end(), myLS.diagVal.data() + n * bsize,
				myLS.diagVal.data() + n * bsize + bsize);
	}
}

void BulkConn::CalResAIMc(vector<OCP_DBL>& res, const Bulk& myBulk, const OCP_DBL& dt)
{
    // IMPORTANT!!!
    // in AIMc for IMPEC Bulk, P was updated in each Newton Step, but Pj didn't.
    // So here Pj must be replaced with P + Pcj, otherwise wrong results will arise

	OCP_FUNCNAME;

	const USI np = myBulk.numPhase;
	const USI nc = myBulk.numCom;
	const USI len = nc + 1;
	OCP_USI   bId, eId, uId, bIdb;

	// Accumalation Term
	for (OCP_USI n = 0; n < numBulk; n++) {

		bId = n * len;
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
		bId = iteratorConn[c].BId;
		eId = iteratorConn[c].EId;
		Akd = CONV1 * CONV2 * iteratorConn[c].area;

		for (USI j = 0; j < np; j++) {
			bId_np_j = bId * np + j;
			eId_np_j = eId * np + j;

			OCP_BOOL exbegin = myBulk.phaseExist[bId_np_j];
			OCP_BOOL exend = myBulk.phaseExist[eId_np_j];

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
                upblock_Velocity[c * np + j] = 0;
				upblock_Rho[c * np + j] = 0;
				continue;
			}

			uId = bId;
            OCP_BOOL    exup = exbegin;
			dP = (Pbegin - GRAVITY_FACTOR * rho * myBulk.depth[bId]) -
				(Pend - GRAVITY_FACTOR * rho * myBulk.depth[eId]);
			if (dP < 0) {
				uId = eId;
                exup = exend;
			}
            uId_np_j = uId * np + j;

			upblock_Rho[c * np + j] = rho;
			upblock[c * np + j] = uId;

            if (exup) {
                upblock_Velocity[c * np + j] = Akd * myBulk.kr[uId_np_j] / myBulk.mu[uId_np_j] * dP;
            }
            else { 
                upblock_Velocity[c * np + j] = 0; 
                continue;
            }
			
            // in AIMc
			// tmp = dt * upblock_Velocity[c * np + j] * myBulk.xi[uId_np_j];
            // in FIM
            tmp = dt * Akd * myBulk.xi[uId_np_j] *
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