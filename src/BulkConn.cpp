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

void BulkConn::SetupIsoT(const Grid& myGrid, const Bulk& myBulk)
{
    Setup(myGrid);
    CalAkd(myBulk);
}

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
        bId   = iteratorConn[c].bId;
        eId   = iteratorConn[c].eId;
        areaB = iteratorConn[c].areaB;
        areaE = iteratorConn[c].areaE;
        switch (iteratorConn[c].direction) {
            case 1:
                T1 = myBulk.ntg[bId] * myBulk.rockKx[bId] * areaB;
                T2 = myBulk.ntg[eId] * myBulk.rockKx[eId] * areaE;
                break;
            case 2:
                T1 = myBulk.ntg[bId] * myBulk.rockKy[bId] * areaB;
                T2 = myBulk.ntg[eId] * myBulk.rockKy[eId] * areaE;
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

void BulkConn::PrintConnectionInfo(const Grid& myGrid) const
{
    for (OCP_USI i = 0; i < numBulk; i++) {
        cout << "(" << myGrid.map_Act2All[i] << ")"
             << "\t";

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

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/