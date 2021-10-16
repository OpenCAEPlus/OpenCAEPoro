/*! \file    Well.cpp
 *  \brief   Well class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "Well.hpp"
#include <cmath>

WellOpt::WellOpt(const WellOptParam& Optparam)
{
    if (Optparam.type == "INJ") {
        type = INJ;
    } else if (Optparam.type == "PROD") {
        type = PROD;
    } else {
        ERRORcheck("WRONG Well Type");
        exit(0);
    }

    if (type == INJ) {
        if (Optparam.fluidType == "OIL") {
            fluidType = OIL;
        } else if (Optparam.fluidType == "GAS") {
            fluidType = GAS;
        } else if (Optparam.fluidType == "WATER") {
            fluidType = WATER;
        } else if (Optparam.fluidType == "SOLVENT") {
            fluidType = SOLVENT;
        } else {
            ERRORcheck("WRONG Fluid type");
            exit(0);
        }
    }

    if (Optparam.state == "OPEN") {
        state = OPEN;
    } else if (Optparam.state == "CLOSE") {
        state = CLOSE;
    } else {
        ERRORcheck("WRONG State type");
        exit(0);
    }

    if (Optparam.optMode == "RATE") {
        optMode = RATE_MODE;
    } else if (Optparam.optMode == "ORAT") {
        optMode = ORATE_MODE;
    } else if (Optparam.optMode == "GRAT") {
        optMode = GRATE_MODE;
    } else if (Optparam.optMode == "WRAT") {
        optMode = WRATE_MODE;
    } else if (Optparam.optMode == "BHP") {
        optMode = BHP_MODE;
    } else {
        ERRORcheck("WRONG Well Opt Mode");
        exit(0);
    }

    maxRate = Optparam.maxRate;
    maxBHP  = Optparam.maxBHP;
    minBHP  = Optparam.minBHP;
}

void Well::InputPerfo(const WellParam& well)
{
    numPerf = well.I_perf.size();
    perf.resize(numPerf);
    for (USI p = 0; p < numPerf; p++) {
        perf[p].I          = well.I_perf[p] - 1;
        perf[p].J          = well.J_perf[p] - 1;
        perf[p].K          = well.K_perf[p] - 1;
        perf[p].WI         = well.WI[p];
        perf[p].radius     = well.diameter[p] / 2.0;
        perf[p].kh         = well.kh[p];
        perf[p].skinFactor = well.skinFactor[p];
        if (well.direction[p] == "X" || well.direction[p] == "x") {
            perf[p].direction = X_DIRECTION;
        } else if (well.direction[p] == "Y" || well.direction[p] == "y") {
            perf[p].direction = Y_DIRECTION;
        } else if (well.direction[p] == "Z" || well.direction[p] == "z") {
            perf[p].direction = Z_DIRECTION;
        } else {
            ERRORcheck("Wrong direction of perforations");
            exit(0);
        }
    }
}

void Well::Setup(const Grid& myGrid, const Bulk& myBulk)
{
    // zi
    if (myBulk.blackOil) {
        for (auto& opt : optSet) {

            opt.zi.resize(myBulk.numCom, 0);
            if (opt.type == INJ) {
                // INJ
                switch (myBulk.PVTmode) {
                    case PHASE_W:
                    case PHASE_OW:
                        opt.zi.back() = 1;
                        break;
                    case PHASE_OGW:
                        if (opt.fluidType == GAS)
                            opt.zi[1] = 1;
                        else
                            opt.zi[2] = 1;
                        break;
                    default:
                        ERRORcheck("WRONG Blackoil type!");
                        exit(0);
                }
            } else {
                // PROD
                switch (myBulk.PVTmode) {
                    case PHASE_W:
                        opt.zi.back() = 1;
                        break;
                    case PHASE_OW:
                        if (opt.optMode == ORATE_MODE)
                            opt.zi[0] = 1;
                        else
                            opt.zi[1] = 1;
                        break;
                    case PHASE_OGW:
                        if (opt.optMode == ORATE_MODE)
                            opt.zi[0] = 1;
                        else if (opt.optMode == GRATE_MODE)
                            opt.zi[1] = 1;
                        else
                            opt.zi[2] = 1;
                        break;
                    default:
                        ERRORcheck("WRONG Blackoil type!");
                        exit(0);
                }
            }
        }
    } else if (myBulk.comps) {

    } else {
        ERRORcheck("Wrong Mixture Type !");
        exit(0);
    }

    qi_lbmol.resize(myBulk.numCom);
    // perf
    dG.resize(numPerf, 0);
    ldG = dG;
    for (USI p = 0; p < numPerf; p++) {
        perf[p].state = OPEN;
        OCP_USI Idg =
            perf[p].K * myGrid.nx * myGrid.ny + perf[p].J * myGrid.nx + perf[p].I;
        if (!myGrid.activeMap_G2B[Idg].GetAct()) {
            ERRORcheck("Perforation is in inactive bulk !");
            exit(0);
        }
        perf[p].location   = myGrid.activeMap_G2B[Idg].GetId();
        perf[p].depth      = myBulk.depth[perf[p].location];
        perf[p].multiplier = 1;
        perf[p].qi_lbmol.resize(myBulk.numCom);
        perf[p].transj.resize(myBulk.numPhase);
    }
    if (depth < 0) depth = perf[0].depth;

    CalWI_Peaceman_Vertical(myBulk);
    cout << "Well::Setup" << endl;
}

void Well::Init(const Bulk& myBulk) { BHP = myBulk.P[perf[0].location]; }

OCP_DBL Well::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_DBL cfl = 0;
    OCP_DBL tmp = 0;
    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;
        tmp       = fabs(perf[p].qt_ft3) * dt;
        tmp /= myBulk.rockVp[k];
        if (cfl < tmp) cfl = tmp;
    }
    return cfl;
}

void Well::CalWI_Peaceman_Vertical(const Bulk& myBulk)
{
    // this fomular needs to be carefully checked !
    // especially the dz

    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].WI > 0) {
            break;
        } else {
            OCP_USI Idb = perf[p].location;
            OCP_DBL dx  = myBulk.dx[Idb];
            OCP_DBL dy  = myBulk.dy[Idb];
            OCP_DBL dz  = myBulk.dz[Idb] * myBulk.ntg[Idb];
            switch (perf[p].direction) {
                case X_DIRECTION:
                    break;
                case Y_DIRECTION:
                    break;
                case Z_DIRECTION: {
                    OCP_DBL kxky  = myBulk.rockKx[Idb] * myBulk.rockKy[Idb];
                    OCP_DBL kx_ky = myBulk.rockKx[Idb] / myBulk.rockKy[Idb];
                    assert(kx_ky > 0);

                    OCP_DBL ro =
                        0.28 *
                        pow((dx * dx * pow(1 / kx_ky, 0.5) + dy * dy * pow(kx_ky, 0.5)),
                            0.5);
                    ro /= (pow(kx_ky, 0.25) + pow(1 / kx_ky, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].WI = (2 * PI) * (dz * pow(kxky, 0.5)) /
                                     (log(ro / perf[p].radius) + perf[p].skinFactor);
                    } else {
                        perf[p].WI = (2 * PI) * perf[p].kh /
                                     (log(ro / perf[p].radius) + perf[p].skinFactor);
                    }
                    break;
                }
                default:
                    ERRORcheck("Wrong direction of perforations!");
                    exit(0);
            }
        }
    }
}

void Well::AssembleMat_INJ_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver,
                                 const OCP_DBL& dt) const
{
    USI     nc  = myBulk.numCom;
    OCP_USI wId = mySolver.dim;
    // important !
    mySolver.dim++;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;

        OCP_DBL Vfi_zi = 0;
        for (USI i = 0; i < nc; i++) {
            Vfi_zi += myBulk.vfi[k * nc + i] * opt.zi[i];
        }

        USI pvtnum = myBulk.PVTNUM[k];
        perf[p].xi =
            myBulk.flashCal[pvtnum]->XiPhase(myBulk.P[k], myBulk.T, &opt.zi[0]);
        OCP_DBL valw = dt * perf[p].xi * perf[p].transj[0];
        OCP_DBL bw   = valw * dG[p];
        OCP_DBL valb = valw * Vfi_zi;
        OCP_DBL bb   = valb * dG[p];

        // Bulk to Well

        // diag
        USI ptr = mySolver.diagPtr[k];
        mySolver.val[k][ptr] += valb;
        // off diag
        mySolver.colId[k].push_back(wId);
        mySolver.val[k].push_back(-valb);
        // b
        mySolver.b[k] += bb;

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // diag
                mySolver.diagVal[wId] += valw;
                // off diag
                mySolver.colId[wId].push_back(k);
                mySolver.val[wId].push_back(-valw);
                // b
                mySolver.b[wId] -= bw;
                break;
            case BHP_MODE:
                mySolver.colId[wId].push_back(k);
                mySolver.val[wId].push_back(0);
                break;
            default:
                ERRORcheck("Wrong Well Opt mode");
                exit(0);
        }
    }

    // Well Self
    assert(mySolver.val[wId].size() == numPerf);
    // the order of perforation is not necessarily in order
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            // diag
            mySolver.colId[wId].push_back(wId);
            mySolver.diagPtr[wId] = numPerf;
            mySolver.val[wId].push_back(mySolver.diagVal[wId]);
            // b
            mySolver.b[wId] += dt * opt.maxRate;
            break;
        case BHP_MODE:
            // diag
            mySolver.colId[wId].push_back(wId);
            mySolver.diagPtr[wId] = numPerf;
            mySolver.val[wId].push_back(dt);
            // b
            mySolver.b[wId] += dt * opt.maxBHP;
            // u   initial value
            mySolver.u[wId] = opt.maxBHP;
            break;
        default:
            ERRORcheck("Wrong Well Opt mode in function");
            exit(0);
    }
}

void Well::AssembleMat_PROD_BLK_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver,
                                      const OCP_DBL& dt) const
{
    USI     np  = myBulk.numPhase;
    USI     nc  = myBulk.numCom;
    OCP_USI wId = mySolver.dim;
    // important !
    mySolver.dim++;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;

        OCP_DBL valb = 0;
        OCP_DBL bb   = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw   = 0;

        for (USI j = 0; j < np; j++) {
            if (!myBulk.phaseExist[k * np + j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < nc; i++) {
                tempb += myBulk.vfi[k * nc + i] * myBulk.cij[k * np * nc + j * nc + i];
                tempw += opt.zi[i] * myBulk.cij[k * np * nc + j * nc + i];
            }
            OCP_DBL trans = dt * perf[p].transj[j] * myBulk.xi[k * np + j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = dG[p] - myBulk.Pc[k * np + j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        // diag
        USI ptr = mySolver.diagPtr[k];
        mySolver.val[k][ptr] += valb;
        // off diag
        mySolver.colId[k].push_back(wId);
        mySolver.val[k].push_back(-valb);
        // b
        mySolver.b[k] += bb;

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // diag  !!! attention! sign is -
                mySolver.diagVal[wId] -= valw;
                // off diag
                mySolver.colId[wId].push_back(k);
                mySolver.val[wId].push_back(valw);
                // b
                mySolver.b[wId] += bw;
                break;
            case BHP_MODE:
                // off diag
                mySolver.colId[wId].push_back(k);
                mySolver.val[wId].push_back(0);
                break;
            default:
                ERRORcheck("Wrong Well Opt mode");
                exit(0);
        }
    }

    // Well Self
    assert(mySolver.val[wId].size() == numPerf);
    // the order of perforation is not necessarily in order
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            // diag
            mySolver.colId[wId].push_back(wId);
            mySolver.diagPtr[wId] = numPerf;
            mySolver.val[wId].push_back(mySolver.diagVal[wId]);
            // b
            mySolver.b[wId] += dt * opt.maxRate;
            break;
        case BHP_MODE:
            // diag
            mySolver.colId[wId].push_back(wId);
            mySolver.diagPtr[wId] = numPerf;
            mySolver.val[wId].push_back(dt);
            // b
            mySolver.b[wId] += dt * opt.minBHP;
            // u   initial value
            mySolver.u[wId] = opt.minBHP;
            break;
        default:
            ERRORcheck("Wrong Well Opt mode");
            exit(0);
    }
}

void Well::SmoothdG()
{
    for (USI p = 0; p < numPerf; p++) {
        dG[p] = (ldG[p] + dG[p]) / 2; // seems better
                                      // dG[p] = ldG[p] + 0.618 * (dG[p] - ldG[p]);
    }
}

void Well::CaldG(const Bulk& myBulk)
{
    if (opt.type == INJ)
        CalInjdG(myBulk);
    else
        CalProddG(myBulk);
}

void Well::CalInjdG(const Bulk& myBulk)
{
    OCP_DBL         maxlen  = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil((perf[0].depth - depth) / maxlen);
                seg_len = (perf[0].depth - depth) / seg_num;
            } else {
                seg_num = ceil((perf[p].depth - perf[p - 1].depth) / maxlen);
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp -=
                    myBulk.flashCal[pvtnum]->RhoPhase(Ptmp, myBulk.T, opt.zi.data()) *
                    GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    } else if (depth >= perf[numPerf - 1].depth) {
        // Well is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil((depth - perf[numPerf - 1].depth) / maxlen);
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            } else {
                seg_num = ceil((perf[p + 1].depth - perf[p].depth) / maxlen);
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp +=
                    myBulk.flashCal[pvtnum]->RhoPhase(Ptmp, myBulk.T, opt.zi.data()) *
                    GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}

void Well::CalProddG(const Bulk& myBulk)
{
    USI             np      = myBulk.numPhase;
    USI             nc      = myBulk.numCom;
    OCP_DBL         maxlen  = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> tmpNi(nc, 0);
    vector<OCP_DBL> dGperf(numPerf, 0);
    OCP_DBL         qtacc  = 0;
    OCP_DBL         rhoacc = 0;
    OCP_DBL         rhotmp = 0;

    if (depth <= perf.front().depth) {
        // Well is higher

        // check qi_lbmol   ----   test
        if (perf[numPerf - 1].state == CLOSE) {
            for (OCP_INT p = numPerf - 2; p >= 0; p--) {
                if (perf[p].state == OPEN) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (OCP_INT p = numPerf - 1; p >= 0; p--) {

            if (p == 0) {
                seg_num = ceil((perf[0].depth - depth) / maxlen);
                seg_len = (perf[0].depth - depth) / seg_num;
            } else {
                seg_num = ceil((perf[p].depth - perf[p - 1].depth) / maxlen);
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }

            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            tmpNi.assign(nc, 0);
            for (OCP_INT p1 = numPerf - 1; p1 >= p; p1--) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }
            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->v[j] / seg_num;
                        rhoacc += myBulk.flashCal[pvtnum]->v[j] * rhotmp *
                                  GRAVITY_FACTOR / seg_num;
                    }
                }
                Ptmp -= rhoacc / qtacc * seg_len;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    } else if (depth >= perf.back().depth) {
        // Well is lower

        // check qi_lbmol   ----   test
        if (perf[0].state == CLOSE) {
            for (USI p = 1; p <= numPerf; p++) {
                if (perf[p].state == OPEN) {
                    for (USI i = 0; i < nc; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil((depth - perf[numPerf - 1].depth) / maxlen);
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            } else {
                seg_num = ceil((perf[p + 1].depth - perf[p].depth) / maxlen);
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }

            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            tmpNi.assign(nc, 0);
            for (OCP_INT p1 = numPerf - 1; p1 >= p; p1--) {
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }
            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
                for (USI j = 0; j < np; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->v[j] / seg_num;
                        rhoacc += myBulk.flashCal[pvtnum]->v[j] * rhotmp *
                                  GRAVITY_FACTOR / seg_num;
                    }
                }
                Ptmp += rhoacc / qtacc * seg_len;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }

    } else {
        ERRORcheck("Wrong Well position");
        exit(0);
    }
}

void Well::CalTrans(const Bulk& myBulk)
{
    USI np = myBulk.numPhase;

    if (opt.type == INJ) {
        for (USI p = 0; p < numPerf; p++) {
            perf[p].transj.assign(np, 0);
            OCP_USI k    = perf[p].location;
            OCP_DBL temp = CONV1 * CONV2 * perf[p].WI * perf[p].multiplier;

            // single phase
            for (USI j = 0; j < np; j++) {
                OCP_USI id = k * np + j;
                if (myBulk.phaseExist[id]) {
                    perf[p].transj[0] += myBulk.kr[id] / myBulk.mu[id];
                }
            }
            perf[p].transj[0] *= temp;
        }
    } else {
        for (USI p = 0; p < numPerf; p++) {
            perf[p].transj.assign(np, 0);
            OCP_USI k    = perf[p].location;
            OCP_DBL temp = CONV1 * CONV2 * perf[p].WI * perf[p].multiplier;

            // multi phase
            for (USI j = 0; j < np; j++) {
                OCP_USI id = k * np + j;
                if (myBulk.phaseExist[id]) {
                    perf[p].transj[j] = temp * myBulk.kr[id] / myBulk.mu[id];
                }
            }
        }
    }
}

void Well::CalFlux(const Bulk& myBulk, const bool flag)
{
    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;

    if (opt.type == INJ) {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P  = BHP + dG[p];
            OCP_USI k  = perf[p].location;
            OCP_DBL dP = perf[p].P - myBulk.P[k];
            dP *= -1.0;

            perf[p].qt_ft3 = perf[p].transj[0] * dP;

            USI pvtnum = myBulk.PVTNUM[k];
            if (flag)
                perf[p].xi =
                    myBulk.flashCal[pvtnum]->XiPhase(myBulk.P[k], myBulk.T, &opt.zi[0]);
            OCP_DBL xi = perf[p].xi;
            for (USI i = 0; i < nc; i++) {
                perf[p].qi_lbmol[i] = perf[p].qt_ft3 * xi * opt.zi[i];
            }
        }
    } else {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P      = BHP + dG[p];
            OCP_USI k      = perf[p].location;
            perf[p].qt_ft3 = 0;
            perf[p].qi_lbmol.assign(nc, 0);

            for (USI j = 0; j < np; j++) {
                OCP_USI id = k * np + j;
                if (myBulk.phaseExist[id]) {
                    OCP_DBL dP = myBulk.Pj[id] - perf[p].P;
                    perf[p].qt_ft3 += perf[p].transj[j] * dP;
                    // cout << p << " p[" << j << "] = " << myBulk.Pj[id] << endl;
                    // cout << p << " perf = " << perf[p].p << endl;

                    OCP_DBL xi = myBulk.xi[id];
                    OCP_DBL xij;
                    for (USI i = 0; i < nc; i++) {
                        xij = myBulk.cij[id * nc + i];
                        perf[p].qi_lbmol[i] += perf[p].transj[j] * dP * xi * xij;
                    }
                }
            }
        }
    }
}

void Well::MassConserve(Bulk& myBulk, const OCP_DBL& dt) const
{
    USI nc = myBulk.numCom;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;
        for (USI i = 0; i < nc; i++) {
            myBulk.Ni[k * nc + i] -= perf[p].qi_lbmol[i] * dt;
        }
    }
}

OCP_DBL Well::CalInjRate_Blk(const Bulk& myBulk)
{
    USI     np = myBulk.numPhase;
    USI     nc = myBulk.numCom;
    OCP_DBL qj = 0;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = opt.maxBHP + dG[p];
        OCP_USI k     = perf[p].location;

        USI     pvtnum = myBulk.PVTNUM[k];
        OCP_DBL xi =
            myBulk.flashCal[pvtnum]->XiPhase(myBulk.P[k], myBulk.T, &opt.zi[0]);
        OCP_DBL dP = Pperf - myBulk.P[k];
        qj += perf[p].transj[0] * xi * dP;
    }
    return qj;
}

OCP_DBL Well::CalProdRate_Blk(const Bulk& myBulk)
{
    USI     np = myBulk.numPhase;
    USI     nc = myBulk.numCom;
    OCP_DBL qj = 0;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = opt.minBHP + dG[p];
        OCP_USI k     = perf[p].location;

        for (USI j = 0; j < np; j++) {
            OCP_USI id = k * np + j;
            if (myBulk.phaseExist[id]) {
                OCP_DBL temp = 0;
                for (USI i = 0; i < nc; i++) {
                    temp += opt.zi[i] * myBulk.cij[id * nc + i];
                }
                OCP_DBL xi = myBulk.xi[id];
                OCP_DBL dP = myBulk.Pj[id] - Pperf;
                qj += perf[p].transj[j] * xi * dP * temp;
            }
        }
    }
    return qj;
}

void Well::CalInjQi_Blk(const Bulk& myBulk, const OCP_DBL& dt)
{
    USI     nc = myBulk.numCom;
    OCP_DBL qj = 0;

    for (USI p = 0; p < numPerf; p++) {

        perf[p].P = BHP + dG[p];

        // OCP_USI k = perf[p].location;
        // OCP_DBL xi = perf[p].xi;
        // OCP_DBL dP = perf[p].p - myBulk.p[k];
        // qj += perf[p].transj[0] * xi * dP;

        for (USI i = 0; i < nc; i++) qj += perf[p].qi_lbmol[i];
    }
    if (opt.fluidType == WATER) {
        WWIR = -qj;
        WWIT += WWIR * dt;
    } else {
        WGIR = -qj;
        WGIT += WGIR * dt;
    }
}

void Well::CalProdQi_Blk(const Bulk& myBulk, const OCP_DBL& dt)
{
    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;

    qi_lbmol.assign(nc, 0);

    for (USI p = 0; p < numPerf; p++) {

        // perf[p].qi_lbmol.assign(nc, 0);
        // perf[p].p = BHP + dG[p];
        // int k = perf[p].location;

        // for (int j = 0; j < np; j++) {
        //	int id = k * np + j;
        //	if (myBulk.phaseExist[id]) {
        //		OCP_DBL xi = myBulk.xi[id];
        //		OCP_DBL dP = myBulk.Pj[id] - perf[p].p;
        //		OCP_DBL xij;
        //		for (int i = 0; i < nc; i++) {
        //			xij = myBulk.cij[id * nc + i];
        //			perf[p].qi_lbmol[i] += perf[p].transj[j] * xi * xij * dP;
        //		}
        //	}
        //}
        for (USI i = 0; i < nc; i++) {
            qi_lbmol[i] += perf[p].qi_lbmol[i];
        }
    }

    for (USI i = 0; i < nc; i++) {
        if (myBulk.index2Phase[i] == OIL) {
            WOPR = qi_lbmol[i];
            WOPT += WOPR * dt;
        } else if (myBulk.index2Phase[i] == GAS) {
            WGPR = qi_lbmol[i];
            WGPT += WGPR * dt;
        } else if (myBulk.index2Phase[i] == WATER) {
            WWPR = qi_lbmol[i];
            WWPT += WWPR * dt;
        }
    }
}

void Well::CheckOptMode(const Bulk& myBulk)
{
    if (opt.type == INJ) {
        if (CalInjRate_Blk(myBulk) > opt.maxRate) {
            opt.optMode = RATE_MODE;
        } else {
            opt.optMode = BHP_MODE;
        }
    } else {
        if (CalProdRate_Blk(myBulk) > opt.maxRate) {
            opt.optMode = RATE_MODE;
        } else {
            opt.optMode = BHP_MODE;
        }
    }
}

OCP_INT Well::CheckP(const Bulk& myBulk)
{
    // 0 : all correct
    // 1 : negative P
    // 2 : outlimited P
    // 3 : crossflow happens

    if (BHP < 0) return 1;
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == OPEN && perf[p].P < 0) {
#ifdef _DEBUG
            cout << "WARNING: Well " << name << " Perf[" << p << "].P = " << perf[p].P
                 << endl;
#endif // _DEBUG
            return 1;
        }
    }

    if (opt.type == INJ) {
        if (opt.optMode != BHP_MODE && BHP > opt.maxBHP) {
            opt.optMode = BHP_MODE;
            return 2;
        }
    } else {
        if (opt.optMode != BHP_MODE && BHP < opt.minBHP) {
            opt.optMode = BHP_MODE;
            return 2;
        }
    }
    OCP_INT flag = 0;
    flag         = CheckCrossFlow(myBulk);

    return flag;
}

OCP_INT Well::CheckCrossFlow(const Bulk& myBulk)
{
    OCP_USI k;
    bool    flagC = true;

    if (opt.type == PROD) {
        USI np = myBulk.numPhase;
        for (USI p = 0; p < numPerf; p++) {
            k            = perf[p].location;
            OCP_DBL minP = myBulk.P[k];
            // THINK MORE !!!
            // for (USI j = 0; j < np; j++) {
            //	if (myBulk.phaseExist[k * np + j])
            //		minP = minP < myBulk.Pj[k * np + j] ? minP : myBulk.Pj[k * np + j];
            //}
            if (perf[p].state == OPEN && minP < perf[p].P) {
                perf[p].state      = CLOSE;
                perf[p].multiplier = 0;
                flagC              = false;
            } else if (perf[p].state == CLOSE && minP > perf[p].P) {
                perf[p].state      = OPEN;
                perf[p].multiplier = 1;
            }
        }
    } else {
        for (USI p = 0; p < numPerf; p++) {
            k = perf[p].location;
            if (perf[p].state == OPEN && myBulk.P[k] > perf[p].P) {
                perf[p].state      = CLOSE;
                perf[p].multiplier = 0;
                flagC              = false;
            } else if (perf[p].state == CLOSE && myBulk.P[k] < perf[p].P) {
                perf[p].state      = OPEN;
                perf[p].multiplier = 1;
            }
        }
    }

    bool flag = false;
    // check well --  if all perf are closed, open the depthest perf
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == OPEN) {
            flag = true;
            break;
        }
    }
    if (!flag) {
        return 1;
        // open the depthest perf
        /*perf.back().state = OPEN;
        perf.back().multiplier = 1;*/
    }

    if (!flagC) {
        dG = ldG;
        CalTrans(myBulk);
        CalFlux(myBulk);
        CaldG(myBulk);
        SmoothdG();
        // CheckOptMode(myBulk);
        return 3;
    }

    return 0;
}

void Well::ShowPerfStatus() const
{
    cout << "----------------------------" << endl;
    cout << name << ":    " << opt.optMode << endl;
    for (USI p = 0; p < numPerf; p++) {
        cout << "Perf[" << p << "].State = " << perf[p].state << endl;
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