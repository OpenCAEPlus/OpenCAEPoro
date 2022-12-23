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


void Well::InputPerfo(const WellParam& well)
{
    OCP_FUNCNAME;

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
            OCP_ABORT("Wrong direction of perforations!");
        }
    }
}

void Well::Setup(const Grid& myGrid, const Bulk& myBulk, const vector<SolventINJ>& sols)
{
    OCP_FUNCNAME;

    numCom = myBulk.numCom;
    numPhase = myBulk.numPhase;
    flashCal = myBulk.flashCal;

    qi_lbmol.resize(numCom);
    prodWeight.resize(numCom);
    prodRate.resize(numPhase);

    for (auto& opt : optSet) {
        if (!opt.state) continue;

        opt.Tinj = myBulk.RTemp;
        flashCal[0]->SetupWellOpt(opt, sols, Psurf, Tsurf);
    }
    
    // Perf
    USI pp = 0;
    for (USI p = 0; p < numPerf; p++) {
        OCP_USI pId =
            perf[p].K * myGrid.nx * myGrid.ny + perf[p].J * myGrid.nx + perf[p].I;
        if (myGrid.map_All2Flu[pId].IsAct()) {

            perf[pp]            = perf[p];
            perf[pp].state      = OPEN;
            perf[pp].location   = myGrid.map_All2Act[pId].GetId();
            perf[pp].depth      = myBulk.depth[perf[pp].location];
            perf[pp].multiplier = 1;
            perf[pp].qi_lbmol.resize(numCom);
            perf[pp].transj.resize(numPhase);
            perf[pp].qj_ft3.resize(numPhase);
            pp++;
        } else {
            OCP_WARNING("Perforation is in non-fluid myBulk!");
        }
    }
    numPerf = pp;
    perf.resize(numPerf);
    // dG
    dG.resize(numPerf, 0);
    ldG = dG;

    if (depth < 0) depth = perf[0].depth;

    CalWI_Peaceman_Vertical(myBulk);
    // test
    // ShowPerfStatus(myBulk);
}

void Well::InitBHP(const Bulk& myBulk) { BHP = myBulk.P[perf[0].location]; }

void Well::CalWI_Peaceman_Vertical(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    // this fomular needs to be carefully checked !
    // especially the dz

    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].WI > 0) {
            break;
        } else {
            const OCP_USI Idb = perf[p].location;
            const OCP_DBL dx  = myBulk.dx[Idb];
            const OCP_DBL dy  = myBulk.dy[Idb];
            const OCP_DBL dz  = myBulk.dz[Idb] * myBulk.ntg[Idb];
            OCP_DBL ro  = 0;
            switch (perf[p].direction) {
                case X_DIRECTION:
                    {
                        const OCP_DBL kykz  = myBulk.rockKy[Idb] * myBulk.rockKz[Idb];
                        const OCP_DBL ky_kz = myBulk.rockKy[Idb] / myBulk.rockKz[Idb];
                        assert(kykz > 0);
                        ro = 0.28 * pow((dy * dy * pow(1 / ky_kz, 0.5) +
                                         dz * dz * pow(ky_kz, 0.5)),
                                        0.5);
                        ro /= (pow(ky_kz, 0.25) + pow(1 / ky_kz, 0.25));

                        if (perf[p].kh < 0) {
                            perf[p].kh = (dx * pow(kykz, 0.5));
                        }
                        break;
                    }
                case Y_DIRECTION:
                    {
                        const OCP_DBL kzkx  = myBulk.rockKz[Idb] * myBulk.rockKx[Idb];
                        const OCP_DBL kz_kx = myBulk.rockKz[Idb] / myBulk.rockKx[Idb];
                        assert(kzkx > 0);
                        ro = 0.28 * pow((dz * dz * pow(1 / kz_kx, 0.5) +
                                         dx * dx * pow(kz_kx, 0.5)),
                                        0.5);
                        ro /= (pow(kz_kx, 0.25) + pow(1 / kz_kx, 0.25));

                        if (perf[p].kh < 0) {
                            perf[p].kh = (dy * pow(kzkx, 0.5));
                        }
                        break;
                    }
                case Z_DIRECTION:
                    {
                        const OCP_DBL kxky  = myBulk.rockKx[Idb] * myBulk.rockKy[Idb];
                        const OCP_DBL kx_ky = myBulk.rockKx[Idb] / myBulk.rockKy[Idb];
                        assert(kxky > 0);
                        ro = 0.28 * pow((dx * dx * pow(1 / kx_ky, 0.5) +
                                         dy * dy * pow(kx_ky, 0.5)),
                                        0.5);
                        ro /= (pow(kx_ky, 0.25) + pow(1 / kx_ky, 0.25));

                        if (perf[p].kh < 0) {
                            perf[p].kh = (dz * pow(kxky, 0.5));
                        }
                        break;
                    }
                default:
                    OCP_ABORT("Wrong direction of perforations!");
            }
            perf[p].WI = CONV2 * (2 * PI) * perf[p].kh /
                         (log(ro / perf[p].radius) + perf[p].skinFactor);
        }
    }
}

void Well::CalTrans(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    if (opt.type == INJ) {
        for (USI p = 0; p < numPerf; p++) {
            perf[p].transINJ = 0;
            OCP_USI k        = perf[p].location;
            OCP_DBL temp     = CONV1 * perf[p].WI * perf[p].multiplier;

            // single phase
            for (USI j = 0; j < numPhase; j++) {
                perf[p].transj[j] = 0;
                OCP_USI id        = k * numPhase + j;
                if (myBulk.phaseExist[id]) {
                    perf[p].transj[j] = temp * myBulk.kr[id] / myBulk.mu[id];
                    perf[p].transINJ += perf[p].transj[j];
                }
            }
        }
    } else {
        for (USI p = 0; p < numPerf; p++) {
            OCP_USI k    = perf[p].location;
            OCP_DBL temp = CONV1 * perf[p].WI * perf[p].multiplier;

            // multi phase
            for (USI j = 0; j < numPhase; j++) {
                perf[p].transj[j] = 0;
                OCP_USI id        = k * numPhase + j;
                if (myBulk.phaseExist[id]) {
                    perf[p].transj[j] = temp * myBulk.kr[id] / myBulk.mu[id];
                }
            }
        }
    }
}

void Well::CalFlux(const Bulk& myBulk, const OCP_BOOL ReCalXi)
{
    OCP_FUNCNAME;

    // cout << name << endl;
    fill(qi_lbmol.begin(), qi_lbmol.end(), 0.0);

    if (opt.type == INJ) {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P  = BHP + dG[p];
            OCP_USI k  = perf[p].location;
            OCP_DBL dP = perf[p].P - myBulk.P[k];
            dP *= -1.0;

            perf[p].qt_ft3 = perf[p].transINJ * dP;

            if (ReCalXi) {
                USI pvtnum = myBulk.PVTNUM[k];
                perf[p].xi =
                    myBulk.flashCal[pvtnum]->XiPhase(perf[p].P, opt.Tinj, &opt.injZi[0], opt.injProdPhase);
            }
            for (USI i = 0; i < numCom; i++) {
                perf[p].qi_lbmol[i] = perf[p].qt_ft3 * perf[p].xi * opt.injZi[i];
                qi_lbmol[i] += perf[p].qi_lbmol[i];
            }
        }
    } else {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P      = BHP + dG[p];
            OCP_USI k      = perf[p].location;
            perf[p].qt_ft3 = 0;
            fill(perf[p].qi_lbmol.begin(), perf[p].qi_lbmol.end(), 0.0);
            fill(perf[p].qj_ft3.begin(), perf[p].qj_ft3.end(), 0.0);

            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = k * numPhase + j;
                if (myBulk.phaseExist[id]) {
                    OCP_DBL dP = myBulk.Pj[id] - perf[p].P;

                    perf[p].qj_ft3[j] = perf[p].transj[j] * dP;
                    perf[p].qt_ft3 += perf[p].qj_ft3[j];

                    OCP_DBL xi = myBulk.xi[id];
                    OCP_DBL xij;
                    for (USI i = 0; i < numCom; i++) {
                        xij = myBulk.xij[id * numCom + i];
                        perf[p].qi_lbmol[i] += perf[p].qj_ft3[j] * xi * xij;
                    }
                }
            }
            for (USI i = 0; i < numCom; i++) qi_lbmol[i] += perf[p].qi_lbmol[i];
        }
    }
}

/// Pressure in injection well equals maximum ones in injection well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL Well::CalInjRate(const Bulk& myBulk, const OCP_BOOL& maxBHP)
{
    OCP_FUNCNAME;

    OCP_DBL qj    = 0;
    OCP_DBL Pwell = maxBHP ? opt.maxBHP : BHP;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        USI pvtnum = myBulk.PVTNUM[k];
        OCP_DBL xi =
            myBulk.flashCal[pvtnum]->XiPhase(Pperf, opt.Tinj, &opt.injZi[0], opt.injProdPhase);

        OCP_DBL dP = Pperf - myBulk.P[k];
        qj += perf[p].transINJ * xi * dP;
    }
    return qj;
}


/// Pressure in production well equals minial ones in production well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL Well::CalProdRate(const Bulk& myBulk, const OCP_BOOL& minBHP)
{
    OCP_FUNCNAME;

    OCP_DBL   qj    = 0;
    OCP_DBL   Pwell = minBHP ? opt.minBHP : BHP;

    vector<OCP_DBL> tmpQi_lbmol(numCom, 0);
    vector<OCP_DBL> tmpQj(numPhase, 0);

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        for (USI j = 0; j < numPhase; j++) {
            OCP_USI id = k * numPhase + j;
            if (myBulk.phaseExist[id]) {
                OCP_DBL dP = myBulk.Pj[id] - Pperf;
                OCP_DBL temp = perf[p].transj[j] * myBulk.xi[id] * dP;
                for (USI i = 0; i < numCom; i++) {
                    tmpQi_lbmol[i] += myBulk.xij[id * numCom + i] * temp;
                }
            }
        }
    }
    flashCal[0]->CalProdRate(Psurf, Tsurf, &tmpQi_lbmol[0], tmpQj);
    for (USI j = 0; j < numPhase; j++) {
        qj += tmpQj[j] * opt.prodPhaseWeight[j];
    }

    return qj;
}

void Well::CalInjQi(const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    OCP_DBL   qj = 0;

    for (USI i = 0; i < numCom; i++) {
        qj += qi_lbmol[i];
    }
    if (opt.fluidType == "WAT") {
        WWIR = -qj;
        WWIT += WWIR * dt;
    } else {
        WGIR = -qj / opt.factorINJ; // Mscf or lbmol -> Mscf
        WGIT += WGIR * dt;
    }
}

void Well::CalProdQj(const Bulk& myBulk, const OCP_DBL& dt)
{

    flashCal[0]->CalProdRate(Psurf, Tsurf, &qi_lbmol[0], prodRate);
    USI iter = 0;
    if (myBulk.oil) WOPR = prodRate[iter++];
    if (myBulk.gas) WGPR = prodRate[iter++];
    if (myBulk.water) WWPR = prodRate[iter++];

    WOPT += WOPR * dt;
    WGPT += WGPR * dt;
    WWPT += WWPR * dt;
}


/// It calculates pressure difference between perforations iteratively.
/// This function can be used in both black oil model and compositional model.
/// stability of this method shoule be tested.
void Well::CaldG(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    if (opt.type == INJ)
        CalInjdG(myBulk);
    else
        CalProddG01(myBulk);
}

void Well::CalInjdG(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const OCP_DBL   maxlen  = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp -=
                    myBulk.flashCal[pvtnum]->RhoPhase(Ptmp, 0, opt.Tinj, opt.injZi.data(), opt.injProdPhase) *
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
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp +=
                    myBulk.flashCal[pvtnum]->RhoPhase(Ptmp, 0, opt.Tinj, opt.injZi.data(), opt.injProdPhase) *
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

void Well::CalProddG01(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const OCP_DBL   maxlen  = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(numCom, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] +=
                        perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * numCom + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(numCom, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < numCom; i++) {
                    tmpNi[i] = myBulk.Ni[n * numCom + i];
                }
            }

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j];
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * GRAVITY_FACTOR * seg_len;
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
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] +=
                        perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * numCom + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(numCom, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < numCom; i++) {
                    tmpNi[i] = myBulk.Ni[n * numCom + i];
                }
            }

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j];
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Ptmp - Pperf;
        }
        dG[numPerf - 1] = dGperf[numPerf - 1];
        for (OCP_INT p = numPerf - 2; p >= 0; p--) {
            dG[p] = dG[p + 1] + dGperf[p];
        }
    }
}


void Well::CalProddG02(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const OCP_DBL   maxlen = 10;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> dGperf(numPerf, 0);
    vector<OCP_DBL> tmpNi(numCom, 0);
    OCP_DBL         rhotmp, qtacc, rhoacc;

    if (depth <= perf.front().depth) {
        // Well is higher
        for (OCP_INT p = numPerf - 1; p >= 0; p--) {
            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;


            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] +=
                        (perf[p].transj[j] > 0) * myBulk.xi[id] * myBulk.xij[id * numCom + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(numCom, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < numCom; i++) {
                    tmpNi[i] = myBulk.Ni[n * numCom + i];
                }
            }

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j];
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp;
                    }
                }
                Ptmp -= rhoacc / qtacc * GRAVITY_FACTOR * seg_len;
            }
            dGperf[p] = Pperf - Ptmp;
        }
        dG[0] = dGperf[0];
        for (USI p = 1; p < numPerf; p++) {
            dG[p] = dG[p - 1] + dGperf[p];
        }
    }
    else if (depth >= perf[numPerf - 1].depth) {
        // Well is lower
        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            }
            else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }
            OCP_USI n = perf[p].location;
            perf[p].P = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] +=
                        (perf[p].transj[j] > 0) * myBulk.xi[id] * myBulk.xij[id * numCom + k];
                }
            }
            OCP_DBL tmpSum = Dnorm1(numCom, &tmpNi[0]);
            if (tmpSum < TINY) {
                for (USI i = 0; i < numCom; i++) {
                    tmpNi[i] = myBulk.Ni[n * numCom + i];
                }
            }

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                qtacc = rhoacc = 0;
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j];
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp;
                    }
                }
                Ptmp += rhoacc / qtacc * GRAVITY_FACTOR * seg_len;
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
    OCP_FUNCNAME;

    const OCP_DBL   maxlen  = 5;
    USI             seg_num = 0;
    OCP_DBL         seg_len = 0;
    vector<OCP_DBL> tmpNi(numCom, 0);
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
                    for (USI i = 0; i < numCom; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (OCP_INT p = numPerf - 1; p >= 0; p--) {

            if (p == 0) {
                seg_num = ceil(fabs((perf[0].depth - depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[0].depth - depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p].depth - perf[p - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p].depth - perf[p - 1].depth) / seg_num;
            }

            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            // tmpNi.assign(numCom, 0);
            // for (OCP_INT p1 = numPerf - 1; p1 >= p; p1--) {
            //    for (USI i = 0; i < numCom; i++) {
            //        tmpNi[i] += perf[p1].qi_lbmol[i];
            //    }
            //}

            for (USI i = 0; i < numCom; i++) {
                tmpNi[i] += perf[p].qi_lbmol[i];
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j] / seg_num;
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp *
                                  GRAVITY_FACTOR / seg_num;
#ifdef DEBUG
                        if (rhotmp <= 0 || !isfinite(rhotmp)) {
                            OCP_ABORT("Wrong rho " + to_string(rhotmp));
                        }
#endif // DEBUG
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
                    for (USI i = 0; i < numCom; i++) {
                        perf[numPerf - 1].qi_lbmol[i] = perf[p].qi_lbmol[i];
                    }
                    break;
                }
            }
        }

        for (USI p = 0; p < numPerf; p++) {
            if (p == numPerf - 1) {
                seg_num = ceil(fabs((depth - perf[numPerf - 1].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (depth - perf[numPerf - 1].depth) / seg_num;
            } else {
                seg_num = ceil(fabs((perf[p + 1].depth - perf[p].depth) / maxlen));
                if (seg_num == 0) continue;
                seg_len = (perf[p + 1].depth - perf[p].depth) / seg_num;
            }

            OCP_USI n     = perf[p].location;
            perf[p].P     = BHP + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (OCP_INT p1 = numPerf - 1; p1 - p >= 0; p1--) {
                for (USI i = 0; i < numCom; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T[n], tmpNi.data());
                for (USI j = 0; j < numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->vj[j] / seg_num;
                        rhoacc += myBulk.flashCal[pvtnum]->vj[j] * rhotmp *
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
        OCP_ABORT("Wrong well position!");
    }
}


void Well::CalProdWeight(const Bulk& myBulk) const
{

    if (opt.type == PROD) {
        // in some cases, qi_lbmol may be zero, so use other methods
        OCP_DBL  qt   = 0;
        OCP_BOOL flag = OCP_TRUE;
        for (USI i = 0; i < myBulk.numCom; i++) {
            qt += qi_lbmol[i];
            if (qi_lbmol[i] < 0) flag = OCP_FALSE;
        }
        if (qt > TINY && flag) {
            flashCal[0]->CalProdWeight(Psurf, Tsurf, &qi_lbmol[0],
                                              opt.prodPhaseWeight, prodWeight);
        } else {
            vector<OCP_DBL> tmpNi(numCom, 0);
            for (USI p = 0; p < numPerf; p++) {
                OCP_USI n = perf[p].location;

                for (USI j = 0; j < numPhase; j++) {
                    OCP_USI id = n * numPhase + j;
                    if (!myBulk.phaseExist[id]) continue;
                    for (USI k = 0; k < numCom; k++) {
                        tmpNi[k] +=
                            perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * numCom + k];
                    }
                }
            }
            qt = Dnorm1(numCom, &tmpNi[0]);
            flashCal[0]->CalProdWeight(Psurf, Tsurf, &tmpNi[0],
                                              opt.prodPhaseWeight, prodWeight);
        }
    }
}


void Well::CalReInjFluid(const Bulk& myBulk, vector<OCP_DBL>& myZi)
{
    CalTrans(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;

        for (USI j = 0; j < numPhase; j++) {
            OCP_USI id = n * numPhase + j;
            if (!myBulk.phaseExist[id]) continue;
            for (USI k = 0; k < numCom; k++) {
                myZi[k] += perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * numCom + k];
            }
        }
    }
}

void Well::SetBHP()
{
    if (opt.type == PROD && opt.optMode == BHP_MODE) {
        BHP = opt.minBHP;
    } else if (opt.type == INJ && opt.optMode == BHP_MODE) {
        BHP = opt.maxBHP;
    }
}

/// It's just a test now to make dG more stable.
void Well::SmoothdG()
{
    OCP_FUNCNAME;

    for (USI p = 0; p < numPerf; p++) {
        dG[p] = (ldG[p] + dG[p]) / 2; // seems better
                                      // dG[p] = ldG[p] + 0.618 * (dG[p] - ldG[p]);
    }
}

/// Constant well pressure would be applied if flow rate is too large.
/// Constant flow rate would be applied if well pressure is outranged.
void Well::CheckOptMode(const Bulk& myBulk)
{
    OCP_FUNCNAME;
    if (opt.initOptMode == BHP_MODE) {
        if (opt.type == INJ) {
            OCP_DBL q = CalInjRate(myBulk, OCP_TRUE);
            // for INJ well, maxRate has been switch to lbmols
            OCP_DBL tarRate = opt.maxRate;
            if (opt.reInj) {
                if (opt.reInjPhase == GAS)
                    tarRate = WGIR;
                else if (opt.reInjPhase == WATER)
                    tarRate = WWIR;
            }
            if (q > tarRate) {
                opt.optMode = RATE_MODE;
            } else {
                opt.optMode = BHP_MODE;
                BHP         = opt.maxBHP;
            }
        } else {
            opt.optMode = BHP_MODE;
            BHP = opt.minBHP;
        }
    } else {
        if (opt.type == INJ) {
            OCP_DBL q = CalInjRate(myBulk, OCP_TRUE);
            // for INJ well, maxRate has been switch to lbmols
            OCP_DBL tarRate = opt.maxRate;
            if (opt.reInj) {
                if (opt.reInjPhase == GAS)
                    tarRate = WGIR;
                else if (opt.reInjPhase == WATER)
                    tarRate = WWIR;
            }

            if (q > tarRate) {
                opt.optMode = opt.initOptMode;
            } else {
                opt.optMode = BHP_MODE;
                BHP         = opt.maxBHP;
            }
        } else {
            OCP_DBL q = CalProdRate(myBulk, OCP_TRUE);
            // cout << q << endl;
            if (q > opt.maxRate) {
                opt.optMode = opt.initOptMode;
            } else {
                opt.optMode = BHP_MODE;
                BHP         = opt.minBHP;
            }
        }
    }
}

OCP_INT Well::CheckP(const Bulk& myBulk)
{
    OCP_FUNCNAME;
    // 0 : all correct
    // 1 : negative P
    // 2 : outlimited P
    // 3 : crossflow happens

    if (BHP < 0) {
        cout << "### WARNING: Negative BHP " << BHP << endl;
        return 1;
    }
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == OPEN && perf[p].P < 0) {
#ifdef DEBUG
            cout << "### WARNING: Well " << name << " Perf[" << p
                 << "].P = " << perf[p].P << endl;
#endif // DEBUG
            cout << "### WARNING: Negative perforation P " << perf[p].P << endl;
            return 1;
        }
    }

    if (opt.type == INJ) {
        if (opt.optMode != BHP_MODE && BHP > opt.maxBHP) {
#if _DEBUG
            cout << "### WARNING: Well " << name << " switch to BHPMode" << endl;
#endif
            opt.optMode = BHP_MODE;
            BHP         = opt.maxBHP;
            return 2;
        }
    } else {
        if (opt.optMode != BHP_MODE && BHP < opt.minBHP) {
#if _DEBUG
            cout << "### WARNING: Well " << name << " switch to BHPMode" << endl;
#endif
            opt.optMode = BHP_MODE;
            BHP         = opt.minBHP;
            return 2;
        }
    }

    return CheckCrossFlow(myBulk);
}

OCP_INT Well::CheckCrossFlow(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    OCP_USI  k;
    OCP_BOOL flagC = OCP_TRUE;

    if (opt.type == PROD) {
        for (USI p = 0; p < numPerf; p++) {
            k            = perf[p].location;
            OCP_DBL minP = myBulk.P[k];
            if (perf[p].state == OPEN && minP < perf[p].P) {
                cout << std::left << std::setw(12) << name << "  "
                     << "Well P = " << perf[p].P << ", "
                     << "Bulk P = " << minP << endl;
                perf[p].state      = CLOSE;
                perf[p].multiplier = 0;
                flagC              = OCP_FALSE;
                break;
            } else if (perf[p].state == CLOSE && minP > perf[p].P) {
                perf[p].state      = OPEN;
                perf[p].multiplier = 1;
            }
        }
    } else {
        for (USI p = 0; p < numPerf; p++) {
            k = perf[p].location;
            if (perf[p].state == OPEN && myBulk.P[k] > perf[p].P) {
                cout << std::left << std::setw(12) << name << "  "
                     << "Well P = " << perf[p].P << ", "
                     << "Bulk P = " << myBulk.P[k] << endl;
                perf[p].state      = CLOSE;
                perf[p].multiplier = 0;
                flagC              = OCP_FALSE;
                break;
            } else if (perf[p].state == CLOSE && myBulk.P[k] < perf[p].P) {
                perf[p].state      = OPEN;
                perf[p].multiplier = 1;
            }
        }
    }

    OCP_BOOL flag = OCP_FALSE;
    // check well --  if all perf are closed, open the depthest perf
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == OPEN) {
            flag = OCP_TRUE;
            break;
        }
    }

    if (!flag) {
        // open the deepest perf
        perf.back().state      = OPEN;
        perf.back().multiplier = 1;
        cout << "### WARNING: All perfs of " << name
             << " are closed! Open the last perf!\n";
    }

    if (!flagC) {
        // if crossflow happens, then corresponding perforation will be closed,
        // then the multiplier of perforation will be set to zero, so trans of well
        // should be recalculated!
        //
        // dG = ldG;
        CalTrans(myBulk);
        // CalFlux(myBulk);
        // CaldG(myBulk);
        // SmoothdG();
        // CheckOptMode(myBulk);
        return 3;
    }

    return 0;
}

void Well::ShowPerfStatus(const Bulk& myBulk) const
{
    OCP_FUNCNAME;

    cout << fixed;
    cout << "----------------------------" << endl;
    cout << name << ":    " << opt.optMode << "   " << setprecision(3) << BHP << endl;
    for (USI p = 0; p < numPerf; p++) {
        vector<OCP_DBL> Qitmp(perf[p].qi_lbmol);
        // OCP_DBL         qt = Dnorm1(myBulk.numCom, &Qitmp[0]);
        OCP_USI n = perf[p].location;
        cout << setw(3) << p << "   " << perf[p].state << "   " << setw(6)
             << perf[p].location << "  " << setw(2) << perf[p].I + 1 << "  " << setw(2)
             << perf[p].J + 1 << "  " << setw(2) << perf[p].K + 1 << "  " << setw(10)
             << setprecision(6) << perf[p].WI << "  "               // ccf
             << setprecision(3) << perf[p].radius << "  "           // ccf
             << setw(8) << setprecision(4) << perf[p].kh << "  "    // kh
             << setw(8) << setprecision(2) << perf[p].depth << "  " // depth
             << setprecision(3) << perf[p].P << "  "                // Pp
             << setw(10) << setprecision(3) << myBulk.P[n] << "   " // Pb
             << setw(6) << setprecision(3) << dG[p] << "   "        // dG
             << setw(8) << perf[p].qi_lbmol[myBulk.numCom - 1] << "   " << setw(6)
             << setprecision(6) << myBulk.S[n * myBulk.numPhase + 0] << "   " << setw(6)
             << setprecision(6) << myBulk.S[n * myBulk.numPhase + 1] << "   " << setw(6)
             << setprecision(6) << myBulk.S[n * myBulk.numPhase + 2] << endl;
    }
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////

void Well::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    if (opt.type == PROD) {
        for (USI p = 0; p < numPerf; p++) {
            if (perf[p].state == OPEN) {
                OCP_USI k = perf[p].location;

                for (USI j = 0; j < numPhase; j++) {
                    myBulk.cfl[k * numPhase + j] += fabs(perf[p].qj_ft3[j]) * dt;
                }
            }
        }
    }
}

void Well::MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;
        for (USI i = 0; i < numCom; i++) {
            myBulk.Ni[k * numCom + i] -= perf[p].qi_lbmol[i] * dt;
        }
    }
}

void Well::AssembleMatINJ_IMPEC(const Bulk&    myBulk,
                                LinearSystem&  myLS,
                                const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.AddDim(1) - 1;
    myLS.NewDiag(wId, 0.0);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;

        OCP_DBL Vfi_zi = 0;
        for (USI i = 0; i < numCom; i++) {
            Vfi_zi += myBulk.vfi[k * numCom + i] * opt.injZi[i];
        }

        OCP_DBL valw = dt * perf[p].xi * perf[p].transINJ;
        OCP_DBL bw   = valw * dG[p];
        OCP_DBL valb = valw * Vfi_zi;
        OCP_DBL bb   = valb * dG[p];

        // Bulk to Well
        myLS.AddDiag(k, valb);
        myLS.NewOffDiag(k, wId, -valb);

        myLS.AddRhs(k, bb);

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                myLS.AddDiag(wId, valw);
                myLS.NewOffDiag(wId, k, -valw);
                myLS.AddRhs(wId, -bw);
                break;
            case BHP_MODE:
                myLS.NewOffDiag(wId, k, 0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            myLS.AddRhs(wId, dt * opt.maxRate);
            break;
        case BHP_MODE:
            myLS.AddDiag(wId, dt);
            myLS.AddRhs(wId, dt * opt.maxBHP);
            myLS.AssignGuess(wId, opt.maxBHP);     
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void Well::AssembleMatPROD_IMPEC(const Bulk&    myBulk,
                                 LinearSystem&  myLS,
                                 const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.AddDim(1) - 1;
    myLS.NewDiag(wId, 0.0);

    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;

        OCP_DBL valb = 0;
        OCP_DBL bb   = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw   = 0;

        for (USI j = 0; j < numPhase; j++) {
            if (!myBulk.phaseExist[n * numPhase + j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < numCom; i++) {
                tempb += myBulk.vfi[n * numCom + i] * myBulk.xij[n * numPhase * numCom + j * numCom + i];
                tempw += prodWeight[i] * myBulk.xij[n * numPhase * numCom + j * numCom + i];
            }
            OCP_DBL trans = dt * perf[p].transj[j] * myBulk.xi[n * numPhase + j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = dG[p] - myBulk.Pc[n * numPhase + j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        myLS.AddDiag(n, valb);
        myLS.NewOffDiag(n, wId, -valb);
        myLS.AddRhs(n, bb);

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                myLS.AddDiag(wId, -valw);
                myLS.NewOffDiag(wId, n, valw);
                myLS.AddRhs(wId, bw);
                break;
            case BHP_MODE:
                myLS.NewOffDiag(wId, n, 0.0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            myLS.AddRhs(wId, dt * opt.maxRate);
            break;
        case BHP_MODE:
            myLS.AddDiag(wId, dt);
            myLS.AddRhs(wId, dt * opt.minBHP);
            myLS.AssignGuess(wId, opt.minBHP);
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void Well::AssembleMatReinjection_IMPEC(const Bulk&         myBulk,
                                        LinearSystem&       myLS,
                                        const OCP_DBL&      dt,
                                        const vector<Well>& allWell,
                                        const vector<USI>&  injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (USI w = 0; w < injId.size(); w++) {
        if (allWell[w].IsOpen() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wOId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        const OCP_DBL factor = allWell[injId[0]].opt.reInjFactor * dt;
        // cout << "Factor(assemble):   " << allWell[injId[0]].opt.factor << endl;
        const OCP_USI prodId = wOId + myBulk.numBulk;
        OCP_USI       n, bId;
        OCP_DBL       tmp, valb;
        OCP_DBL       valw = 0;
        OCP_DBL       rhsw = 0;
        OCP_USI       tar;
        USI           tarsize;
        for (USI p = 0; p < numPerf; p++) {
            n    = perf[p].location;
            valb = 0;

            for (USI j = 0; j < numPhase; j++) {
                bId = n * numPhase + j;
                if (myBulk.phaseExist[bId]) {
                    tmp = perf[p].transj[j] * myBulk.xi[bId];
                    valb += tmp;
                    rhsw += tmp * (myBulk.Pc[bId] - dG[p]);
                }
            }
            valb *= factor;
            // Insert bulk val into equations in ascending order
            for (USI t = 0; t < tlen; t++) {
                tar     = tarId[t];
                tarsize = myLS.colId[tar].size();
                USI i;
                for (i = 0; i < tarsize; i++) {
                    if (n < myLS.colId[tar][i]) {
                        // insert
                        // attention that bulk id is less than well id
                        myLS.colId[tar].insert(myLS.colId[tar].begin() + i, n);
                        myLS.val[tar].insert(myLS.val[tar].begin() + i, -valb);
                        myLS.diagPtr[tar]++;
                        break;
                    }
                }
            }
            valw += valb;
        }
        rhsw *= factor;
        // insert prod well var and rhs into equations in ascending order
        for (USI t = 0; t < tlen; t++) {
            tar     = tarId[t];
            tarsize = myLS.colId[tar].size();
            myLS.b[tar] += rhsw; // rhsw
            USI i;
            for (i = 0; i < tarsize; i++) {
                if (prodId < myLS.colId[tar][i]) {
                    // insert
                    myLS.colId[tar].insert(myLS.colId[tar].begin() + i, prodId);
                    myLS.val[tar].insert(myLS.val[tar].begin() + i, valw);
                    if (i <= myLS.diagPtr[tar]) myLS.diagPtr[tar]++;
                    break;
                }
            }
            if (i == tarsize) {
                // pushback
                myLS.colId[tar].push_back(prodId);
                myLS.val[tar].push_back(valw);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////
// FIM
/////////////////////////////////////////////////////////////////////

void Well::AssembleMatINJ_FIM(const Bulk&    myBulk,
                              LinearSystem&  myLS,
                              const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol   = numCom + 1;
    const USI ncol2  = numPhase * numCom + numPhase;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    OCP_DBL mu, muP;
    OCP_DBL dP;
    OCP_DBL transIJ;

    OCP_USI n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = myBulk.P[n] - BHP - dG[p];

        for (USI j = 0; j < numPhase; j++) {
            n_np_j = n * numPhase + j;
            if (!myBulk.phaseExist[n_np_j]) continue;

            mu  = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < numPhase; k++) {
                    dQdXsB[(i + 1) * ncol2 + k] +=
                        CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                        opt.injZi[i] * myBulk.dKr_dS[n_np_j * numPhase + k] * dP / mu;
                }
                // dQ / dxij
                for (USI k = 0; k < numCom; k++) {
                    dQdXsB[(i + 1) * ncol2 + numPhase + j * numCom + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * numCom + k];
                }
            }
        }
        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * bsize2], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * bsize2], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.maxBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);

    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);
}

void Well::AssembleMatPROD_FIM(const Bulk&    myBulk,
                               LinearSystem&  myLS,
                               const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol   = numCom + 1;
    const USI ncol2  = numPhase * numCom + numPhase;
    const USI bsize  = ncol * ncol;
    const USI bsize2 = ncol * ncol2;

    OCP_DBL xij, xi, mu, muP, xiP;
    OCP_DBL dP;
    OCP_DBL transIJ;
    OCP_DBL tmp;

    OCP_USI n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);

    // Set Prod Weight   
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        for (USI j = 0; j < numPhase; j++) {
            n_np_j = n * numPhase + j;
            if (!myBulk.phaseExist[n_np_j]) continue;

            dP  = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi  = myBulk.xi[n_np_j];
            mu  = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];
            xiP = myBulk.xiP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                xij = myBulk.xij[n_np_j * numCom + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < numPhase; k++) {
                    tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi * xij *
                          myBulk.dKr_dS[n_np_j * numPhase + k];
                    // capillary pressure
                    tmp += transIJ * myBulk.dPcj_dS[n_np_j * numPhase + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dCij
                for (USI k = 0; k < numCom; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                          (myBulk.xix[n_np_j * numCom + k] -
                           xi / mu * myBulk.mux[n_np_j * numCom + k]);
                    dQdXsB[(i + 1) * ncol2 + numPhase + j * numCom + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + numPhase + j * numCom + i] +=
                    perf[p].transj[j] * xi * dP;
            }
        }
        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * bsize2], 1,
                bmat.data());

        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * prodWeight[i];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * bsize2], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.minBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);
    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);

}

void Well::AssembleMatReinjection_FIM(const Bulk&         myBulk,
                                      LinearSystem&       myLS,
                                      const OCP_DBL&      dt,
                                      const vector<Well>& allWell,
                                      const vector<USI>&  injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (USI w = 0; w < injId.size(); w++) {
        if (allWell[w].IsOpen() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wOId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        OCP_USI       tar;
        USI           tarsize;
        const OCP_DBL factor = allWell[injId[0]].opt.reInjFactor;
        const OCP_USI prodId = wOId + myBulk.numBulk;

        // cout << "Factor(assemble):    " << factor << endl;

        const USI ncol   = numCom + 1;
        const USI ncol2  = numPhase * numCom + numPhase;
        const USI bsize  = ncol * ncol;
        const USI bsize2 = ncol * ncol2;

        OCP_DBL xij, xi, mu, muP, xiP, dP, transIJ, tmp;
        OCP_USI n_np_j;

        vector<OCP_DBL> bmat(bsize, 0);
        vector<OCP_DBL> bmat2(bsize, 0);
        vector<OCP_DBL> tmpMat(bsize, 0);
        vector<OCP_DBL> dQdXpB(bsize, 0);
        vector<OCP_DBL> dQdXpW(bsize, 0);
        vector<OCP_DBL> dQdXsB(bsize2, 0);

        for (USI p = 0; p < numPerf; p++) {
            OCP_USI n = perf[p].location;
            fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
            fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
            fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

            for (USI j = 0; j < numPhase; j++) {
                n_np_j = n * numPhase + j;
                if (!myBulk.phaseExist[n_np_j]) continue;

                dP  = myBulk.Pj[n_np_j] - BHP - dG[p];
                xi  = myBulk.xi[n_np_j];
                mu  = myBulk.mu[n_np_j];
                muP = myBulk.muP[n_np_j];
                xiP = myBulk.xiP[n_np_j];

                for (USI i = 0; i < numCom; i++) {
                    xij = myBulk.xij[n_np_j * numCom + i];
                    // dQ / dP
                    transIJ = perf[p].transj[j] * xi * xij;
                    dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                              dP * perf[p].transj[j] * xij * xiP;
                    dQdXpW[(i + 1) * ncol] += -transIJ;

                    // dQ / dS
                    for (USI k = 0; k < numPhase; k++) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                              xij * myBulk.dKr_dS[n * numPhase * numPhase + j * numPhase + k];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n * numPhase * numPhase + j * numPhase + k];
                        dQdXsB[(i + 1) * ncol2 + k] += tmp;
                    }
                    // dQ / dCij
                    for (USI k = 0; k < numCom; k++) {
                        tmp = dP * perf[p].transj[j] * xij *
                              (myBulk.xix[n_np_j * numCom + k] -
                               xi / mu * myBulk.mux[n_np_j * numCom + k]);
                        if (k == i) {
                            tmp += perf[p].transj[j] * xi * dP;
                        }
                        dQdXsB[(i + 1) * ncol2 + numPhase + j * numCom + k] += tmp;
                    }
                }
            }

            // for Prod Well, be careful!
            for (USI i = 0; i < numCom; i++) {
                // tmpMat[0] -= dQdXpW[(i + 1) * ncol] * factor;
                tmpMat[0] += dQdXpW[(i + 1) * ncol] * factor;
            }

            // for perf(bulk) of Prod Well
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * bsize2],
                    1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < numCom; i++) {
                // becareful '-' before factor
                // Daxpy(ncol, -factor, bmat.data() + (i + 1) * ncol, bmat2.data());
                Daxpy(ncol, factor, bmat.data() + (i + 1) * ncol, bmat2.data());
            }

            // Insert bulk val into equations in ascending order
            for (USI t = 0; t < tlen; t++) {
                tar     = tarId[t];
                tarsize = myLS.colId[tar].size();
                USI i;
                for (i = 0; i < tarsize; i++) {
                    if (n < myLS.colId[tar][i]) {
                        // insert
                        // attention that bulk id is less than well id
                        myLS.colId[tar].insert(myLS.colId[tar].begin() + i, n);
                        myLS.val[tar].insert(myLS.val[tar].begin() + i * bsize,
                                             bmat.begin(), bmat.end());
                        myLS.diagPtr[tar]++;
                        break;
                    }
                }
            }
        }
        // insert prod well var into equations in ascending order
        for (USI t = 0; t < tlen; t++) {
            tar     = tarId[t];
            tarsize = myLS.colId[tar].size();
            USI i;
            for (i = 0; i < tarsize; i++) {
                if (prodId < myLS.colId[tar][i]) {
                    // insert
                    myLS.colId[tar].insert(myLS.colId[tar].begin() + i, prodId);
                    myLS.val[tar].insert(myLS.val[tar].begin() + i * bsize,
                                         tmpMat.begin(), tmpMat.end());
                    if (i <= myLS.diagPtr[tar]) myLS.diagPtr[tar]++;
                    break;
                }
            }
            if (i == tarsize) {
                // insert into end
                myLS.colId[tar].push_back(prodId);
                myLS.val[tar].insert(myLS.val[tar].end(), tmpMat.begin(), tmpMat.end());
            }
        }
    }
}

void Well::CalResFIM(OCPRes& resFIM, const Bulk& myBulk, const OCP_DBL& dt,
    const OCP_USI& wId, const vector<Well>& allWell) const
{
    OCP_FUNCNAME;

    // Well to Bulk
    const USI len = numCom + 1;
    OCP_USI   k;

    for (USI p = 0; p < numPerf; p++) {
        k = perf[p].location;
        for (USI i = 0; i < numCom; i++) {
            resFIM.res[k * len + 1 + i] += perf[p].qi_lbmol[i] * dt;
        }
    }

    // Well Self
    OCP_USI bId = (myBulk.numBulk + wId) * len;
    if (opt.type == INJ) {
        // Injection
        switch (opt.optMode) {
            case BHP_MODE:
                BHP             = opt.maxBHP;
                resFIM.res[bId] = BHP - opt.maxBHP;
                break;
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                resFIM.res[bId] = opt.maxRate;
                for (USI i = 0; i < numCom; i++) {
                    resFIM.res[bId] += qi_lbmol[i];
                }
                if (opt.reInj) {
                    for (auto& w : opt.connWell) {
                        OCP_DBL tmp = 0;
                        for (USI i = 0; i < numCom; i++) {
                            tmp += allWell[w].qi_lbmol[i];
                        }
                        tmp *= opt.reInjFactor;
                        resFIM.res[bId] += tmp;
                    }
                }
                resFIM.maxWellRelRes_mol =
                    max(resFIM.maxWellRelRes_mol, fabs(resFIM.res[bId] / opt.maxRate));
                break;
            default:
                OCP_ABORT("Wrong well opt mode!");
                break;
        }
    } else {
        // Production
        switch (opt.optMode) {
            case BHP_MODE:
                BHP             = opt.minBHP;
                resFIM.res[bId] = BHP - opt.minBHP;
                break;
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                CalProdWeight(myBulk);
                resFIM.res[bId] = -opt.maxRate;
                for (USI i = 0; i < numCom; i++) {
                    resFIM.res[bId] += qi_lbmol[i] * prodWeight[i];
                }
                resFIM.maxWellRelRes_mol =
                    max(resFIM.maxWellRelRes_mol, fabs(resFIM.res[bId] / opt.maxRate));
                break;
            default:
                OCP_ABORT("Wrong well opt mode!");
                break;
        }
    }
}


/////////////////////////////////////////////////////////////////////
// FIM(new)
/////////////////////////////////////////////////////////////////////

void Well::AssembleMatINJ_FIM_new(const Bulk&    myBulk,
                                  LinearSystem&  myLS,
                                  const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol    = numCom + 1;
    const USI ncol2   = numPhase * numCom + numPhase;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    OCP_DBL mu, muP;
    OCP_DBL dP;
    OCP_DBL transIJ;

    OCP_USI n_np_j;

    vector<OCP_DBL>  bmat(bsize, 0);
    vector<OCP_DBL>  bmat2(bsize, 0);
    vector<OCP_DBL>  dQdXpB(bsize, 0);
    vector<OCP_DBL>  dQdXpW(bsize, 0);
    vector<OCP_DBL>  dQdXsB(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(numPhase, OCP_FALSE);
    vector<OCP_BOOL> phasedS_B(numPhase, OCP_FALSE);
    vector<USI>      pVnumComB(numPhase, 0);
    USI              ncolB;

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = myBulk.P[n] - BHP - dG[p];

        USI jxB = 0;
        ncolB   = 0;

        for (USI j = 0; j < numPhase; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * numPhase + j];
            phasedS_B[j]   = myBulk.pSderExist[n * numPhase + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * numPhase + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < numPhase; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * numPhase + j;
            mu     = myBulk.mu[n_np_j];
            muP    = myBulk.muP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < numPhase; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.injZi[i] * myBulk.dKr_dS[n_np_j * numPhase + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * numCom + k];
                }
            }
            jxB += pVnumComB[j];
        }
        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * lendSdP], 1,
                bmat.data());

        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * lendSdP], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.maxBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);
    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);
}

void Well::AssembleMatPROD_FIM_new(const Bulk&    myBulk,
                                   LinearSystem&  myLS,
                                   const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol    = numCom + 1;
    const USI ncol2   = numPhase * numCom + numPhase;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    OCP_DBL xij, xi, mu, muP, xiP;
    OCP_DBL dP;
    OCP_DBL transIJ;
    OCP_DBL tmp;

    OCP_USI n_np_j;

    vector<OCP_DBL>  bmat(bsize, 0);
    vector<OCP_DBL>  bmat2(bsize, 0);
    vector<OCP_DBL>  dQdXpB(bsize, 0);
    vector<OCP_DBL>  dQdXpW(bsize, 0);
    vector<OCP_DBL>  dQdXsB(bsize2, 0);
    vector<OCP_BOOL> phaseExistB(numPhase, OCP_FALSE);
    vector<OCP_BOOL> phasedS_B(numPhase, OCP_FALSE);
    vector<USI>      pVnumComB(numPhase, 0);
    USI              ncolB;

    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        USI jxB = 0;
        ncolB   = 0;
        for (USI j = 0; j < numPhase; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * numPhase + j];
            phasedS_B[j]   = myBulk.pSderExist[n * numPhase + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * numPhase + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < numPhase; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * numPhase + j;
            dP     = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi     = myBulk.xi[n_np_j];
            mu     = myBulk.mu[n_np_j];
            muP    = myBulk.muP[n_np_j];
            xiP    = myBulk.xiP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                xij = myBulk.xij[n_np_j * numCom + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < numPhase; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                              xij * myBulk.dKr_dS[n_np_j * numPhase + j1];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n_np_j * numPhase + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                          (myBulk.xix[n_np_j * numCom + k] -
                           xi / mu * myBulk.mux[n_np_j * numCom + k]);
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] += perf[p].transj[j] * xi * dP;
                ;
            }
            jxB += pVnumComB[j];
        }
        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * prodWeight[i];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * lendSdP], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.minBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);
    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);
}

void Well::AssembleMatINJ_FIM_new_n(const Bulk&    myBulk,
                                    LinearSystem&  myLS,
                                    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol    = numCom + 1;
    const USI ncol2   = numPhase * numCom + numPhase;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    OCP_DBL mu, muP;
    OCP_DBL dP;
    OCP_DBL transIJ;

    OCP_USI n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);
    vector<char>    phaseExistB(numPhase, OCP_FALSE);
    vector<USI>     pVnumComB(numPhase, 0);
    vector<char>    phasedS_B(numPhase, OCP_FALSE);
    USI             ncolB;

    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = myBulk.P[n] - BHP - dG[p];

        const USI npB = myBulk.phaseNum[n];
        USI       jxB = 0;
        ncolB         = 0;
        for (USI j = 0; j < numPhase; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * numPhase + j];
            phasedS_B[j]   = myBulk.pSderExist[n * numPhase + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * numPhase + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < numPhase; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * numPhase + j;
            mu     = myBulk.mu[n_np_j];
            muP    = myBulk.muP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.injZi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < numPhase; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.injZi[i] * myBulk.dKr_dS[n_np_j * numPhase + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * numCom + k];
                }
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(), &myBulk.res_n[myBulk.resIndex[n]],
                    1.0, &myLS.b[n * ncol]);
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * lendSdP], 1,
                bmat.data());
        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * lendSdP], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.maxBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);
    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);
}

void Well::AssembleMatPROD_FIM_new_n(const Bulk&    myBulk,
                                     LinearSystem&  myLS,
                                     const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI ncol    = numCom + 1;
    const USI ncol2   = numPhase * numCom + numPhase;
    const USI bsize   = ncol * ncol;
    const USI bsize2  = ncol * ncol2;
    const USI lendSdP = myBulk.maxLendSdP;

    OCP_DBL xij, xi, mu, muP, xiP;
    OCP_DBL dP;
    OCP_DBL transIJ;
    OCP_DBL tmp;

    OCP_USI n_np_j;

    vector<OCP_DBL> bmat(bsize, 0);
    vector<OCP_DBL> bmat2(bsize, 0);
    vector<OCP_DBL> dQdXpB(bsize, 0);
    vector<OCP_DBL> dQdXpW(bsize, 0);
    vector<OCP_DBL> dQdXsB(bsize2, 0);
    vector<char>    phaseExistB(numPhase, OCP_FALSE);
    vector<char>    phasedS_B(numPhase, OCP_FALSE);
    vector<USI>     pVnumComB(numPhase, 0);
    USI             ncolB;

    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        const USI npB = myBulk.phaseNum[n];
        USI       jxB = 0;
        ncolB         = 0;
        for (USI j = 0; j < numPhase; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * numPhase + j];
            phasedS_B[j]   = myBulk.pSderExist[n * numPhase + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * numPhase + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < numPhase; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * numPhase + j;
            dP     = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi     = myBulk.xi[n_np_j];
            mu     = myBulk.mu[n_np_j];
            muP    = myBulk.muP[n_np_j];
            xiP    = myBulk.xiP[n_np_j];

            for (USI i = 0; i < numCom; i++) {
                xij = myBulk.xij[n_np_j * numCom + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < numPhase; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                              xij * myBulk.dKr_dS[n_np_j * numPhase + j1];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n_np_j * numPhase + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                          (myBulk.xix[n_np_j * numCom + k] -
                           xi / mu * myBulk.mux[n_np_j * numCom + k]);
                    tmp -= transIJ * dP / myBulk.nj[n_np_j];
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] +=
                        perf[p].transj[j] * xi * dP / myBulk.nj[n_np_j];
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(), &myBulk.res_n[myBulk.resIndex[n]],
                    1.0, &myLS.b[n * ncol]);
        }

        // Bulk to Well
        bmat = dQdXpB;
        DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * lendSdP], 1,
                bmat.data());

        Dscalar(bsize, dt, bmat.data());
        // Add
        USI ptr = myLS.diagPtr[n];
        for (USI i = 0; i < bsize; i++) {
            myLS.val[n][ptr * bsize + i] += bmat[i];
        }
        // Insert
        bmat = dQdXpW;
        Dscalar(bsize, dt, bmat.data());
        myLS.val[n].insert(myLS.val[n].end(), bmat.begin(), bmat.end());
        myLS.colId[n].push_back(wId);

        // Well
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    bmat[0] += dQdXpW[(i + 1) * ncol] * prodWeight[i];
                    bmat[(i + 1) * ncol + i + 1] = 1;
                }
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }

                // OffDiag
                bmat = dQdXpB;
                DaABpbC(ncol, ncol, ncolB, 1, dQdXsB.data(),
                        &myBulk.dSec_dPri[n * lendSdP], 1, bmat.data());
                fill(bmat2.begin(), bmat2.end(), 0.0);
                for (USI i = 0; i < numCom; i++) {
                    Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < numCom + 1; i++) {
                    bmat[i * ncol + i] = 1;
                }
                // Add
                for (USI i = 0; i < bsize; i++) {
                    myLS.diagVal[wId * bsize + i] += bmat[i];
                }
                // OffDiag
                fill(bmat.begin(), bmat.end(), 0.0);
                // Insert
                myLS.val[wId].insert(myLS.val[wId].end(), bmat.begin(), bmat.end());
                myLS.colId[wId].push_back(n);
                // Solution
                // myLS.u[wId * ncol] = opt.minBHP - BHP;
                break;

            default:
                OCP_ABORT("Wrong Well Opt mode!");
                break;
        }
    }
    assert(myLS.val[wId].size() == numPerf * bsize);
    // Well self
    myLS.colId[wId].push_back(wId);
    myLS.diagPtr[wId] = numPerf;
    myLS.val[wId].insert(myLS.val[wId].end(), myLS.diagVal.data() + wId * bsize,
                         myLS.diagVal.data() + wId * bsize + bsize);
}

void Well::SetPolyhedronWell(const Grid& myGrid, OCPpolyhedron& mypol)
{
    // set a virtual point
    mypol.numPoints = numPerf + 1;
    mypol.Points.resize(mypol.numPoints);

    OCP_USI k;
    Point3D tmpP;

    for (USI p = 0; p < numPerf; p++) {
        tmpP.Reset();
        k = myGrid.map_Act2All[perf[p].location];
        for (USI i = 0; i < myGrid.polyhedronGrid[k].numPoints; i++) {
            tmpP += myGrid.polyhedronGrid[k].Points[i];
        }
        tmpP /= myGrid.polyhedronGrid[k].numPoints;
        mypol.Points[p + 1] = tmpP;
    }

    // Set virtual perf
    mypol.Points[0] = mypol.Points[1];
    mypol.Points[0].z /= 2;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/