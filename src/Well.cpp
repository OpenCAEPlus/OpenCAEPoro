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

void Well::Setup(const Grid& gd, const Bulk& bk, const vector<SolventINJ>& sols)
{
    OCP_FUNCNAME;

    numCom   = bk.numCom;
    numPhase = bk.numPhase;
    flashCal = bk.flashCal;

    qi_lbmol.resize(numCom);
    prodWeight.resize(numCom);
    prodRate.resize(numPhase);

    for (auto& opt : optSet) {
        if (!opt.state) continue;
        if (!bk.ifThermal) {
            opt.injTemp = bk.rsTemp;
        }
        flashCal[0]->SetupWellOpt(opt, sols, Psurf, Tsurf);
    }

    // Perf
    USI pp = 0;
    for (USI p = 0; p < numPerf; p++) {
        OCP_USI pId = perf[p].K * gd.nx * gd.ny + perf[p].J * gd.nx + perf[p].I;
        if (gd.map_All2Flu[pId].IsAct()) {

            perf[pp]            = perf[p];
            perf[pp].state      = OPEN;
            perf[pp].location   = gd.map_All2Act[pId].GetId();
            perf[pp].depth      = bk.depth[perf[pp].location];
            perf[pp].multiplier = 1;
            perf[pp].qi_lbmol.resize(numCom);
            perf[pp].transj.resize(numPhase);
            perf[pp].qj_ft3.resize(numPhase);
            pp++;
        } else {
            OCP_WARNING("Perforation is in non-fluid Bulk!");
        }
    }
    numPerf = pp;
    perf.resize(numPerf);
    // dG
    dG.resize(numPerf, 0);
    ldG = dG;

    if (depth < 0) depth = perf[0].depth;

    CalWI_Peaceman(bk);
    // test
    // ShowPerfStatus(bk);
}

void Well::CalWI_Peaceman(const Bulk& myBulk)
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
            OCP_DBL       ro  = 0;
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
            if (ifUseUnweight) {
                perf[p].transINJ = perf[p].WI;
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
            perf[p].P  = bhp + dG[p];
            OCP_USI k  = perf[p].location;
            OCP_DBL dP = myBulk.P[k] - perf[p].P;

            perf[p].qt_ft3 = perf[p].transINJ * dP;

            if (ReCalXi) {
                USI pvtnum = myBulk.PVTNUM[k];
                perf[p].xi = myBulk.flashCal[pvtnum]->XiPhase(
                    perf[p].P, opt.injTemp, &opt.injZi[0], opt.injProdPhase);
            }
            for (USI i = 0; i < numCom; i++) {
                perf[p].qi_lbmol[i] = perf[p].qt_ft3 * perf[p].xi * opt.injZi[i];
                qi_lbmol[i] += perf[p].qi_lbmol[i];
            }
        }
    } else {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P      = bhp + dG[p];
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
OCP_DBL Well::CalInjRateMaxBHP(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    OCP_DBL qj    = 0;
    OCP_DBL Pwell = opt.maxBHP;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        USI     pvtnum = myBulk.PVTNUM[k];
        OCP_DBL xi = myBulk.flashCal[pvtnum]->XiPhase(Pperf, opt.injTemp, &opt.injZi[0],
                                                      opt.injProdPhase);

        OCP_DBL dP = Pperf - myBulk.P[k];
        qj += perf[p].transINJ * xi * dP;
    }
    return qj;
}

/// Pressure in production well equals minial ones in production well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL Well::CalProdRateMinBHP(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    OCP_DBL qj    = 0;
    OCP_DBL Pwell = opt.minBHP;

    vector<OCP_DBL> tmpQi_lbmol(numCom, 0);
    vector<OCP_DBL> tmpQj(numPhase, 0);

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        for (USI j = 0; j < numPhase; j++) {
            OCP_USI id = k * numPhase + j;
            if (myBulk.phaseExist[id]) {
                OCP_DBL dP   = myBulk.Pj[id] - Pperf;
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

void Well::CalInjQj(const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    OCP_DBL qj = 0;

    for (USI i = 0; i < numCom; i++) {
        qj += qi_lbmol[i];
    }
    if (opt.fluidType == "WAT") {
        WWIR = -qj / opt.factorINJ;
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp -= myBulk.flashCal[pvtnum]->RhoPhase(
                            Ptmp, 0, opt.injTemp, opt.injZi.data(), opt.injProdPhase) *
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            USI pvtnum = myBulk.PVTNUM[n];
            for (USI i = 0; i < seg_num; i++) {
                Ptmp += myBulk.flashCal[pvtnum]->RhoPhase(
                            Ptmp, 0, opt.injTemp, opt.injZi.data(), opt.injProdPhase) *
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

// Use transj
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                const OCP_USI n_np_j = n * numPhase + j;
                if (!myBulk.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] += (myBulk.P[n] - perf[p].P) * perf[p].transj[j] *
                                myBulk.xi[n_np_j] * myBulk.xij[n_np_j * numCom + k];
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            // fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                const OCP_USI n_np_j = n * numPhase + j;
                if (!myBulk.phaseExist[n_np_j]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] += (myBulk.P[n] - perf[p].P) * perf[p].transj[j] *
                                myBulk.xi[n_np_j] * myBulk.xij[n_np_j * numCom + k];
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

// Use bulk
void Well::CalProddG02(const Bulk& myBulk)
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * myBulk.xi[id] *
                                myBulk.xij[id * numCom + k];
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
            perf[p].P     = bhp + dG[p];
            OCP_DBL Pperf = perf[p].P;
            OCP_DBL Ptmp  = Pperf;

            fill(tmpNi.begin(), tmpNi.end(), 0.0);
            for (USI j = 0; j < numPhase; j++) {
                OCP_USI id = n * numPhase + j;
                if (!myBulk.phaseExist[id]) continue;
                for (USI k = 0; k < numCom; k++) {
                    tmpNi[k] += (perf[p].transj[j] > 0) * myBulk.xi[id] *
                                myBulk.xij[id * numCom + k];
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

// Use qi_lbmol
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
            perf[p].P     = bhp + dG[p];
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
            perf[p].P     = bhp + dG[p];
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
            flashCal[0]->CalProdWeight(Psurf, Tsurf, &qi_lbmol[0], opt.prodPhaseWeight,
                                       prodWeight);
        } else {
            vector<OCP_DBL> tmpNi(numCom, 0);
            for (USI p = 0; p < numPerf; p++) {
                OCP_USI n = perf[p].location;

                for (USI j = 0; j < numPhase; j++) {
                    OCP_USI id = n * numPhase + j;
                    if (!myBulk.phaseExist[id]) continue;
                    for (USI k = 0; k < numCom; k++) {
                        tmpNi[k] += perf[p].transj[j] * myBulk.xi[id] *
                                    myBulk.xij[id * numCom + k];
                    }
                }
            }
            qt = Dnorm1(numCom, &tmpNi[0]);
            flashCal[0]->CalProdWeight(Psurf, Tsurf, &tmpNi[0], opt.prodPhaseWeight,
                                       prodWeight);
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
                myZi[k] +=
                    perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * numCom + k];
            }
        }
    }
}

void Well::CorrectBHP()
{
    if (opt.type == PROD && opt.optMode == BHP_MODE) {
        bhp = opt.minBHP;
    } else if (opt.type == INJ && opt.optMode == BHP_MODE) {
        bhp = opt.maxBHP;
    }
}

/// Constant well pressure would be applied if flow rate is too large.
/// Constant flow rate would be applied if well pressure is outranged.
void Well::CheckOptMode(const Bulk& myBulk)
{
    OCP_FUNCNAME;
    if (opt.initOptMode == BHP_MODE) {
        if (opt.type == INJ) {
            OCP_DBL q = CalInjRateMaxBHP(myBulk);
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
                bhp         = opt.maxBHP;
            }
        } else {
            opt.optMode = BHP_MODE;
            bhp         = opt.minBHP;
        }
    } else {
        if (opt.type == INJ) {
            OCP_DBL q = CalInjRateMaxBHP(myBulk);
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
                bhp         = opt.maxBHP;
            }
        } else {
            OCP_DBL q = CalProdRateMinBHP(myBulk);
            // cout << q << endl;
            if (q > opt.maxRate) {
                opt.optMode = opt.initOptMode;
            } else {
                opt.optMode = BHP_MODE;
                bhp         = opt.minBHP;
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

    if (bhp < 0) {
        cout << "### WARNING: Well " << name << " BHP = " << bhp << endl;
        return WELL_NEGATIVE_PRESSURE;
    }
    for (USI p = 0; p < numPerf; p++) {
        if (perf[p].state == OPEN && perf[p].P < 0) {
#ifdef DEBUG
            cout << "### WARNING: Well " << name << " Perf[" << p
                 << "].P = " << perf[p].P << endl;
#endif // DEBUG
            return WELL_NEGATIVE_PRESSURE;
        }
    }

    if (opt.type == INJ) {
        if (opt.optMode != BHP_MODE && bhp > opt.maxBHP) {
#if _DEBUG
            cout << "### WARNING: Well " << name << " switch to BHPMode" << endl;
#endif
            opt.optMode = BHP_MODE;
            bhp         = opt.maxBHP;
            return WELL_SWITCH_TO_BHPMODE;
        }
    } else {
        if (opt.optMode != BHP_MODE && bhp < opt.minBHP) {
#if _DEBUG
            cout << "### WARNING: Well " << name << " switch to BHPMode" << endl;
#endif
            opt.optMode = BHP_MODE;
            bhp         = opt.minBHP;
            return WELL_SWITCH_TO_BHPMODE;
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
        // the multiplier of perforation will be set to zero, so trans of well
        // should be recalculated!
        //
        // dG = ldG;
        CalTrans(myBulk);
        // CalFlux(myBulk);
        // CaldG(myBulk);
        // CheckOptMode(myBulk);
        return WELL_CROSSFLOW;
    }

    return WELL_SUCCESS;
}

void Well::ShowPerfStatus(const Bulk& myBulk) const
{
    OCP_FUNCNAME;

    cout << fixed;
    cout << "----------------------------" << endl;
    cout << name << ":    " << opt.optMode << "   " << setprecision(3) << bhp << endl;
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

void Well::AssembleMatReinjection_IMPEC(const Bulk&         myBulk,
                                        LinearSystem&       myLS,
                                        const OCP_DBL&      dt,
                                        const vector<Well>& allWell,
                                        const vector<USI>&  injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (auto& w : injId) {
        if (allWell[w].IsOpen() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wOId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        // All inj well has the same factor
        const OCP_DBL factor = allWell[injId[0]].opt.reInjFactor * dt;
        const OCP_USI prodId = wOId + myBulk.numBulk;
        OCP_USI       n, bId;
        OCP_DBL       tmp, valb;
        OCP_DBL       valw = 0;
        OCP_DBL       rhsw = 0;
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
            for (USI t = 0; t < tlen; t++) {
                myLS.NewOffDiag(tarId[t], n, -valb);
            }
            valw += valb;
        }
        rhsw *= factor;
        // rhs and prod well
        for (USI t = 0; t < tlen; t++) {
            myLS.AddRhs(tarId[t], rhsw);
            myLS.NewOffDiag(tarId[t], prodId, valw);
        }
    }
}

void Well::AssembleMatReinjection_FIM(const Bulk&         myBulk,
                                      LinearSystem&       myLS,
                                      const OCP_DBL&      dt,
                                      const vector<Well>& allWell,
                                      const vector<USI>&  injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (auto& w : injId) {
        if (allWell[w].IsOpen() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wOId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        // All inj well has the same factor
        const OCP_DBL factor = allWell[injId[0]].opt.reInjFactor;
        const OCP_USI prodId = wOId + myBulk.numBulk;

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
            const OCP_USI n = perf[p].location;
            fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
            fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
            fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

            for (USI j = 0; j < numPhase; j++) {
                n_np_j = n * numPhase + j;
                if (!myBulk.phaseExist[n_np_j]) continue;

                dP  = myBulk.Pj[n_np_j] - bhp - dG[p];
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
                              xij * myBulk.dKr_dS[n_np_j * numPhase + k];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n_np_j * numPhase + k];
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

            // for Prod Well
            for (USI i = 0; i < numCom; i++) {
                tmpMat[0] += dQdXpW[(i + 1) * ncol] * factor;
            }

            // for Perf(bulk) of Prod Well
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * bsize2],
                    1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < numCom; i++) {
                // becareful '-' before factor
                Daxpy(ncol, -factor, bmat.data() + (i + 1) * ncol, bmat2.data());
            }

            // Insert bulk val into equations in ascending order
            for (USI t = 0; t < tlen; t++) {
                myLS.NewOffDiag(tarId[t], n, bmat2);
            }
        }
        // prod well
        for (USI t = 0; t < tlen; t++) {
            myLS.NewOffDiag(tarId[t], prodId, tmpMat);
        }
    }
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