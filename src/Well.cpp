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
        OCP_ABORT("Wrong well type!");
    }

    if (type == INJ) {
        fluidType = Optparam.fluidType;
        if (fluidType == "WAT" || fluidType == "WATER") {
            fluidType = "WAT";
        }
    }

    if (Optparam.state == "OPEN") {
        state = OPEN;
    } else if (Optparam.state == "CLOSE") {
        state = CLOSE;
    } else {
        OCP_ABORT("Wrong state type!");
    }

    if (Optparam.optMode == "RATE") {
        optMode = RATE_MODE;
    } else if (Optparam.optMode == "ORAT") {
        optMode = ORATE_MODE;
    } else if (Optparam.optMode == "GRAT") {
        optMode = GRATE_MODE;
    } else if (Optparam.optMode == "WRAT") {
        optMode = WRATE_MODE;
    } else if (Optparam.optMode == "LRAT") {
        optMode = LRATE_MODE;
    } else if (Optparam.optMode == "BHP") {
        optMode = BHP_MODE;
    } else {
        OCP_ABORT("Wrong well option mode!");
    }

    initOptMode = optMode;
    maxRate     = Optparam.maxRate;
    maxBHP      = Optparam.maxBHP;
    minBHP      = Optparam.minBHP;
}


OCP_BOOL WellOpt::operator !=(const WellOpt& Opt) const
{
    if (this->type != Opt.type)                 return OCP_TRUE;
    if (this->state != Opt.state)               return OCP_TRUE;
    if (this->optMode != Opt.optMode)           return OCP_TRUE;
    if (this->initOptMode != Opt.initOptMode)   return OCP_TRUE;
    if (this->maxRate != Opt.maxRate)           return OCP_TRUE;
    if (this->maxBHP != Opt.maxBHP)             return OCP_TRUE;
    if (this->minBHP != Opt.minBHP)             return OCP_TRUE;
    for (USI i = 0; i < zi.size(); i++) {
        if (fabs(zi[i] - Opt.zi[i]) > TINY)     return OCP_TRUE;
    }
    return OCP_FALSE;
}


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
    qi_lbmol.resize(myBulk.numCom);
    prodWeight.resize(myBulk.numCom);
    factor.resize(3); // oil, gas, liquid
    Mtype = myBulk.flashCal[0]->GetType();
    // zi
    if (myBulk.blackOil) {
        for (auto& opt : optSet) {

            if (!opt.state) continue;

            opt.zi.resize(myBulk.numCom, 0);
            if (opt.type == INJ) {
                // INJ
                switch (myBulk.PVTmode) {
                    case PHASE_W:
                    case PHASE_OW:
                        opt.zi.back() = 1;
                        break;
                    case PHASE_ODGW:
                    case PHASE_DOGW:
                        if (opt.fluidType == "GAS")
                            opt.zi[1] = 1;
                        else
                            opt.zi[2] = 1;
                        break;
                    default:
                        OCP_ABORT("Wrong blackoil type!");
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
                        else if (opt.optMode == WRATE_MODE)
                            opt.zi[1] = 1;
                        else // LRATE_MODE
                            opt.zi[0] = opt.zi[1] = 1;
                        break;
                    case PHASE_DOGW:
                    case PHASE_ODGW:                   
                        if (opt.optMode == ORATE_MODE)
                            opt.zi[0] = 1;
                        else if (opt.optMode == GRATE_MODE)
                            opt.zi[1] = 1;
                        else if (opt.optMode == WRATE_MODE)
                            opt.zi[2] = 1;
                        else if (opt.optMode == LRATE_MODE)
                            opt.zi[2] = opt.zi[0] = 1;
                        break;
                    default:
                        OCP_ABORT("Wrong blackoil type!");
                }
            }
        }
    } else if (myBulk.comps) {

        USI len = sols.size();

        for (auto& opt : optSet) {

            if (!opt.state) continue;

            if (opt.type == INJ) {
                // INJ Well
                if (opt.fluidType == "WAT") {
                    opt.zi.resize(myBulk.numCom, 0);
                    opt.zi.back() = 1;
                } else {
                    for (USI i = 0; i < len; i++) {
                        if (opt.fluidType == sols[i].name) {
                            opt.zi = sols[i].data;
                            opt.zi.resize(myBulk.numCom);
                            // Convert volume units Mscf/stb to molar units lbmoles for
                            // injfluid Use flash in Bulk in surface condition
                            opt.xiINJ = myBulk.flashCal[0]->XiPhase(
                                PRESSURE_STD, TEMPERATURE_STD, &opt.zi[0]);
                            opt.maxRate *=
                                (opt.xiINJ *
                                 1000); // lbmol / ft3 -> lbmol / Mscf for gas
                            break;
                        }
                        if (i == len - 1) {
                            OCP_ABORT("Wrong FluidType!");
                        }
                    }
                }
            }
            // else {
            //     // PROD Well use EoS
            // }
        }
    } else {
        OCP_ABORT("Wrong mixture type!");
    }

    // perf
    USI pp = 0;
    for (USI p = 0; p < numPerf; p++) {
        OCP_USI Idg =
            perf[p].K * myGrid.nx * myGrid.ny + perf[p].J * myGrid.nx + perf[p].I;
        if (myGrid.activeMap_G2B[Idg].IsAct()) {

            perf[pp]            = perf[p];
            perf[pp].state      = OPEN;
            perf[pp].location   = myGrid.activeMap_G2B[Idg].GetId();
            perf[pp].depth      = myBulk.depth[perf[pp].location];
            perf[pp].multiplier = 1;
            perf[pp].qi_lbmol.resize(myBulk.numCom);
            perf[pp].transj.resize(myBulk.numPhase);
            perf[pp].qj_ft3.resize(myBulk.numPhase);
            pp++;
        } else {
            OCP_WARNING("Perforation is in inactive bulk!");
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
            OCP_USI Idb = perf[p].location;
            OCP_DBL dx  = myBulk.dx[Idb];
            OCP_DBL dy  = myBulk.dy[Idb];
            OCP_DBL dz  = myBulk.dz[Idb] * myBulk.ntg[Idb];
            OCP_DBL ro  = 0;
            switch (perf[p].direction) {
                case X_DIRECTION: {
                    OCP_DBL kykz  = myBulk.rockKy[Idb] * myBulk.rockKz[Idb];
                    OCP_DBL ky_kz = myBulk.rockKy[Idb] / myBulk.rockKz[Idb];
                    assert(kykz > 0);
                    ro =
                        0.28 *
                        pow((dy * dy * pow(1 / ky_kz, 0.5) + dz * dz * pow(ky_kz, 0.5)),
                            0.5);
                    ro /= (pow(ky_kz, 0.25) + pow(1 / ky_kz, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].kh = (dx * pow(kykz, 0.5));
                    }
                    break;
                }
                case Y_DIRECTION: {
                    OCP_DBL kzkx  = myBulk.rockKz[Idb] * myBulk.rockKx[Idb];
                    OCP_DBL kz_kx = myBulk.rockKz[Idb] / myBulk.rockKx[Idb];
                    assert(kzkx > 0);
                    ro =
                        0.28 *
                        pow((dz * dz * pow(1 / kz_kx, 0.5) + dx * dx * pow(kz_kx, 0.5)),
                            0.5);
                    ro /= (pow(kz_kx, 0.25) + pow(1 / kz_kx, 0.25));

                    if (perf[p].kh < 0) {
                        perf[p].kh = (dy * pow(kzkx, 0.5));
                    }
                    break;
                }
                case Z_DIRECTION: {
                    OCP_DBL kxky  = myBulk.rockKx[Idb] * myBulk.rockKy[Idb];
                    OCP_DBL kx_ky = myBulk.rockKx[Idb] / myBulk.rockKy[Idb];
                    assert(kxky > 0);
                    ro =
                        0.28 *
                        pow((dx * dx * pow(1 / kx_ky, 0.5) + dy * dy * pow(kx_ky, 0.5)),
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

    const USI np = myBulk.numPhase;

    if (opt.type == INJ) {
        for (USI p = 0; p < numPerf; p++) {
            perf[p].transINJ = 0;
            OCP_USI k        = perf[p].location;
            OCP_DBL temp     = CONV1 * perf[p].WI * perf[p].multiplier;

            // single phase
            for (USI j = 0; j < np; j++) {
                perf[p].transj[j] = 0;
                OCP_USI id        = k * np + j;
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
            for (USI j = 0; j < np; j++) {
                perf[p].transj[j] = 0;
                OCP_USI id        = k * np + j;
                if (myBulk.phaseExist[id]) {
                    perf[p].transj[j] = temp * myBulk.kr[id] / myBulk.mu[id];
                }
            }
        }
    }
    // check Trans
    // if (opt.type == PROD) {

    //    USI count = 0;
    //    for (USI p = 0; p < numPerf; p++) {
    //        if (perf[p].transj[0] != 0) {
    //            count++;
    //            break;
    //        }
    //    }
    //    if (count == 0) {
    //        cout << name << endl;
    //        for (USI p = 0; p < numPerf; p++) {
    //            OCP_USI k = perf[p].location;
    //            cout << "perf " << p << "  " << perf[p].multiplier << "   " <<
    //                myBulk.S[k * np] << "   " << myBulk.kr[k * np] << "   " <<
    //                myBulk.S[k * np + 1] << "   " << myBulk.kr[k * np + 1] << "   " <<
    //                myBulk.S[k * np + 2] << "   " << myBulk.kr[k * np + 2] << "   " <<
    //                endl;
    //        }

    //    }
    //}
}

void Well::CalFlux(const Bulk& myBulk, const OCP_BOOL flag)
{
    OCP_FUNCNAME;

    // cout << name << endl;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    fill(qi_lbmol.begin(), qi_lbmol.end(), 0.0);

    if (opt.type == INJ) {

        for (USI p = 0; p < numPerf; p++) {
            perf[p].P  = BHP + dG[p];
            OCP_USI k  = perf[p].location;
            OCP_DBL dP = perf[p].P - myBulk.P[k];
            dP *= -1.0;

            perf[p].qt_ft3 = perf[p].transINJ * dP;

            if (flag) {
                USI pvtnum = myBulk.PVTNUM[k];
                perf[p].xi =
                    myBulk.flashCal[pvtnum]->XiPhase(myBulk.P[k], myBulk.T, &opt.zi[0]);
            }
            // cout << "xi:  " << setprecision(6) << perf[p].xi << endl;
            OCP_DBL xi = perf[p].xi;
            for (USI i = 0; i < nc; i++) {
                perf[p].qi_lbmol[i] = perf[p].qt_ft3 * xi * opt.zi[i];
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

            for (USI j = 0; j < np; j++) {
                OCP_USI id = k * np + j;
                if (myBulk.phaseExist[id]) {
                    OCP_DBL dP = myBulk.Pj[id] - perf[p].P;

                    perf[p].qj_ft3[j] = perf[p].transj[j] * dP;
                    perf[p].qt_ft3 += perf[p].qj_ft3[j];
                    // cout << p << " p[" << j << "] = " << myBulk.Pj[id] << endl;
                    // cout << p << " perf = " << perf[p].p << endl;

                    OCP_DBL xi = myBulk.xi[id];
                    OCP_DBL xij;
                    for (USI i = 0; i < nc; i++) {
                        xij = myBulk.xij[id * nc + i];
                        perf[p].qi_lbmol[i] += perf[p].qj_ft3[j] * xi * xij;
                    }
                }
            }

            // check if perf[p].qi_lbmol is zero vector
            // OCP_BOOL flag = OCP_FALSE;
            // for (USI i = 0; i < nc; i++) {
            //    if (perf[p].qi_lbmol[i] != 0) {
            //        flag = OCP_TRUE;
            //        break;
            //    }
            //}
            // if (!flag) {
            //    OCP_ABORT("Qi is all zero!");
            //}

            for (USI i = 0; i < nc; i++) qi_lbmol[i] += perf[p].qi_lbmol[i];
        }
    }
    // test
    // cout << name << "----" << endl;
    // vector<OCP_DBL> tmpNiP(qi_lbmol);
    // OCP_DBL qt = Dnorm1(nc, &qi_lbmol[0]);
    // Dscalar(myBulk.numCom, 1 / qt * 100, &tmpNiP[0]);
    // cout << qt << "   ";
    // PrintDX(myBulk.numCom, &tmpNiP[0]);
}

/// Pressure in injection well equals maximum ones in injection well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL Well::CalInjRate(const Bulk& myBulk, const OCP_BOOL& maxBHP)
{
    OCP_FUNCNAME;

    OCP_DBL qj = 0;
    OCP_DBL Pwell = maxBHP ? opt.maxBHP : BHP;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        OCP_DBL dP = Pperf - myBulk.P[k];
        qj += perf[p].transINJ * perf[p].xi * dP;
    }
    return qj;
}

/// Pressure in production well equals minial ones in production well,
/// which is input by users. this function is used to check if operation mode of
/// well shoubld be swtched.
OCP_DBL Well::CalProdRate(const Bulk& myBulk, const OCP_BOOL& minBHP)
{
    OCP_FUNCNAME;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    OCP_DBL   qj = 0;
    OCP_DBL Pwell = minBHP ? opt.minBHP : BHP;

    for (USI p = 0; p < numPerf; p++) {

        OCP_DBL Pperf = Pwell + dG[p];
        OCP_USI k     = perf[p].location;

        for (USI j = 0; j < np; j++) {
            OCP_USI id = k * np + j;
            if (myBulk.phaseExist[id]) {
                OCP_DBL temp = 0;
                for (USI i = 0; i < nc; i++) {
                    temp += prodWeight[i] * myBulk.xij[id * nc + i];
                }
                OCP_DBL xi = myBulk.xi[id];
                OCP_DBL dP = myBulk.Pj[id] - Pperf;
                qj += perf[p].transj[j] * xi * dP * temp;
            }
        }
    }
    return qj;
}

void Well::CalInjQi(const Bulk& myBulk, const OCP_DBL& dt)
{
    OCP_FUNCNAME;

    const USI nc = myBulk.numCom;
    OCP_DBL   qj = 0;

    for (USI i = 0; i < nc; i++) {
        qj += qi_lbmol[i];
    }
    if (opt.fluidType == "WAT") {
        WWIR = -qj;
        WWIT += WWIR * dt;
    } else {
        if (Mtype == BLKOIL) {
            WGIR = -qj;
        } else {
            // EoS PVTW
            WGIR = -qj / (opt.xiINJ * 1000); // lbmol -> Mscf
        }
        WGIT += WGIR * dt;
    }
}

void Well::CalProdQj(const Bulk& myBulk, const OCP_DBL& dt)
{
    if (Mtype == BLKOIL) {
        CalProdQiBO(myBulk);
    } else {
        // EoS PVTW
        CalProdQjCOMP(myBulk);
    }
    WOPT += WOPR * dt;
    WGPT += WGPR * dt;
    WWPT += WWPR * dt;
}

void Well::CalProdQjCOMP(const Bulk& myBulk)
{
    USI     nc = myBulk.numCom;
    OCP_DBL qt = Dnorm1(nc, &qi_lbmol[0]);
    WOPR       = qt * factor[0];
    WGPR       = qt * factor[1];
    WWPR       = qi_lbmol[nc - 1];
}

void Well::CalProdQiBO(const Bulk& myBulk)
{
    OCP_FUNCNAME;

    const USI nc = myBulk.numCom;

    for (USI i = 0; i < nc; i++) {
        if (myBulk.index2Phase[i] == OIL) {
            WOPR = qi_lbmol[i];
        } else if (myBulk.index2Phase[i] == GAS) {
            WGPR = qi_lbmol[i];
        } else if (myBulk.index2Phase[i] == WATER) {
            WWPR = qi_lbmol[i];
        }
    }
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
        CalProddG(myBulk);

    // if (name == "PROD17") {
    //    for (USI p = 0; p < numPerf; p++) {
    //        cout << perf[p].state << "   ";
    //        cout << setprecision(2);
    //        cout << dG[p] << "   ";
    //        OCP_USI n = perf[p].location * myBulk.numPhase;
    //        cout << setprecision(6);
    //        cout << myBulk.kr[n + 0] << "   ";
    //        cout << myBulk.kr[n + 1] << "   ";
    //        cout << myBulk.kr[n + 2] << "   ";
    //        cout << myBulk.S[n + 0] << "   ";
    //        cout << myBulk.S[n + 1] << "   ";
    //        cout << myBulk.S[n + 2] << "   ";
    //        cout << setprecision(2);
    //        for (USI i = 0; i < myBulk.numCom; i++) {
    //            cout << perf[p].qi_lbmol[i] << "   ";
    //        }
    //        cout << endl;
    //    }
    //}
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
    OCP_FUNCNAME;

    const USI       np      = myBulk.numPhase;
    const USI       nc      = myBulk.numCom;
    const OCP_DBL   maxlen  = 5;
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
            // tmpNi.assign(nc, 0);
            // for (OCP_INT p1 = numPerf - 1; p1 >= p; p1--) {
            //    for (USI i = 0; i < nc; i++) {
            //        tmpNi[i] += perf[p1].qi_lbmol[i];
            //    }
            //}

            for (USI i = 0; i < nc; i++) {
                tmpNi[i] += perf[p].qi_lbmol[i];
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T, tmpNi.data(), 0, 0, 0);
                for (USI j = 0; j < myBulk.numPhase; j++) {
                    if (myBulk.flashCal[pvtnum]->phaseExist[j]) {
                        rhotmp = myBulk.flashCal[pvtnum]->rho[j];
                        qtacc += myBulk.flashCal[pvtnum]->v[j] / seg_num;
                        rhoacc += myBulk.flashCal[pvtnum]->v[j] * rhotmp *
                                  GRAVITY_FACTOR / seg_num;
#ifdef _DEBUG
                        if (rhotmp <= 0 || !isfinite(rhotmp)) {
                            OCP_ABORT("Wrong rho " + to_string(rhotmp));
                        }
#endif // _DEBUG
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
                for (USI i = 0; i < nc; i++) {
                    tmpNi[i] += perf[p1].qi_lbmol[i];
                }
            }

            // check tmpNi
            for (auto& v : tmpNi) {
                v = fabs(v);
            }

            for (USI k = 0; k < seg_num; k++) {
                myBulk.flashCal[pvtnum]->Flash(Ptmp, myBulk.T, tmpNi.data(), 0, 0, 0);
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
        OCP_ABORT("Wrong well position!");
    }
}

void Well::CalProdWeight(const Bulk& myBulk) const
{
    if (opt.type == PROD) {
        if (Mtype == BLKOIL || opt.optMode == WRATE_MODE) {
            prodWeight = opt.zi;
        } else if (Mtype == EOS_PVTW) {
            // use qi_lbmol, so it must be calculated before (Well::CalFlux())
            // when reset the step, qi_lbmol may be calculated
            USI pid = 0; // phase Id for oil
            if (opt.optMode == GRATE_MODE) {
                pid = 1; // phase Id for gas
            } else if (opt.optMode == LRATE_MODE) {
                pid = 2; // phase Id for liquid
            }
            // in some cases, qi_lbmol may be zero, so use other methods
            OCP_DBL qt = 0;
            for (USI i = 0; i < myBulk.numCom; i++) {
                qt += qi_lbmol[i];
            }
            if (fabs(qt) > TINY && qt > 0) {

                /*cout << name << endl;
                vector<OCP_DBL> tmpNiP(qi_lbmol);
                Dscalar(myBulk.numCom, 1 / qt * 100, &tmpNiP[0]);
                PrintDX(myBulk.numCom, &tmpNiP[0]);*/

                myBulk.flashCal[0]->Flash(PRESSURE_STD, TEMPERATURE_STD, &qi_lbmol[0],
                                          0, 0, 0);
            } else {
                USI             np = myBulk.numPhase;
                USI             nc = myBulk.numCom;
                vector<OCP_DBL> tmpNi(nc, 0);
                for (USI p = 0; p < numPerf; p++) {
                    OCP_USI n = perf[p].location;

                    for (USI j = 0; j < np; j++) {
                        OCP_USI id = n * np + j;
                        if (!myBulk.phaseExist[id]) continue;
                        for (USI k = 0; k < nc; k++) {
                            tmpNi[k] += perf[p].transj[j] * myBulk.xi[id] *
                                        myBulk.xij[id * nc + k];
                        }
                    }
                }
                qt = Dnorm1(nc, &tmpNi[0]);

                /*cout << name << endl;
                vector<OCP_DBL> tmpNiP(tmpNi);
                Dscalar(nc, 1 / qt, &tmpNiP[0]);
                PrintDX(nc, &tmpNiP[0]);*/

                myBulk.flashCal[0]->Flash(PRESSURE_STD, TEMPERATURE_STD, &tmpNi[0], 0,
                                          0, 0);
            }
            // test
            // USI nc = myBulk.numCom;
            // USI np = myBulk.numPhase;
            // cout << "OIL : ";
            // for (USI i = 0; i < nc; i++) {
            //    cout << myBulk.flashCal[0]->xij[i] * 100 << "   ";
            //}
            // cout << endl << "GAS : ";
            // for (USI i = 0; i < nc; i++) {
            //    cout << myBulk.flashCal[0]->xij[nc + i] * 100 << "   ";
            //}
            // cout << endl;

            fill(factor.begin(), factor.end(), 0);
            if (myBulk.flashCal[0]->phaseExist[0]) { // improve
                OCP_DBL nuo = myBulk.flashCal[0]->xi[0] * myBulk.flashCal[0]->v[0] /
                              qt; // mole fraction of oil phase
                factor[0] = nuo / (myBulk.flashCal[0]->xi[0] *
                                   CONV1); // lbmol / ft3 -> lbmol / stb for oil
            }
            if (myBulk.flashCal[0]->phaseExist[1]) {
                OCP_DBL nug = myBulk.flashCal[0]->xi[1] * myBulk.flashCal[0]->v[1] /
                              qt; // mole fraction of gas phase
                factor[1] = nug / (myBulk.flashCal[0]->xi[1] *
                                   1000); // lbmol / ft3 -> lbmol / Mscf  for gas
            }
            factor[2] = factor[0];
            if (myBulk.flashCal[0]->phaseExist[2]) {
                OCP_DBL nuw = myBulk.flashCal[0]->xi[2] * myBulk.flashCal[0]->v[2] /
                              qt; // mole fraction of water phase
                factor[2] += nuw;
            }

            if (factor[pid] < 1E-12 || !isfinite(factor[pid])) {
                OCP_ABORT("Wrong Condition!");
            }
            fill(prodWeight.begin(), prodWeight.end(), factor[pid]);

            // cout << "Factor  " << factor[0] << "   " << factor[1] << endl;
        }
    }
}

void Well::CalReInjFluid(const Bulk& myBulk, vector<OCP_DBL>& myZi)
{
    CalTrans(myBulk);

    USI np = myBulk.numPhase;
    USI nc = myBulk.numCom;
    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;

        for (USI j = 0; j < np; j++) {
            OCP_USI id = n * np + j;
            if (!myBulk.phaseExist[id]) continue;
            for (USI k = 0; k < nc; k++) {
                myZi[k] += perf[p].transj[j] * myBulk.xi[id] * myBulk.xij[id * nc + k];
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
                if (opt.injPhase == GAS)
                    tarRate = WGIR;
                else if (opt.injPhase == WATER)
                    tarRate = WWIR;
            }
            if (q > tarRate) {
                opt.optMode = RATE_MODE;
            }
            else {
                opt.optMode = BHP_MODE;
                BHP = opt.maxBHP;
            }
        }
        else {
            opt.optMode = BHP_MODE;
            if (opt.type == INJ)
                BHP = opt.maxBHP;
            else
                BHP = opt.minBHP;
        }
        //else {
        //    OCP_DBL q = CalProdRate(myBulk, OCP_FALSE);
        //    // cout << q << endl;
        //    if (q > opt.maxRate) {
        //        opt.optMode = opt.initOptMode;
        //    }
        //    else {
        //        opt.optMode = BHP_MODE;
        //        BHP = opt.minBHP;
        //    }
        //}
    } else {
        if (opt.type == INJ) {
            OCP_DBL q = CalInjRate(myBulk, OCP_TRUE);
            // for INJ well, maxRate has been switch to lbmols
            OCP_DBL tarRate = opt.maxRate;
            if (opt.reInj) {
                if (opt.injPhase == GAS)
                    tarRate = WGIR;
                else if (opt.injPhase == WATER)
                    tarRate = WWIR;
            }

            if (q > tarRate) {
                opt.optMode = opt.initOptMode;
            }
            else {
                opt.optMode = BHP_MODE;
                BHP = opt.maxBHP;
            }
        }
        else {
            OCP_DBL q = CalProdRate(myBulk, OCP_TRUE);
            // cout << q << endl;
            if (q > opt.maxRate) {
                opt.optMode = opt.initOptMode;
            }
            else {
                opt.optMode = BHP_MODE;
                BHP = opt.minBHP;
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
#ifdef _DEBUG
            cout << "### WARNING: Well " << name << " Perf[" << p << "].P = " << perf[p].P
                 << endl;
#endif // _DEBUG
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

    OCP_USI k;
    OCP_BOOL    flagC = OCP_TRUE;

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

void Well::AllocateMat(LinearSystem& myLS) const
{
    OCP_FUNCNAME;

    for (USI p = 0; p < numPerf; p++) {
        myLS.rowCapacity[perf[p].location]++;
    }
}

void Well::SetupWellBulk(Bulk& myBulk) const
{
    // Attention that a bulk can only be penetrated by one well now!
    for (USI p = 0; p < numPerf; p++) {
        myBulk.wellBulkId.push_back(perf[p].location);
    }
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
        OCP_USI         n  = perf[p].location;
        cout << setw(3) << p << "   " << perf[p].state << "   " << setw(6)
            << perf[p].location << "  " << setw(2) << perf[p].I + 1 << "  " << setw(2)
            << perf[p].J + 1 << "  " << setw(2) << perf[p].K + 1 << "  " << setw(10)
            << setprecision(6) << perf[p].WI << "  "               // ccf
            << setprecision(3) << perf[p].radius << "  "               // ccf
            << setw(8) << setprecision(4) << perf[p].kh << "  "    // kh
            << setw(8) << setprecision(2) << perf[p].depth << "  " // depth
            << setprecision(3) << perf[p].P << "  "                // Pp
            << setw(10) << setprecision(3) << myBulk.P[n] << "   " // Pb
            << setw(6) << setprecision(3) << dG[p] << "   "        // dG
            << setw(8) << perf[p].qi_lbmol[myBulk.numCom - 1] << "   " << setw(6) << setprecision(6)
            << myBulk.S[n * myBulk.numPhase + 0] << "   " << setw(6) << setprecision(6)
            << myBulk.S[n * myBulk.numPhase + 1] << "   " << setw(6) << setprecision(6)
            << myBulk.S[n * myBulk.numPhase + 2] << endl;
    }
}

/////////////////////////////////////////////////////////////////////
// IMPEC
/////////////////////////////////////////////////////////////////////


void Well::CalCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    if (opt.type == PROD) {
        const USI np = myBulk.numPhase;
        for (USI p = 0; p < numPerf; p++) {
            if (perf[p].state == OPEN) {
                OCP_USI k = perf[p].location;

                for (USI j = 0; j < np; j++) {
                    myBulk.cfl[k * np + j] += fabs(perf[p].qj_ft3[j]) * dt;
                }
            }
        }
    }
}

void Well::MassConserveIMPEC(Bulk& myBulk, const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI nc = myBulk.numCom;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;
        for (USI i = 0; i < nc; i++) {
            myBulk.Ni[k * nc + i] -= perf[p].qi_lbmol[i] * dt;
        }
    }
}

void Well::AssembleMatINJ_IMPEC(const Bulk& myBulk, LinearSystem& myLS,
                                const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI     nc  = myBulk.numCom;
    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI k = perf[p].location;

        OCP_DBL Vfi_zi = 0;
        for (USI i = 0; i < nc; i++) {
            Vfi_zi += myBulk.vfi[k * nc + i] * opt.zi[i];
        }

        OCP_DBL valw = dt * perf[p].xi * perf[p].transINJ;
        OCP_DBL bw   = valw * dG[p];
        OCP_DBL valb = valw * Vfi_zi;
        OCP_DBL bb   = valb * dG[p];

        // Bulk to Well

        // diag
        USI ptr = myLS.diagPtr[k];
        myLS.val[k][ptr] += valb;
        // off diag
        myLS.colId[k].push_back(wId);
        myLS.val[k].push_back(-valb);
        // b
        myLS.b[k] += bb;

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // diag
                myLS.diagVal[wId] += valw;
                // off diag
                myLS.colId[wId].push_back(k);
                myLS.val[wId].push_back(-valw);
                // b
                myLS.b[wId] -= bw;
                break;
            case BHP_MODE:
                myLS.colId[wId].push_back(k);
                myLS.val[wId].push_back(0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    assert(myLS.val[wId].size() == numPerf);
    // the order of perforation is not necessarily in order
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            // diag
            myLS.colId[wId].push_back(wId);
            myLS.diagPtr[wId] = numPerf;
            myLS.val[wId].push_back(myLS.diagVal[wId]);
            // b
            myLS.b[wId] += dt * opt.maxRate;
            break;
        case BHP_MODE:
            // diag
            myLS.colId[wId].push_back(wId);
            myLS.diagPtr[wId] = numPerf;
            myLS.val[wId].push_back(dt);
            // b
            myLS.b[wId] += dt * opt.maxBHP;
            // u   initial value
            myLS.u[wId] = opt.maxBHP;
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void Well::AssembleMatPROD_IMPEC(const Bulk& myBulk, LinearSystem& myLS,
                                 const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const USI     np  = myBulk.numPhase;
    const USI     nc  = myBulk.numCom;
    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;

        OCP_DBL valb = 0;
        OCP_DBL bb   = 0;
        OCP_DBL valw = 0;
        OCP_DBL bw   = 0;

        for (USI j = 0; j < np; j++) {
            if (!myBulk.phaseExist[n * np + j]) continue;

            OCP_DBL tempb = 0;
            OCP_DBL tempw = 0;

            for (USI i = 0; i < nc; i++) {
                tempb += myBulk.vfi[n * nc + i] * myBulk.xij[n * np * nc + j * nc + i];
                tempw += prodWeight[i] * myBulk.xij[n * np * nc + j * nc + i];
            }
            OCP_DBL trans = dt * perf[p].transj[j] * myBulk.xi[n * np + j];
            valb += tempb * trans;
            valw += tempw * trans;

            OCP_DBL dP = dG[p] - myBulk.Pc[n * np + j];
            bb += tempb * trans * dP;
            bw += tempw * trans * dP;
        }

        // Bulk to Well
        // diag
        USI ptr = myLS.diagPtr[n];
        myLS.val[n][ptr] += valb;
        // off diag
        myLS.colId[n].push_back(wId);
        myLS.val[n].push_back(-valb);
        // b
        myLS.b[n] += bb;

        // Well to Bulk
        switch (opt.optMode) {
            case RATE_MODE:
            case ORATE_MODE:
            case GRATE_MODE:
            case WRATE_MODE:
            case LRATE_MODE:
                // diag  !!! attention! sign is -
                myLS.diagVal[wId] -= valw;
                // off diag
                myLS.colId[wId].push_back(n);
                myLS.val[wId].push_back(valw);
                // b
                myLS.b[wId] += bw;
                break;
            case BHP_MODE:
                // off diag
                myLS.colId[wId].push_back(n);
                myLS.val[wId].push_back(0);
                break;
            default:
                OCP_ABORT("Wrong well option mode!");
        }
    }

    // Well Self
    assert(myLS.val[wId].size() == numPerf);
    // the order of perforation is not necessarily in order
    switch (opt.optMode) {
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            // diag
            myLS.colId[wId].push_back(wId);
            myLS.diagPtr[wId] = numPerf;
            myLS.val[wId].push_back(myLS.diagVal[wId]);
            // b
            myLS.b[wId] += dt * opt.maxRate;
            break;
        case BHP_MODE:
            // diag
            myLS.colId[wId].push_back(wId);
            myLS.diagPtr[wId] = numPerf;
            myLS.val[wId].push_back(dt);
            // b
            myLS.b[wId] += dt * opt.minBHP;
            // u   initial value
            myLS.u[wId] = opt.minBHP;
            break;
        default:
            OCP_ABORT("Wrong well option mode!");
    }
}

void Well::AssembleMatReinjection_IMPEC(const Bulk& myBulk, LinearSystem& myLS,
                                        const OCP_DBL& dt, const vector<Well>& allWell,
                                        const vector<USI>& injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (USI w = 0; w < injId.size(); w++) {
        if (allWell[w].WellState() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wEId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        const OCP_DBL factor = allWell[injId[0]].opt.factor * dt;
        // cout << "Factor(assemble):   " << allWell[injId[0]].opt.factor << endl;
        const OCP_USI prodId = wEId + myBulk.numBulk;
        USI           np     = myBulk.numPhase;
        OCP_USI       n, bId;
        OCP_DBL       tmp, valb;
        OCP_DBL       valw = 0;
        OCP_DBL       rhsw = 0;
        OCP_USI       tar;
        USI           tarsize;
        for (USI p = 0; p < numPerf; p++) {
            n    = perf[p].location;
            valb = 0;

            for (USI j = 0; j < np; j++) {
                bId = n * np + j;
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

void Well::AssembleMatINJ_FIM(const Bulk& myBulk, LinearSystem& myLS,
                              const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np     = myBulk.numPhase;
    const USI nc     = myBulk.numCom;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
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

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!myBulk.phaseExist[n_np_j]) continue;

            mu  = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.zi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    dQdXsB[(i + 1) * ncol2 + k] +=
                        CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                        opt.zi[i] * myBulk.dKr_dS[n_np_j * np + k] * dP / mu;
                }
                // dQ / dxij
                for (USI k = 0; k < nc; k++) {
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * nc + k];
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
                for (USI i = 0; i < nc; i++) {
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
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc + 1; i++) {
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

void Well::AssembleMatPROD_FIM(const Bulk& myBulk, LinearSystem& myLS,
                               const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np     = myBulk.numPhase;
    const USI nc     = myBulk.numCom;
    const USI ncol   = nc + 1;
    const USI ncol2  = np * nc + np;
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

        for (USI j = 0; j < np; j++) {
            n_np_j = n * np + j;
            if (!myBulk.phaseExist[n_np_j]) continue;

            dP  = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi  = myBulk.xi[n_np_j];
            mu  = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];
            xiP = myBulk.xiP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // if (dQdXpW[ncol] == 0) {
                //    cout << name << " perf " << p << "  " << j << endl;
                //}

                // dQ / dS
                for (USI k = 0; k < np; k++) {
                    tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi * xij *
                          myBulk.dKr_dS[n_np_j * np + k];
                    // capillary pressure
                    tmp += transIJ * myBulk.dPcj_dS[n_np_j * np + k];
                    dQdXsB[(i + 1) * ncol2 + k] += tmp;
                }
                // dQ / dCij
                for (USI k = 0; k < nc; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                          (myBulk.xix[n_np_j * nc + k] -
                           xi / mu * myBulk.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                }
                dQdXsB[(i + 1) * ncol2 + np + j * nc + i] += perf[p].transj[j] * xi * dP;
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
                for (USI i = 0; i < nc; i++) {
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
                for (USI i = 0; i < nc; i++) {
                    Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                          bmat2.data());
                }
                myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
                myLS.colId[wId].push_back(n);
                break;

            case BHP_MODE:
                // Diag
                fill(bmat.begin(), bmat.end(), 0.0);
                for (USI i = 0; i < nc + 1; i++) {
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

    // test
    // for (USI i = 0; i < bsize; i++) {
    //    cout << myLS.diagVal[wId * bsize + i] << endl;
    //}
}

void Well::AssembleMatReinjection_FIM(const Bulk& myBulk, LinearSystem& myLS,
                                      const OCP_DBL& dt, const vector<Well>& allWell,
                                      const vector<USI>& injId) const
{
    // find Open injection well under Rate control
    vector<OCP_USI> tarId;
    for (USI w = 0; w < injId.size(); w++) {
        if (allWell[w].WellState() && allWell[w].opt.optMode != BHP_MODE)
            tarId.push_back(allWell[w].wEId + myBulk.numBulk);
    }

    USI tlen = tarId.size();
    if (tlen > 0) {
        OCP_USI       tar;
        USI           tarsize;
        const OCP_DBL factor = allWell[injId[0]].opt.factor;
        const OCP_USI prodId = wEId + myBulk.numBulk;

        // cout << "Factor(assemble):    " << factor << endl;

        const USI np     = myBulk.numPhase;
        const USI nc     = myBulk.numCom;
        const USI ncol   = nc + 1;
        const USI ncol2  = np * nc + np;
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

            for (USI j = 0; j < np; j++) {
                n_np_j = n * np + j;
                if (!myBulk.phaseExist[n_np_j]) continue;

                dP  = myBulk.Pj[n_np_j] - BHP - dG[p];
                xi  = myBulk.xi[n_np_j];
                mu  = myBulk.mu[n_np_j];
                muP = myBulk.muP[n_np_j];
                xiP = myBulk.xiP[n_np_j];

                for (USI i = 0; i < nc; i++) {
                    xij = myBulk.xij[n_np_j * nc + i];
                    // dQ / dP
                    transIJ = perf[p].transj[j] * xi * xij;
                    dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu) +
                                              dP * perf[p].transj[j] * xij * xiP;
                    dQdXpW[(i + 1) * ncol] += -transIJ;

                    // dQ / dS
                    for (USI k = 0; k < np; k++) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi *
                              xij * myBulk.dKr_dS[n * np * np + j * np + k];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n * np * np + j * np + k];
                        dQdXsB[(i + 1) * ncol2 + k] += tmp;
                    }
                    // dQ / dCij
                    for (USI k = 0; k < nc; k++) {
                        tmp = dP * perf[p].transj[j] * xij *
                              (myBulk.xix[n_np_j * nc + k] -
                               xi / mu * myBulk.mux[n_np_j * nc + k]);
                        if (k == i) {
                            tmp += perf[p].transj[j] * xi * dP;
                        }
                        dQdXsB[(i + 1) * ncol2 + np + j * nc + k] += tmp;
                    }
                }
            }

            // for Prod Well, be careful!
            for (USI i = 0; i < nc; i++) {
                // tmpMat[0] -= dQdXpW[(i + 1) * ncol] * factor;
                tmpMat[0] += dQdXpW[(i + 1) * ncol] * factor;
            }

            // for perf(bulk) of Prod Well
            bmat = dQdXpB;
            DaABpbC(ncol, ncol, ncol2, 1, dQdXsB.data(), &myBulk.dSec_dPri[n * bsize2],
                    1, bmat.data());
            fill(bmat2.begin(), bmat2.end(), 0.0);
            for (USI i = 0; i < nc; i++) {
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

void Well::CalResFIM(ResFIM& resFIM, const Bulk& myBulk, const OCP_DBL& dt,
                     const OCP_USI& wId, const vector<Well>& allWell) const
{
    OCP_FUNCNAME;

    // Well to Bulk
    const USI nc  = myBulk.numCom;
    const USI len = nc + 1;
    OCP_USI   k;

    for (USI p = 0; p < numPerf; p++) {
        k = perf[p].location;
        for (USI i = 0; i < nc; i++) {
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
                for (USI i = 0; i < nc; i++) {
                    resFIM.res[bId] += qi_lbmol[i];
                }
                if (opt.reInj) {
                    for (auto& w : opt.connWell) {
                        OCP_DBL tmp = 0;
                        for (USI i = 0; i < nc; i++) {
                            tmp += allWell[w].qi_lbmol[i];
                            // tmp += opt.factor * allWell[w].qi_lbmol[i];
                            // resFIM.res[bId] += opt.factor * allWell[w].qi_lbmol[i];
                        }
                        tmp *= opt.factor;
                        resFIM.res[bId] += tmp;
                        // cout << "Temp(INJ):    " << tmp / opt.xiINJ / 1000 << endl;
                    }
                    // cout << "Factor(res)    " << opt.factor << endl;
                }
                 //cout << name << "   " << resFIM.res[bId] << "   " << opt.maxRate << "   " <<
                 //    fabs(resFIM.res[bId] / opt.maxRate) << endl;
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
                resFIM.res[bId] = -opt.maxRate;
                for (USI i = 0; i < nc; i++) {
                    resFIM.res[bId] += qi_lbmol[i] * prodWeight[i];
                }
                // cout << "Temp(Prod):   " << tmp << endl;
                 //cout << name << "   " << resFIM.res[bId] << "   " << opt.maxRate << "   " 
                 //    << fabs(resFIM.res[bId] / opt.maxRate) << endl;
                resFIM.maxWellRelRes_mol =
                    max(resFIM.maxWellRelRes_mol, fabs(resFIM.res[bId] / opt.maxRate));
                break;
            default:
                OCP_ABORT("Wrong well opt mode!");
                break;
        }
    }
}

void Well::ShowRes(const OCP_USI& wId, const vector<OCP_DBL>& res, const Bulk& myBulk) const
{
    // Well to Bulk
    const USI nc = myBulk.numCom;
    const USI len = nc + 1;

    OCP_USI bId = (myBulk.numBulk + wId) * len;
    cout << name << "   " << res[bId] << "   ";
    // Well Self
    if (opt.type == INJ) {
        // Injection
        switch (opt.optMode) {
        case BHP_MODE:
            cout << "BHPMode" << endl;
            break;
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            cout << opt.maxRate << "   " << fabs(res[bId] / opt.maxRate) << endl;
            break;
        default:
            OCP_ABORT("Wrong well opt mode!");
            break;
        }
    }
    else {
        // Production
        switch (opt.optMode) {
        case BHP_MODE:
            cout << "BHPMode" << endl;
            break;
        case RATE_MODE:
        case ORATE_MODE:
        case GRATE_MODE:
        case WRATE_MODE:
        case LRATE_MODE:
            cout << opt.maxRate << "   " << fabs(res[bId] / opt.maxRate) << endl;
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


void Well::AssembleMatINJ_FIM_new(const Bulk& myBulk, LinearSystem& myLS,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
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
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    USI             ncolB;


    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = myBulk.P[n] - BHP - dG[p];

        USI jxB = 0;
        ncolB = 0;

        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * np + j];
            phasedS_B[j] = myBulk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            mu = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.zi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.zi[i] * myBulk.dKr_dS[n_np_j * np + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij               
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * nc + k];
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
            for (USI i = 0; i < nc; i++) {
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
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
            myLS.colId[wId].push_back(n);
            break;

        case BHP_MODE:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc + 1; i++) {
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


void Well::AssembleMatPROD_FIM_new(const Bulk& myBulk, LinearSystem& myLS,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
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
    vector<OCP_BOOL>    phaseExistB(np, OCP_FALSE);
    vector<OCP_BOOL>    phasedS_B(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    USI             ncolB;


    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        USI jxB = 0;
        ncolB = 0;  
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * np + j];
            phasedS_B[j] = myBulk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;


        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            dP = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi = myBulk.xi[n_np_j];
            mu = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];
            xiP = myBulk.xiP[n_np_j];


            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi * xij *
                            myBulk.dKr_dS[n_np_j * np + j1];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n_np_j * np + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (myBulk.xix[n_np_j * nc + k] -
                            xi / mu * myBulk.mux[n_np_j * nc + k]);
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] += perf[p].transj[j] * xi * dP;;
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
            for (USI i = 0; i < nc; i++) {
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
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
            myLS.colId[wId].push_back(n);
            break;

        case BHP_MODE:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc + 1; i++) {
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


void Well::AssembleMatINJ_FIM_new_n(const Bulk& myBulk, LinearSystem& myLS,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
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
    vector<char>    phaseExistB(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    vector<char>    phasedS_B(np, OCP_FALSE);
    USI             ncolB;


    for (USI p = 0; p < numPerf; p++) {
        const OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        dP = myBulk.P[n] - BHP - dG[p];

        const USI npB = myBulk.phaseNum[n] + 1;
        USI jxB = 0;
        ncolB = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * np + j];
            phasedS_B[j] = myBulk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            mu = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];

            for (USI i = 0; i < nc; i++) {
                // dQ / dP
                transIJ = perf[p].transj[j] * perf[p].xi * opt.zi[i];
                dQdXpB[(i + 1) * ncol] += transIJ * (1 - dP * muP / mu);
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        dQdXsB[(i + 1) * ncolB + j1B] +=
                            CONV1 * perf[p].WI * perf[p].multiplier * perf[p].xi *
                            opt.zi[i] * myBulk.dKr_dS[n_np_j * np + j1] * dP / mu;
                        j1B++;
                    }
                }

                // dQ / dxij               
                for (USI k = 0; k < pVnumComB[j]; k++) {
                    dQdXsB[(i + 1) * ncolB + jxB + k] +=
                        -transIJ * dP / mu * myBulk.mux[n_np_j * nc + k];
                }
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(),
                &myBulk.res_n[myBulk.resIndex[n]], 1.0, &myLS.b[n * ncol]);
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
            for (USI i = 0; i < nc; i++) {
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
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, 1.0, bmat.data() + (i + 1) * ncol, bmat2.data());
            }
            myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
            myLS.colId[wId].push_back(n);
            break;

        case BHP_MODE:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc + 1; i++) {
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


void Well::AssembleMatPROD_FIM_new_n(const Bulk& myBulk, LinearSystem& myLS,
    const OCP_DBL& dt) const
{
    OCP_FUNCNAME;

    const OCP_USI wId = myLS.dim;
    // important !
    myLS.dim++;

    const USI np = myBulk.numPhase;
    const USI nc = myBulk.numCom;
    const USI ncol = nc + 1;
    const USI ncol2 = np * nc + np;
    const USI bsize = ncol * ncol;
    const USI bsize2 = ncol * ncol2;
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
    vector<char>    phaseExistB(np, OCP_FALSE);
    vector<char>    phasedS_B(np, OCP_FALSE);
    vector<USI>     pVnumComB(np, 0);
    USI             ncolB;


    // Set Prod Weight
    if (opt.optMode != BHP_MODE) CalProdWeight(myBulk);

    for (USI p = 0; p < numPerf; p++) {
        OCP_USI n = perf[p].location;
        fill(dQdXpB.begin(), dQdXpB.end(), 0.0);
        fill(dQdXpW.begin(), dQdXpW.end(), 0.0);
        fill(dQdXsB.begin(), dQdXsB.end(), 0.0);

        const USI npB = myBulk.phaseNum[n] + 1;
        USI jxB = 0;
        ncolB = 0;
        for (USI j = 0; j < np; j++) {
            phaseExistB[j] = myBulk.phaseExist[n * np + j];
            phasedS_B[j] = myBulk.pSderExist[n * np + j];
            if (phasedS_B[j]) jxB++;
            pVnumComB[j] = myBulk.pVnumCom[n * np + j];
            ncolB += pVnumComB[j];
        }
        ncolB += jxB;

        for (USI j = 0; j < np; j++) {

            if (!phaseExistB[j]) {
                jxB += pVnumComB[j];
                continue;
            }

            n_np_j = n * np + j;
            dP = myBulk.Pj[n_np_j] - BHP - dG[p];
            xi = myBulk.xi[n_np_j];
            mu = myBulk.mu[n_np_j];
            muP = myBulk.muP[n_np_j];
            xiP = myBulk.xiP[n_np_j];


            for (USI i = 0; i < nc; i++) {
                xij = myBulk.xij[n_np_j * nc + i];
                // dQ / dP
                transIJ = perf[p].transj[j] * xi * xij;
                dQdXpB[(i + 1) * ncol] +=
                    transIJ * (1 - dP * muP / mu) + dP * perf[p].transj[j] * xij * xiP;
                dQdXpW[(i + 1) * ncol] += -transIJ;

                // dQ / dS
                USI j1B = 0;
                for (USI j1 = 0; j1 < np; j1++) {
                    if (phasedS_B[j1]) {
                        tmp = CONV1 * perf[p].WI * perf[p].multiplier * dP / mu * xi * xij *
                            myBulk.dKr_dS[n_np_j * np + j1];
                        // capillary pressure
                        tmp += transIJ * myBulk.dPcj_dS[n_np_j * np + j1];
                        dQdXsB[(i + 1) * ncolB + j1B] += tmp;
                        j1B++;
                    }
                }

                for (USI k = 0; k < pVnumComB[j]; k++) {
                    tmp = dP * perf[p].transj[j] * xij *
                        (myBulk.xix[n_np_j * nc + k] -
                            xi / mu * myBulk.mux[n_np_j * nc + k]);
                    tmp -= transIJ * dP / myBulk.nj[n_np_j];
                    dQdXsB[(i + 1) * ncolB + jxB + k] += tmp;
                }
                // WARNING !!!
                if (i < pVnumComB[j])
                    dQdXsB[(i + 1) * ncolB + jxB + i] += perf[p].transj[j] * xi * dP / myBulk.nj[n_np_j];
            }
            jxB += pVnumComB[j];
        }

        // Assemble rhs
        // Well
        if (npB > 2) {
            DaAxpby(ncol, ncolB, -1.0, dQdXsB.data(),
                &myBulk.res_n[myBulk.resIndex[n]], 1.0, &myLS.b[n * ncol]);
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
            for (USI i = 0; i < nc; i++) {
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
            for (USI i = 0; i < nc; i++) {
                Daxpy(ncol, prodWeight[i], bmat.data() + (i + 1) * ncol,
                    bmat2.data());
            }
            myLS.val[wId].insert(myLS.val[wId].end(), bmat2.begin(), bmat2.end());
            myLS.colId[wId].push_back(n);
            break;

        case BHP_MODE:
            // Diag
            fill(bmat.begin(), bmat.end(), 0.0);
            for (USI i = 0; i < nc + 1; i++) {
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
        k = myGrid.activeMap_B2G[perf[p].location];
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