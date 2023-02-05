/*! \file    MixtureThermal_K.cpp
 *  \brief   MixtureThermal_K class declaration
 *  \author  Shizhe Li
 *  \date    Nov/10/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

// OpenCAEPoro header files
#include "MixtureThermal.hpp"
#include "OCPTable.hpp"

MixtureThermal_K01::MixtureThermal_K01(const ParamReservoir& param, const USI& tarId)
{
    mixtureType = THERMAL;
    numPhase    = 2;
    numCom      = param.numCom; // water is included
    Allocate();

    if (param.comsParam.molden.activity)
        xi_ref = param.comsParam.molden.data[tarId];
    else
        OCP_ABORT("ACF hasn't been input!");
    if (param.comsParam.avisc.activity || param.comsParam.bvisc.activity) {
        if (param.comsParam.avisc.activity)
            avisc = param.comsParam.avisc.data[tarId];
        else
            avisc.resize(numCom, 0);
        if (param.comsParam.bvisc.activity)
            bvisc = param.comsParam.bvisc.data[tarId];
        else
            bvisc.resize(numCom, 0);
        useViscTab = OCP_FALSE;
    } else {
        if (param.comsParam.viscTab.data.size() <= tarId) {
            OCP_ABORT("VISCTAB hasn't been input for " + to_string(tarId + 1) +
                      " th Region!");
        }
        useViscTab = OCP_TRUE;
        visc.Setup(param.comsParam.viscTab.data[tarId]);
        // unit convert: F -> R
        for (auto& v : visc.GetCol(0)) {
            v += CONV5;
        }
    }

    if (param.comsParam.cp.activity)
        cp = param.comsParam.cp.data[tarId];
    else
        cp.resize(numCom, 0);

    if (param.comsParam.ct1.activity)
        ct1 = param.comsParam.ct1.data[tarId];
    else
        ct1.resize(numCom, 0);

    if (param.comsParam.ct2.activity)
        ct2 = param.comsParam.ct2.data[tarId];
    else
        ct2.resize(numCom, 0);

    if (param.comsParam.cpt.activity)
        cpt = param.comsParam.cpt.data[tarId];
    else
        cpt.resize(numCom, 0);

    if (param.comsParam.cpl1.activity)
        cpl1 = param.comsParam.cpl1.data[tarId];
    else
        cpl1.resize(numCom, 0);

    if (param.comsParam.cpl2.activity)
        cpl2 = param.comsParam.cpl2.data[tarId];
    else
        cpl2.resize(numCom, 0);

    if (param.comsParam.cpl3.activity)
        cpl3 = param.comsParam.cpl3.data[tarId];
    else
        cpl3.resize(numCom, 0);

    if (param.comsParam.cpl4.activity)
        cpl4 = param.comsParam.cpl4.data[tarId];
    else
        cpl4.resize(numCom, 0);

    if (param.comsParam.cpg1.activity)
        cpg1 = param.comsParam.cpg1.data[tarId];
    else
        cpg1.resize(numCom, 0);

    if (param.comsParam.cpg2.activity)
        cpg2 = param.comsParam.cpg2.data[tarId];
    else
        cpg2.resize(numCom, 0);

    if (param.comsParam.cpg3.activity)
        cpg3 = param.comsParam.cpg3.data[tarId];
    else
        cpg3.resize(numCom, 0);

    if (param.comsParam.cpg4.activity)
        cpg4 = param.comsParam.cpg4.data[tarId];
    else
        cpg4.resize(numCom, 0);

    if (param.comsParam.hvapr.activity)
        hvapr = param.comsParam.hvapr.data[tarId];
    else
        hvapr.resize(numCom, 0);

    if (param.comsParam.hvr.activity)
        hvr = param.comsParam.hvr.data[tarId];
    else
        hvr.resize(numCom, 0);

    if (param.comsParam.ev.activity)
        ev = param.comsParam.ev.data[tarId];
    else
        ev.resize(numCom, 0);

    if (param.comsParam.Tc.activity)
        Tcrit = param.comsParam.Tc.data[tarId];
    else
        OCP_ABORT("TCRIT hasn't been input!");

    if (param.comsParam.MW.activity)
        MWc = param.comsParam.MW.data[tarId];
    else
        OCP_ABORT("MW hasn't been input!");

    Tref = param.comsParam.Tref[tarId] + CONV5;
    Pref = param.comsParam.Pref[tarId];

    dXsdXp.resize((numCom + 2) * (numPhase + numPhase * numCom));
    MWp.resize(numPhase);

    data.resize(3, 0);
    cdata.resize(3, 0);

    // Init
    phaseExist[0] = true;
    phaseExist[1] = true;

    MWp[0] = MWc[0];
    MWp[1] = MWc[1];

    xij[0 * 2 + 0] = 1;
    xij[0 * 2 + 1] = 0;
    xij[1 * 2 + 0] = 0;
    xij[1 * 2 + 1] = 1;

    // d mu / dP
    fill(rhox.begin(), rhox.end(), 0.0);
    fill(xix.begin(), xix.end(), 0.0);
    fill(muP.begin(), muP.end(), 0.0);
    fill(mux.begin(), mux.end(), 0.0);
}

void MixtureThermal_K01::Flash(const OCP_DBL& Pin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Niin)
{
    FlashIMPEC(Pin, Tin, Niin, 0, 0, 0);
}

void MixtureThermal_K01::InitFlashIMPEC(const OCP_DBL& Pin,
                                        const OCP_DBL& Pbbin,
                                        const OCP_DBL& Tin,
                                        const OCP_DBL* Sjin,
                                        const OCP_DBL& Vpore,
                                        const OCP_DBL* Ziin,
                                        const OCP_USI& bId)
{
    // assign value
    P                = Pin;
    T                = Tin + CONV5;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    S[1]             = Sjin[1];
    S[0]             = 1 - S[1];

    // phase molar density
    xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] / 2 * pow(dT, 2) + cpt[0] * dP * dT);
    xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] / 2 * pow(dT, 2) + cpt[1] * dP * dT);

    // components moles
    Ni[0] = Vpore * S[0] * xi[0];
    Ni[1] = Vpore * S[1] * xi[1];
}

void MixtureThermal_K01::InitFlashFIM(const OCP_DBL& Pin,
                                      const OCP_DBL& Pbbin,
                                      const OCP_DBL& Tin,
                                      const OCP_DBL* Sjin,
                                      const OCP_DBL& Vpore,
                                      const OCP_DBL* Ziin,
                                      const OCP_USI& bId)
{
    InitFlashIMPEC(Pin, Pbbin, Tin, Sjin, Vpore, Ziin, bId);

    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, bId);
}

void MixtureThermal_K01::FlashIMPEC(const OCP_DBL& Pin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Niin,
                                    const USI&     lastNP,
                                    const OCP_DBL* xijin,
                                    const OCP_USI& bId)
{
    // assign Value
    P                = Pin;
    T                = Tin + CONV5;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    Ni[0]            = Niin[0];
    Ni[1]            = Niin[1];

    // phase viscosity
    if (useViscTab) {
        visc.Eval_All(0, T, data, cdata);
        mu[0] = data[1];
        mu[1] = data[2];
    } else {
        mu[0] = avisc[0] * exp(bvisc[0] / T);
        mu[1] = avisc[1] * exp(bvisc[1] / T);
    }

    // phase molar density
    xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] / 2 * pow(dT, 2) + cpt[0] * dP * dT);
    xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] / 2 * pow(dT, 2) + cpt[1] * dP * dT);

    // phase mass density
    rho[0] = MWp[0] * xi[0];
    rho[1] = MWp[1] * xi[1];

    // phase volume
    vj[0] = Ni[0] / xi[0];
    vj[1] = Ni[1] / xi[1];
    // total fluid volume
    vf = vj[0] + vj[1];

    // phase saturation
    S[0] = vj[0] / vf;
    S[1] = vj[1] / vf;

    // d Vf / d Ni
    vfi[0] = 1 / xi[0];
    vfi[1] = 1 / xi[1];

    // d Vf / dP
    const OCP_DBL xiop = xi[0] * (cp[0] + cpt[0] * dT);
    const OCP_DBL xiwp = xi[1] * (cp[1] + cpt[1] * dT);
    vfP                = -(vj[0] * xiop / xi[0] + vj[1] * xiwp / xi[1]);

    // d Vf / dT
    const OCP_DBL xioT = xi[0] * (-ct1[0] - ct2[0] * dT + cpt[0] * dP);
    const OCP_DBL xiwT = xi[1] * (-ct1[1] - ct2[1] * dT + cpt[1] * dP);
    vfT                = -(vj[0] * xioT / xi[0] + vj[1] * xiwT / xi[1]);
}

void MixtureThermal_K01::FlashFIM(const OCP_DBL& Pin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Niin,
                                  const OCP_DBL* Sjin,
                                  const USI&     lastNP,
                                  const OCP_DBL* xijin,
                                  const OCP_USI& bId)
{
    // Assign value
    P                = Pin;
    T                = Tin + CONV5;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    Ni[0]            = Niin[0];
    Ni[1]            = Niin[1];
    Nt               = Ni[0] + Ni[1];

    // phase viscosity
    if (useViscTab) {
        visc.Eval_All(0, T, data, cdata);
        mu[0] = data[1];
        mu[1] = data[2];
        // d mu / dT
        muT[0] = cdata[1];
        muT[1] = cdata[2];

    } else {
        mu[0] = avisc[0] * exp(bvisc[0] / T);
        mu[1] = avisc[1] * exp(bvisc[1] / T);
        // d mu / dT
        muT[0] = mu[0] * (-bvisc[0] / (T * T));
        muT[1] = mu[1] * (-bvisc[1] / (T * T));
    }

    // phase molar density
    xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] * pow(dT, 2) / 2 + cpt[0] * dP * dT);
    xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] * pow(dT, 2) / 2 + cpt[1] * dP * dT);
    // d xi / dP
    xiP[0] = xi[0] * (cp[0] + cpt[0] * dT);
    xiP[1] = xi[1] * (cp[1] + cpt[1] * dT);
    // d xi / dT
    xiT[0] = xi[0] * (-ct1[0] - ct2[0] * dT + cpt[0] * dP);
    xiT[1] = xi[1] * (-ct1[1] - ct2[1] * dT + cpt[1] * dP);

    // phase mass density
    rho[0] = MWp[0] * xi[0];
    rho[1] = MWp[1] * xi[1];
    // d rho / d P
    rhoP[0] = MWp[0] * xiP[0];
    rhoP[1] = MWp[1] * xiP[1];
    // d rho / d T
    rhoT[0] = MWp[0] * xiT[0];
    rhoT[1] = MWp[1] * xiT[1];

    // phase volume
    vj[0] = Ni[0] / xi[0];
    vj[1] = Ni[1] / xi[1];
    // total volume
    vf = vj[0] + vj[1];

    // phase saturation
    S[0] = vj[0] / vf;
    S[1] = vj[1] / vf;

    // d vf/ d Ni
    vfi[0] = 1 / xi[0];
    vfi[1] = 1 / xi[1];
    // d vf / d P
    vfP = -(vj[0] * xiP[0] / xi[0] + vj[1] * xiP[1] / xi[1]);
    // d vf / d T
    vfT = -(vj[0] * xiT[0] / xi[0] + vj[1] * xiT[1] / xi[1]);

    // Derivative of secondary vars with respect to primary vars
    dXsdXp[0] = (-vj[0] * xiP[0] / xi[0] - S[0] * vfP) / vf; // dSo / dP
    dXsdXp[1] = (1 / xi[0] - S[0] * vfi[0]) / vf;            // dSo / dNo
    dXsdXp[2] = -S[0] * vfi[1] / vf;                         // dSo / dNw
    dXsdXp[3] = (-vj[0] * xiT[0] / xi[0] - S[0] * vfT) / vf; // dSo / dT

    dXsdXp[4] = -dXsdXp[0]; // dSw / dP
    dXsdXp[5] = -dXsdXp[1]; // dSw / dNo
    dXsdXp[6] = -dXsdXp[2]; // dSw / dNw
    dXsdXp[7] = -dXsdXp[3]; // dSw / dT

    CalEnthalpy();
}

void MixtureThermal_K01::CalEnthalpy()
{
    fill(H.begin(), H.end(), 0.0);
    fill(HT.begin(), HT.end(), 0.0);

    const OCP_DBL dT = T - Tref;

    if (liquid_based || simple_hvap) {
        for (USI j = 0; j < numPhase; j++) {
            for (USI i = 0; i < numCom; i++) {

                Hx[j * numCom + i] =
                    (cpl1[i] * dT + 1.0 / 2 * cpl2[i] * pow(dT, 2) +
                     1.0 / 3 * cpl3[i] * pow(dT, 3) + 1.0 / 4 * cpl4[i] * pow(dT, 4));

                H[j] += xij[j * numCom + i] * Hx[j * numCom + i];

                HT[j] +=
                    xij[j * numCom + i] * (cpl1[i] + cpl2[i] * dT +
                                           cpl3[i] * pow(dT, 2) + cpl4[i] * pow(dT, 3));
            }
        }
    } else if (gas_based) {
        for (USI j = 0; j < numPhase; j++) {
            for (USI i = 0; i < numCom; i++) {
                Hx[j * numCom + i] = cpg1[i] * dT + 1.0 / 2 * cpg2[i] * pow(dT, 2) +
                                     1.0 / 3 * cpg3[i] * pow(dT, 3) +
                                     1.0 / 4 * cpg4[i] * pow(dT, 4);
                H[j] += xij[j * numCom + i] * Hx[j * numCom + i];
                HT[j] +=
                    xij[j * numCom + i] * (cpg1[i] + cpg2[i] * dT +
                                           cpg3[i] * pow(dT, 2) + cpg4[i] * pow(dT, 3));

                if (T < Tcrit[i]) {
                    Hx[j * numCom + i] -= hvr[i] * pow((Tcrit[i] - T), ev[i]);

                    H[j] -= xij[j * numCom + i] * hvr[i] * pow((Tcrit[i] - T), ev[i]);

                    HT[j] += xij[j * numCom + i] * hvr[i] * ev[i] *
                             pow((Tcrit[i] - T), ev[i] - 1);
                }
            }
        }
    } else {
        OCP_ABORT("WRONG Type !");
    }

    // Internal energy per unit volume of fluid

    // Uf, d Uf / d T, d Uf / d P
    Uf  = -P / (GRAVITY_FACTOR * CONV6);
    UfP = -1 / (GRAVITY_FACTOR * CONV6);
    UfT = 0;

    for (USI j = 0; j < numPhase; j++) {
        // Uf
        Uf += S[j] * xi[j] * H[j];
        // dUf / dP
        UfP += -(vj[j] * xiP[j] / xi[j] + S[j] * vfP) / vf * xi[j] * H[j];
        UfP += xiP[j] * S[j] * H[j];
        // dUf / dT
        UfT += -(vj[j] * xiT[j] / xi[j] + S[j] * vfT) / vf * xi[j] * H[j];
        UfT += (xiT[j] * H[j] + HT[j] * xi[j]) * S[j];
    }

    // d Uf / d Ni
    Ufi[0] = dXsdXp[1] * xi[0] * H[0] + dXsdXp[5] * xi[1] * H[1];
    Ufi[1] = dXsdXp[2] * xi[0] * H[0] + dXsdXp[6] * xi[1] * H[1];
}

OCP_DBL MixtureThermal_K01::CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin)
{
    T                  = Tin + CONV5;
    const OCP_DBL   dT = T - Tref;
    vector<OCP_DBL> zi{0, 1};

    OCP_DBL Hw = 0;
    if (liquid_based) {
        for (USI i = 0; i < numCom; i++) {
            Hw += zi[i] *
                  (cpl1[i] * dT + 1.0 / 2 * cpl2[i] * pow(dT, 2) +
                   1.0 / 3 * cpl3[i] * pow(dT, 3) + 1.0 / 4 * cpl4[i] * pow(dT, 4));
        }
    } else if (gas_based) {
        for (USI i = 0; i < numCom; i++) {
            Hw += zi[i] *
                  (cpg1[i] * dT + 1.0 / 2 * cpg2[i] * pow(dT, 2) +
                   1.0 / 3 * cpg3[i] * pow(dT, 3) + 1.0 / 4 * cpg4[i] * pow(dT, 4));
            if (T < Tcrit[i]) {
                Hw -= zi[i] * hvr[i] * pow((Tcrit[i] - T), ev[i]);
            }
        }
    }
    return Hw;
}

OCP_DBL MixtureThermal_K01::XiPhase(const OCP_DBL& Pin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Ziin,
                                    const USI&     tarPhase)
{
    P                = Pin;
    T                = Tin + CONV5;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;

    if (tarPhase == WATER) {
        // inj fluid is water
        OCP_DBL xiw = xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT -
                                      ct2[1] * pow(dT, 2) / 2 + cpt[1] * dP * dT);
        return xiw;
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

OCP_DBL
MixtureThermal_K01::RhoPhase(const OCP_DBL& Pin,
                             const OCP_DBL& Pbb,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Ziin,
                             const USI&     tarPhase)
{
    P                = Pin;
    T                = Tin + CONV5;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;

    if (tarPhase == OIL) {
        const OCP_DBL xio = xi_ref[0] * exp(cp[0] * dP - ct1[0] * dT -
                                            ct2[0] * pow(dT, 2) / 2 + cpt[0] * dP * dT);
        return MWp[0] * xio;
    } else if (tarPhase == WATER) {
        // inj fluid is water
        const OCP_DBL xiw = xi_ref[1] * exp(cp[1] * dP - ct1[1] * dT -
                                            ct2[1] * pow(dT, 2) / 2 + cpt[1] * dP * dT);
        return MWp[1] * xiw;
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

void MixtureThermal_K01::CalProdWeight(const OCP_DBL&         Pin,
                                       const OCP_DBL&         Tin,
                                       const OCP_DBL*         Niin,
                                       const vector<OCP_DBL>& prodPhase,
                                       vector<OCP_DBL>&       prodWeight)
{
    P                = Pin;
    T                = Tin + 460;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    Ni[0]            = Niin[0];
    Ni[1]            = Niin[1];
    Nt               = Ni[0] + Ni[1];

    // phase molar density
    xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] / 2 * pow(dT, 2) + cpt[0] * dP * dT);
    xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] / 2 * pow(dT, 2) + cpt[1] * dP * dT);

    // phase volume
    vj[0] = Ni[0] / xi[0];
    vj[1] = Ni[1] / xi[1];

    OCP_DBL         qt = Nt;
    vector<OCP_DBL> factor(numPhase, 0);
    factor[0] = vj[0] / qt / CONV1; // stb / lbmol
    factor[1] = vj[1] / qt / CONV1; // stb  / lbmol

    OCP_DBL tmp = 0;
    for (USI i = 0; i < 2; i++) {
        tmp += factor[i] * prodPhase[i];
    }
    if (tmp < 1E-12 || !isfinite(tmp)) {
        OCP_ABORT("Wrong Condition!");
    }
    fill(prodWeight.begin(), prodWeight.end(), tmp);
}

void MixtureThermal_K01::CalProdRate(const OCP_DBL&   Pin,
                                     const OCP_DBL&   Tin,
                                     const OCP_DBL*   Niin,
                                     vector<OCP_DBL>& prodRate)
{
    // assign Value
    P                = Pin;
    T                = Tin + 460;
    const OCP_DBL dP = P - Pref;
    const OCP_DBL dT = T - Tref;
    Ni[0]            = Niin[0];
    Ni[1]            = Niin[1];

    // phase molar density
    xi[0] = xi_ref[0] *
            exp(cp[0] * dP - ct1[0] * dT - ct2[0] / 2 * pow(dT, 2) + cpt[0] * dP * dT);
    xi[1] = xi_ref[1] *
            exp(cp[1] * dP - ct1[1] * dT - ct2[1] / 2 * pow(dT, 2) + cpt[1] * dP * dT);

    // phase volume
    vj[0] = Ni[0] / xi[0];
    vj[1] = Ni[1] / xi[1];

    prodRate[0] = vj[0] / CONV1; // stb
    prodRate[1] = vj[1] / CONV1; // stb
}

void MixtureThermal_K01::SetupWellOpt(WellOpt&                  opt,
                                      const vector<SolventINJ>& sols,
                                      const OCP_DBL&            Psurf,
                                      const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string    fluidName = opt.InjFluidType();
        vector<OCP_DBL> tmpZi(numCom, 0);
        if (fluidName == "WAT") {
            tmpZi.back() = 1;
            opt.SetInjProdPhase(WATER);
            // lbmol / ft3 -> lbmol  / bbl  for
            // injfluid Use flash in Bulk in surface condition
            OCP_DBL tmp = CONV1 * XiPhase(Psurf, Tsurf, &tmpZi[0], WATER);
            opt.SetInjFactor(tmp);
        } else {
            OCP_ABORT("WRONG Fluid Type!");
        }
        opt.SetInjZi(tmpZi);
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(numPhase, 0);
        switch (opt.OptMode()) {
            case ORATE_MODE:
                tmpWght[0] = 1;
                break;
            case WRATE_MODE:
                tmpWght[1] = 1;
                break;
            case LRATE_MODE:
                tmpWght[0] = 1;
                tmpWght[1] = 1;
                break;
            default:
                OCP_ABORT("WRONG optMode");
                break;
        }
        opt.SetProdPhaseWeight(tmpWght);
    } else {
        OCP_ABORT("Wrong Well Type!");
    }
}
/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           NOV/10/2022      Create file                          */
/*----------------------------------------------------------------------------*/
