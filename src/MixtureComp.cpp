/*! \file    MixtureComp.cpp
 *  \brief   MixtureComp class definition for compositional models
 *  \author  Shizhe Li
 *  \date    Jan/05/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureComp.hpp"

COMP::COMP(const vector<string>& comp)
{
    name   = comp[0];
    Pc     = stod(comp[1]);
    Tc     = stod(comp[2]);
    acf    = stod(comp[3]);
    MW     = stod(comp[4]);
    VcMW   = stod(comp[5]);
    Vc     = MW * VcMW;
    OmegaA = stod(comp[6]);
    OmegaB = stod(comp[7]);
    Vshift = stod(comp[8]);
}

MixtureComp::MixtureComp(const ComponentParam& param, const USI& tarId)
{
    mixtureType = EOS_PVTW;
    // if Water don't exist?
    // for Mixture class
    // Now, only one case is considered: oil, gas, water could exist
    numPhase = param.numPhase + 1;
    numCom   = param.numCom + 1;
    Allocate();
    vjp.resize(numPhase, 0);
    vji.resize(numPhase);
    for (auto& v : vji) {
        v.resize(numCom, 0);
    }
    rhoN.resize(numPhase * numCom);
    xiN.resize(numPhase * numCom);
    muN.resize(numPhase * numCom);

    // for MixtureComp class
    NC    = param.numCom;
    NPmax = param.numPhase;

    zi.resize(NC);

    Cname = param.Cname;
    if (param.Tc.activity)
        Tc = param.Tc.data[tarId];
    else
        OCP_ABORT("TCRIT hasn't been input!");
    if (param.Pc.activity)
        Pc = param.Pc.data[tarId];
    else
        OCP_ABORT("PCRIT hasn't been input!");

    if (param.Vc.activity)
        Vc = param.Vc.data[tarId];
    else if (param.Zc.activity) {
        Zc = param.Zc.data[tarId];
        Vc.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vc[i] = 10.73159 * Zc[i] * Tc[i] / Pc[i];
        }
    } else
        OCP_ABORT("VCRIT or ZCRIT hasn't been input!");

    if (param.MW.activity)
        MWC = param.MW.data[tarId];
    else
        OCP_ABORT("MW hasn't been input!");
    if (param.Acf.activity)
        Acf = param.Acf.data[tarId];
    else
        OCP_ABORT("ACF hasn't been input!");
    if (param.OmegaA.activity)
        OmegaA = param.OmegaA.data[tarId];
    else
        OmegaA.resize(NC, 0.457235529);
    if (param.OmegaB.activity)
        OmegaB = param.OmegaB.data[tarId];
    else
        OmegaB.resize(NC, 0.077796074);

    if (param.Vshift.activity) {
        Vshift = param.Vshift.data[tarId];
        for (USI i = 0; i < NC; i++)
            Vshift[i] *= (GAS_CONSTANT * OmegaB[i] * Tc[i] / Pc[i]);
    } else
        Vshift.resize(NC, 0);

    if (param.Vcvis.activity)
        Vcvis = param.Vcvis.data[tarId];
    else if (param.Zcvis.activity) {
        Zcvis = param.Zcvis.data[tarId];
        Vcvis.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vcvis[i] = GAS_CONSTANT * Zcvis[i] * Tc[i] / Pc[i];
        }
    } else
        Vcvis = Vc;

    LBCcoef = param.LBCcoef;
    for (auto& lbc : LBCcoef) {
        lbc *= 10;
    }

    InputMiscibleParam(param, tarId);

    CallId();

    USI len = NC * NC;
    BIC.resize(len, 0);

    if (param.BIC[tarId].size() != len) {
        USI iter = 0;
        for (USI i = 1; i < NC; i++) {
            for (USI j = 0; j < i; j++) {
                BIC[i * NC + j] = param.BIC[tarId][iter];
                BIC[j * NC + i] = BIC[i * NC + j];
                iter++;
            }
        }
    } else {
        BIC = param.BIC[tarId];
    }

    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NC; j++) {
            cout << setw(10) << BIC[i * NC + j] << "   ";
        }
        cout << endl;
    }

    EoSctrl.SSMsta.maxIt = stoi(param.SSMparamSTA[0]);
    EoSctrl.SSMsta.tol   = stod(param.SSMparamSTA[1]);
    EoSctrl.SSMsta.eYt   = stod(param.SSMparamSTA[2]);
    EoSctrl.SSMsta.tol2  = EoSctrl.SSMsta.tol * EoSctrl.SSMsta.tol;

    EoSctrl.NRsta.maxIt = stoi(param.NRparamSTA[0]);
    EoSctrl.NRsta.tol   = stod(param.NRparamSTA[1]);
    EoSctrl.NRsta.tol2  = EoSctrl.NRsta.tol * EoSctrl.NRsta.tol;

    EoSctrl.SSMsp.maxIt = stoi(param.SSMparamSP[0]);
    EoSctrl.SSMsp.tol   = stod(param.SSMparamSP[1]);
    EoSctrl.SSMsp.tol2  = EoSctrl.SSMsp.tol * EoSctrl.SSMsp.tol;

    EoSctrl.NRsp.maxIt = stoi(param.NRparamSP[0]);
    EoSctrl.NRsp.tol   = stod(param.NRparamSP[1]);
    EoSctrl.NRsp.tol2  = EoSctrl.NRsp.tol * EoSctrl.NRsp.tol;

    EoSctrl.RR.maxIt = stoi(param.RRparam[0]);
    EoSctrl.RR.tol   = stod(param.RRparam[1]);
    EoSctrl.RR.tol2  = EoSctrl.RR.tol * EoSctrl.RR.tol;

    AllocateEoS();
    AllocatePhase();
    AllocateMethod();
    AllocateOthers();
}

void MixtureComp::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    ftype = 0;
    lNP   = 0;
    InitPTN(Pin, Tin + CONV5, Niin);
    CalFlash();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt                        = Nh + Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    mu[Wpid]    = data[3];
    xi[Wpid]    = 1 / (CONV1 * bw);
    rho[Wpid]   = std_RhoW / bw;
    vj[Wpid]    = CONV1 * Ni[Wcid] * bw;
    vf += vj[Wpid];

    // Calculate Sj
    CalSaturation();
}

void MixtureComp::InitFlashIMPEC(const OCP_DBL& Pin,
                                 const OCP_DBL& Pbbin,
                                 const OCP_DBL& Tin,
                                 const OCP_DBL* Sjin,
                                 const OCP_DBL& Vpore,
                                 const OCP_DBL* Ziin,
                                 const OCP_USI& bId)
{
    // Attention: zi[numCom - 1] = 0 here, that's Zw = 0;
    SetBulkId(bId);
    ftype = 0;
    lNP   = 0;

    InitPTZ(Pin, Tin + CONV5, Ziin);
    PhaseEquilibrium();
    // Attention Nt = 1 now
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();

    // Calulate Nt, water is exclued
    OCP_DBL Sw = Sjin[numPhase - 1];
    Nh         = Vpore * (1 - Sw) / vf;

    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nh, &vj[0]);
    vf *= Nh;
    // CalVfiVfp_full01();
    CalVfiVfp_full02();
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nh;
    }

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;
    mu[Wpid]    = data[3];
    xi[Wpid]    = 1 / (CONV1 * bw);
    rho[Wpid]   = std_RhoW / bw;
    Ni[Wcid]    = Vpore * Sw * xi[Wpid];
    Nt          = Nh + Ni[Wcid];
    vj[Wpid]    = CONV1 * Ni[Wcid] * bw;
    vf += vj[Wpid];
    vfi[Wcid] = CONV1 * bw;
    vfP += CONV1 * Ni[Wcid] * bwp;

    // Calculate Sj
    CalSaturation();

    CalSkipForNextStep();
}

void MixtureComp::InitFlashFIM(const OCP_DBL& Pin,
                               const OCP_DBL& Pbbin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Sjin,
                               const OCP_DBL& Vpore,
                               const OCP_DBL* Ziin,
                               const OCP_USI& bId)
{
    // Attention: zi[numCom - 1] = 0 here, that's Zw = 0;
    SetBulkId(bId);
    ftype = 0;
    lNP   = 0;

    InitPTZ(Pin, Tin + CONV5, Ziin);
    PhaseEquilibrium();
    // Attention Nt = 1 now
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();

    // Calulate Nt, water is exclued
    OCP_DBL Sw = Sjin[numPhase - 1];
    Nh         = Vpore * (1 - Sw) / vf;

    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nh, &vj[0]);
    vf *= Nh;
    // CalVfiVfp_full01();
    CalVfiVfp_full02();
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nh;
    }

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0       = data[0];
    OCP_DBL bw0       = data[1];
    OCP_DBL cbw       = data[2];
    OCP_DBL bw        = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp       = -cbw * bw0;
    mu[Wpid]          = data[3];
    xi[Wpid]          = 1 / (CONV1 * bw);
    rho[Wpid]         = std_RhoW / bw;
    muP[Wpid]         = cdata[3];
    xiP[Wpid]         = -bwp / (bw * bw * CONV1);
    rhoP[Wpid]        = CONV1 * xiP[Wpid] * std_RhoW;
    Ni[Wcid]          = Vpore * Sw * xi[Wpid];
    nj[Wpid]          = Ni[Wcid];
    Nt                = Nh + Ni[Wcid];
    vj[Wpid]          = CONV1 * Ni[Wcid] * bw;
    vfi[Wcid]         = CONV1 * bw;
    const OCP_DBL vwp = CONV1 * Ni[Wcid] * bwp;
    vf += vj[Wpid];
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

#ifdef OCP_OLD_FIM
    CaldXsdXpAPI02();
#else
    CaldXsdXpAPI02p();
#endif // OCP_OLD_FIM
    CalXiPNX_partial();
    CalRhoPX_partial();
    CalMuPX_partial();

    // Calculate pSderExist and pVnumCom
    fill(pSderExist.begin(), pSderExist.end(), OCP_FALSE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0.0);
    if (phaseExist[0]) {
        pSderExist[0] = OCP_TRUE;
        pVnumCom[0]   = NC;
    }
    if (phaseExist[1]) {
        pSderExist[1] = OCP_TRUE;
        pVnumCom[1]   = NC;
    }
    pSderExist[2] = OCP_TRUE;

    CalSkipForNextStep();
}

void MixtureComp::InitFlashFIMn(const OCP_DBL& Pin,
                                const OCP_DBL& Pbbin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Sjin,
                                const OCP_DBL& Vpore,
                                const OCP_DBL* Ziin,
                                const OCP_USI& bId)
{
    // Attention: zi[numCom - 1] = 0 here, that's Zw = 0;
    SetBulkId(bId);
    ftype = 0;
    lNP   = 0;

    InitPTZ(Pin, Tin + CONV5, Ziin);
    PhaseEquilibrium();
    // Attention Nt = 1 now
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();

    // Calulate Nt, water is exclued
    OCP_DBL Sw = Sjin[numPhase - 1];
    Nh         = Vpore * (1 - Sw) / vf;

    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nh, &vj[0]);
    vf *= Nh;
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nh;
    }

    CalVjpVfpVfn_partial();
    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0       = data[0];
    OCP_DBL bw0       = data[1];
    OCP_DBL cbw       = data[2];
    OCP_DBL bw        = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp       = -cbw * bw0;
    mu[Wpid]          = data[3];
    xi[Wpid]          = 1 / (CONV1 * bw);
    rho[Wpid]         = std_RhoW / bw;
    muP[Wpid]         = cdata[3];
    xiP[Wpid]         = -bwp / (bw * bw * CONV1);
    rhoP[Wpid]        = CONV1 * xiP[Wpid] * std_RhoW;
    Ni[Wcid]          = Vpore * Sw * xi[Wpid];
    nj[Wpid]          = Ni[Wcid];
    Nt                = Nh + Ni[Wcid];
    vj[Wpid]          = CONV1 * Ni[Wcid] * bw;
    vfi[Wcid]         = CONV1 * bw;
    const OCP_DBL vwp = CONV1 * Ni[Wcid] * bwp;
    vf += vj[Wpid];
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();
    CaldXsdXpAPI03();
    CalXiPn_partial();
    CalRhoPn_partial();
    CalMuPn_partial();
    CalVfiVfp_full03();

    // Calculate pSderExist and pVnumCom
    fill(pSderExist.begin(), pSderExist.end(), OCP_FALSE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0.0);
    if (phaseExist[0]) {
        pSderExist[0] = OCP_TRUE;
        pVnumCom[0]   = NC;
    }
    if (phaseExist[1]) {
        pSderExist[1] = OCP_TRUE;
        pVnumCom[1]   = NC;
    }
    pSderExist[2] = OCP_TRUE;

    CalSkipForNextStep();
}

void MixtureComp::FlashIMPEC(const OCP_DBL& Pin,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Niin,
                             const USI&     lastNP,
                             const OCP_DBL* xijin,
                             const OCP_USI& bId)
{
    SetBulkId(bId);
    // Hydroncarbon phase, if lNp = 0, then strict stability analysis will be used
    lNP = lastNP > 0 ? lastNP - 1 : 0;
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeIMPEC();
    CalFlash();
    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    // CalVfiVfp_full01();
    CalVfiVfp_full02();

    // Water Properties
    const USI Wpid            = numPhase - 1;
    const USI Wcid            = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt                        = Nh + Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;
    mu[Wpid]    = data[3];
    xi[Wpid]    = 1 / (CONV1 * bw);
    rho[Wpid]   = std_RhoW / bw;
    vj[Wpid]    = CONV1 * Ni[Wcid] * bw;
    vf += vj[Wpid];
    vfi[Wcid] = CONV1 * bw;
    vfP += CONV1 * Ni[Wcid] * bwp;

    // Calculate Sj
    CalSaturation();

    CalSkipForNextStep();
}

void MixtureComp::FlashFIM(const OCP_DBL& Pin,
                           const OCP_DBL& Tin,
                           const OCP_DBL* Niin,
                           const OCP_DBL* Sjin,
                           const USI&     lastNP,
                           const OCP_DBL* xijin,
                           const OCP_USI& bId)
{

    SetBulkId(bId);

    // Hydroncarbon phase, if lNp = 0, then strict stability analysis will be used
    lNP = lastNP > 0 ? lastNP - 1 : 0;
    if (lNP == 2) {
        for (USI i = 0; i < NC; i++) {
            lKs[i] = xijin[i] / xijin[i + numCom];
        }
    }

    InitPTN(Pin, Tin + CONV5, Niin);
    CalFtypeFIM(Sjin);
    CalFlash();

    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    CalVfiVfp_full02();

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid]                  = Ni[Wcid];
    Nt                        = Nh + Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0       = data[0];
    OCP_DBL bw0       = data[1];
    OCP_DBL cbw       = data[2];
    OCP_DBL bw        = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp       = -cbw * bw0;
    mu[Wpid]          = data[3];
    xi[Wpid]          = 1 / (CONV1 * bw);
    rho[Wpid]         = std_RhoW / bw;
    muP[Wpid]         = cdata[3];
    xiP[Wpid]         = -bwp / (bw * bw * CONV1);
    rhoP[Wpid]        = CONV1 * xiP[Wpid] * std_RhoW;
    vj[Wpid]          = CONV1 * Ni[Wcid] * bw;
    vfi[Wcid]         = CONV1 * bw;
    const OCP_DBL vwp = CONV1 * Ni[Wcid] * bwp;
    vf += vj[Wpid];
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

#ifdef OCP_OLD_FIM
    CaldXsdXpAPI02();
#else
    CaldXsdXpAPI02p();
#endif // OCP_OLD_FIM
    CalXiPNX_partial();
    CalRhoPX_partial();
    CalMuPX_partial();

    // Calculate pSderExist and pVnumCom
    fill(pSderExist.begin(), pSderExist.end(), OCP_FALSE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0);
    if (phaseExist[0]) {
        pSderExist[0] = OCP_TRUE;
        pVnumCom[0]   = NC;
    }
    if (phaseExist[1]) {
        pSderExist[1] = OCP_TRUE;
        pVnumCom[1]   = NC;
    }
    pSderExist[2] = OCP_TRUE;

    CalSkipForNextStep();
}

void MixtureComp::FlashFIMn(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Niin,
                            const OCP_DBL* Sjin,
                            const OCP_DBL* xijin,
                            const OCP_DBL* njin,
                            const USI*     phaseExistin,
                            const USI&     lastNP,
                            const OCP_USI& bId)
{
    SetBulkId(bId);
    inputNP = 0;
    for (USI j = 0; j < numPhase; j++) {
        if (phaseExistin[j] == 1) inputNP++;
    }
    inputNP--;

    if (inputNP == 1 || OCP_TRUE) {
        // Hydroncarbon phase, if lNp = 0, then strict stability analysis will be used
        lNP = lastNP > 0 ? lastNP - 1 : 0;
        if (lNP == 2) {
            for (USI i = 0; i < NC; i++) {
                lKs[i] = xijin[i] / xijin[i + numCom];
            }
        }
        InitPTN(Pin, Tin + CONV5, Niin);
        CalFtypeFIM(Sjin);
        CalFlash();
    } else {
        //! Becareful if NP > 2 (temp)
        NP = inputNP;
        InitPTN(Pin, Tin + CONV5, Niin);
        CalAiBi();
        for (USI j = 0; j < NP; j++) {
            phaseExist[j] = phaseExistin[j];
            S[j]          = Sjin[j];
            nu[j]         = njin[j];
            nj[j]         = njin[j];
            Dcopy(NC, &x[j][0], &xijin[j * numCom]);
            Dcopy(NC, &xij[j * numCom], &xijin[j * numCom]);
        }
        CalFugPhiAll();
        S[numPhase - 1] = Sjin[numPhase - 1];
        phaseLabel[0]   = OIL;
        phaseLabel[1]   = GAS;
        CalMW();
        CalVfXiRho();
        CalViscosity();
        CopyPhase();
        CalSurfaceTension();
    }

    CalVjpVfpVfn_partial();

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = OCP_TRUE;
    xij[Wpid * numCom + Wcid] = 1.0;
    nj[Wpid]                  = Ni[Wcid];
    Nt                        = Nh + Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0       = data[0];
    OCP_DBL bw0       = data[1];
    OCP_DBL cbw       = data[2];
    OCP_DBL bw        = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp       = -cbw * bw0;
    mu[Wpid]          = data[3];
    xi[Wpid]          = 1 / (CONV1 * bw);
    rho[Wpid]         = std_RhoW / bw;
    muP[Wpid]         = cdata[3];
    xiP[Wpid]         = -bwp / (bw * bw * CONV1);
    rhoP[Wpid]        = CONV1 * xiP[Wpid] * std_RhoW;
    vj[Wpid]          = CONV1 * Ni[Wcid] * bw;
    vfi[Wcid]         = CONV1 * bw;
    const OCP_DBL vwp = CONV1 * Ni[Wcid] * bwp;
    vf += vj[Wpid];
    vfP += vwp;
    vji[numPhase - 1][numCom - 1] = vfi[Wcid];
    vjp[numPhase - 1]             = vwp;

    CalSaturation();

    // calculate diff
    if (inputNP == NP && NP == 2 && OCP_FALSE) {
        cout << scientific << setprecision(3);
        cout << "--------- " << NP << " --------- " << inputNP << endl;
        for (USI j = 0; j < NP; j++) {
            USI j1 = phaseLabel[j];
            cout << setw(10) << S[j1] - Sjin[j1] << "   " << setw(10)
                 << (nj[j1] - njin[j1]) / njin[j1] << "   ";
            for (USI i = 0; i < NC; i++) {
                cout << setw(10) << xij[j1 * numCom + i] - xijin[j1 * numCom + i]
                     << "   ";
            }
            cout << endl;
        }
        cout << S[2] - Sjin[2] << endl;
    }

    CaldXsdXpAPI03();
    CalXiPn_partial();
    CalRhoPn_partial();
    CalMuPn_partial();
    CalVfiVfp_full03();

    // Calculate pSderExist and pVnumCom
    fill(pSderExist.begin(), pSderExist.end(), OCP_FALSE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0.0);
    if (phaseExist[0]) {
        pSderExist[0] = OCP_TRUE;
        pVnumCom[0]   = NC;
    }
    if (phaseExist[1]) {
        pSderExist[1] = OCP_TRUE;
        pVnumCom[1]   = NC;
    }
    pSderExist[2] = OCP_TRUE;

    CalSkipForNextStep();
}

void MixtureComp::CalFlash()
{
    PhaseEquilibrium();
    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nh, &nu[0]);
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();
}

OCP_DBL
MixtureComp::XiPhase(const OCP_DBL& Pin,
                     const OCP_DBL& Tin,
                     const OCP_DBL* Ziin,
                     const USI&     tarPhase)
{
    // assume that only single phase exists here
    if (tarPhase == WATER) {
        // water phase
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0   = data[0];
        OCP_DBL bw0   = data[1];
        OCP_DBL cbw   = data[2];
        OCP_DBL bw    = bw0 * (1 - cbw * (Pin - Pw0));
        OCP_DBL xitmp = 1 / (CONV1 * bw);
        return xitmp;
    } else {
        // hydrocarbon phase
        InitPTZ(Pin, Tin + CONV5, Ziin); // P, T has been Set !!!
        NP = 1;
        CalAiBi();
        CalAjBj(Aj[0], Bj[0], zi);
        SolEoS(Zj[0], Aj[0], Bj[0]);
        OCP_DBL vtmp = Zj[0] * GAS_CONSTANT * T / P;
        for (USI i = 0; i < NC; i++) {
            vtmp -= zi[i] * Vshift[i];
        }
        OCP_DBL xitmp = 1 / vtmp;
        return xitmp;
    }
}

OCP_DBL
MixtureComp::RhoPhase(const OCP_DBL& Pin,
                      const OCP_DBL& Pbb,
                      const OCP_DBL& Tin,
                      const OCP_DBL* Ziin,
                      const USI&     tarPhase)
{

    // assume that only single phase exists here
    if (tarPhase == WATER) {
        // water phase
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0 = data[0];
        OCP_DBL bw0 = data[1];
        OCP_DBL cbw = data[2];
        OCP_DBL bw  = (bw0 * (1 - cbw * (Pin - Pw0)));

        return std_RhoW / bw;
    } else {
        // hydrocarbon phase
        OCP_DBL xitmp = XiPhase(Pin, Tin, Ziin, tarPhase);
        NP            = 1;
        x[0]          = zi;
        CalMW();
        return MW[0] * xitmp;
    }
}

void MixtureComp::SetupWellOpt(WellOpt&                  opt,
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
        } else {
            // inj phase is gas
            opt.SetInjProdPhase(GAS);
            const USI len = sols.size();
            for (USI i = 0; i < len; i++) {
                if (fluidName == sols[i].name) {
                    tmpZi = sols[i].data;
                    tmpZi.resize(numCom);
                    // Convert volume units Mscf/stb to molar units lbmoles for
                    // injfluid Use flash in Bulk in surface condition
                    OCP_DBL tmp = 1000 * XiPhase(Psurf, Tsurf, &tmpZi[0], GAS);
                    opt.SetInjFactor(tmp);
                    break;
                }
                if (i == len - 1) {
                    OCP_ABORT("Wrong FluidType!");
                }
            }
        }
        opt.SetInjZi(tmpZi);
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(numPhase, 0);
        switch (opt.OptMode()) {
            case BHP_MODE:
                break;
            case ORATE_MODE:
                tmpWght[0] = 1;
                break;
            case GRATE_MODE:
                tmpWght[1] = 1;
                break;
            case WRATE_MODE:
                tmpWght[2] = 1;
                break;
            case LRATE_MODE:
                tmpWght[0] = 1;
                tmpWght[2] = 1;
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

void MixtureComp::CalProdWeight(const OCP_DBL&         Pin,
                                const OCP_DBL&         Tin,
                                const OCP_DBL*         Niin,
                                const vector<OCP_DBL>& prodPhase,
                                vector<OCP_DBL>&       prodWeight)
{
    Flash(Pin, Tin, Niin);

    OCP_DBL         qt = Nt;
    vector<OCP_DBL> factor(numPhase, 0);

    factor[0] = vj[0] / qt / CONV1; // stb / lbmol
    factor[1] = vj[1] / qt / 1000;  // Mscf / lbmol
    factor[2] = xi[2] * vj[2] / qt; // stb  / stb

    OCP_DBL tmp = 0;
    for (USI i = 0; i < 3; i++) {
        tmp += factor[i] * prodPhase[i];
    }
    if (tmp < 1E-12 || !isfinite(tmp)) {
        OCP_ABORT("Wrong Condition!");
    }
    fill(prodWeight.begin(), prodWeight.end(), tmp);
}

void MixtureComp::CalProdRate(const OCP_DBL&   Pin,
                              const OCP_DBL&   Tin,
                              const OCP_DBL*   Niin,
                              vector<OCP_DBL>& prodRate)
{
    Flash(Pin, Tin, Niin);

    prodRate[0] = vj[0] / CONV1; // stb
    prodRate[1] = vj[1] / 1000;  // Mscf
    prodRate[2] = vj[2] * xi[2]; // stb
}

void MixtureComp::CallId()
{
    lId = 0;
    for (USI i = 1; i < NC; i++) {
        if (MWC[i] < MWC[lId]) lId = i;
    }
}

void MixtureComp::AllocateEoS()
{
    // Allocate Memoery for EoS variables
    Ai.resize(NC);
    Bi.resize(NC);
    Aj.resize(NPmax);
    Bj.resize(NPmax);
    Zj.resize(NPmax);
    Ztmp.resize(3);
    delta1P2 = delta1 + delta2;
    delta1M2 = delta1 - delta2;
    delta1T2 = delta1 * delta2;
}

void MixtureComp::SolEoS(OCP_DBL& ZjT, const OCP_DBL& AjT, const OCP_DBL& BjT) const
{
    const OCP_DBL aj = AjT;
    const OCP_DBL bj = BjT;

    const OCP_DBL a = (delta1 + delta2 - 1) * bj - 1;
    const OCP_DBL b =
        (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1));
    const OCP_DBL c = -(aj * bj + delta1 * delta2 * bj * bj * (bj + 1));

    USI flag = CubicRoot(a, b, c, OCP_TRUE); // True with NT
    if (flag == 1) {
        ZjT = Ztmp[0];
    } else {
        OCP_DBL zj1 = Ztmp[0];
        OCP_DBL zj2 = Ztmp[2];
        OCP_DBL dG  = (zj2 - zj1) + log((zj1 - bj) / (zj2 - bj)) -
                     aj / (bj * (delta2 - delta1)) *
                         log((zj1 + delta1 * bj) * (zj2 + delta2 * bj) /
                             ((zj1 + delta2 * bj) * (zj2 + delta1 * bj)));
        if (dG > 0)
            ZjT = zj1;
        else
            ZjT = zj2;
    }
}

void MixtureComp::CalAiBi()
{
    // Calculate Ai, Bi
    OCP_DBL acf, mwi;
    OCP_DBL Pri, Tri;

    for (USI i = 0; i < NC; i++) {
        acf = Acf[i];
        // PR
        if (acf <= 0.49) {
            mwi = 0.37464 + 1.54226 * acf - 0.26992 * pow(acf, 2);
        } else {
            mwi = 0.379642 + 1.48503 * acf - 0.164423 * pow(acf, 2) +
                  0.016667 * pow(acf, 3);
        }

        Pri   = P / Pc[i];
        Tri   = T / Tc[i];
        Ai[i] = OmegaA[i] * Pri / pow(Tri, 2) * pow((1 + mwi * (1 - sqrt(Tri))), 2);
        Bi[i] = OmegaB[i] * Pri / Tri;
    }
}

void MixtureComp::CalAjBj(OCP_DBL& AjT, OCP_DBL& BjT, const vector<OCP_DBL>& xj) const
{
    AjT = 0;
    BjT = 0;

    for (USI i1 = 0; i1 < NC; i1++) {
        BjT += Bi[i1] * xj[i1];
        AjT += xj[i1] * xj[i1] * Ai[i1] * (1 - BIC[i1 * NC + i1]);

        for (USI i2 = 0; i2 < i1; i2++) {
            AjT +=
                2 * xj[i1] * xj[i2] * sqrt(Ai[i1] * Ai[i2]) * (1 - BIC[i1 * NC + i2]);
        }
    }
}

void MixtureComp::CalAjBj(OCP_DBL& AjT, OCP_DBL& BjT, const OCP_DBL* xj) const
{
    AjT = 0;
    BjT = 0;

    for (USI i1 = 0; i1 < NC; i1++) {
        BjT += Bi[i1] * xj[i1];
        AjT += xj[i1] * xj[i1] * Ai[i1] * (1 - BIC[i1 * NC + i1]);

        for (USI i2 = 0; i2 < i1; i2++) {
            AjT +=
                2 * xj[i1] * xj[i2] * sqrt(Ai[i1] * Ai[i2]) * (1 - BIC[i1 * NC + i2]);
        }
    }
}

void MixtureComp::AllocatePhase()
{
    // Allocate Memoery for Phase variables
    vC.resize(NPmax);
    nu.resize(NPmax);
    xiC.resize(NPmax);
    rhoC.resize(NPmax);
    MW.resize(NPmax);
    phaseLabel.resize(NPmax);
    x.resize(NPmax);
    phi.resize(NPmax);
    fug.resize(NPmax);
    n.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        x[j].resize(NC);
        phi[j].resize(NC);
        fug[j].resize(NC);
        n[j].resize(NC);
    }
    ln = n;
}

void MixtureComp::CalFugPhi(vector<OCP_DBL>&       phiT,
                            vector<OCP_DBL>&       fugT,
                            const vector<OCP_DBL>& xj)
{
    OCP_DBL aj, bj, zj;
    CalAjBj(aj, bj, xj);
    SolEoS(zj, aj, bj);

    const OCP_DBL m1 = delta1;
    const OCP_DBL m2 = delta2;
    // const OCP_DBL m1Mm2 = delta1M2;

    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {
        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
        }
        phiT[i] = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) -
                      aj / (m1 - m2) / bj * (tmp / aj - Bi[i] / bj) *
                          log((zj + m1 * bj) / (zj + m2 * bj)));
        fugT[i] = phiT[i] * xj[i] * P;
    }

    //   OCP_DBL       tmp01 = log(zj - bj);
    //   OCP_DBL       tmp02 = aj / (m1Mm2) / bj * log((zj + m1 * bj) / (zj + m2 * bj));

    // for (USI i = 0; i < NC; i++) {
    //	tmp = 0;
    //	for (int k = 0; k < NC; k++) {
    //		tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
    //	}
    //       phiT[i] = exp(Bi[i] / bj * (zj - 1) - tmp01 - tmp02 * (tmp / aj - Bi[i] /
    //       bj));
    //	fugT[i] = phiT[i] * xj[i] * P;
    //}

    Asta = aj;
    Bsta = bj;
    Zsta = zj;
}

void MixtureComp::CalFugPhi(OCP_DBL* phiT, OCP_DBL* fugT, const OCP_DBL* xj)
{
    OCP_DBL aj, bj, zj;
    CalAjBj(aj, bj, xj);
    SolEoS(zj, aj, bj);

    const OCP_DBL m1    = delta1;
    const OCP_DBL m2    = delta2;
    const OCP_DBL m1Mm2 = delta1M2;

    // OCP_DBL tmp01 = log(zj - bj);
    // OCP_DBL tmp02 = aj / (m1Mm2) / bj * log((zj + m1 * bj) / (zj + m2 * bj));

    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {
        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
        }
        phiT[i] = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) -
                      aj / (m1Mm2) / bj * (tmp / aj - Bi[i] / bj) *
                          log((zj + m1 * bj) / (zj + m2 * bj)));
        fugT[i] = phiT[i] * xj[i] * P;
    }
    // for (USI i = 0; i < NC; i++) {
    //    tmp = 0;
    //    for (int k = 0; k < NC; k++) {
    //        tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
    //    }
    //    phiT[i] = exp(Bi[i] / bj * (zj - 1) - tmp01 - tmp02 * (tmp / aj - Bi[i] /
    //    bj)); fugT[i] = phiT[i] * xj[i] * P;
    //}

    Asta = aj;
    Bsta = bj;
    Zsta = zj;
}

void MixtureComp::CalFugPhi(OCP_DBL* fugT, const OCP_DBL* xj)
{
    OCP_DBL aj, bj, zj;
    CalAjBj(aj, bj, xj);
    SolEoS(zj, aj, bj);

    const OCP_DBL m1    = delta1;
    const OCP_DBL m2    = delta2;
    const OCP_DBL m1Mm2 = delta1M2;

    // OCP_DBL tmp01 = log(zj - bj);
    // OCP_DBL tmp02 = aj / (m1Mm2) / bj * log((zj + m1 * bj) / (zj + m2 * bj));

    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {
        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
        }
        tmp     = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) -
                      aj / (m1Mm2) / bj * (tmp / aj - Bi[i] / bj) *
                          log((zj + m1 * bj) / (zj + m2 * bj)));
        fugT[i] = tmp * xj[i] * P;
    }

    Asta = aj;
    Bsta = bj;
    Zsta = zj;
}

void MixtureComp::CalFugPhiAll()
{
    OCP_DBL       tmp; // , tmp01, tmp02;
    const OCP_DBL m1    = delta1;
    const OCP_DBL m2    = delta2;
    const OCP_DBL m1Mm2 = delta1M2;

    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj   = x[j];
        vector<OCP_DBL>&       phiT = phi[j];
        vector<OCP_DBL>&       fugT = fug[j];
        OCP_DBL&               aj   = Aj[j];
        OCP_DBL&               bj   = Bj[j];
        OCP_DBL&               zj   = Zj[j];

        CalAjBj(aj, bj, xj);
        SolEoS(zj, aj, bj);

        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
            }
            phiT[i] = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) -
                          aj / (m1Mm2) / bj * (tmp / aj - Bi[i] / bj) *
                              log((zj + m1 * bj) / (zj + m2 * bj)));
            fugT[i] = phiT[i] * xj[i] * P;
        }

        // tmp01 = log(zj - bj);
        // tmp02 = aj / (m1Mm2) / bj * log((zj + m1 * bj) / (zj + m2 * bj));
        // for (USI i = 0; i < NC; i++) {
        //	tmp = 0;
        //	for (int k = 0; k < NC; k++) {
        //		tmp += 2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
        //	}
        //          phiT[i] =
        //              exp(Bi[i] / bj * (zj - 1) - tmp01 - tmp02 * (tmp / aj - Bi[i] /
        //              bj));
        //	fugT[i] = phiT[i] * xj[i] * P;
        //}
    }
}

void MixtureComp::CalMW()
{
    // Calculate Molecular Weight of phase
    for (USI j = 0; j < NP; j++) {
        MW[j] = 0;
        for (USI i = 0; i < NC; i++) {
            MW[j] += x[j][i] * MWC[i];
        }
    }
}

void MixtureComp::CalVfXiRho()
{
    // Attention that nu is moles of phase instead of molar fraction of phase now
    vf = 0;
    OCP_DBL tmp;
    for (USI j = 0; j < NP; j++) {

        vector<OCP_DBL>& xj = x[j];
        tmp                 = Zj[j] * GAS_CONSTANT * T / P;
        for (USI i = 0; i < NC; i++) {
            tmp -= xj[i] * Vshift[i];
        }

        vC[j] = tmp * nu[j];
        vf += vC[j];
        xiC[j]  = 1 / tmp;
        rhoC[j] = MW[j] * xiC[j];
    }
}

void MixtureComp::CalSaturation()
{
    for (USI j = 0; j < numPhase; j++) {
        S[j] = 0;
        if (phaseExist[j]) {
            S[j] = vj[j] / vf;
        }
    }
}

USI MixtureComp::FindMWmax()
{
    // find the phase id with the highest Molecular Weight
    USI     tmpId = 0;
    OCP_DBL tmpMW = MW[0];
    for (USI j = 1; j < NP; j++) {
        if (tmpMW < MW[j]) {
            tmpMW = MW[j];
            tmpId = j;
        }
    }
    return tmpId;
}

void MixtureComp::x2n()
{
    // Total moles are supposed to be 1
    for (USI j = 0; j < NP; j++) {

        // nu[j] = fabs(nu[j]);
        for (USI i = 0; i < NC; i++) {
            n[j][i] = nu[j] * x[j][i];
        }
    }
}

void MixtureComp::PrintX()
{
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            cout << x[j][i] << "   ";
        }
        cout << endl;
    }
    cout << "----------------------" << endl;
}

void MixtureComp::AllocateMethod()
{

    Kw.resize(4);
    for (USI i = 0; i < 4; i++) {
        Kw[i].resize(NC);
    }
    Ks.resize(NPmax - 1);
    for (USI i = 0; i < NPmax - 1; i++) {
        Ks[i].resize(NC);
    }
    phiSta.resize(NC);
    fugSta.resize(NC);

    Y.resize(NC);
    di.resize(NC);
    resSTA.resize(NC);

    JmatSTA.resize(NC * NC);
    Ax.resize(NC);
    Bx.resize(NC);
    Zx.resize(NC);
    lKs.resize(NC);

    resRR.resize(NPmax - 1);
    resSP.resize(static_cast<size_t>(NC) * NPmax);
    JmatSP.resize(static_cast<size_t>(NC * NC) * NPmax * NPmax);
    fugX.resize(NPmax);
    fugN.resize(NPmax);
    Zn.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        fugX[j].resize(NC * NC);
        fugN[j].resize(NC * NC);
        Zn[j].resize(NC);
    }
    An.resize(NC);
    Bn.resize(NC);
    pivot.resize(static_cast<size_t>(NC + 1) * NPmax, 1);
    lJmatWork = NC * (NPmax - 1);
    JmatWork.resize(lJmatWork);

    pivot.resize(NPmax * static_cast<size_t>(NC) + numPhase, 1);
}

void MixtureComp::CalKwilson()
{
    for (USI i = 0; i < NC; i++) {
        Kw[0][i] = (Pc[i] / P) * exp(5.373 * (1 + Acf[i]) * (1 - Tc[i] / T));
        Kw[1][i] = 1 / Kw[0][i];
        Kw[2][i] = pow(Kw[0][i], 1.0 / 3);
        Kw[3][i] = pow(Kw[1][i], 1.0 / 3);
    }
}

void MixtureComp::PhaseEquilibrium()
{
    // Attention: sum of components' moles equals 1
    switch (ftype) {
        case 0:
            // flash from single phase
            flagSkip = OCP_TRUE;
            NP       = 1;
            nu[0]    = 1;
            x[0]     = zi;
            CalAiBi();
            CalAjBj(Aj[0], Bj[0], x[0]);
            SolEoS(Zj[0], Aj[0], Bj[0]);
            CalKwilson();
            while (!PhaseStable()) {
                NP++;
                PhaseSplit();
                if (NP == NPmax || NP == 1) break;
            }
            if (NP > 1) {
                flagSkip = OCP_FALSE;
            }
            // record error
            if (NP == 1) {
                if (EoSctrl.NRsta.conflag)
                    ePEC = EoSctrl.NRsta.realTol;
                else if (EoSctrl.SSMsta.conflag)
                    ePEC = EoSctrl.SSMsta.realTol;
                else
                    ePEC = 1E8;
            } else {
                if (EoSctrl.NRsp.conflag)
                    ePEC = EoSctrl.NRsp.realTol;
                else if (EoSctrl.SSMsp.conflag)
                    ePEC = EoSctrl.SSMsp.realTol;
                else
                    ePEC = 1E8;
            }

            break;
        case 1:
            // Skip Phase Stability analysis, only single phase exists
            flagSkip = OCP_TRUE;
            NP       = 1;
            nu[0]    = 1;
            x[0]     = zi;
            CalAiBi();
            CalAjBj(Aj[0], Bj[0], x[0]);
            SolEoS(Zj[0], Aj[0], Bj[0]);
            // record error
            ePEC = 0.0;
            break;

        case 2:
            // Skip Phase Stability analysis, two phases exist
            flagSkip = OCP_FALSE;
            NP       = 2;
            Yt       = 1.01;
            CalAiBi();
            CalKwilson();
            PhaseSplit();

            if (EoSctrl.NRsp.conflag)
                ePEC = EoSctrl.NRsp.realTol;
            else if (EoSctrl.SSMsp.conflag)
                ePEC = EoSctrl.SSMsp.realTol;
            else
                ePEC = 1E8;
            break;

        default:
            OCP_ABORT("Wrong flash type!");
            break;
    }
}

OCP_BOOL MixtureComp::PhaseStable()
{
    if (NP == 1) {
        testPId = 0;
    } else {
        CalMW();
        testPId = FindMWmax();
    }

    EoSctrl.SSMsta.conflag = OCP_FALSE;
    EoSctrl.SSMsta.curIt   = 0;
    EoSctrl.NRsta.conflag  = OCP_FALSE;
    EoSctrl.NRsta.curIt    = 0;

    // Test if a phase is stable, if stable return OCP_TRUE, else return OCP_FALSE
    OCP_BOOL flag;
    USI      tmpNP = NP;

    if (lNP == 0) {
        // strict stability ananlysis
        flag = StableSSM(testPId);
    } else {
        flag = StableSSM01(testPId);
        if (!flag) tmpNP++;

        if (tmpNP != lNP) {
            flag     = StableSSM(testPId);
            flagSkip = OCP_FALSE;
        }
    }
    itersSSMSTA += EoSctrl.SSMsta.curIt;
    itersNRSTA += EoSctrl.NRsta.curIt;
    countsSSMSTA++;
    countsNRSTA++;
    return flag;
}

OCP_BOOL MixtureComp::StableSSM01(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    OCP_DBL Stol  = EoSctrl.SSMsta.tol2;
    USI     maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL Ktol  = EoSctrl.SSMsta.Ktol;
    OCP_DBL dYtol = EoSctrl.SSMsta.dYtol;
    // OCP_DBL& Sk = EoSctrl.SSMsta.curSk;
    OCP_DBL Se, Sk, dY;

    OCP_BOOL flag, Tsol; // Tsol, trivial solution
    USI      iter, k;

    const vector<OCP_DBL>& xj = x[Id];
    CalFugPhi(phi[Id], fug[Id], xj);
    const vector<OCP_DBL>& fugId = fug[Id];

    vector<OCP_DBL>& ks = Ks[0];
    for (k = 0; k < 2; k++) {

        ks   = Kw[k];
        iter = 0;
        flag = OCP_FALSE;
        Tsol = OCP_FALSE;
        while (OCP_TRUE) {
            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = xj[i] * ks[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);

            if (iter > 0) {
                dY = 0;
                for (USI i = 0; i < NC; i++) {
                    dY += pow((Y[i] - di[i]), 2);
                }
                if (dY < dYtol) {
                    // converges
                    flag = OCP_TRUE;
                    break;
                }
            }

            CalFugPhi(&fugSta[0], &Y[0]);
            Se = 0;
            Sk = 0;
            for (USI i = 0; i < NC; i++) {
                di[i] = fugId[i] / (fugSta[i] * Yt);
                ks[i] *= di[i];
                Se += pow(di[i] - 1, 2);
                // Sk += pow(ks[i] - 1, 2);
                Sk += pow(log(ks[i]), 2);
            }

            iter++;
            if (Se < Stol) {
                flag = OCP_TRUE;
                break;
            }
            if (Sk < Ktol) {
                // Sk < Ktol -> trivial solution
                flag = OCP_TRUE;
                Tsol = OCP_TRUE;
                break;
            }
            if (iter > maxIt) {
                break;
            }

            // Record last Y with di
            di = Y;
        }

        if (!Tsol) {
            // flag = StableNR(Id);
        }

        // cout << "Yt = " << setprecision(8) << scientific << Yt << "   " << setw(2)
        //     << "Sk = " << setprecision(3) << scientific << Sk << "   " << setw(2)
        //     << iter << "  ";
        EoSctrl.SSMsta.curIt += iter;

        if (flag && Yt > 1 - 0.1 && Sk > 1) {
            // close to phase boundary, or more than 1 phase, So don't skip at next step
            flagSkip = OCP_FALSE;
        }
        if (flag && Yt > 1 + eYt) {
            EoSctrl.SSMsta.conflag = OCP_TRUE;
            EoSctrl.SSMsta.realTol = sqrt(Se);
            return OCP_FALSE;
        }
    }
    EoSctrl.SSMsta.realTol = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL MixtureComp::StableSSM(const USI& Id)
{
    // if unsatble, return OCP_FALSE
    // if stable, return OCP_TRUE

    const vector<OCP_DBL>& xj = x[Id];
    CalFugPhi(phi[Id], fug[Id], xj);
    const vector<OCP_DBL>& fugId = fug[Id];

    for (USI i = 0; i < NC; i++) {
        di[i] = phi[Id][i] * xj[i];
    }

    OCP_DBL  Stol  = EoSctrl.SSMsta.tol2;
    USI      maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL  eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL  Se;
    OCP_BOOL flag;
    USI      iter;
    USI      k;

    // cout << "SSMBEGINS" << endl;

    for (k = 0; k < 4; k++) {

        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Y[i] = xj[i] / Kw[k][i];
            Yt += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);
        // CalFugPhi(phiSta, fugSta, Y);
        CalFugPhi(&phiSta[0], &fugSta[0], &Y[0]);

        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
        }

        flag = OCP_TRUE;
        iter = 0;

        // cout << "ssmbegins" << endl;

        while (Se > Stol) {

            // cout << setprecision(6) << scientific << Se;
            // cout << "         ";
            // cout << setprecision(12) << scientific << Yt;
            // cout << "         " << iter << endl;
            // PrintDX(NC, &xj[0]);
            // PrintDX(NC, &Y[0]);

            Yt = 0;
            for (USI i = 0; i < NC; i++) {
                Y[i] = di[i] / phiSta[i];
                Yt += Y[i];
            }
            Dscalar(NC, 1 / Yt, &Y[0]);
            // CalFugPhi(phiSta, fugSta, Y);
            CalFugPhi(&phiSta[0], &fugSta[0], &Y[0]);
            Se = 0;
            for (USI i = 0; i < NC; i++) {
                Se += pow(log(fugSta[i] / fugId[i] * Yt), 2);
            }

            iter++;
            if (iter > maxIt) {
                flag = OCP_FALSE;
                break;
            }
        }
        EoSctrl.SSMsta.curIt += iter;
        flag = StableNR(Id);
        // here, a relaxation is needed, on the one hand it can prevent the influence
        // of rounding error, on the other hand, if Yt is too close to 1, phase
        // splitting calculation may get into trouble and single phase is indentified
        // finally
        if (flag && Yt > 1 + eYt) {
            // cout << "Yt = " << setprecision(12) << scientific << Yt << "    " <<
            // setw(3)
            //     << EoSctrl.SSMsta.curIt << "    " << setw(3) << EoSctrl.NRsta.curIt
            //     << "   " << flag << "   " << k << "   " << 2 << "   ";
            EoSctrl.SSMsta.conflag = OCP_TRUE;
            EoSctrl.SSMsta.realTol = sqrt(Se);
            return OCP_FALSE;
        }
    }
    // cout << "Yt = " << setprecision(12) << scientific << Yt << "    " << setw(3)
    //     << EoSctrl.SSMsta.curIt << "    " << setw(3) << EoSctrl.NRsta.curIt << "   "
    //     << flag << "   " << k << "   " << 1 << "   ";
    /*if (!flag) {
        OCP_WARNING("SSM not converged in Stability Analysis");
    }*/
    EoSctrl.SSMsta.realTol = sqrt(Se);
    return OCP_TRUE;
}

OCP_BOOL MixtureComp::StableNR(const USI& Id)
{

    for (USI i = 0; i < NC; i++) {
        resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
    }

    USI     maxIt = EoSctrl.NRsta.maxIt;
    OCP_DBL Stol  = EoSctrl.NRsta.tol;
    OCP_DBL Se    = Dnorm2(NC, &resSTA[0]);
    OCP_DBL alpha = 1;
    USI     iter  = 0;
    // OCP_DBL Se0   = Se;

    while (Se > Stol) {

        CalFugXSTA();
        AssembleJmatSTA();
        // LUSolve(1, NC, &JmatSTA[0], &resSTA[0], &pivot[0]);
        SYSSolve(1, &uplo, NC, &JmatSTA[0], &resSTA[0], &pivot[0], &JmatWork[0],
                 lJmatWork);
        Dscalar(NC, Yt, &Y[0]);
        Daxpy(NC, alpha, &resSTA[0], &Y[0]);
        Yt = 0;
        for (USI i = 0; i < NC; i++) {
            Yt += Y[i];
        }
        Dscalar(NC, 1 / Yt, &Y[0]);

        CalFugPhi(&fugSta[0], &Y[0]);
        for (USI i = 0; i < NC; i++) {
            resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
        }
        Se = Dnorm2(NC, &resSTA[0]);
        iter++;
        if (iter > maxIt) {
            EoSctrl.NRsta.curIt += iter;
            EoSctrl.NRsta.conflag = OCP_FALSE;
            EoSctrl.NRsta.realTol = Se;
            return OCP_FALSE;
        }
    }
    EoSctrl.NRsta.curIt += iter;
    EoSctrl.NRsta.conflag = OCP_TRUE;
    EoSctrl.NRsta.realTol = Se;
    return OCP_TRUE;
}

void MixtureComp::CalFugXSTA()
{
    // Y sums to be 1 now, it's actually the mole fraction of spliting phase
    // for stability analysis
    vector<OCP_DBL>& fugx = fugX[0];
    OCP_DBL          aj   = Asta;
    OCP_DBL          bj   = Bsta;
    OCP_DBL          zj   = Zsta;
    OCP_DBL          tmp  = 0;

    const OCP_DBL m1    = delta1;
    const OCP_DBL m2    = delta2;
    const OCP_DBL m1Pm2 = delta1P2;
    const OCP_DBL m1Mm2 = delta1M2;
    const OCP_DBL m1Tm2 = delta1T2;

    Bx = Bi;
    for (USI i = 0; i < NC; i++) {
        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += Y[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
        }
        Ax[i] = 2 * tmp;
        Zx[i] = ((bj - zj) * Ax[i] + ((aj + m1Tm2 * (3 * bj * bj + 2 * bj)) +
                                      ((m1Pm2) * (2 * bj + 1) - 2 * m1Tm2 * bj) * zj -
                                      (m1Pm2 - 1) * zj * zj) *
                                         Bx[i]) /
                (3 * zj * zj + 2 * ((m1Pm2 - 1) * bj - 1) * zj +
                 (aj + m1Tm2 * bj * bj - (m1Pm2)*bj * (bj + 1)));
    }

    OCP_DBL C, E, G;
    OCP_DBL Cxk, Dxk, Exk, Gxk;
    OCP_DBL aik;
    G = (zj + m1 * bj) / (zj + m2 * bj);

    for (USI i = 0; i < NC; i++) {

        C = Y[i] * P / (zj - bj);
        // C = 1 / (zj - bj);
        // D = Bx[i] * (zj - 1) / bj;
        E = -aj / ((m1Mm2)*bj) * (Ax[i] / aj - Bx[i] / bj);

        for (USI k = 0; k < NC; k++) {
            aik = (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);

            // Cxk = -Y[i] * (Zx[k] - Bx[k]) / ((zj - bj) * (zj - bj));
            Cxk = ((zj - bj) * delta(i, k) - Y[i] * (Zx[k] - Bx[k])) * P /
                  ((zj - bj) * (zj - bj));
            Dxk = Bx[i] / bj * (Zx[k] - Bx[k] * (zj - 1) / bj);
            /*Exk = (Ax[k] * bj - aj * Bx[k]) / (bj * bj) * (Ax[i] / aj - Bx[i] / bj) +
               aj / bj * (2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) / aj -
                Ax[k] * Ax[i] / (aj * aj) + Bx[i] * Bx[k] / (bj * bj));*/
            Exk = (2 * (aj / bj * Bx[k] * Bx[i] + bj * aik) - Ax[i] * Bx[k] -
                   Ax[k] * Bi[i]) /
                  (bj * bj);
            Exk /= -(m1Mm2);
            Gxk = (m1Mm2) / (zj + m2 * bj) / (zj + m2 * bj) * (zj * Bx[k] - Zx[k] * bj);
            fugx[i * NC + k] = 1 / C * Cxk + Dxk + Exk * log(G) + E / G * Gxk;
        }
    }
    // cout << "PhiX" << endl;
    // for (USI i = 0; i < NC; i++) {
    //	for (USI k = 0; k < NC; k++) {
    //		cout << fugx[i * NC + k] << "      ";
    //	}
    //	cout << endl;
    //}

#ifdef OCP_NANCHECK
    for (USI j = 0; j < NP; j++) {
        if (!CheckNan(fugX[j].size(), &fugX[j][0])) {
            OCP_ABORT("INF or NAN in fugX !");
        }
    }
#endif // NANCHECK
}

void MixtureComp::AssembleJmatSTA()
{
    vector<OCP_DBL>& fugx = fugX[0];
    fill(JmatSTA.begin(), JmatSTA.end(), 0.0);
    OCP_DBL tmp;
    for (USI i = 0; i < NC; i++) {

        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += Y[k] * fugx[i * NC + k];
        }

        for (USI j = 0; j < NC; j++) {
            // Symmetric
            // JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + delta(i, j) / Y[i]) / Yt;
            JmatSTA[i * NC + j] = (fugx[i * NC + j] - tmp + 1) / Yt;
        }
    }
    // cout << "JmatSTA" << endl;
    // for (USI i = 0; i < NC; i++) {
    //	for (USI k = 0; k < NC; k++) {
    //		cout << JmatSTA[k * NC + i] << "      ";
    //	}
    //	cout << endl;
    //}
}

void MixtureComp::PhaseSplit()
{
    EoSctrl.SSMsp.conflag = OCP_FALSE;
    EoSctrl.SSMsp.curIt   = 0;
    EoSctrl.NRsp.conflag  = OCP_FALSE;
    EoSctrl.NRsp.curIt    = 0;
    EoSctrl.RR.curIt      = 0;

    // cout << "begin" << endl;
    SplitSSM(OCP_FALSE);
    SplitNR();
    while (!EoSctrl.NRsp.conflag) {
        SplitSSM(OCP_TRUE);
        SplitNR();
        if (!CheckSplit()) break;
        if (EoSctrl.SSMsp.conflag) break;
    }
    CheckSplit();

    itersSSMSP += EoSctrl.SSMsp.curIt;
    itersNRSP += EoSctrl.NRsp.curIt;
    itersRR += EoSctrl.RR.curIt;
    countsSSMSP++;
    countsNRSP++;
    countsRR++;

    // cout << scientific << setprecision(8);
    // for (USI i = 0; i < NC; i++) {
    //     cout << x[0][i] / x[1][i] << "   ";
    // }
    // cout << endl;
    // cout << "Yt = " << scientific << setprecision(12) << Yt << "   "
    //     << setw(3) << "SSMtol = " << setprecision(6) << sqrt(EoSctrl.SSMsp.realTol)
    //     << "   "
    //     << setw(3) << EoSctrl.SSMsp.curIt << "   "
    //     << setw(3) << "NRtol = " << setprecision(6) << EoSctrl.NRsp.realTol << "   "
    //     << setw(3) << EoSctrl.NRsp.curIt << "   "
    //     << setw(2) << lNP << "   " << setw(2) << NP << "   "
    //     << (lNP == NP ? "N" : "Y") << "   ";
}

OCP_BOOL MixtureComp::CheckSplit()
{
    if (NP == 2) {

        OCP_DBL eX    = 0;
        OCP_DBL nuMax = max(nu[0], nu[1]);

        for (USI i = 0; i < NC; i++) {
            eX += (x[0][i] - x[1][i]) * (x[0][i] - x[1][i]);
        }

        if (OCP_TRUE) {
            // Calculate Gibbs Energy

            CalFugPhi(phiSta, fugSta, zi);
            GibbsEnergyB = 0;
            GibbsEnergyE = 0;
            for (USI i = 0; i < NC; i++) {
                GibbsEnergyB += zi[i] * log(fugSta[i]);
                GibbsEnergyE += (n[0][i] * log(fug[0][i]) + n[1][i] * log(fug[1][i]));
            }

            cout << scientific << setprecision(6);
            // cout << GibbsEnergyB << "   " << GibbsEnergyE << endl;
            if (GibbsEnergyE > GibbsEnergyB) {
                cout << ftype << "   ";
                cout << GibbsEnergyB << "   ";
                cout << GibbsEnergyE << "   ";
                cout << nuMax << "   ";
                cout << eX << "   ";
                cout << EoSctrl.NRsp.conflag << "   ";
                cout << bulkId << endl;
            }
        }

        if (nuMax < 1 && EoSctrl.NRsp.conflag && isfinite(eX)) {
            // accept this result
        } else {
            if (!isfinite(eX) || 1 - nuMax < 1E-3) {
                NP    = 1;
                x[0]  = zi;
                nu[0] = 1;
                CalAjBj(Aj[0], Bj[0], x[0]);
                SolEoS(Zj[0], Aj[0], Bj[0]);

                EoSctrl.SSMsta.conflag = OCP_FALSE;
                EoSctrl.NRsta.conflag  = OCP_FALSE;
                return OCP_FALSE;
            }
        }
    }
    return OCP_TRUE;
}

void MixtureComp::SplitSSM(const OCP_BOOL& flag)
{
    if (NP == 2) {
        SplitSSM2(flag);
    } else {
        SplitSSM3(flag);
    }
}

void MixtureComp::SplitSSM2(const OCP_BOOL& flag)
{
    // NP = 2 in this case
    // Ks is very IMPORTANT!
    // flag = OCP_TRUE : Restart SSM
    // flag = OCP_FALSE : New SSM
    EoSctrl.SSMsp.conflag = OCP_TRUE;
    OCP_DBL Se            = 1;
    OCP_DBL Stol          = EoSctrl.SSMsp.tol2;
    USI     maxIt         = EoSctrl.SSMsp.maxIt;

    if (!flag) {
        if (lNP == 2) {
            Ks[NP - 2] = lKs;
        } else {
            if (Yt < 1.1 || OCP_TRUE) {
                Ks[NP - 2] = Kw[0];
            } else {
                for (USI i = 0; i < NC; i++) {
                    // Ks[NP - 2][i] = phi[testPId][i] / phiSta[i];
                    Ks[NP - 2][i] = Y[i] / x[testPId][i];
                }
            }
        }
    } else {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    if (Yt < 1.1) {
        Stol *= 1E-1;
        maxIt *= 2;
    }

    USI iter = 0;
    while (Se > Stol) {

        RachfordRice2();
        UpdateXRR();
        CalFugPhiAll();
        Se = 0;
        for (USI i = 0; i < NC; i++) {
            Se += pow(fug[1][i] / fug[0][i] - 1, 2);
            Ks[0][i] = phi[1][i] / phi[0][i];
        }

        iter++;
        if (iter > maxIt) {
            // OCP_WARNING("SSM not converged in Phase Spliting!");
            EoSctrl.SSMsp.conflag = OCP_FALSE;
            break;
        }
    }

    EoSctrl.SSMsp.realTol = sqrt(Se);
    EoSctrl.SSMsp.curIt += iter;
}

void MixtureComp::SplitSSM3(const OCP_BOOL& flag) {}

void MixtureComp::RachfordRice2() ///< Used when NP = 2
{
    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;

    USI           iter  = 0;
    const OCP_DBL tol   = EoSctrl.RR.tol;
    const OCP_DBL maxIt = EoSctrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J  = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj  = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        } else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            } else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    EoSctrl.RR.curIt += iter;
    nu[1] = 1 - nu[0];

    // cout << scientific << setprecision(6) << nu[0] << "   " << nu[1] << endl;
}

void MixtureComp::RachfordRice2P()
{
    // modified RachfordRice equations
    // less iterations but more divergence --- unstable!

    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    const OCP_DBL numin = 1 / (1 - Kmax);
    const OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    OCP_DBL tmp, rj, J, dnuj, tmpnu;
    OCP_DBL f, df;

    USI           iter  = 0;
    const OCP_DBL tol   = EoSctrl.RR.tol;
    const OCP_DBL maxIt = EoSctrl.RR.maxIt;
    while (OCP_TRUE) {

        rj = 0;
        J  = 0;
        for (USI i = 0; i < NC; i++) {
            tmp = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmp;
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmp * tmp);
        }
        f  = (nu[0] - numin) * (numax - nu[0]);
        df = -2 * nu[0] + (numax + numin);
        J *= f;
        J += df * rj;
        rj *= f;

        if (fabs(rj) < tol || iter > maxIt) break;

        dnuj  = -rj / J;
        tmpnu = nu[0] + dnuj;
        if (tmpnu < numax && tmpnu > numin) {
            nu[0] = tmpnu;
        } else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            } else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }
        iter++;
    }

    EoSctrl.RR.curIt += iter;

    cout << iter << "      " << scientific << setprecision(3) << fabs(rj) << "      "
         << nu[0] << "      " << numin << "      " << numax << endl;

    nu[1] = 1 - nu[0];
}

void MixtureComp::RachfordRice3() ///< Used when NP > 2
{
}

void MixtureComp::UpdateXRR()
{
    OCP_DBL tmp = 0;
    for (USI i = 0; i < NC; i++) {
        tmp = 1;
        for (USI j = 0; j < NP - 1; j++) {
            tmp += nu[j] * (Ks[j][i] - 1);
        }
        x[NP - 1][i] = zi[i] / tmp;
        for (USI j = 0; j < NP - 1; j++) {
            x[j][i] = Ks[j][i] * x[NP - 1][i];
        }
    }
}

void MixtureComp::SplitBFGS()
{
    // Use BFGS to calculate phase splitting
    // JmatSP will store the BFGS mat
    // resSP will store the resSP - lresSP if necessary
    // n will store the n - ln if necessary

    // get initial value, n and ln, resSP and lresSP, H0

    // begin BFGS
}

void MixtureComp::SplitNR()
{
    EoSctrl.NRsp.conflag = OCP_FALSE;
    // for (USI j = 0; j < NP; j++) {
    //     nu[j] = fabs(nu[j]);
    // }

    USI len = NC * (NP - 1);
    x2n();
    CalResSP();
    OCP_DBL eNR0;
    OCP_DBL eNR   = Dnorm2(len, &resSP[0]);
    OCP_DBL NRtol = EoSctrl.NRsp.tol;
    OCP_DBL alpha;

    OCP_DBL en;
    USI     iter = 0;
    eNR0         = eNR;
    while (eNR > NRtol) {

        // eNR0 = eNR;
        ln = n;
        CalFugNAll();
        AssembleJmatSP();

        // LUSolve(1, len, &JmatSP[0], &resSP[0], &pivot[0]);
        // PrintDX(NC, &resSP[0]);

        int info = SYSSolve(1, &uplo, len, &JmatSP[0], &resSP[0], &pivot[0],
                            &JmatWork[0], lJmatWork);
        if (info > 0 && false) {
            for (USI i = 0; i < len; i++) {
                for (USI j = 0; j < len; j++) {
                    cout << scientific << setprecision(9) << JmatSP[i * len + j]
                         << "   ";
                }
                cout << ";\n";
            }
        }
        // PrintDX(NC, &resSP[0]);

        alpha = CalStepNRsp();

        n[NP - 1] = zi;
        for (USI j = 0; j < NP - 1; j++) {
            Daxpy(NC, alpha, &resSP[j * NC], &n[j][0]);
            Daxpy(NC, -1, &n[j][0], &n[NP - 1][0]);

            // nu[j] = Dnorm1(NC, &n[j][0]);
            nu[j] = 0;
            for (USI i = 0; i < NC; i++) {
                nu[j] += n[j][i];
            }

            for (USI i = 0; i < NC; i++) {
                // x[j][i] = n[j][i] / nu[j];
                x[j][i] = fabs(n[j][i] / nu[j]);
            }
        }

        // for (USI i = 0; i < NC; i++) {
        //     n[NP - 1][i] = fabs(n[NP - 1][i]);
        // }
        // nu[NP - 1] = Dnorm1(NC, &n[NP - 1][0]);
        // for (USI i = 0; i < NC; i++) {
        //     x[NP - 1][i] = n[NP - 1][i] / nu[NP - 1];
        // }
        nu[NP - 1] = 0;
        for (USI i = 0; i < NC; i++) {
            nu[NP - 1] += n[NP - 1][i];
        }
        for (USI i = 0; i < NC; i++) {
            x[NP - 1][i] = fabs(n[NP - 1][i] / nu[NP - 1]);
        }

        CalFugPhiAll();
        CalResSP();
        eNR = Dnorm2(len, &resSP[0]);
        iter++;
        if (eNR > eNR0 || iter > EoSctrl.NRsp.maxIt) {
            break;
        }

        // Maybe it should be execute before "eNR > eNR0 || iter > EoSctrl.NRsp.maxIt"
        en = 0;
        for (USI j = 0; j < NP; j++) {
            Daxpy(NC, -1, &n[j][0], &ln[j][0]);
            en += Dnorm2(NC, &ln[j][0]);
        }
        if (en / (NP * NC) < 1E-8) {
            EoSctrl.NRsp.conflag = OCP_TRUE;
            break;
        }
    }
    EoSctrl.NRsp.realTol = eNR;
    if (eNR < NRtol) EoSctrl.NRsp.conflag = OCP_TRUE;
    EoSctrl.NRsp.curIt += iter;

    // cout << iter << "   " << scientific << setprecision(3) << eNR << endl;
}

void MixtureComp::CalResSP()
{
    // So it equals -res
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            resSP[j * NC + i] = log(fug[NP - 1][i] / fug[j][i]);
        }
    }
}

void MixtureComp::CalFugNAll(const OCP_BOOL& Znflag)
{
    OCP_DBL C, E, G;
    OCP_DBL Cnk, Dnk, Enk, Gnk;
    OCP_DBL tmp, aik;

    for (USI j = 0; j < NP; j++) {
        // j th phase
        vector<OCP_DBL>&       fugn = fugN[j];
        const OCP_DBL&         aj   = Aj[j];
        const OCP_DBL&         bj   = Bj[j];
        const OCP_DBL&         zj   = Zj[j];
        const vector<OCP_DBL>& xj   = x[j];
        vector<OCP_DBL>&       Znj  = Zn[j];

        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }
            An[i] = 2 / nu[j] * (tmp - aj);
            Bn[i] = 1 / nu[j] * (Bi[i] - bj);
            if (Znflag) {
                Znj[i] =
                    ((bj - zj) * An[i] +
                     ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                      ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) *
                          zj -
                      (delta1 + delta2 - 1) * zj * zj) *
                         Bn[i]) /
                    (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                     (aj + delta1 * delta2 * bj * bj -
                      (delta1 + delta2) * bj * (bj + 1)));
            }
        }

        G = (zj + delta1 * bj) / (zj + delta2 * bj);

        for (USI i = 0; i < NC; i++) {
            // i th fugacity
            C = xj[i] * P / (zj - bj);
            // D = Bi[i] / bj * (zj - 1);
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
            }
            E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

            for (USI k = 0; k <= i; k++) {
                // k th components

                aik = (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);

                Cnk = P / (zj - bj) / (zj - bj) *
                      ((zj - bj) / nu[j] * (delta(i, k) - xj[i]) -
                       xj[i] * (Znj[k] - Bn[k]));
                Dnk = Bi[i] / bj * (Znj[k] - (Bi[k] - bj) * (zj - 1) / (nu[j] * bj));
                Gnk = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) *
                      (Bn[k] * zj - Znj[k] * bj);
                /*Enk = 1 / ((delta1 - delta2) * bj * bj) * (An[k] * bj - Bn[k] * aj) *
                   (Bi[i] / bj - 2 * tmp / aj)
                    + aj / ((delta1 - delta2) * bj) * (-Bi[i] / (bj * bj) * Bn[k] - 2 /
                   (aj * aj) * (aj * (aik - tmp) / nu[j] - An[k] * tmp));*/
                Enk = -1 / (delta1 - delta2) / (bj * bj) *
                      (2 * (bj * aik / nu[j] + Bn[k] * (Bi[i] * aj / bj - tmp)) -
                       An[k] * Bi[i] - aj * Bi[i] / nu[j]);
                Enk -= E / nu[j];
                fugn[i * NC + k] = 1 / C * Cnk + Dnk + Enk * log(G) + E / G * Gnk;
                fugn[k * NC + i] = fugn[i * NC + k];
                /*cout << "fnn[" << j << "][" << i << "][" << k << "] = " << fugn[i * NC
                + k]; cout << endl;*/
            }
        }
    }
    // PrintFugN();
#ifdef OCP_NANCHECK
    for (USI j = 0; j < NP; j++) {
        if (!CheckNan(fugN[j].size(), &fugN[j][0])) {
            OCP_ABORT("INF or NAN in fugN !");
        }
    }
#endif // NANCHECK
}

void MixtureComp::PrintFugN()
{
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            for (USI k = 0; k < NC; k++) {
                cout << fugN[j][i * NC + k] << "   ";
            }
            cout << endl;
        }
        cout << "---------------" << endl;
    }
}

void MixtureComp::AssembleJmatSP()
{
    // Dim: (NP-1)*NC
    // Attention that fugNj is sysmetric
    fill(JmatSP.begin(), JmatSP.end(), 0);

    OCP_DBL*       Jtmp = &JmatSP[0];
    const OCP_DBL* fugNp;
    const OCP_DBL* fugNj;

    for (USI j = 0; j < NP - 1; j++) {
        fugNp = &fugN[NP - 1][0];
        fugNj = &fugN[j][0];

        for (USI i = 0; i < NC; i++) {
            // ith components
            for (USI k = 0; k < NC; k++) {
                // kth fugacity
                Jtmp[k] = fugNj[k] + fugNp[k];
            }
            Jtmp += NC * (NP - 1);
            fugNp += NC;
            fugNj += NC;
        }
        Jtmp += NC;
    }

    // cout << endl << "Jmat" << endl;
    // for (USI i = 0; i < NC; i++) {
    //     for (USI j = 0; j < NC; j++) {
    //         cout << scientific << setprecision(6) << JmatSP[i * NC + j] << "   ";
    //     }
    //     cout << endl;
    // }
}

OCP_DBL MixtureComp::CalStepNRsp()
{
    OCP_DBL alpha = 1;
    OCP_DBL tmp;

    for (USI j = 0; j < NP - 1; j++) {

        const OCP_DBL* nj = &n[j][0];
        const OCP_DBL* r  = &resSP[j * NC];

        for (USI i = 0; i < NC; i++) {
            tmp = nj[i] + alpha * r[i];
            if (tmp < 0) {
                alpha *= 0.9 * fabs(nj[i] / r[i]);
            }
        }
    }
    return alpha;
}

void MixtureComp::AllocateOthers()
{
    sqrtMWi.resize(NC);
    for (USI i = 0; i < NC; i++) sqrtMWi[i] = sqrt(MWC[i]);
    muC.resize(NPmax);
    muAux.resize(NPmax);
    for (USI i = 0; i < NPmax; i++) {
        muAux[i].resize(5);
    }
    muAux1I.resize(NC);
    for (USI i = 0; i < NC; i++) {
        muAux1I[i] = 5.4402 * pow(Tc[i], 1.0 / 6) / pow(Pc[i], 2.0 / 3);
    }
    fugP.resize(NPmax);
    for (USI j = 0; j < NPmax; j++) {
        fugP[j].resize(NC);
    }
    Zp.resize(NPmax);
    JmatDer.resize(NPmax * NPmax * (NC + 1) * (NC + 1));
    JmatTmp = JmatDer;
    rhsDer.resize(NPmax * (NC + 1) * (NC + 1));
    xixC.resize(NPmax * NC);
    xiPC.resize(NPmax);
    xiNC.resize(NPmax * NC);

    // new
    JmatDer.resize((numPhase + NPmax * NC) * (numPhase + NPmax * NC));
    JmatTmp = JmatDer;
    rhsDer.resize((numPhase + NPmax * NC) * (numCom + 1 + 1));
}

void MixtureComp::IdentifyPhase()
{
    phaseExist[0] = OCP_FALSE;
    phaseExist[1] = OCP_FALSE;
    if (NP == 1) {
        // Critical Temperature Method
        OCP_DBL A = 0;
        OCP_DBL B = 0;
        for (USI i = 0; i < NC; i++) {
            A += x[0][i] * Vc[i] * Tc[i];
            B += x[0][i] * Vc[i];
        }
        OCP_DBL Tc = A / B;
        if (T > Tc) {
            phaseLabel[0] = GAS;
            phaseExist[1] = OCP_TRUE;
        } else {
            phaseLabel[0] = OIL;
            phaseExist[0] = OCP_TRUE;
        }
    } else {
        // Compare MW
        if (MW[0] > MW[1]) {
            phaseLabel[0] = OIL;
            phaseLabel[1] = GAS;
        } else {
            phaseLabel[0] = GAS;
            phaseLabel[1] = OIL;
        }
        phaseExist[0] = OCP_TRUE;
        phaseExist[1] = OCP_TRUE;
    }
}

void MixtureComp::CopyPhase()
{
    // copy vj, x, mu, xi, rho
    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];
        vj[j1]       = vC[j];
        Dcopy(NC, &xij[j1 * numCom], &x[j][0]);
        mu[j1]  = muC[j];
        xi[j1]  = xiC[j];
        rho[j1] = rhoC[j];
        nj[j1]  = nu[j];
    }
}

void MixtureComp::CalViscosity() { CalViscoLBC(); }

void MixtureComp::CalViscoLBC()
{
    OCP_DBL tmp;
    OCP_DBL Tri;
    OCP_DBL xijT;
    OCP_DBL xijP;
    OCP_DBL xijV;

    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj  = x[j];
        vector<OCP_DBL>&       muA = muAux[j];
        fill(muA.begin(), muA.end(), 0.0);
        xijT = 0;
        xijP = 0;
        xijV = 0;

        for (USI i = 0; i < NC; i++) {
            tmp = 5.4402 * pow(Tc[i], 1.0 / 6) / sqrt(MW[j]) / pow(Pc[i], 2.0 / 3);
            // tmp = muAux1I[i] / sqrt(MW[j]);
            Tri = T / Tc[i];
            if (Tri <= 1.5) {
                tmp = 34 * 1E-5 * pow(Tri, 0.94) / tmp;
            } else {
                tmp = 17.78 * 1E-5 * pow((4.58 * Tri - 1.67), 0.625) / tmp;
            }
            muA[0] += xj[i] * sqrt(MWC[i]) * tmp;
            muA[1] += xj[i] * sqrt(MWC[i]);
            xijT += xj[i] * Tc[i];
            xijP += xj[i] * Pc[i];
            xijV += xj[i] * Vcvis[i];
        }
        muA[2] = 5.4402 * pow(xijT, 1.0 / 6) / sqrt(MW[j]) / pow(xijP, 2.0 / 3);
        muA[3] = xiC[j] * xijV;

        if (muA[3] <= 0.18 && OCP_FALSE) {
            muC[j] = muA[0] / muA[1] + 2.05 * 1E-4 * muA[3] / muA[2];
        } else {
            muA[4] = muA[3] * (muA[3] * (muA[3] * (LBCcoef[4] * muA[3] + LBCcoef[3]) +
                                         LBCcoef[2]) +
                               LBCcoef[1]) +
                     LBCcoef[0];
            muC[j] = muA[0] / muA[1] + 1E-4 * (pow(muA[4], 4) - 1) / muA[2];
        }
    }
}

void MixtureComp::CalViscoHZYT() {}

void MixtureComp::CalFugXAll()
{
    for (USI j = 0; j < NP; j++) {

        vector<OCP_DBL>& fugx = fugX[j];
        vector<OCP_DBL>& xj   = x[j];
        OCP_DBL          aj   = Aj[j];
        OCP_DBL          bj   = Bj[j];
        OCP_DBL          zj   = Zj[j];
        OCP_DBL          tmp  = 0;

        Bx = Bi;
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += xj[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
            }
            Ax[i] = 2 * tmp;
            Zx[i] =
                ((bj - zj) * Ax[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bx[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        }

        OCP_DBL C, E, G;
        OCP_DBL Cxk, Dxk, Exk, Gxk;
        OCP_DBL aik;
        G = (zj + delta1 * bj) / (zj + delta2 * bj);

        for (USI i = 0; i < NC; i++) {
            // ith fugacity
            C = xj[i] * P / (zj - bj);
            // C = 1 / (zj - bj);
            // D = Bx[i] * (zj - 1) / bj;
            E = -aj / ((delta1 - delta2) * bj) * (Ax[i] / aj - Bx[i] / bj);

            for (USI k = 0; k < NC; k++) {
                aik = (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);

                // kth components
                Cxk = ((zj - bj) * delta(i, k) - xj[i] * (Zx[k] - Bx[k])) * P /
                      ((zj - bj) * (zj - bj));
                Dxk = Bx[i] / bj * (Zx[k] - Bx[k] * (zj - 1) / bj);
                /*Exk = (Ax[k] * bj - aj * Bx[k]) / (bj * bj) * (Ax[i] / aj - Bx[i] /
                   bj) + aj / bj * (2 * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) / aj
                   - Ax[k] * Ax[i] / (aj * aj) + Bx[i] * Bx[k] / (bj * bj));*/
                Exk = (2 * (aj / bj * Bx[k] * Bx[i] + bj * aik) - Ax[i] * Bx[k] -
                       Ax[k] * Bi[i]) /
                      (bj * bj);
                Exk /= -(delta1 - delta2);
                Gxk = (delta1 - delta2) / (zj + delta2 * bj) / (zj + delta2 * bj) *
                      (zj * Bx[k] - Zx[k] * bj);
                fugx[i * NC + k] = 1 / C * Cxk + Dxk + Exk * log(G) + E / G * Gxk;
            }
        }
    }
#ifdef OCP_NANCHECK
    for (USI j = 0; j < NP; j++) {
        if (!CheckNan(fugX[j].size(), &fugX[j][0])) {
            OCP_ABORT("INF or NAN in fugX !");
        }
    }
#endif // NANCHECK
}

void MixtureComp::CalFugPAll(const OCP_BOOL& Zpflag)
{

    OCP_DBL C, E, G;
    OCP_DBL Cp, Dp, Gp;
    OCP_DBL tmp;

    for (USI j = 0; j < NP; j++) {

        vector<OCP_DBL>& fugp = fugP[j];
        vector<OCP_DBL>& xj   = x[j];
        OCP_DBL&         aj   = Aj[j];
        OCP_DBL&         bj   = Bj[j];
        OCP_DBL&         zj   = Zj[j];

        OCP_DBL Ap = aj / P;
        OCP_DBL Bp = bj / P;
        if (Zpflag) {
            Zp[j] =
                ((bj - zj) * Ap +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bp) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        }

        G  = (zj + delta1 * bj) / (zj + delta2 * bj);
        Gp = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) *
             (Bp * zj - Zp[j] * bj);
        for (USI i = 0; i < NC; i++) {

            C = P / (zj - bj);
            // D = Bi[i] / bj * (zj - 1);

            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }

            E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

            Cp = ((zj - bj) - P * (Zp[j] - Bp)) / ((zj - bj) * (zj - bj));
            Dp = Bi[i] / bj * Zp[j];
            // Ep = 0;

            // Attention that if xj[i] = 0, then fugp[i] = nan
            // but Cp also contains xj[i], so eliminate it first
            fugp[i] = Cp / C + Dp + E / G * Gp;
        }
    }

#ifdef OCP_NANCHECK
    for (USI j = 0; j < NP; j++) {
        if (!CheckNan(fugP[j].size(), &fugP[j][0])) {
            OCP_ABORT("INF or NAN in fugP !");
        }
    }
#endif // NANCHECK
}

void MixtureComp::CalVjpVfpVfn_partial()
{
    // Vfi = 0
    // Vjp, Vfp
    // Vjn(dVj / dnij) in Vji
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    fill(vfi.begin(), vfi.end(), 0.0);
    vfP = 0;

    for (USI j = 0; j < NP; j++) {
        const OCP_DBL&         aj  = Aj[j];
        const OCP_DBL&         bj  = Bj[j];
        const OCP_DBL&         zj  = Zj[j];
        const vector<OCP_DBL>& xj  = x[j];
        vector<OCP_DBL>&       Znj = Zn[j];
        OCP_DBL                tmp;
        const USI              j1 = phaseLabel[j];
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }
            An[i] = 2 / nu[j] * (tmp - aj);
            Bn[i] = 1 / nu[j] * (Bi[i] - bj);
            Znj[i] =
                ((bj - zj) * An[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bn[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
            vji[j1][i] = CgTP * (zj + nu[j] * Znj[i]) - Vshift[i];
        }
        Zp[j] = ((bj - zj) * aj +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     bj) /
                P /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        vjp[j1] = CgTP * nu[j] * (Zp[j] - zj / P);
        vfP += vjp[j1];
    }
}

void MixtureComp::CalXiPn_partial()
{
    // xiN = 0
    // xiP
    // xin(dxi_j / d nij) in xix
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    OCP_DBL tmp;
    fill(xiN.begin(), xiN.end() - numCom, 0.0);

    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];

        tmp = -xi[j1] * xi[j1] * CgTP;
        for (USI i = 0; i < NC; i++) {
            xix[j1 * numCom + i] = tmp * Zn[j][i];
        }
        xiP[j1] = tmp * (Zp[j] - Zj[j] / P);
    }
}

void MixtureComp::CalRhoPn_partial()
{
    // rhoN = 0
    // rhoP
    // rhon(drhoj / d nij) in rhox
    fill(rhoN.begin(), rhoN.end() - numCom, 0.0);
    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1       = phaseLabel[j];
        rhoP[j1] = xiP[j1] * MW[j];
        for (USI m = 0; m < NC; m++) {
            rhox[j1 * numCom + m] =
                xix[j1 * numCom + m] * MW[j] + xi[j1] / nu[j] * (MWC[m] - MW[j]);
        }
    }
}

void MixtureComp::CalMuPn_partial() { CalMuPnLBC_partial(); }

void MixtureComp::CalMuPnLBC_partial()
{
    // muN = 0
    // muP
    // mun in mux

    fill(muN.begin(), muN.end() - numCom, 0.0);

    OCP_DBL val1IJ, val2IJ;
    OCP_DBL der1IJ, der2IJ, der3J, der4J, der6J, der7J, der8J;
    OCP_DBL Tri, tmp;
    OCP_DBL xTj, xPj, xVj;
    OCP_DBL derxTj, derxPj, derMWj;

    // dxij / dP = 0, d MJ / dP = 0
    // Calculate dmuj / dP
    // der2IJ = der3J = der4J = der6J = 0;

    for (USI j = 0; j < NP; j++) {
        const USI              j1     = phaseLabel[j];
        const vector<OCP_DBL>& xj     = x[j];
        const vector<OCP_DBL>& muAuxj = muAux[j];

        xTj = xPj = xVj = 0;
        for (USI i = 0; i < NC; i++) {
            xTj += xj[i] * Tc[i];
            xPj += xj[i] * Pc[i];
            xVj += xj[i] * Vcvis[i];
        }
        der7J = xVj * xiP[j1];
        if (muAuxj[3] <= 0.18 && OCP_FALSE) {
            muP[j1] = (2.05 * 1E-4) * der7J / muAuxj[2];
        } else {
            der8J   = der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
            muP[j1] = (4 * 1E-4) * pow(muAuxj[4], 3) * der8J / muAuxj[2];
        }

        // Calculate dmuj / nkj
        const USI bId = numCom * j1;
        for (USI k = 0; k < NC; k++) {
            derxTj = (Tc[k] - xTj) / nu[j];
            derxPj = (Pc[k] - xPj) / nu[j];
            derMWj = (MWC[k] - MW[j]) / nu[j];
            der3J  = 0;
            der4J  = sqrtMWi[k];
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[j]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[j], -1.5) * derMWj;
                Tri    = T / Tc[i];
                if (Tri <= 1.5) {
                    tmp = 34 * 1E-5 * pow(Tri, 0.94);
                } else {
                    tmp = 17.78 * 1E-5 * pow(4.58 * Tri - 1.67, 0.625);
                }
                val2IJ = tmp / val1IJ;
                der2IJ = -tmp * der1IJ / (val1IJ * val1IJ);
                der3J += sqrtMWi[i] *
                         (xj[i] * der2IJ + (delta(i, k) - xj[i]) * val2IJ / nu[j]);
                der4J -= xj[i] * sqrtMWi[i];
            }
            der4J /= nu[j];
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            der7J = xix[j1 * numCom + k] * xVj + (Vcvis[k] - xVj) * xi[j1] / nu[j];
            if (muAuxj[3] <= 0.18 && OCP_FALSE) {
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }
    }
}

void MixtureComp::CalVjpVfpVfx_partial()
{
    // Wrong Now !

    // Vfi = 0
    // Vjp, Vfp
    // Vjn(dVj / dxij) in Vji
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    fill(vfi.begin(), vfi.end(), 0.0);
    vfP = 0;

    for (USI j = 0; j < NP; j++) {
        const USI              j1  = phaseLabel[j];
        const vector<OCP_DBL>& xj  = x[j];
        const OCP_DBL          aj  = Aj[j];
        const OCP_DBL          bj  = Bj[j];
        const OCP_DBL          zj  = Zj[j];
        OCP_DBL                tmp = 0;

        Bx = Bi;
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += xj[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
            }
            Ax[i] = 2 * tmp;
            Zx[i] =
                ((bj - zj) * Ax[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bx[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));

            vji[j1][i] = CgTP * nu[j] * Zx[i];
        }
        Zp[j] = ((bj - zj) * aj +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     bj) /
                P /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        vjp[j1] = CgTP * nu[j] * (Zp[j] - zj / P);
        vfP += vjp[j1];
    }
}

void MixtureComp::CaldXsdXpAPI04()
{
    // Wrong Now !

    // water not in oil or gas; hydroncarbon not in water
    // inexist phase will bot be stored
    // only hydrocarbon phases, hydrocarbon components, water phase will be stored

    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);
    fill(res.begin(), res.end(), 0.0);
    resPc = 0;

    const USI     ncol = numCom + 1;
    const OCP_DBL vf2  = vf * vf;
    const OCP_DBL vwp  = vjp[numPhase - 1];
    const OCP_DBL vw   = vj[numPhase - 1];

    if (NP == 1) {
        // when NP = 1, vfi = vhi, vfP = vhp + vwp, h = hydrocarbon
        const USI j1  = phaseLabel[0];
        OCP_DBL*  bId = &dXsdXp[0];
        // dS / dP
        bId[0] = ((vfP - vwp) * vf - vfP * vj[j1]) / vf2;
        bId++;
        // dSh / dNm
        for (USI m = 0; m < NC; m++) {
            bId[m] = vfi[m] * vw / vf2;
        }
        bId += NC;
        // dSh / dNw
        bId[0] = -vfi[numCom - 1] * vj[j1] / vf2;
        bId++;
        // dSw / dP, dNm
        for (USI m = 0; m < ncol; m++) {
            bId[m] = -dXsdXp[m];
        }
        bId += ncol + 1;
        // dxij / dNm
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += ncol;
        }
        // res = 0
    } else {
        CalFugXAll();
        CalFugPAll(OCP_FALSE);
        CaldXsdXp04();
    }
}

void MixtureComp::CaldXsdXp04()
{
    // Wrong Now !

    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    OCP_DBL* MTmp = &JmatTmp[0];

    // -dF / dXs
    // -dFn / dXs
    const USI dim = NP + 1 + NP * NC;
    for (USI i = 0; i < NC; i++) {
        // -dFn / dSh
        for (USI m = 0; m < NP; m++) {
            MTmp[m] = vf * xi[phaseLabel[m]] * x[m][i];
        }
        MTmp += NP;
        // -dFn / dSw
        MTmp++;
        // -dFn / dxnm
        for (USI m = 0; m < NP; m++) {
            const USI m1 = phaseLabel[m];
            for (USI n = 0; n < NC; n++) {
                MTmp[n] = vf * S[m1] *
                          (xix[m1 * numCom + n] * x[m][i] + delta(i, n) * xi[m1]);
                for (USI j = 0; j < NP; j++) {
                    MTmp[n] +=
                        vji[m1][n] * x[m][i] * S[phaseLabel[j]] * xi[phaseLabel[j]];
                }
            }
            MTmp += NC;
        }
    }
    // -dFnw / dFs
    // -dFnw / dSw
    MTmp[numPhase - 1] = vf * xi[numPhase - 1];
    MTmp += (NP + 1);
    // -dFnw / dxnm
    for (USI m = 0; m < NP; m++) {
        for (USI n = 0; n < NC; n++) {
            MTmp[n] = vji[phaseLabel[m]][n] * xi[numPhase - 1] * S[numPhase - 1];
        }
        MTmp += NC;
    }

    // -dFx / dXs
    for (USI j = 0; j < NP; j++) {
        // -dFx / dSm
        MTmp += (NP + 1);
        // -dFx / dxnm, m = j
        MTmp += j * NC;
        for (USI n = 0; n < NC; n++) {
            MTmp[n] = -1.0;
        }
        MTmp += (NP - j) * NC;
    }

    // -dFf / dXs
    for (USI j = 0; j < NP - 1; j++) {
        const OCP_DBL* fugXJ  = &fugX[j][0];
        const OCP_DBL* fugXNP = &fugX[NP - 1][0];
        for (USI i = 0; i < NC; i++) {
            // -dFf / dSm
            MTmp += (NP + 1);
            // -dFf / dxnm, m = j or m = np-1
            // m = j
            MTmp += j * NC;
            for (USI n = 0; n < NC; n++) {
                MTmp[n] = -fugXJ[n];
            }
            fugXJ += NC;
            // m = np-1
            MTmp += (NP - 1 - j) * NC;
            for (USI n = 0; n < NC; n++) {
                MTmp[n] = fugXNP[n];
            }
            fugXNP += NC;
            MTmp += NC;
        }
    }

    cout << endl << "dFdS" << endl;
    cout << scientific << setprecision(1);
    for (USI i = 0; i < dim; i++) {
        for (USI j = 0; j < dim; j++) cout << setw(8) << JmatTmp[i * dim + j] << " ";
        cout << endl;
    }

    // Transpose JmatTmp

    for (USI i = 0; i < dim; i++) {
        for (USI j = 0; j < dim; j++) {
            JmatDer[i * dim + j] = JmatTmp[j * dim + i];
        }
    }

    // dF / dXp
    // d fij / dP have been calculated in CalVfiVfp_full01() before
    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    MTmp = &JmatTmp[0];
    OCP_DBL   tmp01;
    OCP_DBL   tmp02;
    USI       j1;
    const USI nrhs = numCom + 1;
    const USI nrow = dim;
    // dFn / dXp
    for (USI i = 0; i < NC; i++) {
        // dFn / dP
        tmp01 = 0;
        tmp02 = 0;
        for (USI j = 0; j < NP; j++) {
            j1 = phaseLabel[j];
            tmp01 += S[j1] * x[j][i] * xiP[j1];
            tmp02 += S[j1] * x[j][i] * xi[j1];
        }
        MTmp[0] = -(vf * tmp01 + vfP * tmp02);
        MTmp += 1;
        // dFn / dNk
        MTmp[i] = 1;
        MTmp += NC;
        // dFn / dNw
        MTmp[0] = -vfi[numCom - 1] * tmp02;
        MTmp++;
    }
    // dFnw / dXp
    // dFnw / dP
    MTmp[0] = -S[numPhase - 1] * (vfP * xi[numPhase - 1] + vf * xiP[numPhase - 1]);
    // dFnw / dNw
    MTmp[numCom] = 1 - vfi[numCom - 1] * S[numPhase - 1] * xi[numPhase - 1];
    MTmp += nrhs;

    // dFx / dXp
    MTmp += NP * nrhs;

    // dFf / dXp
    const vector<OCP_DBL>& fugPNP = fugP[NP - 1];
    for (USI j = 0; j < NP - 1; j++) {
        const vector<OCP_DBL>& fugPj = fugP[j];
        for (USI i = 0; i < NC; i++) {
            // dFf / dP
            MTmp[0] = fugPj[i] - fugPNP[i];
            MTmp += nrhs;
        }
    }

    cout << endl << "dFdp" << endl;
    cout << scientific << setprecision(3);
    for (USI i = 0; i < dim; i++) {
        for (USI j = 0; j < numCom + 1; j++)
            cout << setw(10) << JmatTmp[i * (numCom + 1) + j] << "   ";
        cout << endl;
    }

    // Transpose JmatTmp
    for (USI i = 0; i < nrhs; i++) {
        for (USI j = 0; j < nrow; j++) {
            rhsDer[i * nrow + j] = JmatTmp[j * nrhs + i];
        }
    }

    LUSolve(nrhs, dim, &JmatDer[0], &rhsDer[0], &pivot[0]);

    cout << setprecision(6);
    cout << endl << "dXsdXp" << endl;
    for (USI i = 0; i < nrow; i++) {
        for (USI j = 0; j < nrhs; j++) {
            cout << setw(13) << rhsDer[j * nrow + i] << "   ";
        }
        cout << endl;
    }
    cout << endl;
}

void MixtureComp::CalXiPNX_partial()
{
    // Here!
    // dxi/dP = dxi/dP
    // dxi/dNk = dxi/dNk = 0
    // dxij / dP = dxij / dNk = 0
    // See MixtureComp::CalXiPNX_full01()

    // Calculate xix and xiP
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    OCP_DBL tmp;

    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj = x[j];
        OCP_DBL                aj = Aj[j];
        OCP_DBL                bj = Bj[j];
        OCP_DBL                zj = Zj[j];

        Bx = Bi;
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += xj[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
            }
            Ax[i] = 2 * tmp;
            Zx[i] =
                ((bj - zj) * Ax[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bx[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        }

        tmp = -xiC[j] * xiC[j];
        for (USI i = 0; i < NC; i++) {
            xixC[j * NC + i] = tmp * (CgTP * Zx[i] - Vshift[i]);
        }
        xiPC[j] = tmp * CgTP * (Zp[j] - Zj[j] / P);
    }

    // xiNC
    fill(xiNC.begin(), xiNC.end(), 0.0);

    // Copoy to xiP, xix
    // fill(xiP.begin(), xiP.end(), 0.0);
    // fill(xix.begin(), xix.end(), 0.0);
    fill(xiN.begin(), xiN.end(), 0.0);
    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1      = phaseLabel[j];
        xiP[j1] = xiPC[j];
        Dcopy(NC, &xix[j1 * numCom], &xixC[j * NC]);
    }
}

void MixtureComp::CalRhoPX_partial()
{
    // Here!
    // drho/dP = drho/dP
    // drho/dNk = drho/dNk = 0
    // dxij / dP = dxij / dNk = 0
    // See MixtureComp::CalRhoPNX_full()

    // fill(rhox.begin(), rhox.end() - numCom, 0.0);
    fill(rhoN.begin(), rhoN.end() - numCom, 0.0);
    // fill(rhoP.begin(), rhoP.end() - 1, 0.0);

    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1       = phaseLabel[j];
        rhoP[j1] = xiP[j1] * MW[j];
        for (USI i = 0; i < NC; i++) {
            rhox[j1 * numCom + i] = xix[j1 * numCom + i] * MW[j] + xi[j1] * MWC[i];
        }
    }
}

void MixtureComp::CalMuPX_partial()
{
    // fill(muP.begin(), muP.end() - 1, 0.0);
    fill(muN.begin(), muN.end() - numCom, 0.0);
    // fill(mux.begin(), mux.end() - numCom, 0.0);

    CalMuPXLBC_partial();
}

void MixtureComp::CalMuPXLBC_partial()
{
    // Here!
    // dmu/dP = dmu/dP
    // dmu/dNk = dmu/dNk = 0
    // // dxij / dP = dxij / dNk = 0
    // See MixtureComp::CalMuPXLBC_full01()

    OCP_DBL val1IJ, val2IJ;
    OCP_DBL der1IJ, der2IJ, der3J, der4J, der6J, der7J, der8J;
    OCP_DBL Tri, tmp;
    OCP_DBL xTj, xPj, xVj;
    OCP_DBL derxTj, derxPj, derMWj;

    // dxij / dP = 0, d MJ / dP = 0
    // Calculate dmuj / dP
    // der2IJ = der3J = der4J = der6J = 0;

    for (USI j = 0; j < NP; j++) {
        const USI              j1     = phaseLabel[j];
        const vector<OCP_DBL>& xj     = x[j];
        const vector<OCP_DBL>& muAuxj = muAux[j];

        xTj = xPj = xVj = 0;
        for (USI i = 0; i < NC; i++) {
            xTj += xj[i] * Tc[i];
            xPj += xj[i] * Pc[i];
            xVj += xj[i] * Vcvis[i];
        }
        der7J = xVj * xiP[j1];
        if (muAuxj[3] <= 0.18 && OCP_FALSE) {
            muP[j1] = (2.05 * 1E-4) * der7J / muAuxj[2];
        } else {
            der8J   = der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
            muP[j1] = (4 * 1E-4) * pow(muAuxj[4], 3) * der8J / muAuxj[2];
        }

        // Calculate dmuj / xkj
        const USI bId = numCom * j1;
        for (USI k = 0; k < NC; k++) {
            derxTj = Tc[k];
            derxPj = Pc[k];
            derMWj = MWC[k];
            der3J  = 0;
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[j]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[j], -1.5) * derMWj;
                Tri    = T / Tc[i];
                if (Tri <= 1.5) {
                    tmp = 34 * 1E-5 * pow(Tri, 0.94);
                } else {
                    tmp = 17.78 * 1E-5 * pow(4.58 * Tri - 1.67, 0.625);
                }
                val2IJ = tmp / val1IJ;
                der2IJ = -tmp * der1IJ / (val1IJ * val1IJ);
                der3J +=
                    xj[i] * sqrtMWi[i] * der2IJ + delta(i, k) * sqrtMWi[k] * val2IJ;
            }
            der4J = sqrtMWi[k];
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            der7J = xix[j1 * numCom + k] * xVj + xi[j1] * Vcvis[k];
            if (muAuxj[3] <= 0.18 && OCP_FALSE) {
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }
    }
}

void MixtureComp::CalXiPNX_full01()
{
    // call CalVfiVfp_full01() before, use rhsDer
    // Calculate xiPC, xiNC, xixC
    // Calculate xixC first
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    OCP_DBL tmp;
    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj = x[j];
        OCP_DBL                aj = Aj[j];
        OCP_DBL                bj = Bj[j];
        OCP_DBL                zj = Zj[j];

        Bx = Bi;
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += xj[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
            }
            Ax[i] = 2 * tmp;
            Zx[i] =
                ((bj - zj) * Ax[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bx[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        }

        tmp = -xiC[j] * xiC[j] * CgTP;
        for (USI i = 0; i < NC; i++) {
            xixC[j * NC + i] = tmp * Zx[i];
        }
    }
    // Calculate xiPC and xiNC
    if (NP == 1) {
        // NP = 1
        tmp     = -xiC[0] * xiC[0] * CgTP;
        xiPC[0] = tmp * (Zp[0] - Zj[0] / P);
        for (USI i = 0; i < NC; i++) {
            xiNC[i] = tmp * Zn[0][i];
        }
    } else {
        // NP > 1
        const OCP_DBL* dnkjdNP = &rhsDer[0];
        OCP_DBL        dertmp;

        // Calculate xiNC
        // Now dnkjdNP reach to dnkj / dNi after CalVfiVfp_full01
        for (USI i = 0; i < NC; i++) {
            for (USI j = 0; j < NP; j++) {

                dertmp                     = 0;
                tmp                        = -xiC[j] * xiC[j] * CgTP;
                const vector<OCP_DBL>& Znj = Zn[j];

                for (USI k = 0; k < NC; k++) {
                    dertmp += Znj[k] * dnkjdNP[k];
                }
                xiNC[j * NC + i] = tmp * dertmp;
                dnkjdNP += NC;
            }
        }
        // Calculate xiPC
        // Now dnkjdNP reach to dnkj / dP
        for (USI j = 0; j < NP; j++) {
            tmp                        = -xiC[j] * xiC[j] * CgTP;
            xiPC[j]                    = Zp[j] - Zj[j] / P;
            const vector<OCP_DBL>& Znj = Zn[j];

            // in OCP
            for (USI k = 0; k < NC; k++) {
                xiPC[j] += Znj[k] * dnkjdNP[k];
            }
            dnkjdNP += NC;
            xiPC[j] *= tmp;
        }
    }

    // Copoy to xiP, xix, xiN
    // fill(xiP.begin(), xiP.end(), 0.0);
    // fill(xix.begin(), xix.end(), 0.0);
    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1      = phaseLabel[j];
        xiP[j1] = xiPC[j];
        Dcopy(NC, &xiN[j1 * numCom], &xiNC[j * NC]);
        Dcopy(NC, &xix[j1 * numCom], &xixC[j * NC]);
    }
}

void MixtureComp::CalXiPNX_full02()
{
    // Attention!
    // NP = 1 or NP = 2

    // Call CalVfiVfp_full02() before, use rhsder
    // Calculate xiPC, xiNC, xixC
    // Calculate xixC first
    OCP_DBL CgTP = GAS_CONSTANT * T / P;
    OCP_DBL tmp;
    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj = x[j];
        OCP_DBL                aj = Aj[j];
        OCP_DBL                bj = Bj[j];
        OCP_DBL                zj = Zj[j];

        Bx = Bi;
        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI k = 0; k < NC; k++) {
                tmp += xj[k] * (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
            }
            Ax[i] = 2 * tmp;
            Zx[i] =
                ((bj - zj) * Ax[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bx[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        }

        tmp = -xiC[j] * xiC[j] * CgTP;
        for (USI i = 0; i < NC; i++) {
            xixC[j * NC + i] = tmp * Zx[i];
        }
    }

    // Calculate xiPC and xiNC
    if (NP == 1) {
        // NP = 1
        tmp     = -xiC[0] * xiC[0] * CgTP;
        xiPC[0] = tmp * (Zp[0] - Zj[0] / P);
        for (USI i = 0; i < NC; i++) {
            xiNC[i] = tmp * Zn[0][i];
        }
    } else {
        // NP = 2
        OCP_DBL* nijPN = &rhsDer[0];
        // Calculate xiPC
        xiPC[0] = (Zp[0] - Zj[0] / P);
        xiPC[1] = (Zp[1] - Zj[1] / P);
        for (USI i = 0; i < NC; i++) {
            xiPC[0] += Zn[0][i] * nijPN[i];
            xiPC[1] -= Zn[1][i] * nijPN[i];
        }
        nijPN += NC;
        xiPC[0] *= -xiC[0] * xiC[0] * CgTP;
        xiPC[1] *= -xiC[1] * xiC[1] * CgTP;
        // Calculate xiNC
        OCP_DBL tmp0 = -xiC[0] * xiC[0] * CgTP;
        OCP_DBL tmp1 = -xiC[1] * xiC[1] * CgTP;
        for (USI i = 0; i < NC; i++) {
            xiNC[i]      = 0;
            xiNC[NC + i] = 0;
            for (USI k = 0; k < NC; k++) {
                xiNC[i] += Zn[0][k] * nijPN[k];
                xiNC[NC + i] -= Zn[1][k] * nijPN[k];
            }
            xiNC[NC + i] += Zn[1][i];
            nijPN += NC;
            xiNC[i] *= tmp0;
            xiNC[NC + i] *= tmp1;
        }
    }

    // Copoy to xiP, xix, xiN
    // fill(xiP.begin(), xiP.end(), 0.0);
    // fill(xix.begin(), xix.end(), 0.0);
    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1      = phaseLabel[j];
        xiP[j1] = xiPC[j];
        Dcopy(NC, &xiN[j1 * numCom], &xiNC[j * NC]);
        Dcopy(NC, &xix[j1 * numCom], &xixC[j * NC]);
    }
}

void MixtureComp::CalRhoPNX_full01()
{
    // CaldXsdXp01() shoule be called before
    // Using rhsDer
    // Cal rhoP, rhoX
    // Water has been calculated
    // fill(rhox.begin(), rhox.end() - numCom, 0.0);
    // fill(rhoP.begin(), rhoP.end() - 1, 0.0);

    USI j1, bId;
    for (USI j = 0; j < NP; j++) {
        j1       = phaseLabel[j];
        bId      = j1 * numCom;
        rhoP[j1] = xiP[j1] * MW[j];
        for (USI m = 0; m < NC; m++) {
            rhox[bId + m] = xix[bId + m] * MW[j] + xi[j1] * MWC[m];
        }
    }
    if (NP == 1) {
        OCP_DBL tmp;
        j1  = phaseLabel[0];
        bId = j1 * numCom;
        for (USI m = 0; m < NC; m++) {
            rhoN[bId + m] = xiN[bId + m] * MW[0];
            tmp           = MWC[m] * Nh;
            for (USI i = 0; i < NC; i++) {
                tmp -= MWC[i] * Ni[i];
            }
            rhoN[bId + m] += tmp * xi[j1] / (Nh * Nh);
        }
    }
    if (NP > 1) {
        // correct rhoP
        // use rhsDer, see CaldXsdXp01(), attention that it only contains hydrocarbon
        OCP_DBL        tmp;
        const OCP_DBL* xijP = &rhsDer[NP];
        for (USI j = 0; j < NP; j++) {
            j1  = phaseLabel[j];
            tmp = 0;
            for (USI i = 0; i < NC; i++) {
                tmp += MWC[i] * xijP[i];
            }
            rhoP[j1] += xi[j1] * tmp;
            xijP += NC;
        }
        // Calculate rhoN
        const OCP_DBL* xijN = xijP;
        for (USI m = 0; m < NC; m++) {
            xijN += NP;
            for (USI j = 0; j < NP; j++) {
                j1            = phaseLabel[j];
                bId           = j1 * numCom;
                rhoN[bId + m] = xiN[bId + m] * MW[j];
                tmp           = 0;
                for (USI i = 0; i < NC; i++) {
                    tmp += MWC[i] * xijN[i];
                }
                rhoN[bId + m] += tmp * xi[j1];
                xijN += NC;
            }
        }
    }
}

void MixtureComp::CalRhoPNX_full()
{
    // CaldXsdXpAPI01() or CaldXsdXpAPI02() shoule be called before
    // Using dXsdXp
    // Cal rhoP, rhoX
    // Water has been calculated
    // fill(rhox.begin(), rhox.end() - numCom, 0.0);
    // fill(rhoP.begin(), rhoP.end() - 1, 0.0);

    USI     j1, bId;
    OCP_DBL tmp = 0;
    if (NP == 1) {
        j1       = phaseLabel[0];
        bId      = j1 * numCom;
        rhoP[j1] = xiP[j1] * MW[0];
        for (USI m = 0; m < NC; m++) {
            rhox[bId + m] = xix[bId + m] * MW[0] + xi[j1] * MWC[m];
            rhoN[bId + m] = xiN[bId + m] * MW[0];
            tmp           = MWC[m] * Nh;
            for (USI i = 0; i < NC; i++) {
                tmp -= MWC[i] * Ni[i];
            }
            rhoN[bId + m] += tmp * xi[j1] / (Nh * Nh);
        }
    } else {
        // Using dXsdXp
        fill(rhoP.begin(), rhoP.end() - 1, 0.0);
        fill(rhoN.begin(), rhoN.end() - numCom, 0.0);
        OCP_DBL ncol = numCom + 1;
        for (USI j = 0; j < NP; j++) {
            j1                   = phaseLabel[j];
            bId                  = j1 * numCom;
            const OCP_DBL* xijPN = &dXsdXp[(numPhase + bId) * ncol];

            for (USI i = 0; i < NC; i++) {
                rhoP[j1] += MWC[i] * xijPN[0];
                xijPN++;
                for (USI m = 0; m < NC; m++) {
                    rhoN[bId + m] += MWC[i] * xijPN[m];
                }
                xijPN += numCom;
            }
            rhoP[j1] *= xi[j1];
            rhoP[j1] += xiP[j1] * MW[j];
            for (USI m = 0; m < NC; m++) {
                rhoN[bId + m] *= xi[j1];
                rhoN[bId + m] += xiN[bId + m] * MW[j];
                rhox[bId + m] = xix[bId + m] * MW[j] + xi[j1] * MWC[m];
            }
        }
    }
}

void MixtureComp::CalMuPX_full01()
{
    fill(muP.begin(), muP.end() - 1, 0.0);
    fill(muN.begin(), muN.end() - numCom, 0.0);
    fill(mux.begin(), mux.end() - numCom, 0.0);

    CalMuPXLBC_full01();
}

void MixtureComp::CalMuPXLBC_full01()
{
    // CaldXsdXp01() before
    // use rhsDer
    OCP_DBL val1IJ, val2IJ;
    OCP_DBL der1IJ, der2IJ, der3J, der4J, der6J, der7J, der8J;
    OCP_DBL Tri, tmp;
    OCP_DBL xTj, xPj, xVj;
    OCP_DBL derxTj, derxPj, derMWj, derxVj;

    // Calculate dmuj / dxkj
    for (USI j = 0; j < NP; j++) {
        const vector<OCP_DBL>& xj     = x[j];
        const vector<OCP_DBL>& muAuxj = muAux[j];
        const USI              j1     = phaseLabel[j];
        const USI              bId    = numCom * j1;
        xTj = xPj = xVj = 0;
        for (USI i = 0; i < NC; i++) {
            xTj += xj[i] * Tc[i];
            xPj += xj[i] * Pc[i];
            xVj += xj[i] * Vcvis[i];
        }
        for (USI k = 0; k < NC; k++) {
            derxTj = Tc[k];
            derxPj = Pc[k];
            derMWj = MWC[k];
            derxVj = Vcvis[k];
            der3J  = 0;
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[j]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[j], -1.5) * derMWj;
                Tri    = T / Tc[i];
                if (Tri <= 1.5) {
                    tmp = 34 * 1E-5 * pow(Tri, 0.94);
                } else {
                    tmp = 17.78 * 1E-5 * pow(4.58 * Tri - 1.67, 0.625);
                }
                val2IJ = tmp / val1IJ;
                der2IJ = -tmp * der1IJ / (val1IJ * val1IJ);
                der3J += sqrtMWi[i] * (delta(i, k) * val2IJ + xj[i] * der2IJ);
            }
            der4J = sqrtMWi[k];
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            der7J = xix[bId + k] * xVj + xi[j1] * derxVj;
            if (muAuxj[3] <= 0.18 && OCP_FALSE) {
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                mux[bId + k] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }
    }

    if (NP == 1) {
        // NP = 1, then dxij / dP = 0, d MJ / dP = 0
        // Calculate dmuj / dP
        // der2IJ = der3J = der4J = der6J = 0;
        const vector<OCP_DBL>& xj     = x[0];
        const vector<OCP_DBL>& muAuxj = muAux[0];
        xTj = xPj = xVj = 0;
        for (USI i = 0; i < NC; i++) {
            xTj += xj[i] * Tc[i];
            xPj += xj[i] * Pc[i];
            xVj += xj[i] * Vcvis[i];
        }
        der7J = xVj * xiPC[0];
        if (muAuxj[3] <= 0.18 && OCP_FALSE) {
            muP[phaseLabel[0]] = (2.05 * 1E-4) * der7J / muAuxj[2];
        } else {
            der8J              = der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
            muP[phaseLabel[0]] = (4 * 1E-4) * pow(muAuxj[4], 3) * der8J / muAuxj[2];
        }
    } else {
        // NP > 1
        // Calculate dmuj / dP
        // use rhsDer (after CaldXsdXpAPI01)
        for (USI j = 0; j < NP; j++) {
            const OCP_DBL*         xijP   = &rhsDer[NP + j * NC];
            const vector<OCP_DBL>& xj     = x[j];
            const vector<OCP_DBL>& muAuxj = muAux[j];
            const USI              j1     = phaseLabel[j];
            xTj = xPj = xVj = 0;
            derxTj = derxPj = derMWj = 0;
            for (USI i = 0; i < NC; i++) {
                xTj += xj[i] * Tc[i];
                xPj += xj[i] * Pc[i];
                xVj += xj[i] * Vcvis[i];
                derxTj += xijP[i] * Tc[i];
                derxPj += xijP[i] * Pc[i];
                derMWj += xijP[i] * MWC[i];
            }
            der3J = der4J = der7J = 0;
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[j]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[j], -1.5) * derMWj;
                Tri    = T / Tc[i];
                if (Tri <= 1.5) {
                    tmp = 34 * 1E-5 * pow(Tri, 0.94);
                } else {
                    tmp = 17.78 * 1E-5 * pow(4.58 * Tri - 1.67, 0.625);
                }
                val2IJ = tmp / val1IJ;
                der2IJ = -tmp * der1IJ / (val1IJ * val1IJ);
                der3J += sqrtMWi[i] * (xijP[i] * val2IJ + xj[i] * der2IJ);
                der4J += sqrtMWi[i] * xijP[i];
                der7J += xijP[i] * Vcvis[i];
            }
            der7J *= xi[j1];
            der7J += xiP[j1] * xVj;
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            if (muAuxj[3] <= 0.18 && OCP_FALSE) {
                muP[j1] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                muP[j1] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }
    }
}

void MixtureComp::CalXiRhoMuPN_pfullx()
{
    // Using dXsdXp
    // CalXiPNX_partial(),CalRhoPX_partial(),CalMuPX_partial() have been called
    // s: S,x; p: P,N  Call CaldXsdXpAPI01 or CaldXsdXpAPI02 before

    const USI      ncol  = numCom + 1;
    const OCP_DBL* xijPN = &dXsdXp[numPhase * ncol];

    if (NP == 1) {
        const USI j1  = phaseLabel[0];
        const USI bId = j1 * numCom;
        xijPN += bId * ncol;

        for (USI i = 0; i < NC; i++) {
            xijPN++;
            // dN
            for (USI m = 0; m < NC; m++) {
                xiN[bId + m] += xix[bId + i] * xijPN[m];
                rhoN[bId + m] += rhox[bId + i] * xijPN[m];
                muN[bId + m] += mux[bId + i] * xijPN[m];
            }
            xijPN += numCom;
        }
    } else {
        for (USI j = 0; j < numPhase - 1; j++) {

            if (!phaseExist[j]) { // skip
                xijPN += numCom * ncol;
                continue;
            }

            const USI bId = j * numCom;
            for (USI i = 0; i < NC; i++) {
                // dP
                xiP[j] += xix[bId + i] * xijPN[0];
                rhoP[j] += rhox[bId + i] * xijPN[0];
                muP[j] += mux[bId + i] * xijPN[0];
                xijPN++;
                // xiN
                for (USI m = 0; m < NC; m++) {
                    xiN[bId + m] += xix[bId + i] * xijPN[m];
                    rhoN[bId + m] += rhox[bId + i] * xijPN[m];
                    muN[bId + m] += mux[bId + i] * xijPN[m];
                }
                xijPN += numCom;
            }
            // skip water component
            xijPN += ncol;
        }
    }
}

void MixtureComp::CalXiRhoMuPN_pfullxn(const OCP_BOOL& xflag)
{
    // Using dXsdXp
    // s: S,n; p: P,N  Call CaldXsdXpAPI03, CaldXsdXpAPI02p before
    // CalXiPn_partial(); CalRhoPn_partial(); CalMuPn_partial() have been called

    if (NP == 1 && !xflag) {
        const USI j1  = phaseLabel[0];
        const USI bId = j1 * numCom;

        // dN
        for (USI i = 0; i < NC; i++) {
            xiN[bId + i]  = xix[bId + i];
            rhoN[bId + i] = rhox[bId + i];
            muN[bId + i]  = mux[bId + i];
        }
    } else {
        const USI      ncol  = numCom + 1;
        const OCP_DBL* xijPN = &dXsdXp[(NP + 1) * ncol];

        for (USI j = 0; j < numPhase - 1; j++) {

            if (!phaseExist[j]) // skip
                continue;

            const USI bId = j * numCom;
            for (USI i = 0; i < NC; i++) {
                // dP
                xiP[j] += xix[bId + i] * xijPN[0];
                rhoP[j] += rhox[bId + i] * xijPN[0];
                muP[j] += mux[bId + i] * xijPN[0];
                xijPN++;
                // xiN
                for (USI m = 0; m < NC; m++) {
                    xiN[bId + m] += xix[bId + i] * xijPN[m];
                    rhoN[bId + m] += rhox[bId + i] * xijPN[m];
                    muN[bId + m] += mux[bId + i] * xijPN[m];
                }
                xijPN += numCom;
            }
            // don't skip water component
        }
    }
}

void MixtureComp::CalVfiVfp_full01()
{
    OCP_DBL CgTP = GAS_CONSTANT * T / P;

    if (NP == 1) {
        // NP = 1
        const OCP_DBL&         aj   = Aj[0];
        const OCP_DBL&         bj   = Bj[0];
        const OCP_DBL&         zj   = Zj[0];
        const vector<OCP_DBL>& xj   = x[0];
        vector<OCP_DBL>&       Znij = Zn[0];
        OCP_DBL                tmp;

        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }
            An[i] = 2 / nu[0] * (tmp - aj);
            Bn[i] = 1 / nu[0] * (Bi[i] - bj);
            Znij[i] =
                ((bj - zj) * An[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bn[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
            vfi[i] = CgTP * (zj + nu[0] * Znij[i]);
        }
        Zp[0] = ((bj - zj) * aj +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     bj) /
                P /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        vfP = CgTP * nu[0] * (Zp[0] - zj / P);
    } else {
        // NP > 1
        CalFugNAll();
        CalFugPAll();
        AssembleMatVfiVfp_full01();
        AssembleRhsVfiVfp_full01();
        LUSolve(NC + 1, NC * NP, &JmatDer[0], &rhsDer[0], &pivot[0]);
        // Calculate Vfi
        OCP_DBL* dnkjdNP = &rhsDer[0];
        for (USI i = 0; i < NC; i++) {
            vfi[i] = 0;
            for (USI j = 0; j < NP; j++) {
                for (USI k = 0; k < NC; k++) {
                    vfi[i] += dnkjdNP[j * NC + k] * (Zj[j] + nu[j] * Zn[j][k]);
                }
            }
            vfi[i] *= CgTP;
            dnkjdNP += NP * NC;
        }
        // Calculate Vfp
        vfP = 0;
        for (USI j = 0; j < NP; j++) {
            vfP += nu[j] / P * (Zp[j] * P - Zj[j]);
            for (USI k = 0; k < NC; k++) {
                vfP += dnkjdNP[j * NC + k] * (Zj[j] + nu[j] * Zn[j][k]);
            }
        }
        vfP *= CgTP;
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(vfi.size(), &vfi[0])) {
        OCP_ABORT("INF or NAN in vfi !");
    }
    if (!CheckNan(1, &vfP)) {
        OCP_ABORT("INF or NAN in vfP !");
    }
#endif // NANCHECK
}

void MixtureComp::AssembleMatVfiVfp_full01()
{
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* bId = &JmatDer[0];
    for (USI j = 0; j < NP - 1; j++) {
        // for jth phase
        OCP_DBL* fugNj = &fugN[j][0];
        for (USI i = 0; i < NC; i++) {
            // for ith components
            Dcopy(NC, bId, fugNj);
            bId += (NP - 1 - j) * NC;
            bId[i] = 1.0;
            bId += (1 + j) * NC;
            fugNj += NC;
        }
        bId += NC;
    }
    // NP - 1 phase
    bId            = &JmatDer[(NP - 1) * (NP * NC * NC)];
    OCP_DBL* fugNj = &fugN[NP - 1][0];
    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NP - 1; j++) {
            Dcopy(NC, bId, fugNj);
            bId += NC;
        }
        Dscalar((NP - 1) * NC, -1.0, bId - (NP - 1) * NC);
        fugNj += NC;
        bId[i] = 1.0;
        bId += NC;
    }
}

void MixtureComp::AssembleRhsVfiVfp_full01()
{
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];
    for (USI k = 0; k < NC; k++) {
        // d Nk
        rhstmp[NC * (NP - 1) + k] = 1;
        rhstmp += NP * NC;
    }
    // d P
    for (USI j = 0; j < NP; j++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[j * NC + i] = fugP[NP - 1][i] - fugP[j][i];
        }
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(rhsDer.size(), &rhsDer[0])) {
        OCP_ABORT("INF or NAN in rhsDer !");
    }
#endif // NANCHECK
}

void MixtureComp::CalVfiVfp_full02()
{
    // Attention!
    // NP = 1 or NP = 2

    OCP_DBL CgTP = GAS_CONSTANT * T / P;

    if (NP == 1) {
        // NP = 1
        const OCP_DBL&         aj   = Aj[0];
        const OCP_DBL&         bj   = Bj[0];
        const OCP_DBL&         zj   = Zj[0];
        const vector<OCP_DBL>& xj   = x[0];
        vector<OCP_DBL>&       Znij = Zn[0];
        OCP_DBL                tmp;

        for (USI i = 0; i < NC; i++) {
            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }
            An[i] = 2 / nu[0] * (tmp - aj);
            Bn[i] = 1 / nu[0] * (Bi[i] - bj);
            Znij[i] =
                ((bj - zj) * An[i] +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bn[i]) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
            vfi[i] = CgTP * (zj + nu[0] * Znij[i]) - Vshift[i];
        }
        Zp[0] = ((bj - zj) * aj +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     bj) /
                P /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
        vfP = CgTP * nu[0] * (Zp[0] - zj / P);
    } else {
        // NP = 2,  IF NP > 2  ->  WRONG!
        CalFugNAll();
        CalFugPAll();
        AssembleMatVfiVfp_full02();
        AssembleRhsVfiVfp_full02();
        LUSolve(NC + 1, NC, &JmatDer[0], &rhsDer[0], &pivot[0]);
        // now d nm0 / dP(dNk) has been available
        const OCP_DBL* dnkjdNP = &rhsDer[0];
        // Calculate Vfp
        const USI j0 = phaseLabel[0];
        const USI j1 = phaseLabel[1];
        vjp[j0]      = CgTP * nu[0] * (Zp[0] - Zj[0] / P);
        vjp[j1]      = CgTP * nu[1] * (Zp[1] - Zj[1] / P);
        for (USI k = 0; k < NC; k++) {
            vjp[j0] += (CgTP * (Zj[0] + nu[0] * Zn[0][k]) - Vshift[k]) * dnkjdNP[k];
            vjp[j1] -= (CgTP * (Zj[1] + nu[1] * Zn[1][k]) - Vshift[k]) * dnkjdNP[k];
        }
        vfP = vjp[j0] + vjp[j1];
        dnkjdNP += NC;

        // Calculate Vfi
        for (USI i = 0; i < NC; i++) {
            vfi[i]     = 0;
            vji[j0][i] = 0;
            vji[j1][i] = 0;
            for (USI k = 0; k < NC; k++) {
                vji[j0][i] +=
                    (CgTP * (Zj[0] + nu[0] * Zn[0][k]) - Vshift[k]) * dnkjdNP[k];
                vji[j1][i] += (CgTP * (Zj[1] + nu[1] * Zn[1][k]) - Vshift[k]) *
                              (delta(i, k) - dnkjdNP[k]);
            }
            vfi[i] = vji[j0][i] + vji[j1][i];
            dnkjdNP += NC;
        }
    }
}

void MixtureComp::AssembleMatVfiVfp_full02()
{
    // NP = 2
    fill(JmatDer.begin(), JmatDer.end(), 0.0);
    // Attention 1: JmatDer should be sorted by column
    // Attention 2: d ln fij / d nkj is symetric for each j;
    OCP_DBL* tmpMat = &JmatDer[0];
    for (USI i = 0; i < NC; i++) {
        for (USI m = 0; m < NC; m++) {
            tmpMat[m] = fugN[0][i * NC + m] + fugN[1][i * NC + m];
        }
        tmpMat += NC;
    }
}

void MixtureComp::AssembleRhsVfiVfp_full02()
{
    // NP = 2
    fill(rhsDer.begin(), rhsDer.end(), 0.0);
    OCP_DBL* rhstmp = &rhsDer[0];

    // dP
    for (USI i = 0; i < NC; i++) {
        rhstmp[i] = fugP[1][i] - fugP[0][i];
    }
    rhstmp += NC;

    // dNk
    for (USI k = 0; k < NC; k++) {
        for (USI i = 0; i < NC; i++) {
            rhstmp[i] = fugN[1][k * NC + i]; // d lnfij / d nkj = d lnfkj / d nij
        }
        rhstmp += NC;
    }
}

void MixtureComp::CaldXsdXp01()
{

    // NP > 1 in this function
    // only hydrocarbon was considered
    // if water exists, vf, vfP, and S should be updated
    // CalVfiVfp_full01() and CalXiPNX_full01() should be called before

    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    OCP_DBL* MTmp = &JmatTmp[0];

    // dF / dXs
    // dFn / dXs
    for (USI i = 0; i < NC; i++) {
        // dFn / dSm
        for (USI m = 0; m < NP; m++) {
            MTmp[m] = vf * xiC[m] * x[m][i];
        }
        MTmp += NP;
        // dFn / dxnm
        for (USI m = 0; m < NP; m++) {
            for (USI n = 0; n < NC; n++) {
                MTmp[n] =
                    vf * S[m] * (xixC[m * NC + n] * x[m][i] + delta(i, n) * xiC[m]);
            }
            MTmp += NC;
        }
    }
    // dFx / dXs
    for (USI j = 0; j < NP; j++) {
        // dFx / dSm
        MTmp += NP;
        // dFx / dxnm, m = j
        MTmp += j * NC;
        for (USI n = 0; n < NC; n++) {
            MTmp[n] = 1.0;
        }
        MTmp += (NP - j) * NC;
    }
    // dFf / dXs
    for (USI j = 0; j < NP - 1; j++) {
        const OCP_DBL* fugXJ  = &fugX[j][0];
        const OCP_DBL* fugXNP = &fugX[NP - 1][0];
        for (USI i = 0; i < NC; i++) {
            // dFf / dSm
            MTmp += NP;
            // dFf / dxnm, m = j or m = np-1
            // m = j
            MTmp += j * NC;
            for (USI n = 0; n < NC; n++) {
                MTmp[n] = fugXJ[n];
            }
            fugXJ += NC;
            // m = np-1
            MTmp += (NP - 1 - j) * NC;
            for (USI n = 0; n < NC; n++) {
                MTmp[n] = -fugXNP[n];
            }
            fugXNP += NC;
            MTmp += NC;
        }
    }

    // USI row01 = NP * (NC + 1);
    // USI col01 = row01;
    // cout << "FsXs" << endl;
    // for (USI i = 0; i < row01; i++) {
    //	for (USI j = 0; j < col01; j++) {
    //		cout << JmatTmp[i * row01 + j] << "   ";
    //	}
    //	cout << endl;
    //}

    // Transpose JmatTmp
    USI dim = NP * (NC + 1);
    for (USI i = 0; i < dim; i++) {
        for (USI j = 0; j < dim; j++) {
            JmatDer[i * dim + j] = JmatTmp[j * dim + i];
        }
    }

    // dF / dXp
    // d fij / dP have been calculated in CalVfiVfp_full01() before
    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    MTmp = &JmatTmp[0];
    OCP_DBL tmp01;
    OCP_DBL tmp02;
    // dFn / dXp
    for (USI i = 0; i < NC; i++) {
        // dFn / dP
        tmp01 = 0;
        tmp02 = 0;
        for (USI j = 0; j < NP; j++) {
            tmp01 += S[j] * x[j][i] * xiPC[j];
            tmp02 += S[j] * x[j][i] * xiC[j];
        }
        MTmp[0] = -(vf * tmp01 + vfP * tmp02); // in OCP
        // MTmp[0] = -(vf * tmp01); // in PennSim
        MTmp += 1;
        // dFn / dNk
        // in OCP
        for (USI k = 0; k < NC; k++) {
            tmp01 = 0;
            for (USI j = 0; j < NP; j++) {
                tmp01 += S[j] * x[j][i] * xiNC[j * NC + k];
            }
            MTmp[k] = -(vf * tmp01 + vfi[k] * tmp02 - delta(i, k));
        }
        // in PennSim
        // MTmp[i] = 1.0;
        MTmp += NC;
    }

    // dFx / dXp
    MTmp += NP * (NC + 1);

    // dFf / dXp
    const vector<OCP_DBL>& fugPNP = fugP[NP - 1];
    for (USI j = 0; j < NP - 1; j++) {
        const vector<OCP_DBL>& fugPj = fugP[j];
        for (USI i = 0; i < NC; i++) {
            // dFf / dP
            MTmp[0] = -(fugPj[i] - fugPNP[i]);
            MTmp += (NC + 1);
        }
    }

    // USI row02 = NP * (NC + 1);
    // USI col02 = NC + 1;
    // cout << "FsXp" << endl;
    // for (USI i = 0; i < row02; i++) {
    //	for (USI j = 0; j < col02; j++) {
    //		cout << JmatTmp[i * col02 + j] << "   ";
    //	}
    //	cout << endl;
    //}

    // Transpose JmatTmp
    USI nrhs = NC + 1;
    USI nrow = (NC + 1) * NP;
    for (USI i = 0; i < nrhs; i++) {
        for (USI j = 0; j < nrow; j++) {
            rhsDer[i * nrow + j] = JmatTmp[j * nrhs + i];
        }
    }
    LUSolve(NC + 1, (NC + 1) * NP, &JmatDer[0], &rhsDer[0], &pivot[0]);
    // now in rhsDer: dXs / dP, dXs / dNk
    // print solution
    // cout << "Solution" << endl;
    // for (USI i = 0; i < nrhs; i++) {
    //	for (USI j = 0; j < nrow; j++) {
    //		cout << rhsDer[i * nrow + j] << "    ";
    //	}
    //	cout << endl;
    //}
}

void MixtureComp::CaldXsdXpAPI01()
{
    // Calculate derivates for hydrocarbon phase and components
    // if water exists, vf, vfP, and S should be updated
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);

    OCP_DBL vw  = vj[numPhase - 1];
    OCP_DBL vwp = vjp[numPhase - 1];

    if (NP == 1) {
        USI     bId = (numCom + 1) * phaseLabel[0];
        OCP_DBL tmp = vf * vf;
        // only dSj / dP, dSj / dNk ---- hydrocarbon phase
        dXsdXp[bId] = (vfP * vw - vf * vwp) / tmp;
        for (USI k = 0; k < NC; k++) {
            dXsdXp[bId + k + 1] = (vfi[k] * vw) / tmp;
        }
        // dxij / dNk --- new
        bId = (numCom + 1) * (numPhase + numCom * phaseLabel[0]) + 1;
        for (USI i = 0; i < NC; i++) {
            for (USI k = 0; k < NC; k++) {
                dXsdXp[bId + k] = (delta(i, k) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += (numCom + 1);
        }
    } else {
        // NP > 1
        // Calculate Saturation
        for (USI j = 0; j < NP; j++) {
            S[j] = vC[j] / vf;
        }
        CalFugXAll();
        CaldXsdXp01();
        // copy from rhsDer to dXsdXp
        OCP_DBL*       DTmp;
        const OCP_DBL* STmp;
        USI            nrow  = NP * (NC + 1); // row num of rhsDer
        USI            ncol2 = numCom + 1;    // col num of dXsdXp
        // Saturation
        for (USI j = 0; j < NP; j++) {
            STmp = &rhsDer[j];
            DTmp = &dXsdXp[phaseLabel[j] * ncol2];
            // dS / dP
            DTmp[0] = STmp[0];
            // dS / dNk
            DTmp++;
            for (USI k = 0; k < NC; k++) {
                STmp += nrow;
                DTmp[k] = STmp[0];
            }
            // water is excluded
        }
        // xij ---- mole fraction of component i in phase j
        OCP_DBL* DbId = &dXsdXp[numPhase * ncol2];
        for (USI j = 0; j < NP; j++) {
            DTmp = DbId + phaseLabel[j] * numCom * ncol2;
            for (USI i = 0; i < NC; i++) {
                STmp = &rhsDer[NP + j * NC + i];
                // dxij / dP
                DTmp[0] = STmp[0];
                // dxij / dNk
                DTmp++;
                for (USI k = 0; k < NC; k++) {
                    STmp += nrow;
                    DTmp[k] = STmp[0];
                }
                DTmp += numCom;
                // water is excluded
            }
        }
    }
    // Correct Sj
    CalSaturation();

    USI Wpid = numPhase - 1;
    USI Wcid = numCom - 1;

    // Calculate water derivatives
    // only dSj / dNw, dSw / dP , dSw / dNk needs to be considered
    USI ncol = 1 + numCom;
    // dSj / dNw
    for (USI j = 0; j < NP; j++) {
        dXsdXp[(phaseLabel[j] + 1) * ncol - 1] = -S[phaseLabel[j]] * vfi[Wcid] / vf;
    }
    OCP_DBL* Dtmp = &dXsdXp[Wpid * ncol];
    // dSw / dP
    OCP_DBL vf2 = vf * vf;
    Dtmp[0]     = (vwp * vf - vw * vfP) / vf2;
    Dtmp++;
    // dSw / d Nk
    for (USI k = 0; k < NC; k++) {
        Dtmp[k] = -vw * vfi[k] / vf2;
    }
    // dSw / d Nw
    Dtmp[NC] = (vfi[Wcid] * vf - vw * vfi[Wcid]) / vf2;
}

void MixtureComp::CaldXsdXpAPI02()
{
    // Attention!
    // NP = 1 or NP = 2

    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0);
    USI ncol = numCom + 1;

    OCP_DBL vf2 = vf * vf;
    OCP_DBL vwp = vjp[numPhase - 1];
    OCP_DBL vw  = vj[numPhase - 1];

    if (NP == 1) {
        // when NP = 1, vfi = vhi, vfP = vhp + vwp, h = hydrocarbon
        USI      j0  = phaseLabel[0];
        OCP_DBL* bId = &dXsdXp[j0 * ncol];
        // dS / dP
        bId[0] = ((vfP - vwp) * vf - vfP * vj[j0]) / vf2;
        bId++;
        // dSh / dNm
        for (USI m = 0; m < NC; m++) {
            bId[m] = vfi[m] * vw / vf2;
        }
        bId += NC;
        // dSh / dNw
        bId[0] = -vfi[numCom - 1] * vj[j0] / vf2;
        // dSw / dP, dNm
        bId = &dXsdXp[(numPhase - 1) * ncol];
        for (USI m = 0; m < ncol; m++) {
            bId[m] = -dXsdXp[j0 * ncol + m];
        }
        // dxij / dNm
        bId = &dXsdXp[(numPhase + j0 * numCom) * ncol + 1];
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += ncol;
        }
    } else {
        // NP = 2
        // dSj / dP, dSj / dNm
        OCP_DBL* bId;
        for (USI j = 0; j < 2; j++) {
            const USI j1 = phaseLabel[j];
            bId          = &dXsdXp[j1 * ncol];
            // dSj / dP
            bId[0] = (vjp[j1] * vf - vfP * vj[j1]) / vf2;
            bId++;
            // dSj / dNm
            for (USI m = 0; m < numCom; m++) {
                bId[m] = (vji[j1][m] * vf - vfi[m] * vj[j1]) / vf2;
            }
        }
        // dSw / dP, dSw / dNm
        bId = &dXsdXp[(numPhase - 1) * ncol];
        // dSw / dP
        bId[0] = (vwp * vf - vfP * vw) / vf2;
        bId++;
        // dSw / dNm
        for (USI m = 0; m < numCom; m++) {
            bId[m] = (vji[(numPhase - 1)][m] * vf - vfi[m] * vw) / vf2;
        }
        bId += numCom;

        // dxij / dP, dxij / dNm
        const USI      j0       = phaseLabel[0];
        const USI      j1       = phaseLabel[1];
        const OCP_DBL* dnkjdNP  = &rhsDer[0];
        OCP_DBL*       bId1     = bId + j0 * numCom * ncol;
        OCP_DBL*       bId2     = bId + j1 * numCom * ncol;
        OCP_DBL        njDerSum = 0;
        // dxij / dP
        for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
        for (USI i = 0; i < NC; i++) {
            bId1[0] = (dnkjdNP[i] * nu[0] - njDerSum * Nh * n[0][i]) / (nu[0] * nu[0]);
            bId2[0] = (-dnkjdNP[i] * nu[1] + njDerSum * Nh * n[1][i]) / (nu[1] * nu[1]);
            bId1 += ncol;
            bId2 += ncol;
        }
        dnkjdNP += NC;
        // dxij / dNm
        for (USI m = 0; m < NC; m++) {
            njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
            bId1 = bId + j0 * numCom * ncol + m + 1;
            bId2 = bId + j1 * numCom * ncol + m + 1;
            for (USI i = 0; i < NC; i++) {
                bId1[0] =
                    (dnkjdNP[i] * nu[0] - njDerSum * Nh * n[0][i]) / (nu[0] * nu[0]);
                bId2[0] = ((delta(i, m) - dnkjdNP[i]) * nu[1] -
                           (1 - njDerSum) * Nh * n[1][i]) /
                          (nu[1] * nu[1]);
                bId1 += ncol;
                bId2 += ncol;
            }
            dnkjdNP += NC;
        }
    }
}

void MixtureComp::CaldXsdXpAPI02p()
{
    // Attention!
    // NP = 1 or NP = 2

    // water not in oil or gas; hydroncarbon not in water
    // inexist phase will bot be stored
    // only hydrocarbon phases, hydrocarbon components, water phase will be stored

    // dS / dP
    // S = Sj, xij
    // P = P, Ni
    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0);
    USI ncol = numCom + 1;

    OCP_DBL vf2 = vf * vf;
    OCP_DBL vwp = vjp[numPhase - 1];
    OCP_DBL vw  = vj[numPhase - 1];

    if (NP == 1) {
        // when NP = 1, vfi = vhi, vfP = vhp + vwp, h = hydrocarbon
        const USI j1  = phaseLabel[0];
        OCP_DBL*  bId = &dXsdXp[0];
        // dS / dP
        bId[0] = ((vfP - vwp) * vf - vfP * vj[j1]) / vf2;
        bId++;
        // dSh / dNm
        for (USI m = 0; m < NC; m++) {
            bId[m] = vfi[m] * vw / vf2;
        }
        bId += NC;
        // dSh / dNw
        bId[0] = -vfi[numCom - 1] * vj[j1] / vf2;
        bId++;
        // dSw / dP, dNm
        for (USI m = 0; m < ncol; m++) {
            bId[m] = -dXsdXp[m];
        }
        bId += ncol + 1;
        // dxij / dNm
        for (USI i = 0; i < NC; i++) {
            for (USI m = 0; m < NC; m++) {
                bId[m] = (delta(i, m) * Nh - Ni[i]) / (Nh * Nh);
            }
            bId += ncol;
        }
    } else {
        // NP = 2
        // dSj / dP, dSj / dNm
        OCP_DBL* bId;
        for (USI j = 0; j < 2; j++) {
            const USI j1 = phaseLabel[j];
            bId          = &dXsdXp[j1 * ncol];
            // dSj / dP
            bId[0] = (vjp[j1] * vf - vfP * vj[j1]) / vf2;
            bId++;
            // dSj / dNm
            for (USI m = 0; m < numCom; m++) {
                bId[m] = (vji[j1][m] * vf - vfi[m] * vj[j1]) / vf2;
            }
        }
        // dSw / dP, dSw / dNm
        bId = &dXsdXp[(numPhase - 1) * ncol];
        // dSw / dP
        bId[0] = (vwp * vf - vfP * vw) / vf2;
        bId++;
        // dSw / dNm
        for (USI m = 0; m < numCom; m++) {
            bId[m] = (vji[(numPhase - 1)][m] * vf - vfi[m] * vw) / vf2;
        }
        bId += numCom;

        // dxij / dP, dxij / dNm
        const USI      j0       = phaseLabel[0];
        const USI      j1       = phaseLabel[1];
        const OCP_DBL* dnkjdNP  = &rhsDer[0];
        OCP_DBL*       bId1     = bId + j0 * NC * ncol;
        OCP_DBL*       bId2     = bId + j1 * NC * ncol;
        OCP_DBL        njDerSum = 0;
        // dxij / dP
        for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
        for (USI i = 0; i < NC; i++) {
            bId1[0] = (dnkjdNP[i] * nu[0] - njDerSum * Nh * n[0][i]) / (nu[0] * nu[0]);
            bId2[0] = (-dnkjdNP[i] * nu[1] + njDerSum * Nh * n[1][i]) / (nu[1] * nu[1]);
            bId1 += ncol;
            bId2 += ncol;
        }
        dnkjdNP += NC;
        // dxij / dNm
        for (USI m = 0; m < NC; m++) {
            njDerSum = 0;
            for (USI i = 0; i < NC; i++) njDerSum += dnkjdNP[i];
            bId1 = bId + j0 * NC * ncol + m + 1;
            bId2 = bId + j1 * NC * ncol + m + 1;
            for (USI i = 0; i < NC; i++) {
                bId1[0] =
                    (dnkjdNP[i] * nu[0] - njDerSum * Nh * n[0][i]) / (nu[0] * nu[0]);
                bId2[0] = ((delta(i, m) - dnkjdNP[i]) * nu[1] -
                           (1 - njDerSum) * Nh * n[1][i]) /
                          (nu[1] * nu[1]);
                bId1 += ncol;
                bId2 += ncol;
            }
            dnkjdNP += NC;
        }
    }
}

void MixtureComp::CaldXsdXpAPI03()
{
    // p : P, Ni
    // s : S, nij
    // water not in oil or gas; hydroncarbon not in water
    // inexist phase will bot be stored
    // only hydrocarbon phases, hydrocarbon components, water phase will be stored

    // water is included
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);
    fill(res.begin(), res.end(), 0.0);
    resPc              = 0;
    const USI     ncol = numCom + 1;
    const OCP_DBL vf2  = vf * vf;
    const OCP_DBL vwp  = vjp[numPhase - 1];
    const OCP_DBL vw   = vj[numPhase - 1];

    if (NP == 1) {
        // when NP = 1, vfi = vhi, vfP = vhp + vwp, h = hydrocarbon
        // and vfi = vji
        const USI j0  = phaseLabel[0];
        OCP_DBL*  bId = &dXsdXp[0];
        // dS / dP
        bId[0] = ((vfP - vwp) * vf - vfP * vj[j0]) / vf2;
        bId++;
        // dSh / dNm
        for (USI m = 0; m < NC; m++) {
            bId[m] = vji[j0][m] * vw / vf2;
        }
        bId += NC;
        // dSh / dNw
        bId[0] = -vfi[numCom - 1] * vj[j0] / vf2;
        bId++;
        // dSw / dP, dNm
        for (USI m = 0; m < ncol; m++) {
            bId[m] = -dXsdXp[m];
        }

        // dnij / dNm
        bId = &dXsdXp[2 * ncol + 1];
        for (USI i = 0; i < NC; i++) {
            bId[i] = 1;
            bId += ncol;
        }
        // res = 0

    } else {
        CalFugNAll(OCP_FALSE);
        CalFugPAll(OCP_FALSE);
        CaldXsdXp03();

        // Assemble dXsdXp
        // // inexist phase will bot be stored
        // no d xiw / dP, d xiw / dNk
        // no d xwj / dP, d xwj / dNk

        vector<USI> sIndex(numPhase, 0); // store index
        for (USI j = 0; j < numPhase - 1; j++) {
            if (phaseExist[j]) {
                for (USI j1 = j + 1; j1 < numPhase; j1++) {
                    sIndex[j1]++;
                }
            }
        }

        const OCP_DBL* STmp = &rhsDer[0];
        const USI      nrhs = numCom + 1;
        const USI      wNP  = NP + 1;

        USI bId;
        for (USI i = 0; i < nrhs; i++) {
            // Sh
            for (USI j = 0; j < NP; j++) {
                dXsdXp[sIndex[phaseLabel[j]] * nrhs + i] = STmp[j];
            }
            STmp += NP;
            // Sw
            dXsdXp[NP * nrhs + i] = STmp[0];
            STmp++;
            // nij
            for (USI j = 0; j < NP; j++) {
                bId = wNP + sIndex[phaseLabel[j]] * NC;
                for (USI k = 0; k < NC; k++) {
                    dXsdXp[(bId + k) * nrhs + i] = STmp[k];
                }
                STmp += NC;
            }
        }

        // res
        // Sh
        for (USI j = 0; j < NP; j++) {
            res[sIndex[phaseLabel[j]]] = STmp[j];
        }
        STmp += NP;
        // Sw
        res[NP] = STmp[0];
        STmp++;
        // nij
        for (USI j = 0; j < NP; j++) {
            bId = wNP + sIndex[phaseLabel[j]] * NC;
            for (USI k = 0; k < NC; k++) {
                res[bId + k] = STmp[k];
            }
            STmp += NC;
        }

        OCP_DBL* myres = &res[NP + 1];
        // precalculate a value
        for (USI j = 0; j < NP; j++) {
            const USI j1 = phaseLabel[j];
            for (USI i = 0; i < NC; i++) {
                resPc += vji[j1][i] * myres[i];
            }
            myres += NC;
        }

        // resPc = 0;
        // fill(res.begin(), res.end(), 0.0);

        // cout << "res1" << endl;
        //  const USI nrow = wNP + NP * NC;
        // for (USI i = 0; i < nrow; i++) {
        //     for (USI j = nrhs; j < nrhs + 1; j++) {
        //         cout << setw(12) << rhsDer[j * nrow + i] << "   ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;
        // cout << "resPc = " << resPc << endl;
    }

#ifdef OCP_NANCHECK
    if (!CheckNan(dXsdXp.size(), &dXsdXp[0])) {
        OCP_ABORT("INF or INF in bmat !");
    }
#endif
}

void MixtureComp::CaldXsdXp03()
{
    // p : P, Ni
    // s : S, nij

    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    OCP_DBL*      MTmp = &JmatTmp[0];
    const USI     dim  = NP + 1 + NP * NC;
    const OCP_DBL vf2  = vf * vf;
    // -dF / ds

    // -dFf / ds
    for (USI j = 0; j < NP - 1; j++) {
        const OCP_DBL* fugNJ  = &fugN[j][0];
        const OCP_DBL* fugNNP = &fugN[NP - 1][0];
        for (USI i = 0; i < NC; i++) {
            // -dFf / dS = 0
            MTmp += (NP + 1);
            // -dFf / dnij
            MTmp += j * NC;
            for (USI k = 0; k < NC; k++) {
                MTmp[k] = -fugNJ[k];
            }
            fugNJ += NC;
            // np-1 phase
            MTmp += (NP - 1 - j) * NC;
            for (USI k = 0; k < NC; k++) {
                MTmp[k] = fugNNP[k];
            }
            fugNNP += NC;
            MTmp += NC;
        }
    }
    // -dFn / ds
    for (USI i = 0; i < NC; i++) {
        // -dFn / dS = 0
        MTmp += (NP + 1);
        // -dFn / dnij
        for (USI j = 0; j < NP; j++) {
            MTmp[i] = 1;
            MTmp += NC;
        }
    }
    // -dFs / ds  ---- Hydroncarbon
    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];
        // -dFs / dS
        MTmp[j] = -1;
        MTmp += (NP + 1);
        // -dFs / dnkm
        for (USI m = 0; m < NP; m++) {
            const USI m1 = phaseLabel[m];
            if (m1 == j1) {
                for (USI k = 0; k < NC; k++) {
                    MTmp[k] = vji[j1][k] * (vf - vj[j1]) / vf2;
                }
            } else {
                for (USI k = 0; k < NC; k++) {
                    MTmp[k] = -vji[m1][k] * vj[j1] / vf2;
                }
            }
            MTmp += NC;
        }
    }
    // -dFsw / ds
    // -dFsw / dSw
    MTmp[NP] = -1;
    MTmp += (NP + 1);
    // -dFsw / dnkm
    for (USI m = 0; m < NP; m++) {
        const USI m1 = phaseLabel[m];
        for (USI k = 0; k < NC; k++) {
            MTmp[k] = -vji[m1][k] * vj[numPhase - 1] / vf2;
        }
        MTmp += NC;
    }

    // cout << "dFdS" << endl;
    // cout << scientific << setprecision(1);
    // for (USI i = 0; i < dim; i++) {
    //     for (USI j = 0; j < dim; j++)
    //         cout << JmatTmp[i * dim + j] << " ";
    //     cout << endl;
    //     if (i == dim - 4)
    //         cout << endl;
    // }
    // cout << phaseLabel[0] << "   " << phaseLabel[1] << endl;

    // Transpose JmatTmp
    for (USI i = 0; i < dim; i++) {
        for (USI j = 0; j < dim; j++) {
            JmatDer[i * dim + j] = JmatTmp[j * dim + i];
        }
    }

    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    MTmp = &JmatTmp[0];

    // dF / dp
    const USI ncol2 = numCom + 1;
    // dFf / dp
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            MTmp[0] = fugP[j][i] - fugP[NP - 1][i];
            MTmp += ncol2;
        }
    }
    // dFn / dp
    for (USI i = 0; i < NC; i++) {
        MTmp++;
        MTmp[i] = 1;
        MTmp += numCom;
    }
    // dFs / dp   ---- Hydroncarbon
    for (USI j = 0; j < NP; j++) {
        const USI j1 = phaseLabel[j];
        // dFs / dP
        MTmp[0] = (vfP * vj[j1] - vjp[j1] * vf) / vf2;
        // dFs / dNh
        MTmp += numCom;
        // dFs / dNw
        MTmp[0] = vfi[numCom - 1] * vj[j1] / vf2;
        MTmp++;
    }
    // dFsw / dp
    // dFsw / dP
    MTmp[0] = (vfP * vj[numPhase - 1] - vjp[numPhase - 1] * vf) / vf2;
    // dFsw / dNh
    MTmp += numCom;
    // dFsw / dNw
    MTmp[0] = vfi[numCom - 1] * (vj[numPhase - 1] - vf) / vf2;
    MTmp++;

    // cout << "dFdp" << endl;
    // cout << scientific << setprecision(3);
    // for (USI i = 0; i < dim; i++) {
    //     for (USI j = 0; j < numCom+1; j++)
    //         cout << JmatTmp[i * (numCom + 1) + j] << "   ";
    //     cout << endl;
    // }

    // Transpose JmatTmp
    const USI nrhs = numCom + 1;
    const USI nrow = NP * NC + NP + 1;
    for (USI i = 0; i < nrhs; i++) {
        for (USI j = 0; j < nrow; j++) {
            rhsDer[i * nrow + j] = JmatTmp[j * nrhs + i];
        }
    }

    // Calculate res
    // resF
    OCP_DBL* myRes = &rhsDer[nrhs * nrow];
    for (USI j = 0; j < NP - 1; j++) {
        const OCP_DBL* fugJ  = &fug[j][0];
        const OCP_DBL* fugNP = &fug[NP - 1][0];
        for (USI i = 0; i < NC; i++) {
            myRes[i] = fugJ[i] - fugNP[i];
        }
        myRes += NC;
    }
    // resN
    for (USI i = 0; i < NC; i++) {
        myRes[i] = Ni[i];
        for (USI j = 0; j < NP; j++) {
            myRes[i] -= x[j][i] * nu[j];
        }
    }
    myRes += NC;
    // resS
    for (USI j = 0; j < NP; j++) {
        myRes[j] = S[phaseLabel[j]] - vj[phaseLabel[j]] / vf;
    }
    myRes[NP] = S[numPhase - 1] - vj[numPhase - 1] / vf;

    // cout << scientific << setprecision(4);
    // cout << "res0" << endl;
    // for (USI i = 0; i < nrow; i++) {
    //     for (USI j = nrhs; j < nrhs + 1; j++) {
    //         cout << setw(12) << rhsDer[j * nrow + i] << "   ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    LUSolve(numCom + 1 + 1, nrow, &JmatDer[0], &rhsDer[0], &pivot[0]);
}

void MixtureComp::CalVfiVfp_full03()
{
    // Call CaldXsdXpAPI03() before
    // use dXsdXp: s = Sj,nij; p = P,Ni
    USI j1;
    if (NP == 1) {
        // vfP = vfP; vfi = vji
        j1 = phaseLabel[0];
        for (USI i = 0; i < NC; i++) {
            vfi[i] = vji[j1][i];
        }
    } else {
        const USI      ncol  = numCom + 1;
        const OCP_DBL* nijPN = &dXsdXp[(NP + 1) * ncol];
        for (USI j = 0; j < numPhase - 1; j++) {

            if (!phaseExist[j]) // skip
                continue;

            for (USI i = 0; i < NC; i++) {
                vfP += vji[j][i] * nijPN[0];
                nijPN++;
                for (USI m = 0; m < NC; m++) {
                    vfi[m] += vji[j][i] * nijPN[m];
                }
                nijPN += numCom;
            }
            // don't skip water, no water in dXsdXp here.
        }
    }
}

void MixtureComp::CalKeyDerx()
{
    // Call CaldXsdXpAPI01 or CaldXsdXpAPI02 before
    // Calculate d (xij*xi/mu) / dP or dNk
    // no water component
    // use dXsdXp
    // s = S, xij, p = P, Nk

    fill(keyDer.begin(), keyDer.end(), 0.0);
    const USI ncol = numCom + 1;
    OCP_DBL   tmp;

    const OCP_DBL* sTmp   = &dXsdXp[numPhase * ncol];
    const OCP_DBL* xiNtmp = &xiN[0];
    const OCP_DBL* muNtmp = &muN[0];
    const OCP_DBL* xijtmp = &xij[0];
    OCP_DBL*       dTmp   = &keyDer[0];

    // hydrocarbon phase
    for (USI j = 0; j < numPhase - 1; j++) {
        if (!phaseExist[j]) {
            sTmp += numCom * ncol;
            xiNtmp += numCom;
            muNtmp += numCom;
            xijtmp += numCom;
            continue;
        }
        tmp = xi[j] / (mu[j] * mu[j]);
        for (USI i = 0; i < NC; i++) {
            // dP
            dTmp[0] = (sTmp[0] * xi[j] + xijtmp[i] * xiP[j]) / mu[j] -
                      xijtmp[i] * muP[j] * tmp;
            dTmp++;
            sTmp++;
            // dNk
            for (USI k = 0; k < NC; k++) {
                dTmp[k] = (sTmp[k] * xi[j] + xijtmp[i] * xiNtmp[k]) / mu[j] -
                          xijtmp[i] * muNtmp[k] * tmp;
            }
            dTmp += numCom;
            sTmp += numCom;
        }
        xiNtmp += numCom;
        muNtmp += numCom;
        xijtmp += numCom;
        // skip water components
        sTmp += ncol;
    }

    // water phase
    // dP
    const USI wpid = numPhase - 1;
    keyDer[NP * NC * ncol] =
        (xiP[wpid] * mu[wpid] - muP[wpid] * xi[wpid]) / (mu[wpid] * mu[wpid]);
}

void MixtureComp::CalKeyDern()
{
    // Call CaldXsdXpAPI03 before
    // Calculate d (xij*xi/mu) / dP or dNk
    // no water component
    // use dXsdXp
    // s = S, nij, p = P, Nk
    // xij = nij / nj

    vector<OCP_DBL> njPN(NC + 1, 0);

    fill(keyDer.begin(), keyDer.end(), 0.0);
    const USI      ncol   = numCom + 1;
    const OCP_DBL* sTmp   = &dXsdXp[(NP + 1) * ncol];
    const OCP_DBL* xiNtmp = &xiN[0];
    const OCP_DBL* muNtmp = &muN[0];
    const OCP_DBL* xijtmp = &xij[0];
    OCP_DBL*       dTmp   = &keyDer[0];
    OCP_DBL        tmp01, tmp02;

    if (NP == 1) {
        const USI j = phaseLabel[0];
        xiNtmp += j * numCom;
        muNtmp += j * numCom;
        xijtmp += j * numCom;
        const OCP_DBL tmp = xi[j] / (mu[j] * mu[j]);

        for (USI i = 0; i < NC; i++) {
            // dP
            dTmp[0] = (xijtmp[i] * xiP[j]) / mu[j] - xijtmp[i] * muP[j] * tmp;
            dTmp++;
            // dNk
            for (USI k = 0; k < NC; k++) {
                dTmp[k] = ((delta(i, k) - xijtmp[i]) * xi[j] / nj[j] +
                           xijtmp[i] * xiNtmp[k]) /
                              mu[j] -
                          xijtmp[i] * muNtmp[k] * tmp;
            }
            dTmp += numCom;
        }
    } else {
        // hydrocarbon phase
        for (USI j = 0; j < numPhase - 1; j++) {
            if (!phaseExist[j]) {
                xiNtmp += numCom;
                muNtmp += numCom;
                xijtmp += numCom;
                continue;
            }

            fill(njPN.begin(), njPN.end(), 0.0);
            for (USI i = 0; i < NC; i++) {
                for (USI k = 0; k < NC + 1; k++) {
                    njPN[k] += sTmp[k];
                }
                sTmp += ncol;
            }
            sTmp -= NC * ncol;

            tmp01 = xi[j] / (mu[j] * mu[j]);
            for (USI i = 0; i < NC; i++) {
                tmp02 = (sTmp[0] - njPN[0] * xijtmp[i]) / nj[j];
                // dP
                dTmp[0] = (tmp02 * xi[j] + xijtmp[i] * xiP[j]) / mu[j] -
                          xijtmp[i] * muP[j] * tmp01;
                dTmp++;
                sTmp++;
                // dNk
                for (USI k = 0; k < NC; k++) {
                    tmp02   = (sTmp[k] - njPN[k + 1] * xijtmp[i]) / nj[j];
                    dTmp[k] = (tmp02 * xi[j] + xijtmp[i] * xiNtmp[k]) / mu[j] -
                              xijtmp[i] * muNtmp[k] * tmp01;
                }
                dTmp += numCom;
                sTmp += numCom;
            }
            xiNtmp += numCom;
            muNtmp += numCom;
            xijtmp += numCom;
        }
    }

    // water phase
    // dP
    const USI wpid = numPhase - 1;
    keyDer[NP * NC * ncol] =
        (xiP[wpid] * mu[wpid] - muP[wpid] * xi[wpid]) / (mu[wpid] * mu[wpid]);
}

USI MixtureComp::CubicRoot(const OCP_DBL&  a,
                           const OCP_DBL&  b,
                           const OCP_DBL&  c,
                           const OCP_BOOL& NTflag) const
{

    OCP_DBL Q = (a * a - 3 * b) / 9;
    OCP_DBL R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

    OCP_DBL Q3 = Q * Q * Q;
    OCP_DBL M  = R * R - Q3;

    if (M <= 0) {
        // 3 real roots
        OCP_DBL theta = acos(R / sqrt(Q3));
        Ztmp[0]       = -2 * sqrt(Q) * cos(theta / 3) - a / 3;
        Ztmp[1]       = -2 * sqrt(Q) * cos((theta + 2 * PI) / 3) - a / 3;
        Ztmp[2]       = -2 * sqrt(Q) * cos((theta - 2 * PI) / 3) - a / 3;

        if (NTflag) {
            NTcubicroot(Ztmp[0], a, b, c);
            NTcubicroot(Ztmp[1], a, b, c);
            NTcubicroot(Ztmp[2], a, b, c);
        }

        sort(Ztmp.begin(), Ztmp.end());

        // vector<OCP_DBL> e(3, 0);
        // for (USI i = 0; i < 3; i++) {
        //	e[i] = Ztmp[i] * (Ztmp[i] * (Ztmp[i] + a) + b) + c;
        //}
        // for (USI i = 0; i < 3; i++) {
        //	cout << scientific << e[i] << "\t";
        //}

        return 3;
    } else {
        OCP_DBL tmp1 = -R + sqrt(M);
        OCP_DBL tmp2 = R + sqrt(M);
        OCP_DBL S    = signD(tmp1) * pow(fabs(tmp1), 1.0 / 3);
        OCP_DBL T    = -signD(tmp2) * pow(fabs(tmp2), 1.0 / 3);
        Ztmp[0]      = S + T - a / 3;

        if (NTflag) {
            NTcubicroot(Ztmp[0], a, b, c);
        }

        // vector<OCP_DBL> e(1, 0);
        // for (USI i = 0; i < 1; i++) {
        //	e[i] = Ztmp[i] * (Ztmp[i] * (Ztmp[i] + a) + b) + c;
        //}
        // for (USI i = 0; i < 1; i++) {
        //	cout << scientific << e[i] << "\t";
        //}

        return 1;
    }
}

/// Return the sign of double di
OCP_DBL signD(const OCP_DBL& d)
{
    if (d > 0) {
        return 1.0;
    } else if (d < 0) {
        return -1.0;
    } else {
        return 0.0;
    }
}

OCP_DBL delta(const USI& i, const USI& j)
{
    if (i == j) {
        return 1.0;
    }
    return 0.0;
}

void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c)
{
    OCP_DBL e = root * (root * (root + a) + b) + c;
    OCP_DBL df;
    OCP_DBL iter    = 0;
    OCP_DBL optroot = root;
    OCP_DBL opte    = fabs(e);

    while (fabs(e) > 1E-8) {

        df   = root * (3 * root + 2 * a) + b;
        root = root - e / df;
        iter++;
        if (iter > 10) {
            // std::cout << "WARNING: INEXACT ROOT FOR CUBIC EQUATIONS" << std::endl;
            break;
        }
        e = root * (root * (root + a) + b) + c;
        if (fabs(e) <= opte) {
            opte    = fabs(e);
            optroot = root;
        }
    }
    root = optroot;
}

/////////////////////////////////////////////////////////////////////
// For Output
/////////////////////////////////////////////////////////////////////

void MixtureComp::OutMixtureIters() const
{
    cout << "SSMSTA:     " << setw(12) << itersSSMSTA << setw(15)
         << itersSSMSTA * 1.0 / countsSSMSTA << endl;
    cout << "NRSTA:      " << setw(12) << itersNRSTA << setw(15)
         << itersNRSTA * 1.0 / countsNRSTA << endl;
    cout << "SSMSP:      " << setw(12) << itersSSMSP << setw(15)
         << itersSSMSP * 1.0 / countsSSMSP << endl;
    cout << "NRSP:       " << setw(12) << itersNRSP << setw(15)
         << itersNRSP * 1.0 / countsNRSP << endl;
    cout << "NRRR:       " << setw(12) << itersRR << setw(15)
         << itersRR * 1.0 / countsRR << endl;
}

/////////////////////////////////////////////////////////////////////
// Optional Features
/////////////////////////////////////////////////////////////////////

void MixtureComp::SetupOptionalFeatures(OptionalFeatures& optFeatures,
                                        const OCP_USI&    numBulk)
{
    skipSta = &optFeatures.skipStaAnaly;
    if (skipSta->IfUseSkip()) {
        skipSta->Setup(numBulk, numPhase - 1, numCom - 1);
        AllocateSkip();
    }
    misTerm = &optFeatures.miscible;
    if (misTerm->IfUseMiscible() == ifUseMiscible == OCP_TRUE) {
        misTerm->Setup(numBulk);
    } else if (misTerm->IfUseMiscible() != ifUseMiscible) {
        OCP_ABORT(
            "Keywords MISCIBLE, PARACHOR, MISCSTR aren't been given simultaneously!");
    }
}

/////////////////////////////////////////////////////////////////////
// Accelerate PVT
/////////////////////////////////////////////////////////////////////

void MixtureComp::AllocateSkip()
{
    phiN.resize(NC * NC);
    skipMatSTA.resize(NC * NC);
    eigenSkip.resize(NC);
    eigenWork.resize(2 * NC + 1);
}

void MixtureComp::CalPhiNSkip()
{
    OCP_DBL C, E, G;
    OCP_DBL Cnk, Dnk, Enk, Gnk;
    OCP_DBL tmp, aik;

    // 0 th phase
    const OCP_DBL&         aj  = Aj[0];
    const OCP_DBL&         bj  = Bj[0];
    const OCP_DBL&         zj  = Zj[0];
    const vector<OCP_DBL>& xj  = x[0];
    vector<OCP_DBL>&       Znj = Zn[0];

    for (USI i = 0; i < NC; i++) {
        tmp = 0;
        for (USI m = 0; m < NC; m++) {
            tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
        }
        An[i]  = 2 / nu[0] * (tmp - aj);
        Bn[i]  = 1 / nu[0] * (Bi[i] - bj);
        Znj[i] = ((bj - zj) * An[i] +
                  ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                   ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                   (delta1 + delta2 - 1) * zj * zj) *
                      Bn[i]) /
                 (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                  (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
    }

    G = (zj + delta1 * bj) / (zj + delta2 * bj);

    for (USI i = 0; i < NC; i++) {
        // i th fugacity
        C = 1 / (zj - bj);
        // D = Bi[i] / bj * (zj - 1);
        tmp = 0;
        for (USI k = 0; k < NC; k++) {
            tmp += (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
        }
        E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

        for (USI k = 0; k <= i; k++) {
            // k th components

            aik = (1 - BIC[i * NC + k]) * sqrt(Ai[i] * Ai[k]);

            Cnk = (Bn[k] - Znj[k]) / ((zj - bj) * (zj - bj));
            Dnk = Bi[i] / bj * (Znj[k] - (Bi[k] - bj) * (zj - 1) / (nu[0] * bj));
            Gnk = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) *
                  (Bn[k] * zj - Znj[k] * bj);
            /*Enk = 1 / ((delta1 - delta2) * bj * bj) * (An[k] * bj - Bn[k] * aj) *
               (Bi[i] / bj - 2 * tmp / aj)
                + aj / ((delta1 - delta2) * bj) * (-Bi[i] / (bj * bj) * Bn[k] - 2 /
               (aj * aj) * (aj * (aik - tmp) / nu[j] - An[k] * tmp));*/
            Enk = -1 / (delta1 - delta2) / (bj * bj) *
                  (2 * (bj * aik / nu[0] + Bn[k] * (Bi[i] * aj / bj - tmp)) -
                   An[k] * Bi[i] - aj * Bi[i] / nu[0]);
            Enk -= E / nu[0];
            phiN[i * NC + k] = 1 / C * Cnk + Dnk + Enk * log(G) + E / G * Gnk;
            phiN[k * NC + i] = phiN[i * NC + k];
        }
    }
}

void MixtureComp::AssembleSkipMatSTA()
{
    // Sysmetric Matrix
    // stored by colum
    vector<OCP_DBL>& xj = x[0];

    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j <= i; j++) {
            skipMatSTA[i * NC + j] =
                delta(i, j) + nu[0] * sqrt(xj[i] * xj[j]) * phiN[i * NC + j];
            skipMatSTA[j * NC + i] = skipMatSTA[i * NC + j];
        }
    }
}

void MixtureComp::CalSkipForNextStep()
{
    if (flagSkip && ftype == 0) {
        // 1. Np == 1 is base for Skipping
        // 2. If flagSkip == true, then next stablity analysis is possible to be
        // skipped, it depends on if conditions are met
        // 3. If ftype == 0, then the range should be calculated, which also means last
        // skip is unsatisfied
        CalPhiNSkip();
        AssembleSkipMatSTA();
#ifdef DEBUG
        if (!CheckNan(skipMatSTA.size(), &skipMatSTA[0])) {
            OCP_WARNING("Nan in skipMatSTA!");
        }
#endif // DEBUG

        CalEigenSY(NC, &skipMatSTA[0], &eigenSkip[0], &eigenWork[0], 2 * NC + 1);
        skipSta->AssignValue(bulkId, eigenSkip[0], P, T, zi);
    }
    skipSta->SetFlagSkip(bulkId, flagSkip);
}

/////////////////////////////////////////////////////////////////////
// Miscible
/////////////////////////////////////////////////////////////////////

void MixtureComp::InputMiscibleParam(const ComponentParam& param, const USI& tarId)
{
    ifUseMiscible = param.miscible;
    if (param.Parachor.activity) parachor = param.Parachor.data[tarId];
    if (ifUseMiscible && parachor.empty()) {
        OCP_ABORT("PARACHOR has not been Input!");
    }
}

void MixtureComp::CalSurfaceTension()
{
    // be careful!
    // phase molar densities should be converted into gm-M/cc here
    if (ifUseMiscible) {
        if (NP == 1)
            surTen = 100;
        else {
            const OCP_DBL b0 = xiC[0] * CONV7;
            const OCP_DBL b1 = xiC[1] * CONV7;
            surTen           = 0;
            for (USI i = 0; i < NC; i++)
                surTen += parachor[i] * (b0 * x[0][i] - b1 * x[1][i]);
            surTen = pow(surTen, 4.0);
        }
        misTerm->AssignValue(bulkId, surTen);
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/