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

MixtureComp::MixtureComp(const EoSparam& param, const USI& tar)
{
    // if Water don't exist?
    // for Mixture class
    numPhase = param.numPhase + 1;
    numCom   = param.numComp + 1;
    Allocate();

    // for MixtureComp class
    NC    = param.numComp;
    NPmax = param.numPhase;

    zi.resize(NC);

    // comp.resize(NC);
    // for (USI i = 0; i < NC; i++) {
    //	comp[i] = COMP(param.COM[i]);
    //}

    Cname = param.Cname;
    if (param.Tc.activity)
        Tc = param.Tc.data[tar];
    else
        OCP_ABORT("TCRIT hasn't been input!");
    if (param.Pc.activity)
        Pc = param.Pc.data[tar];
    else
        OCP_ABORT("PCRIT hasn't been input!");

    if (param.Vc.activity)
        Vc = param.Vc.data[tar];
    else if (param.Zc.activity) {
        Zc = param.Zc.data[tar];
        Vc.resize(NC);
        for (USI i = 0; i < NC; i++) {
            Vc[i] = 10.73159 * Zc[i] * Tc[i] / Pc[i];
        }
    } else
        OCP_ABORT("VCRIT or ZCRIT hasn't been input!");

    if (param.MW.activity)
        MWC = param.MW.data[tar];
    else
        OCP_ABORT("MW hasn't been input!");
    if (param.Acf.activity)
        Acf = param.Acf.data[tar];
    else
        OCP_ABORT("ACF hasn't been input!");
    if (param.OmegaA.activity)
        OmegaA = param.OmegaA.data[tar];
    else
        OmegaA.resize(NC, 0.457235529);
    if (param.OmegaB.activity)
        OmegaB = param.OmegaB.data[tar];
    else
        OmegaB.resize(NC, 0.077796074);

    if (param.Vshift.activity) {
        Vshift = param.Vshift.data[tar];
        for (USI i = 0; i < NC; i++)
            Vshift[i] *= (GAS_CONSTANT * OmegaB[i] *  Tc[i] / Pc[i]);
    }      
    else
        Vshift.resize(NC, 0);

    ParachorAct = true;
    if (param.Parachor.activity)
        Parachor = param.Parachor.data[tar];
    else
        ParachorAct = false;

    if (param.Vcvis.activity)
        Vcvis = param.Vcvis.data[tar];
    else if (param.Zcvis.activity) {
        Zcvis = param.Zcvis.data[tar];
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

    miscible = param.miscible;
    if (miscible && !ParachorAct) {
        OCP_ABORT("PARACHOR has not been Input!");
    }
    

    CallId();

    USI len   = NC * NC;
    USI count = 0;
    BIC.resize(len, 0);

    if (param.BIC[tar].size() != len) {
        USI iter = 0;
        for (USI i = 1; i < NC; i++) {
            for (USI j = 0; j < i; j++) {
                BIC[i * NC + j] = param.BIC[tar][iter];
                BIC[j * NC + i] = BIC[i * NC + j];
                iter++;
            }
        }
    } else {
        BIC = param.BIC[tar];
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

void MixtureComp::InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin,
                            const OCP_DBL& Tin, const OCP_DBL* Sjin,
                            const OCP_DBL& Vpore, const OCP_DBL* Ziin)
{
    ///////////////////////////////////////////////////
    // test      
            
    //P = 3978.475;     
    //Ni[0] = 640.533; Ni[1] = 129.459; Ni[2] = 8.631; Ni[3] = 20.138; 
    //Ni[4] = 57.537;  Ni[5] = 43.153;  Ni[6] = 14.384;

    //P = 3980.258;
    //Ni[0] = 630.125; Ni[1] = 136.094; Ni[2] = 9.073; Ni[3] = 21.170;
    //Ni[4] = 60.486;  Ni[5] = 45.365;  Ni[6] = 15.122;

    //T = 620.000;
    //Nt = Dnorm1(NC, &Ni[0]);
    //nu[0] = Nt;
    //setZi();
    //NP = 1;
    //x[0] = zi;
    //CalAiBi();
    //CalAjBj(Aj[0], Bj[0], x[0]);
    //SolEoS(Zj[0], Aj[0], Bj[0]);
    //vf = 0;
    //OCP_DBL tmp;
    //for (USI j = 0; j < NP; j++) {

    //    vector<OCP_DBL>& xj = x[j];
    //    tmp = Zj[j] * GAS_CONSTANT * T / P;
    //    for (USI i = 0; i < NC; i++) {
    //        tmp -= xj[i] * Vshift[i];
    //    }
    //    vC[j] = tmp * nu[j];
    //    vf += vC[j];
    //    xiC[j] = 1 / tmp;
    //}
    //cout << "Done!" << endl;


    // test
    /////////////////////////////////////////////////////////

    ftype = 0;
    lNP   = 0;

    Nt    = 1;
    nu[0] = 1;
    setPT(Pin, Tin);
    setZi(Ziin);
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
    Nt         = Vpore * (1 - Sw) / vf;
    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nt, &nu[0]);
    // correct vj, vf with new Nt
    Dscalar(NPmax, Nt, &v[0]);
    vf *= Nt;
    CalVfiVfp();
    // Calculate Ni
    for (USI i = 0; i < NC; i++) {
        Ni[i] = zi[i] * Nt;
    }

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = true;
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
    Nt += Ni[Wcid];
    v[Wpid] = CONV1 * Ni[Wcid] * bw;
    vf += v[Wpid];
    vfi[Wcid] = CONV1 * bw;
    vfp += CONV1 * Ni[Wcid] * bwp;

    // Calculate Sj
    CalSaturation();
}

void MixtureComp::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin,
                        const USI& Myftype, const USI& lastNP,
    const OCP_DBL* lastKs)
{
    ftype = Myftype;
    lNP   = lastNP;
    if (lNP == 2) {
        Dcopy(NC, &lKs[0], lastKs);
    }
    

    CalFlash(Pin, Tin, Niin);
    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    CalVfiVfp();

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = true;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt += Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;
    mu[Wpid]    = data[3];
    xi[Wpid]    = 1 / (CONV1 * bw);
    rho[Wpid]   = std_RhoW / bw;
    v[Wpid]     = CONV1 * Ni[Wcid] * bw;
    vf += v[Wpid];
    vfi[Wcid] = CONV1 * bw;
    vfp += CONV1 * Ni[Wcid] * bwp;

    // Calculate Sj
    CalSaturation();
}

void MixtureComp::FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
                             const OCP_DBL* Niin, const USI& Myftype, const USI& lastNP,
                             const OCP_DBL* lastKs)
{
    ftype = Myftype;
    lNP   = lastNP;
    if (lNP == 2) {
        Dcopy(NC, &lKs[0], lastKs);
    }

    CalFlash(Pin, Tin, Niin);
    // Calculate derivates for hydrocarbon phase and components
    // d vf / d Ni, d vf / d P
    CalVfiVfp();
    // d xi / dP, d xi / d Ni, d xi / d xij
    CalXiPNX();

    // cout << setprecision(6);
    // for (USI j = 0; j < NP; j++) {
    //	cout << "xi" << endl;
    //	cout << xiC[j] << endl;
    //	cout << "xij" << endl;
    //	for (USI i = 0; i < NC; i++)
    //		cout << x[j][i] << "    ";
    //	cout << endl;
    //	cout << "xiP" << endl;
    //	cout << xiPC[j] << endl;
    //	cout << "xix" << endl;
    //	for (USI i = 0; i < NC; i++)
    //		cout << xixC[j * NC + i] << "    ";
    //	cout << endl;
    //}

    // Water Properties
    USI Wpid                  = numPhase - 1;
    USI Wcid                  = numCom - 1;
    phaseExist[Wpid]          = true;
    xij[Wpid * numCom + Wcid] = 1.0;
    Nt += Ni[Wcid];
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;
    mu[Wpid]    = data[3];
    xi[Wpid]    = 1 / (CONV1 * bw);
    rho[Wpid]   = std_RhoW / bw;
    muP[Wpid]   = cdata[3];
    xiP[Wpid]   = -bwp / (bw * bw * CONV1);
    rhoP[Wpid]  = CONV1 * xiP[Wpid] * std_RhoW;
    v[Wpid]     = CONV1 * Ni[Wcid] * bw;
    vfi[Wcid]   = CONV1 * bw;
    OCP_DBL vwp = CONV1 * Ni[Wcid] * bwp;
    vf += v[Wpid];
    vfp += vwp;

    // hydrocarbon
    CaldXsdXp(v[Wpid], vwp);
    CalRhoPX();
    CalMuPX();
    // Correct Sj
    CalSaturation();

    // Calculate water derivatives
    // only dSj / dNw, dSw / dP , dSw / dNk needs to be considered
    USI ncol = 1 + numCom;
    // dSj / dNw
    for (USI j = 0; j < NP; j++) {
        dXsdXp[(j + 1) * ncol - 1] = -S[j] * CONV1 * bw / vf;
    }
    OCP_DBL* Dtmp = &dXsdXp[Wpid * ncol];
    // dSw / dP
    OCP_DBL vfvf = vf * vf;
    Dtmp[0]      = (vwp * vf - v[Wpid] * vfp) / vfvf;
    Dtmp++;
    // dSw / d Nk
    for (USI k = 0; k < NC; k++) {
        Dtmp[k] = -v[Wpid] * vfi[k] / vfvf;
    }
    // dSw / d Nw
    Dtmp[NC] = (CONV1 * bw * vf - v[Wpid] * vfi[Wcid]) / vfvf;

    // Print dXsdXp
    // USI myrow = (numCom + 1) * numPhase;
    // USI mycol = (numCom + 1);
    // cout << "dXsdXp" << endl;
    // for (USI i = 0; i < myrow; i++) {
    //	for (USI j = 0; j < mycol; j++) {
    //		cout << dXsdXp[i * mycol + j] << "    ";
    //	}
    //	cout << endl;
    //}
}

void MixtureComp::CalFlash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    setPT(Pin, Tin);
    setNi(Niin);

    // P = 3141.918998650599405664;
    // T = 620;
    // Ni[0] = 41602.131481952070316765;
    // Ni[1] = 11440.761357378274624352;
    // Ni[2] = 1955.137836561682888714;
    // Ni[3] = 337.768923412141248264;
    // Ni[4] = 276.965621027262329790;
    // Ni[5] = 93.492433020052502002;
    //
    // P = 3141.891199116627376497;
    // T = 620;
    // Ni[0] = 41602.572555854152597021;
    // Ni[1] = 11440.882654558850845206;
    // Ni[2] = 1955.158565313887720549;
    // Ni[3] = 337.772504504064784214;
    // Ni[4] = 276.968557470729479064;
    // Ni[5] = 93.493424245039747689;

    // P = 3141.891199116627376497;
    // T = 620.000000000000000000;
    // Ni[0] = 41602.572555854152597021;
    // Ni[1] = 11440.882654558850845206;
    // Ni[2] = 1955.158565313887720549;
    // Ni[3] = 337.772504504064784214;
    // Ni[4] = 276.968557470729479064;
    // Ni[5] = 93.493424245039747689;

    // P = 3211.937138268530816276;
    // T = 620.000000000000000000;
    // Ni[0] = 41692.343588740142877214;
    // Ni[1] = 11485.541047459752007853;
    // Ni[2] = 1973.898169547999714268;
    // Ni[3] = 354.912643337034239721;
    // Ni[4] = 289.460231227657970976;
    // Ni[5] = 97.367441114446421579;

    nu[0] = 1;
    nu[1] = 0;
    // Water is excluded
    Nt = Dnorm1(NC, &Ni[0]);
    setZi();
    PhaseEquilibrium();
    // Next, nu represents moles of phase instead of molar fraction of phase
    Dscalar(NP, Nt, &nu[0]);
    CalMW();
    CalVfXiRho();
    CalViscosity();
    CalSurfaceTension();
    IdentifyPhase();
    CopyPhase();
}

OCP_DBL MixtureComp::XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
                             const OCP_DBL* Ziin)
{
    // assume that only single phase exists here
    if (Ziin[numCom - 1] > 1 - 1e-6) {
        // water phase
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0   = data[0];
        OCP_DBL bw0   = data[1];
        OCP_DBL cbw   = data[2];
        OCP_DBL bw    = bw0 * (1 - cbw * (P - Pw0));
        OCP_DBL xitmp = 1 / (CONV1 * bw);
        return xitmp;
    } else {
        // hydrocarbon phase
        setPT(Pin, Tin);
        setZi(Ziin);
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

OCP_DBL MixtureComp::RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
                              const OCP_DBL* Ziin)
{
    OCP_DBL xitmp = XiPhase(Pin, Tin, Ziin);
    OCP_DBL rhotmp;
    // assume that only single phase exists here
    if (Ziin[numCom - 1] > 1 - 1e-6) {
        // water phase
        rhotmp = std_RhoW * (CONV1 * xitmp);
        return rhotmp;
    } else {
        // hydrocarbon phase
        x[0] = zi;
        CalMW();
        rhotmp = MW[0] * xitmp;
        return rhotmp;
    }
}

OCP_DBL MixtureComp::GammaPhaseW(const OCP_DBL& Pin)
{
    PVTW.Eval_All(0, Pin, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = (bw0 * (1 - cbw * (Pin - Pw0)));

    return std_GammaW / bw;
}

OCP_DBL MixtureComp::GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
                                  const OCP_DBL* Ziin)
{
    // assume that only single phase exists here, no matter it's oil or gas
    OCP_DBL rhotmp = RhoPhase(Pin, Tin, Ziin);
    return rhotmp * GRAVITY_FACTOR;
}

void MixtureComp::CallId()
{
    lId = 0;
    for (USI i = 1; i < NC; i++) {
        if (MWC[i] < MWC[lId]) lId = i;
    }
}

void MixtureComp::CalSurfaceTension()
{
    // be careful! 
    // phase molar densities should be converted into gm-M/cc here
    if (miscible) {
        if (NP == 1)  surTen = 100;
        else {
            const OCP_DBL unitF = CONV3 / (CONV4 * 1E3); // lbm/ft3 -> gm-M/cc
            const OCP_DBL b0 = xiC[0] * unitF;
            const OCP_DBL b1 = xiC[1] * unitF;
            surTen = 0;
            for (USI i = 0; i < NC; i++)
                surTen += Parachor[i] * (b0 * x[0][i] - b1 * x[1][i]);
            surTen = pow(surTen, 4.0);
            // cout << surTen << endl;
        }
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

    USI flag = CubicRoot(a, b, c, true); // True with NT
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

void MixtureComp::CalFugPhi(vector<OCP_DBL>& phiT, vector<OCP_DBL>& fugT,
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
        for (int k = 0; k < NC; k++) {
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
        for (int k = 0; k < NC; k++) {
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
        for (int k = 0; k < NC; k++) {
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
    OCP_DBL       tmp;// , tmp01, tmp02;
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
            for (int k = 0; k < NC; k++) {
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
            S[j] = v[j] / vf;
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
    phiN.resize(NC * NC);
    skipMatSTA.resize(NC * NC);
    eigenSkip.resize(NC);
    leigenWork = 2 * NC + 1;
    eigenWork.resize(leigenWork);
    lKs.resize(NC);

    tmpRR.resize(NC);
    resRR.resize(NPmax - 1);
    resSP.resize(NC * NPmax);
    JmatSP.resize(NC * NC * NPmax * NPmax);
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
    pivot.resize((NC + 1) * NPmax, 1);
    lJmatWork = NC * (NPmax - 1);
    JmatWork.resize(lJmatWork);
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
    NP   = 1;
    x[0] = zi;
    CalAiBi();
    CalAjBj(Aj[0], Bj[0], x[0]);
    SolEoS(Zj[0], Aj[0], Bj[0]);

    flagSkip = true;

    // tmpFtype = ftype;
    // ftype = 0;

    switch (ftype) {
        case 0:
            // flash from single phase
            CalKwilson();
            while (!PhaseStable()) {
                NP++;
                PhaseSplit();
                if (NP == NPmax || NP == 1) break;
            }
            break;
        case 1:
            // Skip Phase Stability analysis, only single phase exists
            NP = 1;
            // cout << "NP = 1   ";
            break;

        case 2:
            // Skip Phase Stability analysis, two phases exist
            CalKwilson();
            NP                    = 2;
            EoSctrl.SSMsp.conflag = false;
            EoSctrl.NRsp.conflag  = false;
            Ks[NP - 2]            = Kw[0];
            SplitSSM(true);
            SplitNR();
            while (!EoSctrl.NRsp.conflag) {
                SplitSSM(true);
                SplitNR();
                if (!CheckSplit()) break;
                if (EoSctrl.SSMsp.conflag) break;
            }
            CheckSplit();
            break;

        default:
            OCP_ABORT("Wrong flash type!");
            break;
    }

    if (NP > 1)  flagSkip = false;
    if (NP == 1 && ftype == 0 && flagSkip) {
        CalPhiNSTA();
        AssembleSkipMatSTA();
        MinEigenSY(NC, &skipMatSTA[0], &eigenSkip[0], &eigenWork[0], leigenWork);
        // PrintDX(NC, &eigenSkip[0]);
        // cout << "done!" << endl;
    }
}

bool MixtureComp::PhaseStable()
{
    if (NP == 1) {
        testPId = 0;
    } else {
        CalMW();
        testPId = FindMWmax();
    }

    EoSctrl.SSMsta.curIt = 0;
    EoSctrl.NRsta.curIt  = 0;

    // Test if a phase is stable, if stable return true, else return false
    bool flag;
    USI  tmpNP = NP;

    if (lNP == 0) {
        // strict stability ananlysis
        flag = StableSSM(testPId);
    } 
    else {
        flag = StableSSM01(testPId);
        if (!flag) tmpNP++;

        if (tmpNP != lNP) {
            flag = StableSSM(testPId);
            flagSkip = false;
        }
    }
    SSMSTAiters += EoSctrl.SSMsta.curIt;
    NRSTAiters += EoSctrl.NRsta.curIt;
    //cout << "Yt = " << setprecision(8) << scientific << Yt << "  " << setw(2)
    //    << EoSctrl.SSMsta.curIt << "  " << setw(2) << EoSctrl.NRsta.curIt << "  "
    //    << lNP << "  " << tmpNP << "  " << tmpFtype << "   "
    //    << (lNP == tmpNP ? "N" : "Y") << "  ";
    //if (lNP == 1 && tmpNP == 2 && tmpFtype == 1) {
    //    cout << "AAA" << "  ";
    //}
    return flag;
}

bool MixtureComp::StableSSM01(const USI& Id)
{
    OCP_DBL Stol  = EoSctrl.SSMsta.tol2;
    USI     maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL Ktol  = EoSctrl.SSMsta.Ktol;
    OCP_DBL dYtol = EoSctrl.SSMsta.dYtol;
    // OCP_DBL& Sk = EoSctrl.SSMsta.curSk;
    OCP_DBL Se, Sk, dY;

    bool    flag, Tsol;
    USI     iter, k;

    const vector<OCP_DBL>& xj = x[Id];
    CalFugPhi(phi[Id], fug[Id], xj);
    const vector<OCP_DBL>& fugId = fug[Id];
    vector<OCP_DBL>&       ks    = Ks[0];

    for (k = 0; k < 2; k++) {

        ks   = Kw[k];
        iter = 0;
        flag = false;
        Tsol = false;
        while (true) {
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
                    flag = true;
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

            // if (Yt > 1 + 1E-12 && Yt < 2 && iter > 0) {
            //	cout << setprecision(6) << scientific << Se;
            //	cout << "         ";
            //	cout << setprecision(6) << scientific << Sk;
            //	cout << "         ";
            //	cout << setprecision(6) << scientific << dY;
            //	cout << "         ";
            //	cout << setprecision(12) << scientific << Yt;
            //	cout << "         " << iter << endl;
            //	PrintDX(NC, &xj[0]);
            //	PrintDX(NC, &Y[0]);
            //}

            // if (P < 20 || true) {
            //	cout << setprecision(6) << scientific << Se;
            //	cout << "         ";
            //	cout << setprecision(6) << scientific << Sk;
            //	cout << "         ";
            //	cout << setprecision(6) << scientific << dY;
            //	cout << "         ";
            //	cout << setprecision(12) << scientific << Yt;
            //	cout << "         " << iter << endl;
            //	PrintDX(NC, &xj[0]);
            //	PrintDX(NC, &Y[0]);
            //}

            iter++;
            if (Se < Stol) {
                flag = true;
                break;
            }
            if (Sk < Ktol) {
                // Sk < Ktol -> trivial solution
                flag = true;
                Tsol = true;
                break;
            }
            if (iter > maxIt) {
                break;
            }

            // Record last Y with di
            di = Y;
        }
        // if (!Tsol) {
        // flag = StableNR(Id);
        //}
        //cout << "Yt = " << setprecision(8) << scientific << Yt << "   " << setw(2)
        //    << "Sk = " << setprecision(3) << scientific << Sk << "   " << setw(2)
        //    << iter << "  ";
        EoSctrl.SSMsta.curIt += iter;

        if (flag && Yt > 1 - 0.1 && Sk > 1) {
            flagSkip = false;
        }
        if (flag && Yt > 1 + eYt) {
            return false;
        }
    }
    return true;
}

bool MixtureComp::StableSSM(const USI& Id)
{
    const vector<OCP_DBL>& xj = x[Id];
    CalFugPhi(phi[Id], fug[Id], xj);
    const vector<OCP_DBL>& fugId = fug[Id];

    for (USI i = 0; i < NC; i++) {
        di[i] = phi[Id][i] * xj[i];
    }

    OCP_DBL Stol  = EoSctrl.SSMsta.tol2;
    USI     maxIt = EoSctrl.SSMsta.maxIt;
    OCP_DBL eYt   = EoSctrl.SSMsta.eYt;
    OCP_DBL Se;
    bool    flag;
    USI     iter;
    USI     k;

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

        flag = true;
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
                flag = false;
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
            return false;
        }
    }
    // cout << "Yt = " << setprecision(12) << scientific << Yt << "    " << setw(3)
    //     << EoSctrl.SSMsta.curIt << "    " << setw(3) << EoSctrl.NRsta.curIt << "   "
    //     << flag << "   " << k << "   " << 1 << "   ";
    /*if (!flag) {
        OCP_WARNING("SSM not converged in Stability Analysis");
    }*/
    return true;
}

bool MixtureComp::StableNR(const USI& Id)
{

#ifdef DEBUG
    cout << endl << "Stable NR Begins !" << endl << endl;
#endif // DEBUG

    for (USI i = 0; i < NC; i++) {
        resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
    }

#ifdef DEBUG
    // PrintDX(NC, &resSTA[0]);
#endif // DEBUG

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
        SYSSolve(1, &uplo, NC, &JmatSTA[0], &resSTA[0], &pivot[0], &JmatWork[0], lJmatWork);
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
            // cout << "---------------" << endl;
            // cout << setprecision(6) << scientific << Se0 << endl;
            // cout << setprecision(6) << scientific << Se << endl;
            // cout << setprecision(6) << scientific << Yt << endl;
            // PrintDX(NC, &x[Id][0]);
            // PrintDX(NC, &Y[0]);
            // cout << "---------------" << endl;
            EoSctrl.NRsta.curIt += iter;
            return false;
        }

#ifdef DEBUG
        // PrintDX(NC, &resSTA[0]);
#endif // DEBUG
    }
    // if (iter > 0 && !isfinite(Se)) {
    //	cout << "---------------" << endl;
    //	cout << setprecision(6) << scientific << Se0 << endl;
    //	cout << setprecision(6) << scientific << Se << endl;
    //	cout << setprecision(6) << scientific << Yt << endl;
    //	PrintDX(NC, &x[Id][0]);
    //	PrintDX(NC, &Y[0]);
    //	cout << "done!" << endl;
    //	cout << "---------------" << endl;
    //}
    // if (Yt > 1 + 1E-8 && iter > 0) {
    //	cout << setprecision(6) << scientific << Se0 << endl;
    //	cout << setprecision(6) << scientific << Se << endl;
    //	PrintDX(NC, &x[Id][0]);
    //	PrintDX(NC, &Y[0]);
    //	cout << endl;
    //}

    EoSctrl.NRsta.curIt += iter;
    return true;
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

    OCP_DBL  C, D, E, G;
    OCP_DBL  Cxk, Dxk, Exk, Gxk;
    OCP_DBL  aik;
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
    EoSctrl.SSMsp.conflag = false;
    EoSctrl.NRsp.conflag = false;
    EoSctrl.SSMsp.curIt = 0;
    EoSctrl.NRsp.curIt = 0;

    SplitSSM(false);
    SplitNR();
    while (!EoSctrl.NRsp.conflag) {
        SplitSSM(true);
        SplitNR();
        if (!CheckSplit()) break;
        if (EoSctrl.SSMsp.conflag) break;
    }
    CheckSplit();

    SSMSPiters += EoSctrl.SSMsp.curIt;
    NRSPiters += EoSctrl.NRsp.curIt;

    //cout << scientific << setprecision(8);
    //for (USI i = 0; i < NC; i++) {
    //    cout << x[0][i] / x[1][i] << "   ";
    //}
    //cout << endl;
    //cout << "Yt = " << scientific << setprecision(12) << Yt << "   "
    //    << setw(3) << "SSMtol = " << setprecision(6) << sqrt(EoSctrl.SSMsp.realTol) << "   "
    //    << setw(3) << EoSctrl.SSMsp.curIt << "   "
    //    << setw(3) << "NRtol = " << setprecision(6) << EoSctrl.NRsp.realTol << "   "
    //    << setw(3) << EoSctrl.NRsp.curIt << "   "
    //    << setw(2) << lNP << "   " << setw(2) << NP << "   "
    //    << (lNP == NP ? "N" : "Y") << "   ";
}

bool MixtureComp::CheckSplit()
{
    if (NP == 2) {

        OCP_DBL tmp   = 0;
        OCP_DBL nuMax = max(nu[0], nu[1]);

        for (USI i = 0; i < NC; i++) {
            tmp += (x[0][i] - x[1][i]) * (x[0][i] - x[1][i]);
        }

        if (nuMax < 1 && EoSctrl.NRsp.conflag && isfinite(tmp)) {
            // accept this result
        } else {
            if (!isfinite(tmp) || (1 - nuMax) < 1E-3) {
                // single phase
                // cout << "-------------------------------------------" << endl;
                // cout << fixed << scientific << setprecision(12);
                // cout << Yt << "   " << nu[0] << "   " << nu[1] << "   ";
                // cout << sqrt(EoSctrl.SSMsp.realTol) << "   " << EoSctrl.SSMsp.conflag
                // << "   "; cout << sqrt(EoSctrl.NRsp.realTol) << "   " <<
                // EoSctrl.NRsp.conflag << endl; for (USI i = 0; i < NC; i++) 	cout <<
                //zi[i] << "   "; cout << endl; for (USI j = 0; j < NP; j++) { 	for (USI
                //i = 0; i < NC; i++) { 		cout << x[j][i] << "   ";
                //	}
                //	cout << endl;
                //}
                NP    = 1;
                x[0]  = zi;
                nu[0] = 1;
                CalAjBj(Aj[0], Bj[0], x[0]);
                SolEoS(Zj[0], Aj[0], Bj[0]);
                return false;
            }
        }
    }
    return true;
}

void MixtureComp::SplitSSM(const bool& flag)
{
#ifdef DEBUG
    cout << "SSMSP Begins!" << endl;
#endif

    if (NP == 2) {
        SplitSSM2(flag);
    } else {
        SplitSSM3(flag);
    }
}

void MixtureComp::SplitSSM2(const bool& flag)
{
    // NP = 2 in this case
    // Ks is very IMPORTANT!
    // flag = true : Restart SSM
    // flag = false : New SSM
    EoSctrl.SSMsp.conflag = true;
    OCP_DBL Se            = 1;
    OCP_DBL Stol          = EoSctrl.SSMsp.tol2;
    USI     maxIt         = EoSctrl.SSMsp.maxIt;

    if (!flag) {
        if (lNP == 2) {
            Ks[NP - 2] = lKs;
        }
        else {
            if (Yt < 1.1 || true) {
                Ks[NP - 2] = Kw[0];
            }
            else {
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

#ifdef DEBUG
    PrintX();
#endif // DEBUG

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
#ifdef DEBUG
        PrintX();
#endif // DEBUG

        iter++;
        if (iter > maxIt) {
            // OCP_WARNING("SSM not converged in Phase Spliting!");
            EoSctrl.SSMsp.conflag = false;
            break;
        }
    }

    EoSctrl.SSMsp.realTol = Se;
    EoSctrl.SSMsp.curIt += iter; 
}

void MixtureComp::SplitSSM3(const bool& flag) {}

void MixtureComp::RachfordRice2() ///< Used when NP = 2
{
    const vector<OCP_DBL>& Ktmp = Ks[0];
    OCP_DBL                Kmin = Ktmp[0];
    OCP_DBL                Kmax = Ktmp[0];

    for (USI i = 1; i < NC; i++) {
        if (Ktmp[i] < Kmin) Kmin = Ktmp[i];
        if (Ktmp[i] > Kmax) Kmax = Ktmp[i];
    }

    OCP_DBL numin = 1 / (1 - Kmax);
    OCP_DBL numax = 1 / (1 - Kmin);

    nu[0] = 0.5 * (numin + numax);

    // Solve RR with NR
    fill(tmpRR.begin(), tmpRR.end(), 0);
    OCP_DBL rj, J, dnuj, tmp;

    USI     iter  = 0;
    OCP_DBL RRtol = EoSctrl.RR.tol;
    while (true) {

        rj = 0;
        J  = 0;
        for (USI i = 0; i < NC; i++) {
            tmpRR[i] = 1 + nu[0] * (Ktmp[i] - 1);
            rj += zi[i] * (Ktmp[i] - 1) / tmpRR[i];
            J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmpRR[i] * tmpRR[i]);
        }

        if (fabs(rj) < RRtol || iter > EoSctrl.RR.maxIt) break;

        dnuj = -rj / J;
        tmp  = nu[0] + dnuj;
        if (tmp < numax && tmp > numin) {
            nu[0] = tmp;
        } else {
            if (dnuj > 0) {
                nu[0] = (nu[0] + numax) / 2;
            } else {
                nu[0] = (nu[0] + numin) / 2;
            }
        }

        iter++;
    }

    // if (iter > EoSctrl.RR.maxIt) {
    //	OCP_WARNING("RR2 not converged!");
    //}

    nu[1] = 1 - nu[0];

#ifdef DEBUG
    cout << nu[0] << endl;
#endif // DEBUG
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

void MixtureComp::SplitNR()
{
    EoSctrl.NRsp.conflag = false;
    for (USI j = 0; j < NP; j++) {
        nu[j] = fabs(nu[j]);
    }

    USI len = NC * (NP - 1);
    x2n();
    CalResSP();
    OCP_DBL eNR0;
    OCP_DBL eNR   = Dnorm2(len, &resSP[0]);
    OCP_DBL NRtol = EoSctrl.NRsp.tol;
    OCP_DBL alpha;

#ifdef DEBUG
    cout << "NRSP Begins!\n";
#endif

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

        SYSSolve(1, &uplo, len, &JmatSP[0], &resSP[0], &pivot[0], &JmatWork[0], lJmatWork);
        // PrintDX(NC, &resSP[0]);

        alpha = CalStepNRsp();

        n[NP - 1] = zi;
        for (USI j = 0; j < NP - 1; j++) {
            Daxpy(NC, alpha, &resSP[j * NC], &n[j][0]);
            Daxpy(NC, -1, &n[j][0], &n[NP - 1][0]);

            nu[j] = Dnorm1(NC, &n[j][0]);
            for (USI i = 0; i < NC; i++) {
                x[j][i] = n[j][i] / nu[j];
            }
#ifdef DEBUG
            PrintDX(NC, &x[j][0]);
#endif
        }
        for (USI i = 0; i < NC; i++) {
            n[NP - 1][i] = fabs(n[NP - 1][i]);
        }
        nu[NP - 1] = Dnorm1(NC, &n[NP - 1][0]);
        for (USI i = 0; i < NC; i++) {
            x[NP - 1][i] = n[NP - 1][i] / nu[NP - 1];
        }

#ifdef DEBUG
        PrintDX(NC, &x[NP - 1][0]);
        cout << "---------------------" << endl;
#endif

        CalFugPhiAll();
        CalResSP();
        eNR = Dnorm2(len, &resSP[0]);
        iter++;
        if (eNR > eNR0 || iter > EoSctrl.NRsp.maxIt) {
            break;
        }

        en = 0;
        for (USI j = 0; j < NP; j++) {
            Daxpy(NC, -1, &n[j][0], &ln[j][0]);
            en += Dnorm2(NC, &ln[j][0]);
        }
        if (en / (NP * NC) < 1E-8) {
            EoSctrl.NRsp.conflag = true;
            break;
        }
    }
    EoSctrl.NRsp.realTol = eNR;
    if (eNR < NRtol) EoSctrl.NRsp.conflag = true;
    EoSctrl.NRsp.curIt += iter;
}

void MixtureComp::CalResSP()
{
    // So it equals -res
    for (USI j = 0; j < NP - 1; j++) {
        for (USI i = 0; i < NC; i++) {
            // if zi is too small, resSP[j][i] = 0?
            resSP[j * NC + i] = log(fug[NP - 1][i] / fug[j][i]);
        }
    }
}

void MixtureComp::CalFugNAll()
{
    OCP_DBL C, D, E, G;
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
            Znj[i] =
                ((bj - zj) * An[i] +
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

    //cout << endl << "Jmat" << endl;
    //for (USI i = 0; i < NC; i++) {
    //    for (USI j = 0; j < NC; j++) {
    //        cout << scientific << setprecision(6) << JmatSP[i * NC + j] << "   ";
    //    }
    //    cout << endl;
    //}
}

void MixtureComp::CalPhiNSTA()
{
    OCP_DBL C, D, E, G;
    OCP_DBL Cnk, Dnk, Enk, Gnk;
    OCP_DBL tmp, aik;

	// 0 th phase
	const OCP_DBL& aj = Aj[0];
	const OCP_DBL& bj = Bj[0];
	const OCP_DBL& zj = Zj[0];
	const vector<OCP_DBL>& xj = x[0];
	vector<OCP_DBL>& Znj = Zn[0];

	for (USI i = 0; i < NC; i++) {
		tmp = 0;
		for (USI m = 0; m < NC; m++) {
			tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
		}
		An[i] = 2 / nu[0] * (tmp - aj);
		Bn[i] = 1 / nu[0] * (Bi[i] - bj);
		Znj[i] =
			((bj - zj) * An[i] +
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

    //cout << endl << "phiN" << endl;
    //for (USI i = 0; i < NC; i++) {
    //    for (USI j = 0; j < NC; j++) {
    //        cout << phiN[i * NC + j] << "   ";
    //    }
    //    cout << endl;
    //}
}

void MixtureComp::AssembleSkipMatSTA()
{
    // Sysmetric Matrix
    // stored by colum
    vector<OCP_DBL>& xj = x[0];

    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j <= i; j++) {
            skipMatSTA[i * NC + j] = delta(i, j) + sqrt(xj[i] * xj[j]) * phiN[i * NC + j];
            skipMatSTA[j * NC + i] = skipMatSTA[i * NC + j];
        }
    }

/*    cout << endl << "skipMatSTA" << endl;
    for (USI i = 0; i < NC; i++) {
        for (USI j = 0; j < NC; j++) {
            cout << skipMatSTA[i * NC + j] << "   ";
        }
        cout << endl;
    } */  
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
}

void MixtureComp::IdentifyPhase()
{
    phaseExist[0] = false;
    phaseExist[1] = false;
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
            phaseExist[1] = true;
        } else {
            phaseLabel[0] = OIL;
            phaseExist[0] = true;
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
        phaseExist[0] = true;
        phaseExist[1] = true;
    }
}

void MixtureComp::CopyPhase()
{
    // copy v, x, mu, xi, rho
    for (USI j = 0; j < NP; j++) {
        USI j1 = phaseLabel[j];
        v[j1]  = vC[j];
        Dcopy(NC, &xij[j1 * numCom], &x[j][0]);
        mu[j1]  = muC[j];
        xi[j1]  = xiC[j];
        rho[j1] = rhoC[j];
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

    // test
    //T = 750;
    //P = 4867.594;
    //x[0][0] = 0.657538;
    //x[0][1] = 0.014955;
    //x[0][2] = 0.003126;
    //x[0][3] = 0.008749;
    //x[0][4] = 0.008200;
    //x[0][5] = 0.017913;
    //x[0][6] = 0.035430;
    //x[0][7] = 0.254089;
    //CalAiBi();
    //CalAjBj(Aj[0], Bj[0], x[0]);
    //SolEoS(Zj[0], Aj[0], Bj[0]);
    //CalMW();
    //CalVfXiRho();
    //test

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

        if (muA[3] <= 0.18 && false) {
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

        OCP_DBL  C, D, E, G;
        OCP_DBL  Cxk, Dxk, Exk, Gxk;
        OCP_DBL  aik;
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
}

void MixtureComp::CalFugPAll()
{

    OCP_DBL C, D, E, G;
    OCP_DBL Cp, Dp, Ep, Gp;
    OCP_DBL tmp;

    for (USI j = 0; j < NP; j++) {

        vector<OCP_DBL>& fugp = fugP[j];
        vector<OCP_DBL>& xj   = x[j];
        OCP_DBL&         aj   = Aj[j];
        OCP_DBL&         bj   = Bj[j];
        OCP_DBL&         zj   = Zj[j];

        OCP_DBL Ap = aj / P;
        OCP_DBL Bp = bj / P;
        Zp[j]      = ((bj - zj) * Ap +
                 ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj)) +
                  ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj -
                  (delta1 + delta2 - 1) * zj * zj) *
                     Bp) /
                (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj +
                 (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));

        G  = (zj + delta1 * bj) / (zj + delta2 * bj);
        Gp = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) *
             (Bp * zj - Zp[j] * bj);
        for (USI i = 0; i < NC; i++) {

            C = xj[i] * P / (zj - bj);
            // D = Bi[i] / bj * (zj - 1);

            tmp = 0;
            for (USI m = 0; m < NC; m++) {
                tmp += (1 - BIC[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
            }

            E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

            Cp = xj[i] / ((zj - bj) * (zj - bj)) * ((zj - bj) - P * (Zp[j] - Bp));
            Dp = Bi[i] / bj * Zp[j];
            // Ep = 0;

            fugp[i] = 1 / C * Cp + Dp + E / G * Gp;
        }
    }
}

void MixtureComp::CalXiPNX()
{
    // call CalVfiVfp() before
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

        tmp           = -xiC[j] * xiC[j] * CgTP;
        OCP_DBL* xiXj = &xixC[j * NC];
        for (USI i = 0; i < NC; i++) {
            xiXj[i] = tmp * Zx[i];
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
        // Now dnkjdNP reach to dnkj / dNi after CalVfiVfp
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

    // Copoy to xiP, xix
    fill(xiP.begin(), xiP.end(), 0.0);
    fill(xix.begin(), xix.end(), 0.0);
    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1      = phaseLabel[j];
        xiP[j1] = xiPC[j];
        Dcopy(NC, &xix[j1 * numCom], &xixC[j * NC]);
    }
}

void MixtureComp::CalRhoPX()
{
    // CaldXsdXp() shoule be called before
    // Cal rhoP, rhoX
    // Water has been calculated
    fill(rhox.begin(), rhox.end() - numCom, 0.0);
    fill(rhoP.begin(), rhoP.end() - 1, 0.0);

    USI j1;
    for (USI j = 0; j < NP; j++) {
        j1       = phaseLabel[j];
        rhoP[j1] = xiPC[j] * MW[j];
        for (USI i = 0; i < NC; i++) {
            rhox[j1 * numCom + i] = xixC[j * NC + i] * MW[j] + xiC[j] * MWC[i];
        }
    }
    if (NP > 1) {
        // correct rhoP
        // use rhsDer, see CaldXsdXp(), attention that it only contains hydrocarbon
        OCP_DBL        tmp;
        const OCP_DBL* xijP = &rhsDer[NP];
        for (USI j = 0; j < NP; j++) {
            tmp = 0;
            for (USI i = 0; i < NC; i++) {
                tmp += MWC[i] * xijP[i];
            }
            rhoP[phaseLabel[j]] += xiC[j] * tmp;
            xijP += NC;
        }
    }
}

void MixtureComp::CalMuPX()
{
    fill(muP.begin(), muP.end() - 1, 0.0);
    fill(mux.begin(), mux.end() - numCom, 0.0);

    CalMuPXLBC();
}

void MixtureComp::CalMuPXLBC()
{
    OCP_DBL val1IJ, val2IJ;
    OCP_DBL der1IJ, der2IJ, der3J, der4J, der6J, der7J, der8J;
    OCP_DBL Tri, tmp;
    OCP_DBL xTj, xPj, xVj;
    OCP_DBL derxTj, derxPj, derMWj, derxVj;

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
        if (muAuxj[3] <= 0.18) {
            muP[phaseLabel[0]] = (2.05 * 1E-4) * der7J / muAuxj[2];
        } else {
            der8J              = der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
            muP[phaseLabel[0]] = (4 * 1E-4) * pow(muAuxj[4], 3) * der8J / muAuxj[2];
        }
        // Calculate dmuj / xkj
        const USI bId = numCom * phaseLabel[0];
        for (USI k = 0; k < NC; k++) {
            derxTj = Tc[k];
            derxPj = Pc[k];
            derMWj = MWC[k];
            der3J  = 0;
            for (USI i = 0; i < NC; i++) {
                val1IJ = muAux1I[i] / sqrt(MW[0]);
                der1IJ = -(1 / 2) * muAux1I[i] * pow(MW[0], -1.5) * derMWj;
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
                 pow(xTj, 1.0 / 6) * (0.5 / MW[0] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[0]) * pow(xPj, 2.0 / 3));
            der7J = xixC[k] * xVj + xiC[0] * Vcvis[k];
            if (muAuxj[3] <= 0.18) {
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
    } else {
        // NP > 1
        // Calculate dmuj / dP
        // use rhsDer (after CaldXsdXp)
        for (USI j = 0; j < NP; j++) {
            const OCP_DBL*         xijP   = &rhsDer[NP + j * NC];
            const vector<OCP_DBL>& xj     = x[j];
            const vector<OCP_DBL>& muAuxj = muAux[j];
            const USI              bId    = phaseLabel[j];
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
            der7J *= xiC[j];
            der7J += xiPC[j] * xVj;
            der6J =
                5.4402 *
                (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                 pow(xTj, 1.0 / 6) * (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
            if (muAuxj[3] <= 0.18) {
                muP[bId] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    2.05 * 1E-4 * (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            } else {
                der8J =
                    der7J * (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                muP[bId] =
                    (der3J * muAuxj[1] - muAuxj[0] * der4J) / (muAuxj[1] * muAuxj[1]) +
                    1E-4 *
                        (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                         (pow(muAuxj[4], 4) - 1) * der6J) /
                        (muAuxj[2] * muAuxj[2]);
            }
        }

        // Calculate dmuj / dxkj
        for (USI j = 0; j < NP; j++) {
            const vector<OCP_DBL>& xj     = x[j];
            const vector<OCP_DBL>& muAuxj = muAux[j];
            const USI              bId    = numCom * phaseLabel[j];
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
                der6J = 5.4402 *
                        (1.0 / 6 * pow(xTj, -5.0 / 6) * derxTj -
                         pow(xTj, 1.0 / 6) *
                             (0.5 / MW[j] * derMWj + 2.0 / 3 / xPj * derxPj)) /
                        (sqrt(MW[j]) * pow(xPj, 2.0 / 3));
                der7J = xixC[k] * xVj + xiC[0] * derxVj;
                if (muAuxj[3] <= 0.18) {
                    mux[bId + k] = (der3J * muAuxj[1] - muAuxj[0] * der4J) /
                                       (muAuxj[1] * muAuxj[1]) +
                                   2.05 * 1E-4 *
                                       (der7J * muAuxj[2] - muAuxj[3] * der6J) /
                                       (muAuxj[2] * muAuxj[2]);
                } else {
                    der8J = der7J *
                            (LBCcoef[1] +
                             muAuxj[3] * (2 * LBCcoef[2] +
                                          muAuxj[3] * (3 * LBCcoef[3] +
                                                       muAuxj[3] * 4 * LBCcoef[4])));
                    mux[bId + k] = (der3J * muAuxj[1] - muAuxj[0] * der4J) /
                                       (muAuxj[1] * muAuxj[1]) +
                                   1E-4 *
                                       (4 * pow(muAuxj[4], 3) * der8J * muAuxj[2] -
                                        (pow(muAuxj[4], 4) - 1) * der6J) /
                                       (muAuxj[2] * muAuxj[2]);
                }
            }
        }
    }
}

void MixtureComp::CalVfiVfp()
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
        vfp = CgTP * nu[0] * (Zp[0] - zj / P);
    } else {
        // NP > 1
        CalFugNAll();
        CalFugPAll();
        AssembleMatVfiVfp();
        AssembleRhsVfiVfp();
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
        vfp = 0;
        for (USI j = 0; j < NP; j++) {
            vfp += nu[j] / P * (Zp[j] * P - Zj[j]);
            for (USI k = 0; k < NC; k++) {
                vfp += dnkjdNP[j * NC + k] * (Zj[j] + nu[j] * Zn[j][k]);
            }
        }
        vfp *= CgTP;
    }

    // check
    // for (USI i = 0; i < NC; i++) {
    //	if (!isfinite(vfi[i])) {
    //		cout << "NAN or INF" << endl;
    //	}
    //}
    // if (!isfinite(vfp)) {
    //	cout << "NAN or INF" << endl;
    //}
}

void MixtureComp::AssembleMatVfiVfp()
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

void MixtureComp::AssembleRhsVfiVfp()
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
}

void MixtureComp::CaldXsdXp()
{

    // fugP
    // cout << fixed << setprecision(6);
    // for (USI j = 0; j < NP; j++) {
    //	cout << "fp[" << j << "]" << endl;
    //	for (USI i = 0; i < NC; i++) {
    //		cout << fugP[j][i] << "    ";
    //	}
    //	cout << endl;
    //}
    // cout << "-----------------------------" << endl;

    // NP > 1 in this function
    // only hydrocarbon was considered
    // if water exists, vf, vfp, and S should be updated
    // CalVfiVfp() and CalXiPNX() should be called before
    CalFugXAll();

    // for (USI j = 0; j < NP; j++) {
    //	cout << "fxx[" << j << "]" << endl;
    //	for (USI i = 0; i < NC; i++) {
    //		for (USI k = 0; k < NC; k++) {
    //			cout << fugX[j][i * NC + k] << "    ";
    //		}
    //		cout << endl;
    //	}
    //	cout << endl;
    //}cout << "-----------------------------" << endl;

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
        // dFn / dSm
        MTmp += NP;
        // dFn / dxnm, m = j
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
            // dFn / dSm
            MTmp += NP;
            // dFn / dxnm, m = j or m = np-1
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

    USI row01 = NP * (NC + 1);
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
    // d fij / dP have been calculated in CalVfiVfp() before
    fill(JmatTmp.begin(), JmatTmp.end(), 0.0);
    MTmp = &JmatTmp[0];
    OCP_DBL tmp01;
    OCP_DBL tmp02;
    // dFn / dXs
    for (USI i = 0; i < NC; i++) {
        // dFn / dP
        tmp01 = 0;
        tmp02 = 0;
        for (USI j = 0; j < NP; j++) {
            tmp01 += S[j] * x[j][i] * xiPC[j];
            tmp02 += S[j] * x[j][i] * xiC[j];
        }
        MTmp[0] = -(vf * tmp01 + vfp * tmp02); // in OCP
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

    // dFx / dXs
    MTmp += NP * (NC + 1);

    // dFf / dXs
    const vector<OCP_DBL>& fugPNP = fugP[NP - 1];
    for (USI j = 0; j < NP - 1; j++) {
        const vector<OCP_DBL>& fugPj = fugP[j];
        for (USI i = 0; i < NC; i++) {
            // dFf / dP
            MTmp[0] = -(fugPj[i] - fugPNP[i]);
            MTmp += (NC + 1);
        }
    }

    USI row02 = NP * (NC + 1);
    USI col02 = NC + 1;
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

void MixtureComp::CaldXsdXp(const OCP_DBL& vw, const OCP_DBL& vwp)
{
    // Calculate derivates for hydrocarbon phase and components
    // if water exists, vf, vfp, and S should be updated
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);
    if (NP == 1) {
        USI     bId = (numCom + 1) * phaseLabel[0];
        OCP_DBL tmp = vf * vf;
        // only dSj / dP, dSj / dNk ---- hydrocarbon phase
        dXsdXp[bId] = (vfp * vw - vf * vwp) / tmp;
        for (USI k = 0; k < NC; k++) {
            dXsdXp[bId + k + 1] = (vfi[k] * vw) / tmp;
        }
    } else {
        // NP > 1
        // Calculate Saturation
        for (USI j = 0; j < NP; j++) {
            S[j] = vC[j] / vf;
        }
        CaldXsdXp();
        // copy from rhsDer to dXsdXp
        OCP_DBL*       DTmp;
        const OCP_DBL* STmp;
        USI            nrow  = NP * (NC + 1); // row num of rhsDer
        USI            ncol  = NC + 1;        // col num of rhsDer
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
}

USI MixtureComp::CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c,
                           const bool& NTflag) const
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

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jan/05/2022      Create file                          */
/*----------------------------------------------------------------------------*/