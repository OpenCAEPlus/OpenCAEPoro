/*! \file    MixtureBO3_ODGW.cpp
 *  \brief   Used for the condition where oil, gas, disolve gas, water exist.
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureBO.hpp"

///////////////////////////////////////////////
// BOMixture_ODGW
///////////////////////////////////////////////

BOMixture_ODGW::BOMixture_ODGW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_ODGW;
    numPhase    = 3;
    numCom      = 3;
    BOMixtureInit(rs_param);

    PVTW.Setup(rs_param.PVTW_T.data[i]);
    PVCO.Setup(rs_param.PVCO_T.data[i]);
    PVDG.Setup(rs_param.PVDG_T.data[i]);

    data.resize(6, 0);
    cdata.resize(6, 0);
}

void BOMixture_ODGW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    FlashIMPEC(Pin, Tin, Niin, 0, 0, 0);
}

void BOMixture_ODGW::InitFlashIMPEC(const OCP_DBL& Pin,
                                    const OCP_DBL& Pbbin,
                                    const OCP_DBL& Tin,
                                    const OCP_DBL* Sjin,
                                    const OCP_DBL& Vpore,
                                    const OCP_DBL* Ziin,
                                    const OCP_USI& bId)
{
    for (USI j = 0; j < 3; j++) phaseExist[j] = OCP_FALSE;
    for (USI i = 0; i < 3; i++) Ni[i] = 0;
    fill(xij.begin(), xij.end(), 0.0);

    P    = Pin;
    S[1] = Sjin[1]; // gas saturation
    S[2] = Sjin[2]; // water saturation

    // Water Property
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[2]  = data[3];
    rho[2] = std_RhoW / bw;
    xi[2]  = 1 / (CONV1 * bw);
    Ni[2]  = Vpore * S[2] * xi[2];

    USI phasecae;

    if (1 - S[1] - S[2] < TINY) {
        if (S[1] < TINY)
            phasecae = PHASE_W; // case 1 : water, no oil, no gas
        else
            phasecae = PHASE_GW; // case 2 : water, gas, no oil
    } else if (S[1] < TINY)
        phasecae = PHASE_OW; // case 3 : water, oil, no gas
    else
        phasecae = PHASE_ODGW; // case 4 : water, oil, gas

    switch (phasecae) {
        case PHASE_W:
            {
                // water
                phaseExist[2]  = OCP_TRUE;
                S[0]           = 0;
                S[1]           = 0;
                S[2]           = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // hypothetical Gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg = data[1] * (CONV1 / 1000);

                // total
                vj[0]  = 0;
                vj[1]  = 0;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[2];
                vfP    = CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_GW:
            {
                // water, gas
                phaseExist[1]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

                mu[1]  = data[2];
                xi[1]  = 1 / 1000 / bg;
                rho[1] = std_RhoG / bg;
                Ni[1]  = Vpore * S[1] * xi[1];

                vj[0] = 0;
                vj[1] = 1000 * Ni[1] * bg; // Ni[0] = 0;
                vj[2] = CONV1 * Ni[2] * bw;
                // total
                vf     = vj[1] + vj[2];
                S[0]   = 0;
                S[1]   = vj[1] / vf;
                S[2]   = vj[2] / vf;
                vfP    = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_OW:
            {
                // water, oil, unsaturated
                phaseExist[0] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                OCP_DBL Pbb = Pbbin;
                PVCO.Eval_All(0, Pbb, data, cdata);
                OCP_DBL rs      = data[1];
                OCP_DBL bosat   = data[2];
                OCP_DBL muosat  = data[3];
                OCP_DBL cbosat  = data[4];
                OCP_DBL cmuosat = data[5];
                OCP_DBL bo      = bosat * (1 - cbosat * (P - Pbb));
                OCP_DBL bop     = -bosat * cbosat;
                OCP_DBL dBo_drs = bo / bosat * cdata[2] +
                                  bosat * (cdata[4] * (Pbb - P) + cbosat * cdata[0]);
                dBo_drs /= cdata[1];

                Ni[0]  = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
                Ni[1]  = Ni[0] * rs;
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;
                mu[0]  = muosat * (1 + cmuosat * (P - Pbb));

                xij[0 * 3 + 0] = Ni[0] / (Ni[0] + Ni[1]);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // total
                vj[0]  = CONV1 * Ni[0] * bo;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[0] + vj[2];
                S[0]   = vj[0] / vf;
                S[1]   = 0;
                S[2]   = vj[2] / vf;
                vfP    = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
                vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
                vfi[1] = CONV1 * dBo_drs;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_ODGW:
            {
                phaseExist[0] = OCP_TRUE;
                phaseExist[1] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs     = data[1];
                OCP_DBL bo     = data[2];
                OCP_DBL crs    = cdata[1];
                OCP_DBL cbosat = cdata[2];

                mu[0]  = data[3];
                Ni[0]  = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);
                Ni[1]       = Vpore * S[1] / bg / 1000 + Ni[0] * rs;
                xi[1]       = 1 / CONV1 / data[1];
                rho[1]      = std_RhoG / bg;
                mu[1]       = data[2];

                xij[0 * 3 + 0] = 1 / (1 + rs);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // total
                vj[0] = CONV1 * Ni[0] * bo;
                vj[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
                vj[2] = CONV1 * Ni[2] * bw;

                vf   = vj[0] + vj[1] + vj[2];
                S[0] = vj[0] / vf;
                S[1] = vj[1] / vf;
                S[2] = vj[2] / vf;
                vfP  = CONV1 * Ni[0] * cbosat +
                      1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) +
                      CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
    }
}

void BOMixture_ODGW::InitFlashFIM(const OCP_DBL& Pin,
                                  const OCP_DBL& Pbbin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Sjin,
                                  const OCP_DBL& Vpore,
                                  const OCP_DBL* Ziin,
                                  const OCP_USI& bId)
{
    // Calculate Ni only frist, then call FlashDer()
    for (USI i = 0; i < 3; i++) Ni[i] = 0;

    P    = Pin;
    S[1] = Sjin[1]; // gas saturation
    S[2] = Sjin[2]; // water saturation

    // Water Property
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));

    xi[2] = 1 / (CONV1 * bw);
    Ni[2] = Vpore * S[2] * xi[2];

    USI phasecae;

    if (1 - S[1] - S[2] < TINY) {
        if (S[1] < TINY)
            phasecae = PHASE_W; // case 1 : water, no oil, no gas
        else
            phasecae = PHASE_GW; // case 2 : water, gas, no oil
    } else if (S[1] < TINY)
        phasecae = PHASE_OW; // case 3 : water, oil, no gas
    else
        phasecae = PHASE_ODGW; // case 4 : water, oil, gas

    switch (phasecae) {
        case PHASE_W:
            {
                // water
                break;
            }
        case PHASE_GW:
            {
                // water, gas

                // Ni[1]
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg = data[1] * (CONV1 / 1000);
                xi[1]      = 1 / 1000 / bg;
                Ni[1]      = Vpore * S[1] * xi[1];
                break;
            }
        case PHASE_OW:
            {
                // water, oil, unsaturated
                // Ni[0] and Ni[1]
                OCP_DBL Pbb = Pbbin;
                PVCO.Eval_All(0, Pbb, data, cdata);
                OCP_DBL rs     = data[1];
                OCP_DBL bosat  = data[2];
                OCP_DBL cbosat = data[4];
                OCP_DBL bo     = bosat * (1 - cbosat * (P - Pbb));

                Ni[0] = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
                Ni[1] = Ni[0] * rs;

                break;
            }
        case PHASE_ODGW:
            {
                // oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];
                Ni[0]      = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg = data[1] * (CONV1 / 1000);
                Ni[1]      = Vpore * S[1] / bg / 1000 + Ni[0] * rs;
                break;
            }
    }
    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, 0);
}

void BOMixture_ODGW::FlashIMPEC(const OCP_DBL& Pin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Niin,
                                const USI&     lastNP,
                                const OCP_DBL* xijin,
                                const OCP_USI& bId)
{
    for (USI j = 0; j < 3; j++) phaseExist[j] = OCP_FALSE;
    fill(xij.begin(), xij.end(), 0.0);

    P  = Pin;
    Nt = 0;
    for (USI i = 0; i < numCom; i++) {
        Ni[i] = Niin[i];
        Nt += Ni[i];
    }

    // Water property
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[2]  = data[3];
    xi[2]  = 1 / (CONV1 * bw);
    rho[2] = std_RhoW / bw;

    USI     phasecase;
    OCP_DBL Rs_sat = PVCO.Eval(0, P, 1);

    if (Ni[0] < Nt * TINY) {
        if (Ni[1] <= Ni[0] * Rs_sat)
            phasecase = PHASE_W; // water, no oil, no gas
        else
            phasecase = PHASE_GW; // water, gas, no oil
    } else if (Ni[1] <= Ni[0] * Rs_sat)
        phasecase = PHASE_OW; // water, oil, no gas
    else
        phasecase = PHASE_ODGW; // water, oil ,gas

    switch (phasecase) {
        case PHASE_W:
            {
                // water
                phaseExist[2]  = OCP_TRUE;
                S[0]           = 0;
                S[1]           = 0;
                S[2]           = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // hypothetical Gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg = data[1] * (CONV1 / 1000);

                // total
                vj[0]  = 0;
                vj[1]  = 0;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[2];
                vfP    = CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_GW:
            {
                // water, gas
                phaseExist[1]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

                mu[1]  = data[2];
                xi[1]  = 1 / 1000 / bg;
                rho[1] = std_RhoG / bg;

                // total
                vj[0] = 0;
                vj[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
                vj[2] = CONV1 * Ni[2] * bw;

#ifdef DEBUG
                if (vj[1] <= 0) OCP_ABORT("gas volume <= 0");
#endif // DEBUG

                vf     = vj[1] + vj[2];
                S[0]   = 0;
                S[1]   = vj[1] / vf;
                S[2]   = vj[2] / vf;
                vfP    = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_OW:
            {
                // water, oil
                phaseExist[0]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[0 * 3 + 0] = Ni[0] / (Ni[0] + Ni[1]);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // oil property
                OCP_DBL rs = Ni[1] / Ni[0];
                PVCO.Eval_All(1, rs, data, cdata);
                OCP_DBL pbb     = data[0];
                OCP_DBL bosat   = data[2];
                OCP_DBL muosat  = data[3];
                OCP_DBL cbosat  = data[4];
                OCP_DBL cmuosat = data[5];
                OCP_DBL bo      = bosat * (1 - cbosat * (P - pbb));
                OCP_DBL bop     = -bosat * cbosat;
                OCP_DBL dBo_drs = bo / bosat * cdata[2] +
                                  bosat * (cdata[4] * (pbb - P) + cbosat * cdata[0]);

                mu[0]  = muosat * (1 + cmuosat * (P - pbb));
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

                // total
                vj[0]  = CONV1 * Ni[0] * bo;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[0] + vj[2];
                S[0]   = vj[0] / vf;
                S[1]   = 0;
                S[2]   = vj[2] / vf;
                vfP    = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
                vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
                vfi[1] = CONV1 * dBo_drs;
                vfi[2] = CONV1 * bw;

                break;
            }
        case PHASE_ODGW:
            {
                phaseExist[0] = OCP_TRUE;
                phaseExist[1] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs     = data[1];
                OCP_DBL bo     = data[2];
                OCP_DBL crs    = cdata[1];
                OCP_DBL cbosat = cdata[2];

                mu[0]  = data[3];
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

                mu[1]  = data[2];
                xi[1]  = 1 / CONV1 / data[1];
                rho[1] = std_RhoG / bg;

                // total
                xij[0 * 3 + 0] = 1 / (1 + rs);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                vj[0] = CONV1 * Ni[0] * bo;
                vj[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
                vj[2] = CONV1 * Ni[2] * bw;
                vf    = vj[0] + vj[1] + vj[2];
                S[0]  = vj[0] / vf;
                S[1]  = vj[1] / vf;
                S[2]  = vj[2] / vf;
                vfP   = CONV1 * Ni[0] * cbosat +
                      1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) +
                      CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                break;
            }
    }
}

void BOMixture_ODGW::FlashFIM(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const OCP_DBL* Sjin,
                              const USI&     lastNP,
                              const OCP_DBL* xijin,
                              const OCP_USI& bId)
{
    for (USI j = 0; j < 3; j++) {
        phaseExist[j] = OCP_FALSE;
        rhoP[j]       = 0;
        xiP[j]        = 0;
        muP[j]        = 0;
    }
    fill(xij.begin(), xij.end(), 0.0);
    fill(rhox.begin(), rhox.end(), 0.0);
    fill(xix.begin(), xix.end(), 0.0);
    fill(mux.begin(), mux.end(), 0.0);
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);
    fill(pSderExist.begin(), pSderExist.end(), OCP_FALSE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0);

    P  = Pin;
    Nt = 0;
    for (USI i = 0; i < numCom; i++) {
        Ni[i] = Niin[i];
        Nt += Ni[i];
    }

    // Water property
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[2]  = data[3];
    xi[2]  = 1 / CONV1 / bw;
    rho[2] = std_RhoW / bw;

    muP[2]  = cdata[3];
    xiP[2]  = -bwp / (bw * bw * CONV1);
    rhoP[2] = CONV1 * xiP[2] * std_RhoW;

    USI     phasecase;
    OCP_DBL Rs_sat = PVCO.Eval(0, P, 1);

    if (Ni[0] < Nt * TINY) {
        if (Ni[1] <= Ni[0] * Rs_sat)
            phasecase = PHASE_W; // water, no oil, no gas
        else
            phasecase = PHASE_GW; // water, gas, no oil
    }
    // Beacareful when NO = NG = 0, if it's possible.
    else if (Ni[1] <= Ni[0] * Rs_sat)
        phasecase = PHASE_OW; // water, oil, no gas
    else
        phasecase = PHASE_ODGW; // water, oil ,gas

    switch (phasecase) {
        case PHASE_W:
            {
                // water
                phaseExist[2]  = OCP_TRUE;
                S[0]           = 0;
                S[1]           = 0;
                S[2]           = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // hypothetical Gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg = data[1] * (CONV1 / 1000);

                // total
                vj[0]  = 0;
                vj[1]  = 0;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[2];
                vfP    = CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                // dXsdXp[0] = 0; // dSo / dP
                dXsdXp[1] = CONV1 * bo / vf; // dSo / dNo
                // dXsdXp[2] = 0;   // dSo / dNg;
                // dXsdXp[3] = 0;   // dSo / dNw;

                // dXsdXp[1 * 4 + 0] = 0;   // dSg / dP
                // dXsdXp[1 * 4 + 1] = 0;   // dSg / dNo
                dXsdXp[1 * 4 + 2] = 1000 * bg / vf; // dSg / dNg
                // dXsdXp[1 * 4 + 3] = 0;   // dSg / dNw

                dXsdXp[2 * 4 + 0] = (Ni[2] * bwp * CONV1 - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (CONV1 * bw - S[2] * vfi[2]) / vf; // dSw / dNw

                pSderExist[0] = OCP_TRUE;
                pSderExist[1] = OCP_TRUE;
                pSderExist[2] = OCP_TRUE;
                // pVnumCom[0] = 0; pVnumCom[1] = 0; pVnumCom[2] = 0;

                break;
            }
        case PHASE_GW:
            {
                // water, gas
                phaseExist[1]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // hypothetical Oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs = data[1];
                OCP_DBL bo = data[2];

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

                mu[1]  = data[2];
                xi[1]  = 1 / 1000 / bg;
                rho[1] = std_RhoG / bg;

                muP[1]  = cdata[2];
                xiP[1]  = -cdata[1] / (CONV1 * data[1] * data[1]);
                rhoP[1] = 1000 * std_RhoG * xiP[1];

                // total
                vj[0] = 0;
                vj[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
                vj[2] = CONV1 * Ni[2] * bw;

#ifdef DEBUG
                if (vj[1] <= 0) OCP_ABORT("gas volume <= 0");
#endif // DEBUG

                vf     = vj[1] + vj[2];
                S[0]   = 0;
                S[1]   = vj[1] / vf;
                S[2]   = vj[2] / vf;
                vfP    = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                dXsdXp[0] = 0;               // dSo / dP
                dXsdXp[1] = CONV1 * bo / vf; // dSo / dNo
                dXsdXp[2] = 0;               // dSo / dNg
                dXsdXp[3] = 0;               // dSo / dNw

                dXsdXp[1 * 4 + 0] = (1000 * Ni[1] * cbg - S[1] * vfP) / vf; // dSg / dP
                dXsdXp[1 * 4 + 1] = (-1000 * rs * bg - S[1] * vfi[0]) / vf; // dSg / dNo
                dXsdXp[1 * 4 + 2] = (1000 * bg - S[1] * vfi[1]) / vf;       // dSg / dNg
                dXsdXp[1 * 4 + 3] = -S[1] * vfi[2] / vf;                    // dSg / dNw

                dXsdXp[2 * 4 + 0] = (CONV1 * Ni[2] * bwp - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (CONV1 * bw - S[2] * vfi[2]) / vf; // dSw / dNw

                pSderExist[0] = OCP_TRUE;
                pSderExist[1] = OCP_TRUE;
                pSderExist[2] = OCP_TRUE;
                // pVnumCom[0] = 0; pVnumCom[1] = 0; pVnumCom[2] = 0;

                break;
            }
        case PHASE_OW:
            {
                // water, oil
                phaseExist[0]  = OCP_TRUE;
                phaseExist[2]  = OCP_TRUE;
                xij[0 * 3 + 0] = Ni[0] / (Ni[0] + Ni[1]);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                // oil property
                OCP_DBL rs = Ni[1] / Ni[0];
                PVCO.Eval_All(1, rs, data, cdata);
                OCP_DBL pbb     = data[0];
                OCP_DBL bosat   = data[2];
                OCP_DBL muosat  = data[3];
                OCP_DBL cbosat  = data[4];
                OCP_DBL cmuosat = data[5];
                OCP_DBL bo      = bosat * (1 - cbosat * (P - pbb));
                OCP_DBL bop     = -bosat * cbosat;
                OCP_DBL dBo_drs = (1 - cbosat * (P - pbb)) * cdata[2] +
                                  bosat * (cdata[4] * (pbb - P) + cbosat * cdata[0]);

                mu[0]  = muosat * (1 + cmuosat * (P - pbb));
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

                muP[0] = muosat * cmuosat;
                xiP[0] = -(1 + rs) * bop / (CONV1 * bo * bo);
                rhoP[0] =
                    -(std_RhoO + (1000 / CONV1) * rs * std_RhoG) / (bo * bo) * bop;

                OCP_DBL muo_rs = mu[0] / muosat * cdata[3] +
                                 muosat * (cdata[5] * (P - pbb) - cmuosat * cdata[0]);
                OCP_DBL xio_rs =
                    1 / CONV1 / bo - (1 + rs) * dBo_drs / (CONV1 * bo * bo);
                OCP_DBL rhoo_rs =
                    (1000 / CONV1) * std_RhoG / bo -
                    (std_RhoO + (1000 / CONV1) * rs * std_RhoG) * dBo_drs / (bo * bo);

                // total
                vj[0]  = CONV1 * Ni[0] * bo;
                vj[2]  = CONV1 * Ni[2] * bw;
                vf     = vj[0] + vj[2];
                S[0]   = vj[0] / vf;
                S[1]   = 0;
                S[2]   = vj[2] / vf;
                vfP    = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
                vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
                vfi[1] = CONV1 * dBo_drs;
                vfi[2] = CONV1 * bw;

                dXsdXp[0] =
                    (-CONV1 * Ni[0] * cbosat * bosat - S[0] * vfP) / vf; // dSo / dP
                dXsdXp[1] =
                    (CONV1 * bo - CONV1 * dBo_drs * (Ni[1] / Ni[0]) - S[0] * vfi[0]) /
                    vf;                                             // dSo / dNo
                dXsdXp[2] = (CONV1 * dBo_drs - S[0] * vfi[1]) / vf; // dSo / dNg
                dXsdXp[3] = -S[0] / vf * vfi[2];                    // dSo / dNw

#ifdef OCP_OLD_FIM
                dXsdXp[2 * 4 + 0] = (Ni[2] * bwp * CONV1 - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (CONV1 * bw - S[2] * vfi[2]) / vf; // dSw / dNw

                dXsdXp[3 * 4 + 1] = Ni[1] / pow((Ni[0] + Ni[1]), 2);  // d Xoo / d No
                dXsdXp[3 * 4 + 2] = -Ni[0] / pow((Ni[0] + Ni[1]), 2); // d Xoo / d Ng
                dXsdXp[4 * 4 + 1] = -dXsdXp[3 * 4 + 1];               // d Xgo / d No
                dXsdXp[4 * 4 + 2] = -dXsdXp[3 * 4 + 2];               // d Xgo / d Ng
#else
                dXsdXp[1 * 4 + 0] = (Ni[2] * bwp * CONV1 - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[1 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[1 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[1 * 4 + 3] = (CONV1 * bw - S[2] * vfi[2]) / vf; // dSw / dNw

                dXsdXp[2 * 4 + 1] = Ni[1] / pow((Ni[0] + Ni[1]), 2);  // d Xoo / d No
                dXsdXp[2 * 4 + 2] = -Ni[0] / pow((Ni[0] + Ni[1]), 2); // d Xoo / d Ng
                dXsdXp[3 * 4 + 1] = -dXsdXp[2 * 4 + 1];               // d Xgo / d No
                dXsdXp[3 * 4 + 2] = -dXsdXp[2 * 4 + 2];               // d Xgo / d Ng

#endif // OCP_OLD_FIM

                pSderExist[0] = OCP_TRUE;
                pSderExist[2] = OCP_TRUE;
                pVnumCom[0]   = 2;

                OCP_DBL tmp = xij[0] * xij[0];

                mux[0] = -muo_rs * xij[1] / tmp; // dMuo / dXoo
                mux[1] = muo_rs / xij[0];        // dMuo / dXgo

                xix[0] = -xio_rs * xij[1] / tmp; // dXio / dXoo
                xix[1] = xio_rs / xij[0];        // dXio / dXgo

                rhox[0] = -rhoo_rs * xij[1] / tmp; // dRhoo / dXoo
                rhox[1] = rhoo_rs / xij[0];        // dRhoo / dXgo

                // New! right I think
                // OCP_DBL tmp_new = (1 + rs) * (1 + rs);

                // mux[0] = -muo_rs * tmp_new;         // dMuo / dXoo
                // mux[1] = muo_rs * tmp_new;          // dMuo / dXgo

                // xix[0] = -xio_rs * tmp_new;         // dXio / dXoo
                // xix[1] = xio_rs * tmp_new;          // dXio / dXgo

                // rhox[0] = -rhoo_rs * tmp_new;       // dRhoo / dXoo
                // rhox[1] = rhoo_rs * tmp_new;        // dRhoo / dXgo

                break;
            }
        case PHASE_ODGW:
            {
                phaseExist[0] = OCP_TRUE;
                phaseExist[1] = OCP_TRUE;
                phaseExist[2] = OCP_TRUE;

                // oil property
                PVCO.Eval_All(0, P, data, cdata);
                OCP_DBL rs     = data[1];
                OCP_DBL bo     = data[2];
                OCP_DBL crs    = cdata[1];
                OCP_DBL cbosat = cdata[2];

                mu[0]  = data[3];
                xi[0]  = (1 + rs) / (CONV1 * bo);
                rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

                muP[0] = cdata[3];
                xiP[0] = (crs * bo - (1 + rs) * cbosat) / (bo * bo) / CONV1;
                rhoP[0] =
                    (1000 / CONV1) * std_RhoG * crs / bo -
                    (std_RhoO + (1000 / CONV1) * rs * std_RhoG) * cbosat / (bo * bo);

                // gas property
                PVDG.Eval_All(0, P, data, cdata);
                OCP_DBL bg  = data[1] * (CONV1 / 1000);
                OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

                mu[1]  = data[2];
                xi[1]  = 1 / CONV1 / data[1];
                rho[1] = std_RhoG / bg;

                muP[1]  = cdata[2];
                xiP[1]  = -cdata[1] / (CONV1 * data[1] * data[1]);
                rhoP[1] = 1000 * std_RhoG * xiP[1];

                // total
                xij[0 * 3 + 0] = 1 / (1 + rs);
                xij[0 * 3 + 1] = 1 - xij[0 * 3 + 0];
                xij[1 * 3 + 1] = 1;
                xij[2 * 3 + 2] = 1;

                vj[0] = CONV1 * Ni[0] * bo;
                vj[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
                vj[2] = CONV1 * Ni[2] * bw;
                vf    = vj[0] + vj[1] + vj[2];
                S[0]  = vj[0] / vf;
                S[1]  = vj[1] / vf;
                S[2]  = vj[2] / vf;
                vfP   = CONV1 * Ni[0] * cbosat +
                      1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) +
                      CONV1 * Ni[2] * bwp;
                vfi[0] = CONV1 * bo - 1000 * rs * bg;
                vfi[1] = 1000 * bg;
                vfi[2] = CONV1 * bw;

                dXsdXp[0] = (CONV1 * Ni[0] * cbosat - S[0] * vfP) / vf; // dSo / dP
                dXsdXp[1] = (CONV1 * bo - S[0] * vfi[0]) / vf;          // dSo / dNo
                dXsdXp[2] = -S[0] * vfi[1] / vf;                        // dSo / dNg
                dXsdXp[3] = -S[0] * vfi[2] / vf;                        // dSo / dNw

                dXsdXp[1 * 4 + 0] = (1000 * (Ni[1] - rs * Ni[0]) * cbg -
                                     1000 * Ni[0] * bg * crs - S[1] * vfP) /
                                    vf;                                     // dSg / dP
                dXsdXp[1 * 4 + 1] = (-1000 * rs * bg - S[1] * vfi[0]) / vf; // dSg / dNo
                dXsdXp[1 * 4 + 2] = (1000 * bg - S[1] * vfi[1]) / vf;       // dSg / dNg
                dXsdXp[1 * 4 + 3] = -S[1] * vfi[2] / vf;                    // dSg / dNw

                dXsdXp[2 * 4 + 0] = (Ni[2] * bwp * CONV1 - S[2] * vfP) / vf; // dSw / dP
                dXsdXp[2 * 4 + 1] = -S[2] * vfi[0] / vf;               // dSw / dNo
                dXsdXp[2 * 4 + 2] = -S[2] * vfi[1] / vf;               // dSw / dNg
                dXsdXp[2 * 4 + 3] = (CONV1 * bw - S[2] * vfi[2]) / vf; // dSw / dNw

                dXsdXp[3 * 4 + 0] = -crs / ((1 + rs) * (1 + rs)); // d Xoo / dP
                dXsdXp[4 * 4 + 0] = -dXsdXp[3 * 4 + 0];           // d Xgo / dP

                pSderExist[0] = OCP_TRUE;
                pSderExist[1] = OCP_TRUE;
                pSderExist[2] = OCP_TRUE;
                pVnumCom[0]   = 2;

                break;
            }
    }

#ifdef Debug
    for (auto vj : dXsdXp) {
        if (!isfinite(vj)) {
            OCP_ABORT("Nan in dXsdXp!");
        }
    }
#endif
}

OCP_DBL
BOMixture_ODGW::XiPhase(const OCP_DBL& Pin,
                        const OCP_DBL& Tin,
                        const OCP_DBL* Ziin,
                        const USI&     tarPhase)
{
    if (tarPhase == GAS) {
        OCP_DBL bg = PVDG.Eval(0, Pin, 1);
        // OCP_DBL xig = 1 / (CONV1 * bg);
        OCP_DBL xig = 1 / CONV1 / bg;
        return xig;
    } else if (tarPhase == WATER) {

        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0 = data[0];
        OCP_DBL bw0 = data[1];
        OCP_DBL cbw = data[2];
        OCP_DBL bw  = bw0 * (1 - cbw * (Pin - Pw0));
        OCP_DBL xiw = 1 / CONV1 / bw;
        return xiw;
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

OCP_DBL
BOMixture_ODGW::RhoPhase(const OCP_DBL& Pin,
                         const OCP_DBL& Pbbin,
                         const OCP_DBL& Tin,
                         const OCP_DBL* Ziin,
                         const USI&     tarPhase)
{

    if (tarPhase == OIL) {
        PVCO.Eval_All(0, Pbbin, data, cdata);
        OCP_DBL rs     = data[1];
        OCP_DBL bosat  = data[2];
        OCP_DBL cbosat = data[4];
        OCP_DBL bo     = bosat * (1 - cbosat * (Pin - Pbbin));
        OCP_DBL rhoO   = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;
        return rhoO;
    } else if (tarPhase == GAS) {
        OCP_DBL bg   = PVDG.Eval(0, Pin, 1);
        OCP_DBL rhog = (1000 / CONV1) * std_RhoG / bg;
        return rhog;
    } else if (tarPhase == WATER) {
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0  = data[0];
        OCP_DBL bw0  = data[1];
        OCP_DBL cbw  = data[2];
        OCP_DBL bw   = bw0 * (1 - cbw * (Pin - Pw0));
        OCP_DBL rhow = std_RhoW / bw;
        return rhow;
    } else {
        OCP_ABORT("WRONG tarPhase!");
    }
}

void BOMixture_ODGW::SetupWellOpt(WellOpt&                  opt,
                                  const vector<SolventINJ>& sols,
                                  const OCP_DBL&            Psurf,
                                  const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({0, 0, 1});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        } else if (fluidName == "GAS") {
            vector<OCP_DBL> tmpZi({0, 1, 0});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(GAS);
        } else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }

    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(3, 0);
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
                tmpWght[0] = tmpWght[2] = 1;
                break;
            default:
                OCP_ABORT("WRONG Opt Mode!");
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
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/