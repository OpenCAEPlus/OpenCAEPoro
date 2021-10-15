/*! \file    MixtureBO3.cpp
 *  \brief   MixtureBO3 class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureBO.hpp"

void BOMixture::BOFlash_Sj_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin,
                               const OCP_DBL* Sjin, const OCP_DBL& Vpore)
{
    // initialize    numPhase = 3, numCom = 3
    phaseExist.assign(numPhase, false);
    cij.assign(numPhase * numCom, 0); //
    Ni.assign(numCom, 0);

    P = Pin;
    for (USI j = 0; j < numPhase; j++) S[j] = Sjin[j];

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
        phasecae = PHASE_OGW; // case 4 : water, oil, gas

    switch (phasecae) {
        case PHASE_W: {
            // water
            phaseExist[2] = true;
            S[0]          = 0;
            S[1]          = 0;
            S[2]          = 1;
            cij[8]        = 1;

            // hypothetical Oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs = data[1];
            OCP_DBL bo = data[2];

            // hypothetical Gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg = data[1] * (CONV1 / 1000);

            // total
            v[0]   = 0;
            v[1]   = 0;
            v[2]   = CONV1 * Ni[2] * bw;
            vf     = v[2];
            vfp    = CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_GW: {
            // water, gas
            phaseExist[1] = true;
            phaseExist[2] = true;
            cij[4]        = 1;
            cij[8]        = 1;

            // hypothetical Oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs = data[1];
            OCP_DBL bo = data[2];

            // gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg  = data[1] * (CONV1 / 1000);
            OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

            mu[1]  = data[2];
            xi[1]  = 1 / bg / 1000;
            rho[1] = std_RhoG / bg;
            Ni[1]  = Vpore * S[1] * xi[1];

            v[0] = 0;
            v[1] = 1000 * Ni[1] * bg; // Ni[0] = 0;
            v[2] = CONV1 * Ni[2] * bw;
            // total
            vf     = v[1] + v[2];
            S[0]   = 0;
            S[1]   = v[1] / vf;
            S[2]   = v[2] / vf;
            vfp    = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_OW: {
            // water, oil
            phaseExist[0] = true;
            phaseExist[2] = true;

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

            cij[0] = Ni[0] / (Ni[0] + Ni[1]);
            cij[1] = 1 - cij[0];
            cij[4] = 1;
            cij[8] = 1;

            // total
            v[0]   = CONV1 * Ni[0] * bo;
            v[2]   = CONV1 * Ni[2] * bw;
            vf     = v[0] + v[2];
            S[0]   = v[0] / vf;
            S[1]   = 0;
            S[2]   = v[2] / vf;
            vfp    = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
            vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
            vfi[1] = CONV1 * dBo_drs;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_OGW: {
            phaseExist.assign(3, true);

            // oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs     = data[1];
            OCP_DBL bo     = data[2];
            OCP_DBL crs    = cdata[1];
            OCP_DBL cbosat = cdata[2];

            mu[0]  = data[3];
            Ni[0]  = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
            xi[0]  = (1 + rs) / bo / CONV1;
            rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

            // gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg  = data[1] * (CONV1 / 1000);
            OCP_DBL cbg = cdata[1] * (CONV1 / 1000);
            Ni[1]       = Vpore * S[1] / bg / 1000 + Ni[0] * rs;
            xi[1]       = 1 / data[1] / CONV1;
            rho[1]      = std_RhoG / bg;
            mu[1]       = data[2];

            cij[0] = 1 / (1 + rs);
            cij[1] = 1 - cij[0];
            cij[4] = 1;
            cij[8] = 1;

            // total
            v[0] = CONV1 * Ni[0] * bo;
            v[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
            v[2] = CONV1 * Ni[2] * bw;

            vf   = v[0] + v[1] + v[2];
            S[0] = v[0] / vf;
            S[1] = v[1] / vf;
            S[2] = v[2] / vf;
            vfp  = CONV1 * Ni[0] * cbosat +
                  1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) +
                  CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
    }
}

void BOMixture::BOFlash_Ni_OGW(const OCP_DBL& Pin, const OCP_DBL* Niin)
{
    // initialize    numPhase = 3, numCom = 3
    phaseExist.assign(numPhase, false);
    cij.assign(numPhase * numCom, 0); //

    P          = Pin;
    OCP_DBL NT = 0;
    for (USI i = 0; i < numCom; i++) {
        Ni[i] = Niin[i];
        NT += Ni[i];
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

    if (Ni[0] < NT * TINY) {
        if (Ni[1] < Ni[0] * Rs_sat)
            phasecase = PHASE_W; // water, no oil, no gas
        else
            phasecase = PHASE_GW; // water, gas, no oil
    } else if (Ni[1] < Ni[0] * Rs_sat)
        phasecase = PHASE_OW; // water, oil, no gas
    else
        phasecase = PHASE_OGW; // water, oil ,gas

    switch (phasecase) {
        case PHASE_W: {
            // water
            phaseExist[2] = true;
            S[0]          = 0;
            S[1]          = 0;
            S[2]          = 1;
            cij[8]        = 1;

            // hypothetical Oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs = data[1];
            OCP_DBL bo = data[2];

            // hypothetical Gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg = data[1] * (CONV1 / 1000);

            // total
            v[0]   = 0;
            v[1]   = 0;
            v[2]   = CONV1 * Ni[2] * bw;
            vf     = v[2];
            vfp    = CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_GW: {
            // water, gas
            phaseExist[1] = true;
            phaseExist[2] = true;
            cij[4]        = 1;
            cij[8]        = 1;

            // hypothetical Oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs = data[1];
            OCP_DBL bo = data[2];

            // gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg  = data[1] * (CONV1 / 1000);
            OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

            mu[1]  = data[2];
            xi[1]  = 1 / bg / 1000;
            rho[1] = std_RhoG / bg;

            // total
            v[0] = 0;
            v[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
            v[2] = CONV1 * Ni[2] * bw;
            if (v[1] < 0) v[1] = 0;
            vf     = v[1] + v[2];
            S[0]   = 0;
            S[1]   = v[1] / vf;
            S[2]   = v[2] / vf;
            vfp    = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_OW: {
            // water, oil
            phaseExist[0] = true;
            phaseExist[2] = true;
            cij[0]        = Ni[0] / (Ni[0] + Ni[1]);
            cij[1]        = 1 - cij[0];
            cij[4]        = 1;
            cij[8]        = 1;

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
            v[0]   = CONV1 * Ni[0] * bo;
            v[2]   = CONV1 * Ni[2] * bw;
            vf     = v[0] + v[2];
            S[0]   = v[0] / vf;
            S[1]   = 0;
            S[2]   = v[2] / vf;
            vfp    = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
            vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
            vfi[1] = CONV1 * dBo_drs;
            vfi[2] = CONV1 * bw;

            break;
        }
        case PHASE_OGW: {
            phaseExist.assign(3, true);

            // oil property
            PVCO.Eval_All(0, P, data, cdata);
            OCP_DBL rs     = data[1];
            OCP_DBL bo     = data[2];
            OCP_DBL crs    = cdata[1];
            OCP_DBL cbosat = cdata[2];

            mu[0]  = data[3];
            xi[0]  = (1 + rs) / bo / CONV1;
            rho[0] = (std_RhoO + (1000 / CONV1) * rs * std_RhoG) / bo;

            // gas property
            PVDG.Eval_All(0, P, data, cdata);
            OCP_DBL bg  = data[1] * (CONV1 / 1000);
            OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

            mu[1]  = data[2];
            xi[1]  = 1 / data[1] / CONV1;
            rho[1] = std_RhoG / bg;

            // total
            cij[0] = 1 / (1 + rs);
            cij[1] = 1 - cij[0];
            cij[4] = 1;
            cij[8] = 1;

            v[0] = CONV1 * Ni[0] * bo;
            v[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
            v[2] = CONV1 * Ni[2] * bw;
            vf   = v[0] + v[1] + v[2];
            S[0] = v[0] / vf;
            S[1] = v[1] / vf;
            S[2] = v[2] / vf;
            vfp  = CONV1 * Ni[0] * cbosat +
                  1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) +
                  CONV1 * Ni[2] * bwp;
            vfi[0] = CONV1 * bo - 1000 * rs * bg;
            vfi[1] = 1000 * bg;
            vfi[2] = CONV1 * bw;

            break;
        }
    }
}

OCP_DBL BOMixture::XiPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin)
{
    if (Ziin[1] > 1 - TINY) {
        // inj fluid is gas
        OCP_DBL bg  = PVDG.Eval(0, Pin, 1);
        OCP_DBL xig = (1 / CONV1) / bg;
        return xig;
    } else if (Ziin[2] > 1 - TINY) {
        // inj fluid is water

        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0 = data[0];
        OCP_DBL bw0 = data[1];
        OCP_DBL cbw = data[2];
        OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
        OCP_DBL xiw = (1 / CONV1) / bw;
        return xiw;
    } else {
        ERRORcheck("Wrong Zi!");
        exit(0);
    }
}

OCP_DBL BOMixture::RhoPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin)
{
    if (Ziin[1] > 1 - TINY) {
        // inj fluid is gas
        OCP_DBL bg   = PVDG.Eval(0, Pin, 1);
        OCP_DBL rhog = (1000 / CONV1) * std_RhoG / bg;
        return rhog;
    } else if (Ziin[2] > 1 - TINY) {
        // inj fluid is water

        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0  = data[0];
        OCP_DBL bw0  = data[1];
        OCP_DBL cbw  = data[2];
        OCP_DBL bw   = bw0 * (1 - cbw * (P - Pw0));
        OCP_DBL rhow = std_RhoW / bw;
        return rhow;
    } else {
        ERRORcheck("Wrong Zi!");
        exit(0);
    }
}

OCP_DBL BOMixture::GammaPhaseO_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin)
{

    PVCO.Eval_All(0, Pbbin, data, cdata);
    OCP_DBL rs     = data[1];
    OCP_DBL bosat  = data[2];
    OCP_DBL cbosat = data[4];
    OCP_DBL bo     = bosat * (1 - cbosat * (Pin - Pbbin));
    OCP_DBL gammaO = (std_GammaO + (1000 / CONV1) * rs * std_GammaG) / bo;

    return gammaO;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/