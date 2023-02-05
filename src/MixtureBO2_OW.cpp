/*! \file    MixtureBO2_OW.cpp
 *  \brief   Used for Black Oil model, where only water and oil exist.
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
// BOMixture_OW
///////////////////////////////////////////////

BOMixture_OW::BOMixture_OW(const ParamReservoir& rs_param, const USI& i)
{
    mixtureType = BLKOIL_OW;
    numPhase    = 2;
    numCom      = 2;
    BOMixtureInit(rs_param);

    PVTW.Setup(rs_param.PVTW_T.data[i]);
    PVDO.Setup(rs_param.PVDO_T.data[i]);

    data.resize(5, 0);
    cdata.resize(5, 0);

    phaseExist[0] = OCP_TRUE;
    phaseExist[1] = OCP_TRUE;

    xij[0 * 2 + 0] = 1;
    xij[0 * 2 + 1] = 0;
    xij[1 * 2 + 0] = 0;
    xij[1 * 2 + 1] = 1;
}

void BOMixture_OW::Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
    FlashIMPEC(Pin, Tin, Niin, 0, 0, 0);
}

void BOMixture_OW::InitFlashIMPEC(const OCP_DBL& Pin,
                                  const OCP_DBL& Pbbin,
                                  const OCP_DBL& Tin,
                                  const OCP_DBL* Sjin,
                                  const OCP_DBL& Vpore,
                                  const OCP_DBL* Ziin,
                                  const OCP_USI& bId)
{
    P    = Pin;
    S[1] = Sjin[1];
    // Water Properties
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[1]  = data[3];
    xi[1]  = 1 / (bw * CONV1);
    rho[1] = std_RhoW / bw;

    Ni[1] = Vpore * S[1] * xi[1];

    // Oil Properties
    PVDO.Eval_All(0, P, data, cdata);
    OCP_DBL bo  = data[1];
    OCP_DBL bop = cdata[1];

    mu[0]  = data[2];
    xi[0]  = 1 / (CONV1 * bo);
    rho[0] = std_RhoO / bo;
    Ni[0]  = Vpore * (1 - S[1]) * xi[0];

    vj[0]  = CONV1 * Ni[0] * bo;
    vj[1]  = CONV1 * Ni[1] * bw;
    vf     = vj[0] + vj[1];
    S[0]   = vj[0] / vf;
    S[1]   = vj[1] / vf;
    vfi[0] = CONV1 * bo;
    vfi[1] = CONV1 * bw;
    vfP    = CONV1 * Ni[0] * bop + CONV1 * Ni[1] * bwp;
}

void BOMixture_OW::InitFlashFIM(const OCP_DBL& Pin,
                                const OCP_DBL& Pbbin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Sjin,
                                const OCP_DBL& Vpore,
                                const OCP_DBL* Ziin,
                                const OCP_USI& bId)
{
    P    = Pin;
    S[1] = Sjin[1];
    // Water Properties
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    xi[1]       = 1 / (bw * CONV1);
    Ni[1]       = Vpore * S[1] * xi[1];

    // Oil Properties
    PVDO.Eval_All(0, P, data, cdata);
    OCP_DBL bo = data[1];

    xi[0] = 1 / (CONV1 * bo);
    Ni[0] = Vpore * (1 - S[1]) * xi[0];

    FlashFIM(Pin, Tin, &Ni[0], 0, 0, 0, 0);
}

void BOMixture_OW::FlashIMPEC(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Niin,
                              const USI&     lastNP,
                              const OCP_DBL* xijin,
                              const OCP_USI& bId)
{

    P     = Pin;
    Ni[0] = Niin[0];
    Ni[1] = Niin[1];

    // Water Properties
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[1]  = data[3];
    xi[1]  = 1 / (bw * CONV1);
    rho[1] = std_RhoW / bw;

    // Oil Properties
    PVDO.Eval_All(0, P, data, cdata);
    OCP_DBL bo  = data[1];
    OCP_DBL bop = cdata[1];

    mu[0]  = data[2];
    xi[0]  = 1 / (CONV1 * bo);
    rho[0] = std_RhoO / bo;

    vj[0]  = CONV1 * Ni[0] * bo;
    vj[1]  = CONV1 * Ni[1] * bw;
    vf     = vj[0] + vj[1];
    S[0]   = vj[0] / vf;
    S[1]   = vj[1] / vf;
    vfi[0] = CONV1 * bo;
    vfi[1] = CONV1 * bw;
    vfP    = CONV1 * Ni[0] * bop + CONV1 * Ni[1] * bwp;
}

void BOMixture_OW::FlashFIM(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Niin,
                            const OCP_DBL* Sjin,
                            const USI&     lastNP,
                            const OCP_DBL* xijin,
                            const OCP_USI& bId)
{
    fill(dXsdXp.begin(), dXsdXp.end(), 0.0);
    fill(pSderExist.begin(), pSderExist.end(), OCP_TRUE);
    fill(pVnumCom.begin(), pVnumCom.end(), 0);

    P     = Pin;
    Ni[0] = Niin[0];
    Ni[1] = Niin[1];
    Nt    = Ni[0] + Ni[1];

    // Water Properties
    PVTW.Eval_All(0, P, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = bw0 * (1 - cbw * (P - Pw0));
    OCP_DBL bwp = -cbw * bw0;

    mu[1]  = data[3];
    xi[1]  = 1 / (bw * CONV1);
    rho[1] = std_RhoW / bw;

    muP[1]  = cdata[3];
    xiP[1]  = -bwp / (bw * bw * CONV1);
    rhoP[1] = CONV1 * xiP[1] * std_RhoW;

    // Oil Properties
    PVDO.Eval_All(0, P, data, cdata);
    OCP_DBL bo  = data[1];
    OCP_DBL bop = cdata[1];

    mu[0]  = data[2];
    xi[0]  = 1 / (CONV1 * bo);
    rho[0] = std_RhoO / bo;

    muP[0]  = cdata[2];
    xiP[0]  = -xi[0] * bop / bo;
    rhoP[0] = -rho[0] * bop / bo;

    vj[0]  = CONV1 * Ni[0] * bo;
    vj[1]  = CONV1 * Ni[1] * bw;
    vf     = vj[0] + vj[1];
    S[0]   = vj[0] / vf;
    S[1]   = vj[1] / vf;
    vfi[0] = CONV1 * bo;
    vfi[1] = CONV1 * bw;
    vfP    = CONV1 * Ni[0] * bop + CONV1 * Ni[1] * bwp;

    dXsdXp[0] = (CONV1 * Ni[0] * bop - S[0] * vfP) / vf; // dSo / dP
    dXsdXp[1] = (CONV1 * bo - S[0] * vfi[0]) / vf;       // dSo / dNo
    dXsdXp[2] = -S[0] * vfi[1] / vf;                     // dSo / dNw

    dXsdXp[3] = (CONV1 * Ni[1] * bwp - S[1] * vfP) / vf; // dSw / dP
    dXsdXp[4] = -S[1] * vfi[0] / vf;                     // dSw / dNo
    dXsdXp[5] = (CONV1 * bw - S[1] * vfi[1]) / vf;       // dSw / dNw
}

OCP_DBL BOMixture_OW::XiPhase(const OCP_DBL& Pin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Ziin,
                              const USI&     tarPhase)
{
    if (tarPhase == WATER) {
        // inj fluid is water
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0 = data[0];
        OCP_DBL bw0 = data[1];
        OCP_DBL cbw = data[2];
        OCP_DBL bw  = bw0 * (1 - cbw * (Pin - Pw0));
        OCP_DBL xiw = 1 / (CONV1 * bw);
        return xiw;
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

OCP_DBL BOMixture_OW::RhoPhase(const OCP_DBL& Pin,
                               const OCP_DBL& Pbb,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Ziin,
                               const USI&     tarPhase)
{
    if (tarPhase == OIL) {
        OCP_DBL bo = PVDO.Eval(0, Pin, 1);
        return std_RhoO / bo;
    } else if (tarPhase == WATER) {
        PVTW.Eval_All(0, Pin, data, cdata);
        OCP_DBL Pw0 = data[0];
        OCP_DBL bw0 = data[1];
        OCP_DBL cbw = data[2];
        OCP_DBL bw  = bw0 * (1 - cbw * (Pin - Pw0));
        return std_RhoW / bw;
    } else {
        OCP_ABORT("Wrong tarPhase!");
    }
}

void BOMixture_OW::SetupWellOpt(WellOpt&                  opt,
                                const vector<SolventINJ>& sols,
                                const OCP_DBL&            Psurf,
                                const OCP_DBL&            Tsurf)
{
    const USI wellType = opt.WellType();
    if (wellType == INJ) {
        const string fluidName = opt.InjFluidType();
        opt.SetInjFactor(1.0);

        if (fluidName == "WAT") {
            vector<OCP_DBL> tmpZi({0, 1});
            opt.SetInjZi(tmpZi);
            opt.SetInjProdPhase(WATER);
        } else {
            OCP_ABORT("WRONG Injecting Fluid!");
        }
    } else if (wellType == PROD) {
        vector<OCP_DBL> tmpWght(2, 0);
        switch (opt.OptMode()) {
            case BHP_MODE:
                break;
            case ORATE_MODE:
                tmpWght[0] = 1;
                break;
            case WRATE_MODE:
                tmpWght[1] = 1;
                break;
            case LRATE_MODE:
                tmpWght[0] = tmpWght[1] = 1;
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
/*----------------------------------------------------------------------------*/