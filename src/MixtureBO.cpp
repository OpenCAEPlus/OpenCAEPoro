/*! \file    MixtureBO.cpp
 *  \brief   MixtureBO class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "MixtureBO.hpp"

BOMixture::BOMixture(const ParamReservoir& rs_param, const USI& PVTmode, const USI& i)
{
    mixtureType = BLKOIL;
    mode        = PVTmode;
    numPhase    = rs_param.numPhase;
    numCom      = rs_param.numCom;

    Ni.resize(numCom);
    phaseExist.resize(numPhase);
    v.resize(numPhase);
    S.resize(numPhase);
    xi.resize(numPhase);
    cij.resize(numPhase * numCom);
    rho.resize(numPhase);
    mu.resize(numPhase);
    vfi.resize(numCom);
    // Derivatives for FIM
    rhoP.resize(numPhase);
    xiP.resize(numPhase);
    muP.resize(numPhase);
    rhoC.resize(numPhase * numCom);
    xiC.resize(numPhase * numCom);
    muC.resize(numPhase * numCom);
    dSec_dPri.resize((numCom + 1) * (numPhase + numPhase * numCom));


    if (rs_param.density.activity) {
        std_RhoO = rs_param.density.data[0];
        std_RhoW = rs_param.density.data[1];
        std_RhoG = rs_param.density.data[2];
    } else {
        std_RhoO = (141.5 * RHOW_STD) / (rs_param.gravity.data[0] + 131.5);
        std_RhoW = RHOW_STD * rs_param.gravity.data[1];
        std_RhoG = RHOAIR_STD * rs_param.gravity.data[2];
    }

    std_GammaO = GRAVITY_FACTOR * std_RhoO;
    std_GammaG = GRAVITY_FACTOR * std_RhoG;
    std_GammaW = GRAVITY_FACTOR * std_RhoW;

    if (rs_param.PVTW_T.data.size() != 0) {
        PVTW.Setup(rs_param.PVTW_T.data[i]);
        len = max(len, PVTW.GetCol());
    }
    if (rs_param.PVCO_T.data.size() != 0) {
        PVCO.Setup(rs_param.PVCO_T.data[i]);
        len = max(len, PVCO.GetCol());
    }
    if (rs_param.PVDO_T.data.size() != 0) {
        PVDO.Setup(rs_param.PVDO_T.data[i]);
        len = max(len, PVDO.GetCol());
    }
    if (rs_param.PVDG_T.data.size() != 0) {
        PVDG.Setup(rs_param.PVDG_T.data[i]);
        len = max(len, PVDG.GetCol());
    }
    data.resize(len, 0);
    cdata.resize(len, 0);
}

void BOMixture::Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
                         const OCP_DBL* Sjin, const OCP_DBL& Vpore, const OCP_DBL* Ziin)
{
    switch (mode) {
        case PHASE_W:
            BOFlash_Sj_W(Pin, Sjin, Vpore);
            break;

        case PHASE_OW:
            BOFlash_Sj_OW(Pin, Sjin, Vpore);
            break;

        case PHASE_ODGW:
            BOFlash_Sj_ODGW(Pin, Pbbin, Sjin, Vpore);
            break;

        default:
            OCP_ABORT("Wrong mode specified!");
    }
}

void BOMixture::Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
#ifdef _DEBUG
    // CheckNi(Niin);
#endif // _DEBUG

    switch (mode) {
        case PHASE_W:
            BOFlash_Ni_W(Pin, Niin);
            break;

        case PHASE_OW:
            BOFlash_Ni_OW(Pin, Niin);
            break;

        case PHASE_ODGW:
            BOFlash_Ni_ODGW(Pin, Niin);
            break;

        default:
            OCP_ABORT("Wrong mode specified!");
    }
}

void BOMixture::Flash_Ni_Deriv(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
#ifdef _DEBUG
    // CheckNi(Niin);
#endif // _DEBUG

    switch (mode) {
    case PHASE_W:
        BOFlash_Ni_W_Deriv(Pin, Niin);
        break;

    case PHASE_OW:
        BOFlash_Ni_OW_Deriv(Pin, Niin);
        break;

    case PHASE_ODGW:
        BOFlash_Ni_ODGW_Deriv(Pin, Niin);
        break;

    default:
        OCP_ABORT("Wrong mode specified!");
    }
}

void BOMixture::BOFlash_Ni_W_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin){}

void BOMixture::BOFlash_Ni_OW_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin) {}

void BOMixture::BOFlash_Ni_W(const OCP_DBL& Pin, const OCP_DBL* Niin) {}

void BOMixture::BOFlash_Ni_OW(const OCP_DBL& Pin, const OCP_DBL* Niin) {}

void BOMixture::BOFlash_Sj_W(const OCP_DBL& Pin, const OCP_DBL* Sjin,
                             const OCP_DBL& Vpore)
{
}

void BOMixture::BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin,
                              const OCP_DBL& Vpore)
{
}

OCP_DBL BOMixture::XiPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin)
{
    switch (mode) {
        case PHASE_W:
            OCP_ABORT("Will be added!");
        case PHASE_OW:
            OCP_ABORT("Will be added!");
        case PHASE_ODGW:
            return XiPhase_ODGW(Pin, Ziin);
            break;
        default:
            OCP_ABORT("Not implemented yet!");
    }

    return 0.0; // Should not reach here!
}

OCP_DBL BOMixture::RhoPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin)
{
    switch (mode) {
        case PHASE_W:
            OCP_ABORT("Will be added!");
        case PHASE_OW:
            OCP_ABORT("Will be added!");
        case PHASE_ODGW:
            return RhoPhase_ODGW(Pin, Ziin);
            break;
        default:
            OCP_ABORT("Not implemented yet!");
    }

    return 0.0; // Should not reach here!
}

OCP_DBL BOMixture::GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin)
{
    switch (mode) {
        case PHASE_OW:
            return GammaPhaseO_OW(Pin);
        case PHASE_ODGW:
            return GammaPhaseO_ODGW(Pin, Pbbin);
        default:
            OCP_ABORT("Wrong Mode");
    }
}

OCP_DBL BOMixture::GammaPhaseO_OW(const OCP_DBL& Pin) { return 0; }

OCP_DBL BOMixture::GammaPhaseW(const OCP_DBL& Pin)
{

    PVTW.Eval_All(0, Pin, data, cdata);
    OCP_DBL Pw0 = data[0];
    OCP_DBL bw0 = data[1];
    OCP_DBL cbw = data[2];
    OCP_DBL bw  = (bw0 * (1 - cbw * (Pin - Pw0)));

    return std_GammaW / bw;
}

OCP_DBL BOMixture::GammaPhaseG(const OCP_DBL& Pin)
{
    if (PVDG.IsEmpty()) {
        OCP_ABORT("PVDG is missing!");
    }
    OCP_DBL bg = (CONV1 / 1000) * PVDG.Eval(0, Pin, 1);

    return std_GammaG / bg;
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/