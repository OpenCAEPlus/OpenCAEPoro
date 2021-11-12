/*! \file    FlowUnit.cpp
 *  \brief   FlowUnit class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "FlowUnit.hpp"
#include "UtilError.hpp"

/// Generate the table of Sw vs. Pcwg using SGOF and SWOF.
void FlowUnit::Generate_SWPCWG()
{
    if (SGOF.IsEmpty()) OCP_ABORT("SGOF is missing!");
    if (SWOF.IsEmpty()) OCP_ABORT("SWOF is missing!");

    std::vector<OCP_DBL> Sw(SWOF.GetCol(0));
    std::vector<OCP_DBL> Pcw(SWOF.GetCol(3));
    USI                  n = Sw.size();
    for (USI i = 0; i < n; i++) {
        OCP_DBL Pcg = SGOF.Eval(0, 1 - Sw[i], 3);
        Pcw[i] += Pcg;
    }

    SWPCWG.PushCol(Sw);
    SWPCWG.PushCol(Pcw);
    SWPCWG.SetRowCol();
}

/// TODO: Add Doxygen
void FlowUnit::CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
    switch (mode) {
        case PHASE_W:
            CalKrPc_W(kr_out, pc_out);
            break;
        case PHASE_OW:
            CalKrPc_OW(S_in, kr_out, pc_out);
            break;
        case PHASE_OG:
            CalKrPc_OG(S_in, kr_out, pc_out);
            break;
        case PHASE_ODGW: // should SWOF be used here?
            CalKrPc_ODGW(S_in, kr_out, pc_out);
            break;
        default:
            OCP_ABORT("Wrong mode specified!");
    }
}

/// TODO: Add Doxygen
void FlowUnit::CalKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out)
{
    kr_out[0] = 1;
    pc_out[0] = 0;
}

/// TODO: Add Doxygen
void FlowUnit::CalKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
    OCP_DBL Sw = S_in[1];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw = data[1];
    OCP_DBL kro = data[2];
    OCP_DBL Pcw = -data[3];

    kr_out[0] = kro;
    kr_out[1] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcw;
}

/// TODO: Add Doxygen
void FlowUnit::CalKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
    OCP_DBL Sg = S_in[1];

    // three phase black oil model using stone 2
    SGOF.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg = data[1];
    OCP_DBL kro = data[2];
    OCP_DBL Pcg = data[3];

    kr_out[0] = kro;
    kr_out[1] = krg;
    pc_out[0] = 0;
    pc_out[1] = Pcg;
}

/// TODO: Add Doxygen
void FlowUnit::CalKrPc_ODGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
    OCP_DBL Sg = S_in[1];
    OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw  = data[1];
    OCP_DBL krow = data[2];
    OCP_DBL Pcw  = -data[3];

    SGOF.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg  = data[1];
    OCP_DBL krog = data[2];
    OCP_DBL Pcg  = data[3];

    OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);

    kr_out[0] = kro; kr_out[1] = krg; kr_out[2] = krw;
    pc_out[0] = 0;   pc_out[1] = Pcg; pc_out[2] = Pcw;
}

/// TODO: Add Doxygen
OCP_DBL FlowUnit::CalKro_Stone2(const OCP_DBL& krow, const OCP_DBL& krog,
                                const OCP_DBL& krw, const OCP_DBL& krg) const
{
    // krog : oil relative permeability for a system with oil, gas and connate water
    // krow : oil relative permeability for a system with oil and water only

    OCP_DBL kro =
        kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));
    if (kro < 0) kro = 0;

    return kro;
}

void FlowUnit::CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out, OCP_DBL* dkrdS, OCP_DBL* dPcjdS)
{
    switch (mode)
    {
    case PHASE_OW:
        CalKrPcDeriv_OW(S_in, kr_out, pc_out, dkrdS, dPcjdS);
        break;
    case PHASE_ODGW:
        CalKrPcDeriv_ODGW(S_in, kr_out, pc_out, dkrdS, dPcjdS);
        break;
    default:
        break;
    }
}


void FlowUnit::CalKrPcDeriv_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out, OCP_DBL* dkrdS, OCP_DBL* dPcjdS)
{
    OCP_DBL Sw = S_in[1];
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw = data[1];   OCP_DBL dKrwdSw = cdata[1];
    OCP_DBL krow = data[2];  OCP_DBL dKrowdSw = cdata[2];
    OCP_DBL Pcw = -data[3];  OCP_DBL dPcwdSw = -cdata[3];

    kr_out[0] = krow; kr_out[1] = krw;
    pc_out[0] = 0;    pc_out[1] = Pcw;

    dkrdS[0] = 0; dkrdS[1] = dKrowdSw;
    dkrdS[2] = 0; dkrdS[3] = dKrwdSw;

    dPcjdS[0] = 0; dPcjdS[1] = 0;
    dPcjdS[2] = 0; dPcjdS[3] = dPcwdSw;

}


void FlowUnit::CalKrPcDeriv_ODGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out, OCP_DBL* dkrdS, OCP_DBL* dPcjdS)
{
    OCP_DBL Sg = S_in[1];
    OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw = data[1];   OCP_DBL dKrwdSw = cdata[1];
    OCP_DBL krow = data[2];  OCP_DBL dKrowdSw = cdata[2];
    OCP_DBL Pcw = -data[3];  OCP_DBL dPcwdSw = -cdata[3];

    SGOF.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg = data[1];   OCP_DBL dKrgdSg = cdata[1];
    OCP_DBL krog = data[2];  OCP_DBL dKrogdSg = cdata[2];
    OCP_DBL Pcg = data[3];   OCP_DBL dPcgdSg = cdata[3];

    OCP_DBL dKrodSg, dKrodSw;

    OCP_DBL kro = kro_stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSw, dKrgdSg, dKrogdSg, dKrodSw, dKrodSg);
    //if (kro < 0) {
    //    cout << S_in[0] << "   " << S_in[1] << "   " << S_in[2] << endl;
    //    kro = 0;
    //}

    kr_out[0] = kro; kr_out[1] = krg; kr_out[2] = krw;
    pc_out[0] = 0;   pc_out[1] = Pcg; pc_out[2] = Pcw;

    dkrdS[0] = 0;    dkrdS[1] = dKrodSg;    dkrdS[2] = dKrodSw;
    dkrdS[3] = 0;    dkrdS[4] = dKrgdSg;    dkrdS[5] = 0;
    dkrdS[6] = 0;    dkrdS[7] = 0;          dkrdS[8] = dKrwdSw;

    dPcjdS[0] = 0;   dPcjdS[1] = 0;         dPcjdS[2] = 0;
    dPcjdS[3] = 0;   dPcjdS[4] = dPcgdSg;   dPcjdS[5] = 0;
    dPcjdS[6] = 0;   dPcjdS[7] = 0;         dPcjdS[8] = dPcwdSw;
}


OCP_DBL FlowUnit::kro_stone2Der(OCP_DBL krow, OCP_DBL krog,
    OCP_DBL krw, OCP_DBL krg,
    OCP_DBL dkrwdSw, OCP_DBL dkrowdSw,
    OCP_DBL dkrgdSg, OCP_DBL dkrogdSg,
    OCP_DBL& out_dkrodSw, OCP_DBL& out_dkrodSg) {
    OCP_DBL kro, dkrodSw, dkrodSg;
    kro = kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));

    dkrodSw = kroMax * ((dkrowdSw / kroMax + dkrwdSw) * (krog / kroMax + krg)
        - (dkrwdSw));
    dkrodSg = kroMax * ((krow / kroMax + krw) * (dkrogdSg / kroMax + dkrgdSg) - (dkrgdSg));

    if (kro < 0) {
        kro = 0;
        dkrodSg = 0;
        dkrodSw = 0;
    }
    out_dkrodSw = dkrodSw;
    out_dkrodSg = dkrodSg;
    return kro;
}


/// TODO: Add Doxygen
FlowUnit::FlowUnit(const ParamReservoir& rs_param, const USI& inmode, const USI& i)
{
    mode = inmode;
    if (rs_param.water) {
        SWOF.Setup(rs_param.SWOF_T.data[i]);
        len = max(len, SWOF.GetCol());
    }
    if (rs_param.gas) {
        SGOF.Setup(rs_param.SGOF_T.data[i]);
        len = max(len, SWOF.GetCol());
    }

    data.resize(len, 0);
    cdata.resize(len, 0);

    kroMax = 0;
    if (rs_param.water) {
        kroMax = SWOF.GetCol(2)[0];
    } else if (rs_param.gas) {
        kroMax = SGOF.GetCol(2)[0];
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