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

///////////////////////////////////////////////
// FlowUnit_W
///////////////////////////////////////////////

void FlowUnit_W::CalKrPc(const OCP_DBL* S_in,
                         OCP_DBL*       kr_out,
                         OCP_DBL*       pc_out,
                         const OCP_USI& bId)
{
    kr_out[0] = 1;
    pc_out[0] = 0;
}

void FlowUnit_W::CalKrPcDeriv(const OCP_DBL* S_in,
                              OCP_DBL*       kr_out,
                              OCP_DBL*       pc_out,
                              OCP_DBL*       dkrdS,
                              OCP_DBL*       dPcjdS,
                              const OCP_USI& bId)
{
    kr_out[0] = 1;
    pc_out[0] = 0;
    dkrdS[0]  = 0;
    dPcjdS[0] = 0;
}

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

FlowUnit_OW::FlowUnit_OW(const ParamReservoir& rs_param, const USI& i)
{
    SWOF.Setup(rs_param.SWOF_T.data[i]);
    Swco = SWOF.GetCol(0)[0];

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OW::CalKrPc(const OCP_DBL* S_in,
                          OCP_DBL*       kr_out,
                          OCP_DBL*       pc_out,
                          const OCP_USI& bId)
{
    OCP_DBL Sw = S_in[1];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw  = data[1];
    OCP_DBL kro  = data[2];
    OCP_DBL Pcwo = -data[3];

    kr_out[0] = kro;
    kr_out[1] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcwo;
}

void FlowUnit_OW::CalKrPcDeriv(const OCP_DBL* S_in,
                               OCP_DBL*       kr_out,
                               OCP_DBL*       pc_out,
                               OCP_DBL*       dkrdS,
                               OCP_DBL*       dPcjdS,
                               const OCP_USI& bId)
{
    OCP_DBL Sw = S_in[1];
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw      = data[1];
    OCP_DBL dKrwdSw  = cdata[1];
    OCP_DBL krow     = data[2];
    OCP_DBL dKrowdSw = cdata[2];
    OCP_DBL Pcwo     = -data[3];
    OCP_DBL dPcwdSw  = -cdata[3];

    kr_out[0] = krow;
    kr_out[1] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcwo;

    dkrdS[0] = 0;
    dkrdS[1] = dKrowdSw;
    dkrdS[2] = 0;
    dkrdS[3] = dKrwdSw;

    dPcjdS[0] = 0;
    dPcjdS[1] = 0;
    dPcjdS[2] = 0;
    dPcjdS[3] = dPcwdSw;
}

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

FlowUnit_OG::FlowUnit_OG(const ParamReservoir& rs_param, const USI& i)
{
    SGOF.Setup(rs_param.SGOF_T.data[i]);

    kroMax = SGOF.GetCol(2)[0];

    data.resize(4, 0);
    cdata.resize(4, 0);
}

void FlowUnit_OG::CalKrPc(const OCP_DBL* S_in,
                          OCP_DBL*       kr_out,
                          OCP_DBL*       pc_out,
                          const OCP_USI& bId)
{
    OCP_DBL Sg = S_in[1];

    // three phase black oil model using stone 2
    SGOF.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg  = data[1];
    OCP_DBL kro  = data[2];
    OCP_DBL Pcgo = data[3];

    kr_out[0] = kro;
    kr_out[1] = krg;
    pc_out[0] = 0;
    pc_out[1] = Pcgo;
}

void FlowUnit_OG::CalKrPcDeriv(const OCP_DBL* S_in,
                               OCP_DBL*       kr_out,
                               OCP_DBL*       pc_out,
                               OCP_DBL*       dkrdS,
                               OCP_DBL*       dPcjdS,
                               const OCP_USI& bId)
{
    OCP_ABORT("Not Completed Now!");
}

///////////////////////////////////////////////
// FlowUnit_ODGW
///////////////////////////////////////////////

OCP_DBL FlowUnit_ODGW::CalKro_Stone2(const OCP_DBL& krow,
                                     const OCP_DBL& krog,
                                     const OCP_DBL& krw,
                                     const OCP_DBL& krg) const
{
    // krog : oil relative permeability for a system with oil, gas and connate water
    // krow : oil relative permeability for a system with oil and water only

    OCP_DBL kro =
        kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));
    if (kro < 0) kro = 0;

    return kro;
}

OCP_DBL FlowUnit_ODGW::CalKro_Default(const OCP_DBL& Sg,
                                      const OCP_DBL& Sw,
                                      const OCP_DBL& krog,
                                      const OCP_DBL& krow) const
{
    OCP_DBL tmp = Sg + Sw - Swco;
    if (tmp <= TINY) {
        return kroMax;
    }
    OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    return kro;
}

///////////////////////////////////////////////
// FlowUnit_ODGW01
///////////////////////////////////////////////

FlowUnit_ODGW01::FlowUnit_ODGW01(const ParamReservoir& rs_param, const USI& i)
{
    SWOF.Setup(rs_param.SWOF_T.data[i]);
    SGOF.Setup(rs_param.SGOF_T.data[i]);

    kroMax = SWOF.GetCol(2)[0];
    Swco   = SWOF.GetCol(0)[0];

    data.resize(4, 0);
    cdata.resize(4, 0);

    Generate_SWPCWG();

    Scm.resize(3, 0);
    // oil, set to zero now
    OCP_INT tmprow;
    tmprow = SGOF.GetRowZero(1);
    if (tmprow >= 0) Scm[1] = SGOF.GetCol(0)[tmprow];
    tmprow = SWOF.GetRowZero(1);
    if (tmprow >= 0) Scm[2] = SWOF.GetCol(0)[tmprow];
}

void FlowUnit_ODGW01::CalKrPc(const OCP_DBL* S_in,
                              OCP_DBL*       kr_out,
                              OCP_DBL*       pc_out,
                              const OCP_USI& bId)
{
    const OCP_DBL Sg = S_in[1];
    const OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    const OCP_DBL krw  = data[1];
    const OCP_DBL krow = data[2];
    const OCP_DBL Pcwo = -data[3];

    SGOF.Eval_All(0, Sg, data, cdata);
    const OCP_DBL krg  = data[1];
    const OCP_DBL krog = data[2];
    const OCP_DBL Pcgo = data[3];

    const OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
    // OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);

    kr_out[0] = kro;
    kr_out[1] = krg;
    kr_out[2] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcgo;
    pc_out[2] = Pcwo;
}

void FlowUnit_ODGW01::CalKrPcDeriv(const OCP_DBL* S_in,
                                   OCP_DBL*       kr_out,
                                   OCP_DBL*       pc_out,
                                   OCP_DBL*       dkrdS,
                                   OCP_DBL*       dPcjdS,
                                   const OCP_USI& bId)
{
    OCP_DBL Sg = S_in[1];
    OCP_DBL Sw = S_in[2];

    // three phase black oil model using stone 2
    SWOF.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw      = data[1];
    OCP_DBL dKrwdSw  = cdata[1];
    OCP_DBL krow     = data[2];
    OCP_DBL dKrowdSw = cdata[2];
    OCP_DBL Pcwo     = -data[3];
    OCP_DBL dPcwdSw  = -cdata[3];

    SGOF.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg      = data[1];
    OCP_DBL dKrgdSg  = cdata[1];
    OCP_DBL krog     = data[2];
    OCP_DBL dKrogdSg = cdata[2];
    OCP_DBL Pcgo     = data[3];
    OCP_DBL dPcgdSg  = cdata[3];

    OCP_DBL dKrodSg{0}, dKrodSw{0}, kro{0};

    kro = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSw, dKrgdSg, dKrogdSg,
                           dKrodSw, dKrodSg);
    // if (kro < 0) {
    //    cout << S_in[0] << "   " << S_in[1] << "   " << S_in[2] << endl;
    //    kro = 0;
    //}
    // kro = CalKro_DefaultDer(Sg, Sw, krog, krow, dKrogdSg, dKrowdSw, dKrodSg,
    // dKrodSw);

    kr_out[0] = kro;
    kr_out[1] = krg;
    kr_out[2] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcgo;
    pc_out[2] = Pcwo;

    dkrdS[0] = 0;
    dkrdS[1] = dKrodSg;
    dkrdS[2] = dKrodSw;
    dkrdS[3] = 0;
    dkrdS[4] = dKrgdSg;
    dkrdS[5] = 0;
    dkrdS[6] = 0;
    dkrdS[7] = 0;
    dkrdS[8] = dKrwdSw;

    dPcjdS[0] = 0;
    dPcjdS[1] = 0;
    dPcjdS[2] = 0;
    dPcjdS[3] = 0;
    dPcjdS[4] = dPcgdSg;
    dPcjdS[5] = 0;
    dPcjdS[6] = 0;
    dPcjdS[7] = 0;
    dPcjdS[8] = dPcwdSw;
}

OCP_DBL FlowUnit_ODGW01::CalKro_Stone2Der(OCP_DBL  krow,
                                          OCP_DBL  krog,
                                          OCP_DBL  krw,
                                          OCP_DBL  krg,
                                          OCP_DBL  dkrwdSw,
                                          OCP_DBL  dkrowdSw,
                                          OCP_DBL  dkrgdSg,
                                          OCP_DBL  dkrogdSg,
                                          OCP_DBL& out_dkrodSw,
                                          OCP_DBL& out_dkrodSg) const
{
    OCP_DBL kro, dkrodSw, dkrodSg;
    kro = kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));

    dkrodSw =
        kroMax * ((dkrowdSw / kroMax + dkrwdSw) * (krog / kroMax + krg) - (dkrwdSw));
    dkrodSg =
        kroMax * ((krow / kroMax + krw) * (dkrogdSg / kroMax + dkrgdSg) - (dkrgdSg));

    if (kro < 0) {
        kro     = 0;
        dkrodSg = 0;
        dkrodSw = 0;
    }
    out_dkrodSw = dkrodSw;
    out_dkrodSg = dkrodSg;
    return kro;
}

OCP_DBL FlowUnit_ODGW01::CalKro_DefaultDer(const OCP_DBL& Sg,
                                           const OCP_DBL& Sw,
                                           const OCP_DBL& krog,
                                           const OCP_DBL& krow,
                                           const OCP_DBL& dkrogSg,
                                           const OCP_DBL& dkrowSw,
                                           OCP_DBL&       dkroSg,
                                           OCP_DBL&       dkroSw) const
{
    OCP_DBL tmp = Sg + Sw - Swco;
    if (tmp <= TINY) {
        dkroSg = 0;
        dkroSw = 0;
        return kroMax;
    }
    OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    dkroSg      = (krog + Sg * dkrogSg - kro) / tmp;
    dkroSw      = (krow + (Sw - Swco) * dkrowSw - kro) / tmp;
    return kro;
}

void FlowUnit_ODGW01::Generate_SWPCWG()
{
    if (SGOF.IsEmpty()) OCP_ABORT("SGOF is missing!");
    if (SWOF.IsEmpty()) OCP_ABORT("SWOF is missing!");

    std::vector<OCP_DBL> Sw(SWOF.GetCol(0));
    std::vector<OCP_DBL> Pcow(SWOF.GetCol(3));
    USI                  n = Sw.size();
    for (USI i = 0; i < n; i++) {
        OCP_DBL Pcgo = SGOF.Eval(0, 1 - Sw[i], 3);
        Pcow[i] += Pcgo; // Pcgw
    }

    SWPCGW.PushCol(Sw);
    SWPCGW.PushCol(Pcow);
    SWPCGW.SetRowCol();
}

///////////////////////////////////////////////
// FlowUnit_ODGW01_Miscible
///////////////////////////////////////////////

void FlowUnit_ODGW01_Miscible::SetupScale(const OCP_USI& bId,
                                          OCP_DBL&       Swin,
                                          const OCP_DBL& Pcowin)
{
    if (scaleTerm->IfScale()) {
        const OCP_DBL SwInit = scaleTerm->GetSwInit(bId);
        if (SwInit <= Swco) {
            Swin = Swco;
        } else {
            Swin = SwInit;
            if (Pcowin > 0) {
                const OCP_DBL PcowInit = GetPcowBySw(Swin);
                if (PcowInit > 0) {
                    scaleTerm->AssignScaleValue(
                        bId,
                        (Pcowin / PcowInit * maxPcow - minPcow) / (maxPcow - minPcow));
                }
            }
        }
    }
};

void FlowUnit_ODGW01_Miscible::CalKrPc(const OCP_DBL* S_in,
                                       OCP_DBL*       kr_out,
                                       OCP_DBL*       pc_out,
                                       const OCP_USI& bId)
{
    if (!misTerm->CalFkFp(bId, Fk, Fp)) {
        FlowUnit_ODGW01::CalKrPc(S_in, kr_out, pc_out, bId);
    } else {
        OCP_DBL So = S_in[0];
        OCP_DBL Sg = S_in[1];
        OCP_DBL Sw = S_in[2];

        SWOF.Eval_All(0, Sw, data, cdata);
        const OCP_DBL krw  = data[1];
        OCP_DBL       krow = data[2];
        const OCP_DBL Pcwo = -data[3];

        SGOF.Eval_All(0, Sg, data, cdata);
        OCP_DBL       krg  = data[1];
        OCP_DBL       krog = data[2];
        const OCP_DBL Pcgo = data[3] * Fp;

        OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
        // OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);
        OCP_DBL krgt = SGOF.Eval(0, (1 - Sw), 1);
        OCP_DBL krh  = 0.5 * (krow + krgt);

        // from CMG, see *SIGMA
        kro = Fk * kro + (1 - Fk) * krh * So / (1 - Sw);
        krg = Fk * krg + (1 - Fk) * krh * Sg / (1 - Sw);

        kr_out[0] = kro;
        kr_out[1] = krg;
        kr_out[2] = krw;
        pc_out[0] = 0;
        pc_out[1] = Pcgo;
        pc_out[2] = Pcwo;
    }

    if (scaleTerm->IfScale()) {
        pc_out[2] = -((-pc_out[2] - minPcow) * scaleTerm->GetScaleVal(bId) + minPcow);
    }
}

void FlowUnit_ODGW01_Miscible::CalKrPcDeriv(const OCP_DBL* S_in,
                                            OCP_DBL*       kr_out,
                                            OCP_DBL*       pc_out,
                                            OCP_DBL*       dkrdS,
                                            OCP_DBL*       dPcjdS,
                                            const OCP_USI& bId)
{
    if (!misTerm->CalFkFp(bId, Fk, Fp)) {
        FlowUnit_ODGW01::CalKrPcDeriv(S_in, kr_out, pc_out, dkrdS, dPcjdS, bId);
    } else {
        const OCP_DBL So = S_in[0];
        const OCP_DBL Sg = S_in[1];
        const OCP_DBL Sw = S_in[2];

        SWOF.Eval_All(0, Sw, data, cdata);
        const OCP_DBL krw      = data[1];
        const OCP_DBL dKrwdSw  = cdata[1];
        const OCP_DBL krow     = data[2];
        const OCP_DBL dKrowdSw = cdata[2];
        const OCP_DBL Pcwo     = -data[3];
        const OCP_DBL dPcwdSw  = -cdata[3];

        SGOF.Eval_All(0, Sg, data, cdata);
        OCP_DBL       krg      = data[1];
        const OCP_DBL dKrgdSg  = cdata[1];
        const OCP_DBL krog     = data[2];
        const OCP_DBL dKrogdSg = cdata[2];
        const OCP_DBL Pcgo     = data[3] * Fp;
        const OCP_DBL dPcgdSg  = cdata[3] * Fp;

        OCP_DBL dKrodSg{0}, dKrodSw{0}, kro{0};
        kro = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSw, dKrgdSg,
                               dKrogdSg, dKrodSw, dKrodSg);

        OCP_DBL       dkrgd1_Sw = 0;
        const OCP_DBL krgt      = SGOF.Eval(0, (1 - Sw), 1, dkrgd1_Sw);
        const OCP_DBL krh       = 0.5 * (krow + krgt);

        // from CMG, see *SIGMA
        kro = Fk * kro + (1 - Fk) * krh * So / (1 - Sw);
        krg = Fk * krg + (1 - Fk) * krh * Sg / (1 - Sw);

        kr_out[0] = kro;
        kr_out[1] = krg;
        kr_out[2] = krw;
        pc_out[0] = 0;
        pc_out[1] = Pcgo;
        pc_out[2] = Pcwo;

        OCP_DBL dkrhdSw = 0.5 * (dKrowdSw - dkrgd1_Sw);
        OCP_DBL temp    = (1 - Fk) / (1 - Sw) * (krh / (1 - Sw) + dkrhdSw);

        dkrdS[0] = (1 - Fk) * krh / (1 - Sw);
        dkrdS[1] = Fk * dKrodSg;
        dkrdS[2] = Fk * dKrodSw + temp * So;
        dkrdS[3] = 0;
        dkrdS[4] = Fk * dKrgdSg + (1 - Fk) * krh / (1 - Sw);
        dkrdS[5] = temp * Sg;
        dkrdS[6] = 0;
        dkrdS[7] = 0;
        dkrdS[8] = dKrwdSw;

        dPcjdS[0] = 0;
        dPcjdS[1] = 0;
        dPcjdS[2] = 0;
        dPcjdS[3] = 0;
        dPcjdS[4] = dPcgdSg;
        dPcjdS[5] = 0;
        dPcjdS[6] = 0;
        dPcjdS[7] = 0;
        dPcjdS[8] = dPcwdSw;
    }

    if (scaleTerm->IfScale()) {
        pc_out[2] = -((-pc_out[2] - minPcow) * scaleTerm->GetScaleVal(bId) + minPcow);
        dPcjdS[8] *= scaleTerm->GetScaleVal(bId);
    }
}

///////////////////////////////////////////////
// FlowUnit_ODGW02
///////////////////////////////////////////////

FlowUnit_ODGW02::FlowUnit_ODGW02(const ParamReservoir& rs_param, const USI& i)
{
    SWFN.Setup(rs_param.SWFN_T.data[i]);
    SGFN.Setup(rs_param.SGFN_T.data[i]);
    SOF3.Setup(rs_param.SOF3_T.data[i]);

    kroMax = SOF3.GetCol(1).back();
    Swco   = SWFN.GetCol(0)[0];

    data.resize(3, 0);
    cdata.resize(3, 0);

    Generate_SWPCWG();
}

void FlowUnit_ODGW02::CalKrPc(const OCP_DBL* S_in,
                              OCP_DBL*       kr_out,
                              OCP_DBL*       pc_out,
                              const OCP_USI& bId)
{
    OCP_DBL So = S_in[0];
    OCP_DBL Sg = S_in[1];
    OCP_DBL Sw = S_in[2];

    SWFN.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw  = data[1];
    OCP_DBL Pcwo = -data[2];

    SGFN.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg  = data[1];
    OCP_DBL Pcgo = data[2];

    SOF3.Eval_All(0, So, data, cdata);
    OCP_DBL krow = data[1];
    OCP_DBL krog = data[2];

    // OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);
    OCP_DBL kro = CalKro_Default(Sg, Sw, krog, krow);

    kr_out[0] = kro;
    kr_out[1] = krg;
    kr_out[2] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcgo;
    pc_out[2] = Pcwo;
}

void FlowUnit_ODGW02::CalKrPcDeriv(const OCP_DBL* S_in,
                                   OCP_DBL*       kr_out,
                                   OCP_DBL*       pc_out,
                                   OCP_DBL*       dkrdS,
                                   OCP_DBL*       dPcjdS,
                                   const OCP_USI& bId)
{
    OCP_DBL So = S_in[0];
    OCP_DBL Sg = S_in[1];
    OCP_DBL Sw = S_in[2];

    SWFN.Eval_All(0, Sw, data, cdata);
    OCP_DBL krw      = data[1];
    OCP_DBL dKrwdSw  = cdata[1];
    OCP_DBL Pcwo     = -data[2];
    OCP_DBL dPcwodSw = -cdata[2];

    SGFN.Eval_All(0, Sg, data, cdata);
    OCP_DBL krg      = data[1];
    OCP_DBL dKrgdSg  = cdata[1];
    OCP_DBL Pcgo     = data[2];
    OCP_DBL dPcgodSg = cdata[2];

    SOF3.Eval_All(0, So, data, cdata);
    OCP_DBL krow     = data[1];
    OCP_DBL dKrowdSo = cdata[1];
    OCP_DBL krog     = data[2];
    OCP_DBL dKrogdSo = cdata[2];

    OCP_DBL dKroSo = 0;
    OCP_DBL kro    = CalKro_Stone2Der(krow, krog, krw, krg, dKrwdSw, dKrowdSo, dKrgdSg,
                                      dKrogdSo, dKroSo);
    // OCP_DBL kro = CalKro_DefaultDer(Sg, Sw, krog, krow, dKrogdSo, dKrowdSo, dKroSo);

    kr_out[0] = kro;
    kr_out[1] = krg;
    kr_out[2] = krw;
    pc_out[0] = 0;
    pc_out[1] = Pcgo;
    pc_out[2] = Pcwo;

    dkrdS[0] = dKroSo;
    dkrdS[1] = 0;
    dkrdS[2] = 0;
    dkrdS[3] = 0;
    dkrdS[4] = dKrgdSg;
    dkrdS[5] = 0;
    dkrdS[6] = 0;
    dkrdS[7] = 0;
    dkrdS[8] = dKrwdSw;

    dPcjdS[0] = 0;
    dPcjdS[1] = 0;
    dPcjdS[2] = 0;
    dPcjdS[3] = 0;
    dPcjdS[4] = dPcgodSg;
    dPcjdS[5] = 0;
    dPcjdS[6] = 0;
    dPcjdS[7] = 0;
    dPcjdS[8] = dPcwodSw;
}

OCP_DBL FlowUnit_ODGW02::CalKro_Stone2Der(OCP_DBL  krow,
                                          OCP_DBL  krog,
                                          OCP_DBL  krw,
                                          OCP_DBL  krg,
                                          OCP_DBL  dkrwdSw,
                                          OCP_DBL  dkrowdSo,
                                          OCP_DBL  dkrgdSg,
                                          OCP_DBL  dkrogdSo,
                                          OCP_DBL& out_dkrodSo) const
{
    OCP_DBL kro, dKrodSo;

    kro     = kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));
    dKrodSo = dkrowdSo * (krog / kroMax + krg) + dkrogdSo * (krow / kroMax + krw);

    if (kro < 0) {
        kro     = 0;
        dKrodSo = 0;
    }
    out_dkrodSo = dKrodSo;
    return kro;
}

OCP_DBL FlowUnit_ODGW02::CalKro_DefaultDer(const OCP_DBL& Sg,
                                           const OCP_DBL& Sw,
                                           const OCP_DBL& krog,
                                           const OCP_DBL& krow,
                                           const OCP_DBL& dkrogdSo,
                                           const OCP_DBL& dkrowdSo,
                                           OCP_DBL&       out_dkrodSo) const
{
    OCP_DBL tmp = Sg + Sw - Swco;
    OCP_DBL kro = (Sg * krog + (Sw - Swco) * krow) / tmp;
    out_dkrodSo = (Sg * dkrogdSo + (Sw - Swco) * dkrowdSo) / tmp;

    if (tmp <= TINY) {
        kro         = kroMax;
        out_dkrodSo = 0;
    }
    return kro;
}

void FlowUnit_ODGW02::Generate_SWPCWG()
{
    if (SWFN.IsEmpty()) OCP_ABORT("SWFN is missing!");
    if (SGFN.IsEmpty()) OCP_ABORT("SGFN is missing!");

    std::vector<OCP_DBL> Sw(SWFN.GetCol(0));
    std::vector<OCP_DBL> Pcow(SWFN.GetCol(2));
    USI                  n = Sw.size();
    for (USI i = 0; i < n; i++) {
        OCP_DBL Pcgo = SGFN.Eval(0, 1 - Sw[i], 2);
        Pcow[i] += Pcgo; // pcgw
    }

    SWPCGW.PushCol(Sw);
    SWPCGW.PushCol(Pcow);
    SWPCGW.SetRowCol();
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/