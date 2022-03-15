/*! \file    FlowUnit.hpp
 *  \brief   FlowUnit class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLOWUNIT_HEADER__
#define __FLOWUNIT_HEADER__

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "OCPTable.hpp"
#include "ParamReservoir.hpp"

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure woulbe be calculated here.
class FlowUnit
{
public:
    /// Default constructor.
    FlowUnit() = default;

    /// Pcow = Po - Pw
    virtual OCP_DBL GetPcowBySw(const OCP_DBL& sw) = 0;
    virtual OCP_DBL GetSwByPcow(const OCP_DBL& pcow) = 0;
    /// Pcgo = Pg - Po
    virtual OCP_DBL GetPcgoBySg(const OCP_DBL& sg) = 0;
    virtual OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) = 0;
    /// Pcgw = Pg - Pw
    virtual OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) = 0;

    /// calculate relative permeability and capillary pressure.
    virtual void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) = 0;
    virtual void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) = 0;

};

class FlowUnit_W : public FlowUnit
{
public:
    FlowUnit_W() = default;
    FlowUnit_W(const ParamReservoir& rs_param, const USI& i) {};

    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) override;
    void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) override;


    OCP_DBL GetPcowBySw(const OCP_DBL& sw)  override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow)  override { return 0; }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg)  override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo)  override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw)  override { return 0; }

};


class FlowUnit_OW : public FlowUnit
{
public:
    FlowUnit_OW() = default;
    FlowUnit_OW(const ParamReservoir& rs_param, const USI& i);

    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) override;
    void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) override;
 
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override {
        return SWOF.Eval(0, sw, 3);
    }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override {
        return SWOF.Eval_Inv(3, pcow, 0);
    }

    // useless
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg)  override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo)  override { return 0; }
    virtual OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw)  override { return 0; }

private:
    OCPTable SWOF;   ///< saturation table about water and oil.
    vector<OCP_DBL>
        data; ///< container used to store the results of values of interpolation.
    vector<OCP_DBL>
        cdata; ///< container used to store the results of slopes of interpolation.
    OCP_DBL kroMax;
    OCP_DBL Swco; ///< Saturation of connate water.
};

class FlowUnit_OG : public FlowUnit
{
public:
    FlowUnit_OG() = default;
    FlowUnit_OG(const ParamReservoir & rs_param, const USI & i);

    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) override;
    void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) override;
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg)  override {
        return SGOF.Eval(0, sg, 3);
    }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo)  override {
        return SGOF.Eval(3, pcgo, 0);
    }

    // uesless
    OCP_DBL GetPcowBySw(const OCP_DBL& sw)  override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow)  override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw)  override { return 0; }

private:
    OCPTable SGOF;   ///< saturation table about gas and oil.
    vector<OCP_DBL>
        data; ///< container used to store the results of values of interpolation.
    vector<OCP_DBL>
        cdata; ///< container used to store the results of slopes of interpolation.
    OCP_DBL kroMax;
};


class FlowUnit_ODGW : public FlowUnit
{
public:

    OCP_DBL CalKro_Stone2(const OCP_DBL& krow, const OCP_DBL& krog,
        const OCP_DBL& krw, const OCP_DBL& krg) const;
    OCP_DBL CalKro_Stone2Der(OCP_DBL krow, OCP_DBL krog, OCP_DBL krw, OCP_DBL krg,
        OCP_DBL dkrwdSw, OCP_DBL dkrowdSw, OCP_DBL dkrgdSg,
        OCP_DBL dkrogdSg, OCP_DBL& out_dkrodSw,
        OCP_DBL& out_dkrodSg) const;
    OCP_DBL CalKro_Default(const OCP_DBL& Sg, const OCP_DBL& Sw,
        const OCP_DBL& krog, const OCP_DBL& krow) const;
    OCP_DBL CalKro_DefaultDer(const OCP_DBL& Sg, const OCP_DBL& Sw,
        const OCP_DBL& krog, const OCP_DBL& krow,
        const OCP_DBL& dkrogSg, const OCP_DBL& dkrowSw,
        OCP_DBL& dkroSg, OCP_DBL& dkroSw) const;
    /// oil relative permeability in the presence of connate water only, used in stone2
    OCP_DBL kroMax; 
    OCP_DBL Swco; ///< Saturation of connate water.

};

class FlowUnit_ODGW01 : public FlowUnit_ODGW
{
public:
    FlowUnit_ODGW01() = default;
    FlowUnit_ODGW01(const ParamReservoir& rs_param, const USI& i);

    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) override;
    void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) override;
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override {
        return SWOF.Eval(0, sw, 3);
    }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override {
        return SWOF.Eval_Inv(3, pcow, 0);
    }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg)  override {
        return SGOF.Eval(0, sg, 3);
    }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo)  override {
        return SGOF.Eval(3, pcgo, 0);
    }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override {
        return SWPCWG.Eval_Inv(1, pcgw, 0);
    }
    void Generate_SWPCWG();

private:
    OCPTable SGOF;   ///< saturation table about gas and oil.
    OCPTable SWOF;   ///< saturation table about water and oil.
    OCPTable SWPCWG; ///< auxiliary table: saturation of water vs. capillary
                     ///< pressure between water and gas.
    vector<OCP_DBL>
        data; ///< container used to store the results of values of interpolation.
    vector<OCP_DBL>
        cdata; ///< container used to store the results of slopes of interpolation.
};

class FlowUnit_ODGW02 : public FlowUnit_ODGW
{
public:

    FlowUnit_ODGW02() = default;
    FlowUnit_ODGW02(const ParamReservoir& rs_param, const USI& i) {};

    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out) override;
    void CalKrPcDeriv(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out,
        OCP_DBL* dkrdS, OCP_DBL* dPcjdS) override;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw)  override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow)  override { return 0; }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg)  override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo)  override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw)  override { return 0; }

private:
    OCPTable SWFN;   ///< saturation table about water.
    OCPTable SGFN;   ///< saturation table about gas.
    OCPTable SOF3;   ///< saturation table about oil.
    vector<OCP_DBL>
        data; ///< container used to store the results of values of interpolation.
    vector<OCP_DBL>
        cdata; ///< container used to store the results of slopes of interpolation.
};


#endif /* end if __FLOWUNIT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/