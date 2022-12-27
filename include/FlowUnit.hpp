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

#include <math.h>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "OCPTable.hpp"
#include "ParamReservoir.hpp"
#include "OptionalFeatures.hpp"

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure will be calculated here.
class FlowUnit
{
public:
    /// Default constructor.
    FlowUnit() = default;

    /// Pcow = Po - Pw
    virtual OCP_DBL GetPcowBySw(const OCP_DBL& sw)   = 0;
    virtual OCP_DBL GetSwByPcow(const OCP_DBL& pcow) = 0;
    /// Pcgo = Pg - Po
    virtual OCP_DBL GetPcgoBySg(const OCP_DBL& sg)   = 0;
    virtual OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) = 0;
    /// Pcgw = Pg - Pw
    virtual OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) = 0;
    /// Return the value of Swco
    virtual OCP_DBL GetSwco() const = 0;
    /// Return the value of Scm
    virtual const vector<OCP_DBL>& GetScm() const = 0;

    /// Calculate relative permeability and capillary pressure.
    virtual void CalKrPc(const OCP_DBL* S_in,
                         OCP_DBL*       kr_out,
                         OCP_DBL*       pc_out,
                         const OCP_DBL& MySurTen,
                         OCP_DBL&       MyFk,
                         OCP_DBL&       MyFp) = 0;

    /// Calculate derivatives of relative permeability and capillary pressure.
    virtual void CalKrPcDeriv(const OCP_DBL* S_in,
                              OCP_DBL*       kr_out,
                              OCP_DBL*       pc_out,
                              OCP_DBL*       dkrdS,
                              OCP_DBL*       dPcjdS,
                              const OCP_DBL& MySurTen,
                              OCP_DBL&       MyFk,
                              OCP_DBL&       MyFp) = 0;
};

///////////////////////////////////////////////
// FlowUnit_W
///////////////////////////////////////////////

class FlowUnit_W : public FlowUnit
{
public:
    FlowUnit_W() = default;
    FlowUnit_W(const ParamReservoir& rs_param, const USI& i){};

    void CalKrPc(const OCP_DBL* S_in,
                 OCP_DBL*       kr_out,
                 OCP_DBL*       pc_out,
                 const OCP_DBL& MySurTen,
                 OCP_DBL&       MyFk,
                 OCP_DBL&       MyFp) override;
    void CalKrPcDeriv(const OCP_DBL* S_in,
                      OCP_DBL*       kr_out,
                      OCP_DBL*       pc_out,
                      OCP_DBL*       dkrdS,
                      OCP_DBL*       dPcjdS,
                      const OCP_DBL& MySurTen,
                      OCP_DBL&       MyFk,
                      OCP_DBL&       MyFp) override;

    OCP_DBL GetSwco() const override { return Swco; };

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override { return 0; }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }

protected:
    OCP_DBL Swco;
};

///////////////////////////////////////////////
// FlowUnit_OW
///////////////////////////////////////////////

class FlowUnit_OW : public FlowUnit
{
public:
    FlowUnit_OW() = default;
    FlowUnit_OW(const ParamReservoir& rs_param, const USI& i);

    void CalKrPc(const OCP_DBL* S_in,
                 OCP_DBL*       kr_out,
                 OCP_DBL*       pc_out,
                 const OCP_DBL& MySurTen,
                 OCP_DBL&       MyFk,
                 OCP_DBL&       MyFp) override;
    void CalKrPcDeriv(const OCP_DBL* S_in,
                      OCP_DBL*       kr_out,
                      OCP_DBL*       pc_out,
                      OCP_DBL*       dkrdS,
                      OCP_DBL*       dPcjdS,
                      const OCP_DBL& MySurTen,
                      OCP_DBL&       MyFk,
                      OCP_DBL&       MyFp) override;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return SWOF.Eval(0, sw, 3); }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override
    {
        return SWOF.Eval_Inv(3, pcow, 0);
    }
    OCP_DBL GetSwco() const override { return Swco; };

    // useless
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return 0; }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }

private:
    OCPTable        SWOF;  ///< saturation table about water and oil.
    vector<OCP_DBL> data;  ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata; ///< container to store the slopes of interpolation.
    OCP_DBL         Swco;  ///< Saturation of connate water.
    // OCP_DBL         kroMax;
};

///////////////////////////////////////////////
// FlowUnit_OG
///////////////////////////////////////////////

class FlowUnit_OG : public FlowUnit
{
public:
    FlowUnit_OG() = default;
    FlowUnit_OG(const ParamReservoir& rs_param, const USI& i);

    void    CalKrPc(const OCP_DBL* S_in,
                    OCP_DBL*       kr_out,
                    OCP_DBL*       pc_out,
                    const OCP_DBL& MySurTen,
                    OCP_DBL&       MyFk,
                    OCP_DBL&       MyFp) override;
    void    CalKrPcDeriv(const OCP_DBL* S_in,
                         OCP_DBL*       kr_out,
                         OCP_DBL*       pc_out,
                         OCP_DBL*       dkrdS,
                         OCP_DBL*       dPcjdS,
                         const OCP_DBL& MySurTen,
                         OCP_DBL&       MyFk,
                         OCP_DBL&       MyFp) override;
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGOF.Eval(0, sg, 3); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGOF.Eval(3, pcgo, 0); }

    // useless
    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return 0; }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override { return 0; }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override { return 0; }
    OCP_DBL GetSwco() const override { return 0; };

    const vector<OCP_DBL>& GetScm() const override
    {
        OCP_ABORT("Not Completed in this Condition!");
    }

private:
    OCPTable        SGOF;  ///< saturation table about gas and oil.
    vector<OCP_DBL> data;  ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata; ///< container to store the slopes of interpolation.
    OCP_DBL         kroMax;
};

///////////////////////////////////////////////
// FlowUnit_ODGW
///////////////////////////////////////////////

class FlowUnit_ODGW : public FlowUnit
{
public:
    OCP_DBL CalKro_Stone2(const OCP_DBL& krow,
                          const OCP_DBL& krog,
                          const OCP_DBL& krw,
                          const OCP_DBL& krg) const;

    OCP_DBL CalKro_Default(const OCP_DBL& Sg,
                           const OCP_DBL& Sw,
                           const OCP_DBL& krog,
                           const OCP_DBL& krow) const;

    OCP_DBL                GetSwco() const override { return Swco; };
    const vector<OCP_DBL>& GetScm() const override { return Scm; }

protected:
    /// oil relative permeability in the presence of connate water only, used in stone2
    OCP_DBL         kroMax;
    OCP_DBL         Swco; ///< Saturation of connate water.
    vector<OCP_DBL> Scm;  ///< critical saturation when phase becomes mobile / immobile
};

///////////////////////////////////////////////
// FlowUnit_ODGW01
///////////////////////////////////////////////

class FlowUnit_ODGW01 : public FlowUnit_ODGW
{
public:
    FlowUnit_ODGW01() = default;
    FlowUnit_ODGW01(const ParamReservoir& rs_param, const USI& i);

    virtual void CalKrPc(const OCP_DBL* S_in,
                         OCP_DBL*       kr_out,
                         OCP_DBL*       pc_out,
                         const OCP_DBL& MySurTen,
                         OCP_DBL&       MyFk,
                         OCP_DBL&       MyFp) override;
    virtual void CalKrPcDeriv(const OCP_DBL* S_in,
                              OCP_DBL*       kr_out,
                              OCP_DBL*       pc_out,
                              OCP_DBL*       dkrdS,
                              OCP_DBL*       dPcjdS,
                              const OCP_DBL& MySurTen,
                              OCP_DBL&       MyFk,
                              OCP_DBL&       MyFp) override;

    OCP_DBL CalKro_Stone2Der(OCP_DBL  krow,
                             OCP_DBL  krog,
                             OCP_DBL  krw,
                             OCP_DBL  krg,
                             OCP_DBL  dkrwdSw,
                             OCP_DBL  dkrowdSw,
                             OCP_DBL  dkrgdSg,
                             OCP_DBL  dkrogdSg,
                             OCP_DBL& out_dkrodSw,
                             OCP_DBL& out_dkrodSg) const;
    OCP_DBL CalKro_DefaultDer(const OCP_DBL& Sg,
                              const OCP_DBL& Sw,
                              const OCP_DBL& krog,
                              const OCP_DBL& krow,
                              const OCP_DBL& dkrogSg,
                              const OCP_DBL& dkrowSw,
                              OCP_DBL&       dkroSg,
                              OCP_DBL&       dkroSw) const;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return SWOF.Eval(0, sw, 3); }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override
    {
        return SWOF.Eval_Inv(3, pcow, 0);
    }

    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGOF.Eval(0, sg, 3); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGOF.Eval(3, pcgo, 0); }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override
    {
        return SWPCGW.Eval_Inv(1, pcgw, 0);
    }

    void Generate_SWPCWG();

protected:
    OCPTable        SGOF;   ///< saturation table about gas and oil.
    OCPTable        SWOF;   ///< saturation table about water and oil.
    vector<OCP_DBL> data;   ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata;  ///< container to store the slopes of interpolation.
    OCPTable        SWPCGW; ///< auxiliary table: saturation of water vs. capillary
                            ///< pressure between water and gas.
};

///////////////////////////////////////////////
// FlowUnit_ODGW01_Miscible
///////////////////////////////////////////////

class FlowUnit_ODGW01_Miscible : public FlowUnit_ODGW01
{
public:
    FlowUnit_ODGW01_Miscible(const ParamReservoir& rs_param, const USI& i)
        : FlowUnit_ODGW01(rs_param, i)
    {

        // gas is moveable all the time
        Scm[1] = 0;

        USI len = rs_param.miscstr.surTenRef.size();
        if (len > 0) {
            surTenRef = rs_param.miscstr.surTenRef[0];
            surTenPc  = 1;
            if (len > 2) {
                surTenPc = rs_param.miscstr.surTenRef[2] / surTenRef;
            }
        } else {
            OCP_ABORT("No data in MISCSTR!");
        }
    }

    void CalKrPc(const OCP_DBL* S_in,
                 OCP_DBL*       kr_out,
                 OCP_DBL*       pc_out,
                 const OCP_DBL& MySurTen,
                 OCP_DBL&       MyFk,
                 OCP_DBL&       MyFp) override;
    void CalKrPcDeriv(const OCP_DBL* S_in,
                      OCP_DBL*       kr_out,
                      OCP_DBL*       pc_out,
                      OCP_DBL*       dkrdS,
                      OCP_DBL*       dPcjdS,
                      const OCP_DBL& MySurTen,
                      OCP_DBL&       MyFk,
                      OCP_DBL&       MyFp) override;

private:
    OCP_DBL Fk{0};        ///< The relative permeability interpolation parameter
    OCP_DBL Fp{0};        ///< The capillary pressure interpolation parameter
    OCP_DBL surTenRef{0}; ///< Reference surface tension
    OCP_DBL Fkexp{0.25};  ///< Exponent set used to calculate Fk
    OCP_DBL surTenPc{0};  ///< Maximum surface tension for capillary pressure/surTenRef
    OCP_DBL surTen{0};    ///< Surface tension

    /*
    OCP_DBL kroMis{0};  ///< miscible oil relative permeability
    OCP_DBL krgMis{0};  ///< miscible gas relative permeability
    OCP_DBL PcogMis{0}; ///< miscible gas capillary pressure
    OCP_DBL socrMis{0}; ///< oil critical miscible saturations
    OCP_DBL sgcrMis{0}; ///< gas critical miscible saturations
    */
};

///////////////////////////////////////////////
// FlowUnit_ODGW02
///////////////////////////////////////////////

class FlowUnit_ODGW02 : public FlowUnit_ODGW
{
public:
    FlowUnit_ODGW02() = default;
    FlowUnit_ODGW02(const ParamReservoir& rs_param, const USI& i);

    void    CalKrPc(const OCP_DBL* S_in,
                    OCP_DBL*       kr_out,
                    OCP_DBL*       pc_out,
                    const OCP_DBL& MySurTen,
                    OCP_DBL&       MyFk,
                    OCP_DBL&       MyFp) override;
    void    CalKrPcDeriv(const OCP_DBL* S_in,
                         OCP_DBL*       kr_out,
                         OCP_DBL*       pc_out,
                         OCP_DBL*       dkrdS,
                         OCP_DBL*       dPcjdS,
                         const OCP_DBL& MySurTen,
                         OCP_DBL&       MyFk,
                         OCP_DBL&       MyFp) override;
    OCP_DBL CalKro_Stone2Der(OCP_DBL  krow,
                             OCP_DBL  krog,
                             OCP_DBL  krw,
                             OCP_DBL  krg,
                             OCP_DBL  dkrwdSw,
                             OCP_DBL  dkrowdSo,
                             OCP_DBL  dkrgdSg,
                             OCP_DBL  dkrogdSo,
                             OCP_DBL& out_dkrodSo) const;
    OCP_DBL CalKro_DefaultDer(const OCP_DBL& Sg,
                              const OCP_DBL& Sw,
                              const OCP_DBL& krog,
                              const OCP_DBL& krow,
                              const OCP_DBL& dkrogSo,
                              const OCP_DBL& dkrowSo,
                              OCP_DBL&       dkroSo) const;

    OCP_DBL GetPcowBySw(const OCP_DBL& sw) override { return SWFN.Eval(0, sw, 2); }
    OCP_DBL GetSwByPcow(const OCP_DBL& pcow) override
    {
        return SWFN.Eval_Inv(2, pcow, 0);
    }
    OCP_DBL GetPcgoBySg(const OCP_DBL& sg) override { return SGFN.Eval(0, sg, 2); }
    OCP_DBL GetSgByPcgo(const OCP_DBL& pcgo) override { return SGFN.Eval(2, pcgo, 0); }
    OCP_DBL GetSwByPcgw(const OCP_DBL& pcgw) override
    {
        return SWPCGW.Eval_Inv(1, pcgw, 0);
    }
    void Generate_SWPCWG();

private:
    OCPTable        SWFN;   ///< saturation table about water.
    OCPTable        SGFN;   ///< saturation table about gas.
    OCPTable        SOF3;   ///< saturation table about oil.
    vector<OCP_DBL> data;   ///< container to store the values of interpolation.
    vector<OCP_DBL> cdata;  ///< container to store the slopes of interpolation.
    OCPTable        SWPCGW; ///< auxiliary table: saturation of water vs. capillary
                            ///< pressure between water and gas.
};

#endif /* end if __FLOWUNIT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/05/2022      Format file                          */
/*----------------------------------------------------------------------------*/