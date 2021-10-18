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

    /// TODO: Add Doxygen
    FlowUnit(const ParamReservoir& rs_param, const USI& inmode, const USI& i);

    /// Check whether SGOF is empty.
    bool IsEmpty_SGOF() { return SGOF.IsEmpty(); }

    /// Check whether SWOF is empty.
    bool IsEmpty_SWOF() { return SWOF.IsEmpty(); }

    /// Generate the table of water saturation vs. capillary between water and gas.
    void Generate_SWPCWG();

    /// interpolate the specified monotonically increasing column in SWOF to evaluate
    /// the target column.
    OCP_DBL Eval_SWOF(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SWOF.Eval(j, val, destj);
    }

    /// interpolate the specified monotonically decreasing column in SWOF to evaluate
    /// the target column.
    OCP_DBL EvalInv_SWOF(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SWOF.Eval_Inv(j, val, destj);
    }

    /// interpolate the specified monotonically increasing column in SGOF to evaluate
    /// the target column.
    OCP_DBL Eval_SGOF(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SGOF.Eval(j, val, destj);
    }

    /// interpolate the specified monotonically decreasing column in SGOF to evaluate
    /// the target column.
    OCP_DBL EvalInv_SGOF(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SGOF.Eval_Inv(j, val, destj);
    }

    /// interpolate the specified monotonically increasing column in SWPCWG to evaluate
    /// the target column.
    OCP_DBL Eval_SWPCWG(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SWPCWG.Eval(j, val, destj);
    }

    /// interpolate the specified monotonically decreasing column in SWPCWG to evaluate
    /// the target column.
    OCP_DBL EvalInv_SWPCWG(const USI& j, const OCP_DBL& val, const USI& destj)
    {
        return SWPCWG.Eval_Inv(j, val, destj);
    }

    /// calculate relative permeability and capillary pressure.
    void CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);

    /// calculate relative permeability and capillary pressure for water.
    /// it will be used if water exists or could be exist.
    void CalKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out);

    /// calculate relative permeability and capillary pressure for oil and water.
    /// it will be used if oil and water exist or could be exist.
    void CalKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);

    /// calculate relative permeability and capillary pressure for oil and gas.
    /// it will be used if oil and gas exist or could be exist.
    void CalKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);

    /// calculate relative permeability and capillary pressure for oil, gas and water.
    /// it will be used if oil, gas and water exist or could be exist.
    void CalKrPc_OGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);

    // FlowUnits Model
    /// calculate relative permeability of oil phase with stone2 method.
    OCP_DBL CalKro_Stone2(const OCP_DBL& krow, const OCP_DBL& krog, const OCP_DBL& krw,
                          const OCP_DBL& krg) const;

private:
    USI                mode;   ///< decide which saturation table will be uesd.
    OCP_Table<OCP_DBL> SWOF;   ///< saturation table about water and oil.
    OCP_Table<OCP_DBL> SGOF;   ///< saturation table about gas and oil.
    OCP_Table<OCP_DBL> SWPCWG; ///< auxiliary table: saturation of water vs. capillary
                               ///< pressure between water and gas.
    OCP_DBL
    kroMax; ///< oil relative permeability in the presence of connate water only.

    // Auxiliary parameters for Table interpolation
    USI len{0}; ///< maximum number of columns of tables among all above.
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