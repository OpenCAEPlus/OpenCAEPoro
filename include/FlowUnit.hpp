/*! \file    FlowUnit.hpp
 *  \brief   FlowUnit class declaration
 *  \author  Shizhe Li
 *  \date    Oct/06/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __FLOWUNIT_HEADER__
#define __FLOWUNIT_HEADER__


#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamReservoir.hpp"

/// designed to deal with matters related to saturation table.
/// relative permeability, capillary pressure woulbe be calculated here.
class FlowUnit
{
public:
	FlowUnit() = default;
	FlowUnit(const ParamReservoir& rs_param, const USI& mode, const USI& i);

	/// judge if SGOF is empty.
	bool empty_SGOF() { return SGOF.isempty(); }
	/// judge if SWOF is empty.
	bool empty_SWOF() { return SWOF.isempty(); }

	/// generate the table of saturation of water vs. capillary between water and gas by SGOF and SWOF.
	void generate_SWPCWG();
	/// interpolate the specified monotonically increasing column in SWOF to evaluate the target column.
	OCP_DBL eval_SWOF(const USI& j, const OCP_DBL& val, const USI& destj) { return SWOF.eval(j, val, destj); }
	/// interpolate the specified monotonically decreasing column in SWOF to evaluate the target column.
	OCP_DBL evalinv_SWOF(const USI& j, const OCP_DBL& val, const USI& destj) { return SWOF.eval_inv(j, val, destj); }
	/// interpolate the specified monotonically increasing column in SGOF to evaluate the target column.
	OCP_DBL eval_SGOF(const USI& j, const OCP_DBL& val, const USI& destj) { return SGOF.eval(j, val, destj); }
	/// interpolate the specified monotonically decreasing column in SGOF to evaluate the target column.
	OCP_DBL evalinv_SGOF(const USI& j, const OCP_DBL& val, const USI& destj) { return SGOF.eval_inv(j, val, destj); }
	/// interpolate the specified monotonically increasing column in SWPCWG to evaluate the target column.
	OCP_DBL eval_SWPCWG(const USI& j, const OCP_DBL& val, const USI& destj) { return SWPCWG.eval(j, val, destj); }
	/// interpolate the specified monotonically decreasing column in SWPCWG to evaluate the target column.
	OCP_DBL evalinv_SWPCWG(const USI& j, const OCP_DBL& val, const USI& destj) { return SWPCWG.eval_inv(j, val, destj); }

	/// calculate relative permeability and capillary pressure.
	void calKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	/// calculate relative permeability and capillary pressure for water.
	/// it will be used if water exists or could be exist.
	void calKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out);
	/// calculate relative permeability and capillary pressure for oil and water.
	/// it will be used if oil and water exist or could be exist.
	void calKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	/// calculate relative permeability and capillary pressure for oil and gas.
	/// it will be used if oil and gas exist or could be exist. 
	void calKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	/// calculate relative permeability and capillary pressure for oil, gas and water.
	/// it will be used if oil, gas and water exist or could be exist. 
	void calKrPc_OGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);


	// FlowUnits Model
	/// calculate relative permeability of oil phase with stone2 method.
	OCP_DBL kro_stone2(const OCP_DBL& krow, const OCP_DBL& krog, const OCP_DBL& krw, const OCP_DBL& krg) const;

private:
	USI								Mode;		///< decide which saturation table will be uesd.
	ReservoirTable<OCP_DBL>			SWOF;		///< saturation table about water and oil.
	ReservoirTable<OCP_DBL>			SGOF;		///< saturation table about gas and oil.
	ReservoirTable<OCP_DBL>			SWPCWG;		///< auxiliary table: saturation of water vs. capillary pressure between water and gas.

	OCP_DBL							KroMax;		///< oil relative permeability in the presence of connate water only.

	// Auxiliary parameters for Table interpolation
    USI                             len{0};		///< maximum number of columns of tables among all above.
    vector<OCP_DBL>					data;		///< container used to store the results of values of interpolation.
    vector<OCP_DBL>					cdata;		///< container used to store the results of slopes of interpolation.
};


#endif