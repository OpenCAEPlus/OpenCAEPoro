#pragma once
#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.hpp"
#include "ParamReservoir.hpp"

class FlowUnit
{
public:
	FlowUnit() = default;
	FlowUnit(ParamReservoir& rs_param, int mode, int i);

	bool empty_SGOF() { return SGOF.isempty(); }
	bool empty_SWOF() { return SWOF.isempty(); }

	void generate_SWPCWG();

	OCP_DBL eval_SWOF(int j, OCP_DBL val, int destj) { return SWOF.eval(j, val, destj); }
	OCP_DBL evalinv_SWOF(int j, OCP_DBL val, int destj) { return SWOF.eval_inv(j, val, destj); }
	OCP_DBL eval_SGOF(int j, OCP_DBL val, int destj) { return SGOF.eval(j, val, destj); }
	OCP_DBL evalinv_SGOF(int j, OCP_DBL val, int destj) { return SGOF.eval_inv(j, val, destj); }
	OCP_DBL eval_SWPCWG(int j, OCP_DBL val, int destj) { return SWPCWG.eval(j, val, destj); }
	OCP_DBL evalinv_SWPCWG(int j, OCP_DBL val, int destj) { return SWPCWG.eval_inv(j, val, destj); }


	void calKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	void calKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out);
	void calKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	void calKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);
	void calKrPc_OGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out);


	// FlowUnits Model
	OCP_DBL kro_stone2(OCP_DBL krow, OCP_DBL krog, OCP_DBL krw, OCP_DBL krg);

private:
	int								Mode;
	ReservoirTable<OCP_DBL>			SWOF;
	ReservoirTable<OCP_DBL>			SGOF;
	ReservoirTable<OCP_DBL>			SWPCWG;

	OCP_DBL							KroMax;		// oil relative permeability in the presence of connate water only

	// Auxiliary parameters for Table interpolation
    int                             len{0};
    vector<OCP_DBL>					data;
    vector<OCP_DBL>					cdata;
};


