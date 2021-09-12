#pragma once
#include "ReservoirTable.hxx"
#include "OpenCAEPoro_consts.hpp"

class FlowUnit
{
public:

	bool empty_SGOF() { return SGOF.isempty(); }
	bool empty_SWOF() { return SWOF.isempty(); }

	void generate_SWPCWG();

	int eval_SWOF(int j, double val, int destj) { return SWOF.eval(j, val, destj); }
	int evalinv_SWOF(int j, double val, int destj) { return SWOF.eval_inv(j, val, destj); }
	int eval_SGOF(int j, double val, int destj) { return SGOF.eval(j, val, destj); }
	int evalinv_SGOF(int j, double val, int destj) { return SGOF.eval_inv(j, val, destj); }
	int eval_SWPCWG(int j, double val, int destj) { return SWPCWG.eval(j, val, destj); }
	int evalinv_SWPCWG(int j, double val, int destj) { return SWPCWG.eval_inv(j, val, destj); }


	void calKrPc(const double* S_in, double* kr_out, double* pc_out);
	void calKrPc_W(double* kr_out, double* pc_out);
	void calKrPc_OW(const double* S_in, double* kr_out, double* pc_out);
	void calKrPc_OG(const double* S_in, double* kr_out, double* pc_out);
	void calKrPc_OGW(const double* S_in, double* kr_out, double* pc_out);


	// FlowUnits Model
	double kro_stone2(double krow, double krog, double krw, double krg);

private:
	int								Mode;
	ReservoirTable<double>			SWOF;
	ReservoirTable<double>			SGOF;
	ReservoirTable<double>			SWPCWG;

	double							KroMax;		// oil relative permeability in the presence of connate water only
};


