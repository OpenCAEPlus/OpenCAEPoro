#pragma once
#include "ReservoirTable.hxx"

class FlowUnit
{
public:

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

	double							KroMax;		// oil relative permeability in the presence of connate water only
};



