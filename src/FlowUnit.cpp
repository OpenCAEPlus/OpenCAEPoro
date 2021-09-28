#include "FlowUnit.hpp"
#include "OpenCAEPoro_consts.hpp"

void FlowUnit::generate_SWPCWG()
{
	if (SGOF.isempty()) {
		ERRORcheck("SGOF is missing!");
		exit(0);
	}
	if (SWOF.isempty()) {
		ERRORcheck("SWOF is missing!");
		exit(0);
	}

	std::vector<double>		Sw(SWOF.getCol(0));
	std::vector<double>		Pcw(SWOF.getCol(3));
	int n = Sw.size();
	for (int i = 0; i < n; i++) {
		double Pcg = SGOF.eval(0, 1 - Sw[i], 3);
		Pcw[i] += Pcg;
	}

	SWPCWG.pushCol(Sw);
	SWPCWG.pushCol(Pcw);
	SWPCWG.setRowCol();
}



void FlowUnit::calKrPc(const double* S_in, double* kr_out, double* pc_out)
{
	switch (Mode)
	{
	case PHASE_W:
		calKrPc_W(kr_out, pc_out);
		break;
	case PHASE_OW:
		calKrPc_OW(S_in, kr_out, pc_out);
		break;
	case PHASE_OG:
		calKrPc_OG(S_in, kr_out, pc_out);
		break;
	case PHASE_OGW: // should SWOF be used here?
		calKrPc_OGW(S_in, kr_out, pc_out);
		break;
	default:
		ERRORcheck("Wrong Mode");
		break;
	}
}


void FlowUnit::calKrPc_W(double* kr_out, double* pc_out)
{
	kr_out[0] = 1;
	pc_out[0] = 0;
}


void FlowUnit::calKrPc_OW(const double* S_in, double* kr_out, double* pc_out)
{
	double Sw = S_in[1];

	//three phase black oil model using stone 2
	SWOF.eval_all(0, Sw, data, cdata);
	double krw = data[1];
	double kro = data[2];
	double Pcw = -data[3];

	kr_out[0] = kro;	kr_out[1] = krw;
	pc_out[0] = 0;		pc_out[1] = Pcw;
}

void FlowUnit::calKrPc_OG(const double* S_in, double* kr_out, double* pc_out)
{
	double Sg = S_in[1];

	//three phase black oil model using stone 2
	SGOF.eval_all(0, Sg, data, cdata);
	double krg = data[1];
	double kro = data[2];
	double Pcg = data[3];

	kr_out[0] = kro;	kr_out[1] = krg;
	pc_out[0] = 0;		pc_out[1] = Pcg;
}

void FlowUnit::calKrPc_OGW(const double* S_in, double* kr_out, double* pc_out)
{
	double Sg = S_in[1]; 
	double Sw = S_in[2];

	//three phase black oil model using stone 2
	SWOF.eval_all(0, Sw, data, cdata);
	double krw = data[1];
	double krow = data[2];
	double Pcw = -data[3];

	SGOF.eval_all(0, Sg, data, cdata);
	double krg = data[1];
	double krog = data[2];
	double Pcg = data[3];

	double kro = kro_stone2(krow, krog, krw, krg);

	kr_out[0] = kro;  kr_out[1] = krg;  kr_out[2] = krw;
	pc_out[0] = 0;    pc_out[1] = Pcg;  pc_out[2] = Pcw;
}

double FlowUnit::kro_stone2(double krow, double krog, double krw, double krg)
{
	// krog : oil relative permeability for a system with oil, gas and connate water
	// krow : oil relative permeability for a system with oil and water only

	double kro = KroMax * ((krow / KroMax + krw) * (krog / KroMax + krg) - (krw + krg));
	if (kro < 0)
		kro = 0;

	return kro;
}

FlowUnit::FlowUnit(ParamReservoir& rs_param, int mode, int i)
{
	Mode = mode;
	if (rs_param.WATER) {
		SWOF.setup(rs_param.SWOF_T.data[i]);
        len = max(len, SWOF.getCol());
	}
	if (rs_param.GAS) {
		SGOF.setup(rs_param.SGOF_T.data[i]);
        len = max(len, SWOF.getCol());
	}

	data.resize(len, 0);
    cdata.resize(len, 0);
	
	KroMax = 0;
	if (rs_param.WATER) {
		KroMax = SWOF.getCol(2)[0];
	}
	else if (rs_param.GAS) {
		KroMax = SGOF.getCol(2)[0];
	}
}
