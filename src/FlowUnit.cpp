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

	std::vector<OCP_DBL>		Sw(SWOF.getCol(0));
	std::vector<OCP_DBL>		Pcw(SWOF.getCol(3));
	USI n = Sw.size();
	for (USI i = 0; i < n; i++) {
		OCP_DBL Pcg = SGOF.eval(0, 1 - Sw[i], 3);
		Pcw[i] += Pcg;
	}

	SWPCWG.pushCol(Sw);
	SWPCWG.pushCol(Pcw);
	SWPCWG.setRowCol();
}



void FlowUnit::calKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
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


void FlowUnit::calKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	kr_out[0] = 1;
	pc_out[0] = 0;
}


void FlowUnit::calKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sw = S_in[1];

	//three phase black oil model using stone 2
	SWOF.eval_all(0, Sw, data, cdata);
	OCP_DBL krw = data[1];
	OCP_DBL kro = data[2];
	OCP_DBL Pcw = -data[3];

	kr_out[0] = kro;	kr_out[1] = krw;
	pc_out[0] = 0;		pc_out[1] = Pcw;
}

void FlowUnit::calKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sg = S_in[1];

	//three phase black oil model using stone 2
	SGOF.eval_all(0, Sg, data, cdata);
	OCP_DBL krg = data[1];
	OCP_DBL kro = data[2];
	OCP_DBL Pcg = data[3];

	kr_out[0] = kro;	kr_out[1] = krg;
	pc_out[0] = 0;		pc_out[1] = Pcg;
}

void FlowUnit::calKrPc_OGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sg = S_in[1]; 
	OCP_DBL Sw = S_in[2];

	//three phase black oil model using stone 2
	SWOF.eval_all(0, Sw, data, cdata);
	OCP_DBL krw = data[1];
	OCP_DBL krow = data[2];
	OCP_DBL Pcw = -data[3];

	SGOF.eval_all(0, Sg, data, cdata);
	OCP_DBL krg = data[1];
	OCP_DBL krog = data[2];
	OCP_DBL Pcg = data[3];

	OCP_DBL kro = kro_stone2(krow, krog, krw, krg);

	kr_out[0] = kro;  kr_out[1] = krg;  kr_out[2] = krw;
	pc_out[0] = 0;    pc_out[1] = Pcg;  pc_out[2] = Pcw;
}

OCP_DBL FlowUnit::kro_stone2(const OCP_DBL& krow, const OCP_DBL& krog, const OCP_DBL& krw, const OCP_DBL& krg) const
{
	// krog : oil relative permeability for a system with oil, gas and connate water
	// krow : oil relative permeability for a system with oil and water only

	OCP_DBL kro = KroMax * ((krow / KroMax + krw) * (krog / KroMax + krg) - (krw + krg));
	if (kro < 0)
		kro = 0;

	return kro;
}

FlowUnit::FlowUnit(const ParamReservoir& rs_param, const USI& mode, const USI& i)
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
