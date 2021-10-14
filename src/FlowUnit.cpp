#include "FlowUnit.hpp"
#include "OpenCAEPoro_consts.hpp"

void FlowUnit::Generate_SWPCWG()
{
	if (SGOF.IsEmpty()) {
		ERRORcheck("SGOF is missing!");
		exit(0);
	}
	if (SWOF.IsEmpty()) {
		ERRORcheck("SWOF is missing!");
		exit(0);
	}

	std::vector<OCP_DBL>		Sw(SWOF.GetCol(0));
	std::vector<OCP_DBL>		Pcw(SWOF.GetCol(3));
	USI n = Sw.size();
	for (USI i = 0; i < n; i++) {
		OCP_DBL Pcg = SGOF.Eval(0, 1 - Sw[i], 3);
		Pcw[i] += Pcg;
	}

	SWPCWG.PushCol(Sw);
	SWPCWG.PushCol(Pcw);
	SWPCWG.SetRowCol();
}



void FlowUnit::CalKrPc(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	switch (mode)
	{
	case PHASE_W:
		CalKrPc_W(kr_out, pc_out);
		break;
	case PHASE_OW:
		CalKrPc_OW(S_in, kr_out, pc_out);
		break;
	case PHASE_OG:
		CalKrPc_OG(S_in, kr_out, pc_out);
		break;
	case PHASE_OGW: // should SWOF be used here?
		CalKrPc_OGW(S_in, kr_out, pc_out);
		break;
	default:
		ERRORcheck("Wrong Mode");
		break;
	}
}


void FlowUnit::CalKrPc_W(OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	kr_out[0] = 1;
	pc_out[0] = 0;
}


void FlowUnit::CalKrPc_OW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sw = S_in[1];

	//three phase black oil model using stone 2
	SWOF.Eval_All(0, Sw, data, cdata);
	OCP_DBL krw = data[1];
	OCP_DBL kro = data[2];
	OCP_DBL Pcw = -data[3];

	kr_out[0] = kro;	kr_out[1] = krw;
	pc_out[0] = 0;		pc_out[1] = Pcw;
}

void FlowUnit::CalKrPc_OG(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sg = S_in[1];

	//three phase black oil model using stone 2
	SGOF.Eval_All(0, Sg, data, cdata);
	OCP_DBL krg = data[1];
	OCP_DBL kro = data[2];
	OCP_DBL Pcg = data[3];

	kr_out[0] = kro;	kr_out[1] = krg;
	pc_out[0] = 0;		pc_out[1] = Pcg;
}

void FlowUnit::CalKrPc_OGW(const OCP_DBL* S_in, OCP_DBL* kr_out, OCP_DBL* pc_out)
{
	OCP_DBL Sg = S_in[1]; 
	OCP_DBL Sw = S_in[2];

	//three phase black oil model using stone 2
	SWOF.Eval_All(0, Sw, data, cdata);
	OCP_DBL krw = data[1];
	OCP_DBL krow = data[2];
	OCP_DBL Pcw = -data[3];

	SGOF.Eval_All(0, Sg, data, cdata);
	OCP_DBL krg = data[1];
	OCP_DBL krog = data[2];
	OCP_DBL Pcg = data[3];

	OCP_DBL kro = CalKro_Stone2(krow, krog, krw, krg);

	kr_out[0] = kro;  kr_out[1] = krg;  kr_out[2] = krw;
	pc_out[0] = 0;    pc_out[1] = Pcg;  pc_out[2] = Pcw;
}

OCP_DBL FlowUnit::CalKro_Stone2(const OCP_DBL& krow, const OCP_DBL& krog, const OCP_DBL& krw, const OCP_DBL& krg) const
{
	// krog : oil relative permeability for a system with oil, gas and connate water
	// krow : oil relative permeability for a system with oil and water only

	OCP_DBL kro = kroMax * ((krow / kroMax + krw) * (krog / kroMax + krg) - (krw + krg));
	if (kro < 0)
		kro = 0;

	return kro;
}

FlowUnit::FlowUnit(const ParamReservoir& rs_param, const USI& inmode, const USI& i)
{
	mode = inmode;
	if (rs_param.water) {
		SWOF.Setup(rs_param.SWOF_T.data[i]);
        len = max(len, SWOF.GetCol());
	}
	if (rs_param.gas) {
		SGOF.Setup(rs_param.SGOF_T.data[i]);
        len = max(len, SWOF.GetCol());
	}

	data.resize(len, 0);
    cdata.resize(len, 0);
	
	kroMax = 0;
	if (rs_param.water) {
		kroMax = SWOF.GetCol(2)[0];
	}
	else if (rs_param.gas) {
		kroMax = SGOF.GetCol(2)[0];
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/