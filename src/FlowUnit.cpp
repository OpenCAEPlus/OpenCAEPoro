#include "FlowUnit.hpp"
#include "OpenCAEPoro_consts.hpp"

void FlowUnit::generate_SWPCWG()
{
	if (SGOF.isempty()) {
		ERRORcheck("SGOF is missing!");
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
	case PHASE_OGW:
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
	double So = S_in[0]; 
	double Sw = S_in[1];

	std::vector<double>		dataw(4, 0);
	std::vector<double>		cdataw(4, 0);
	//three phase black oil model using stone 2
	SWOF.eval_all(0, Sw, dataw, cdataw);
	double krw = dataw[1];
	double kro = dataw[2];
	double Pcw = -dataw[3];

	kr_out[0] = kro;	kr_out[1] = krw;
	pc_out[0] = 0;		pc_out[1] = Pcw;
}

void FlowUnit::calKrPc_OG(const double* S_in, double* kr_out, double* pc_out)
{

	double So = S_in[0]; 
	double Sg = S_in[1];

	std::vector<double>		datag(4, 0);
	std::vector<double>		cdatag(4, 0);
	//three phase black oil model using stone 2
	SGOF.eval_all(0, Sg, datag, cdatag);
	double krg = datag[1];
	double kro = datag[2];
	double Pcg = datag[3];

	kr_out[0] = kro;	kr_out[1] = krg;
	pc_out[0] = 0;		pc_out[1] = Pcg;
}

void FlowUnit::calKrPc_OGW(const double* S_in, double* kr_out, double* pc_out)
{
	double So = S_in[0]; 
	double Sg = S_in[1]; 
	double Sw = S_in[2];

	std::vector<double>		data(4, 0);
	std::vector<double>		cdata(4, 0);
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
