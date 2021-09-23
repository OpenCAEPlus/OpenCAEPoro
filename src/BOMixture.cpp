#include "BOMixture.hpp"

BOMixture::BOMixture(ParamReservoir& rs_param, int PVTmode, int i)
{
	MixtureType = BLKOIL;
	Mode = PVTmode;
	Np = rs_param.Np;
	Nc = rs_param.Nc;

	Ni.resize(Nc);
	PhaseExist.resize(Np);
	V.resize(Np);
	S.resize(Np);
	Xi.resize(Np);
	Cij.resize(Np * Nc);
	Rho.resize(Np);
	Mu.resize(Np);
	Vfi.resize(Nc);

	if (rs_param.Density.activity) {
		Std_RhoO = rs_param.Density.data[0];
		Std_RhoW = rs_param.Density.data[1];
		Std_RhoG = rs_param.Density.data[2];
	}
	else {
		Std_RhoO = (141.5 * RHOW_STD) / (rs_param.Gravity.data[0] + 131.5);
		Std_RhoW = RHOW_STD * rs_param.Gravity.data[1];
		Std_RhoG = RHOAIR_STD * rs_param.Gravity.data[2];
	}

	Std_GammaO = GRAVITY_FACTOR * Std_RhoO;
	Std_GammaG = GRAVITY_FACTOR * Std_RhoG;
	Std_GammaW = GRAVITY_FACTOR * Std_RhoW;

	if (rs_param.PVTW_T.data.size() != 0) {
		PVTW.setup(rs_param.PVTW_T.data[i]);
	}
	if (rs_param.PVCO_T.data.size() != 0) {
		PVCO.setup(rs_param.PVCO_T.data[i]);
	}
	if (rs_param.PVDO_T.data.size() != 0) {
		PVDO.setup(rs_param.PVDO_T.data[i]);
	}
	if (rs_param.PVDG_T.data.size() != 0) {
		PVDG.setup(rs_param.PVDG_T.data[i]);
	}
}

void BOMixture::Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin)
{
	switch (Mode)
	{
	case PHASE_W:
		BOFlash_Sj_W(Pin, Sjin, Vpore);
		break;

	case PHASE_OW:
		BOFlash_Sj_OW(Pin, Sjin, Vpore);
		break;

	case PHASE_OGW:
		BOFlash_Sj_OGW(Pin, Pbbin, Sjin, Vpore);
		break;

	default:
		ERRORcheck("Wrong Mode");
		exit(0);
		break;
	}
}

void BOMixture::Flash_Ni(const double Pin, const double Tin, const double* Niin)
{
	switch (Mode)
	{
	case PHASE_W:
		BOFlash_Ni_W(Pin, Niin);
		break;

	case PHASE_OW:
		BOFlash_Ni_OW(Pin, Niin);
		break;

	case PHASE_OGW:
		BOFlash_Ni_OGW(Pin, Niin);
		break;

	default:
		ERRORcheck("Wrong Mode");
		exit(0);
	}
}


void BOMixture::BOFlash_Ni_W(const double Pin, const double* Niin)
{

}

void BOMixture::BOFlash_Ni_OW(const double Pin, const double* Niin)
{

}

void BOMixture::BOFlash_Sj_W(const double Pin, const double* Sjin, double Vpore)
{

}

void BOMixture::BOFlash_Sj_OW(const double Pin, const double* Sjin, double Vpore)
{

}

double BOMixture::xiPhase(double Pin, double T, double* Ziin)
{
	switch (Mode)
	{
	case PHASE_W:
		break;
	case PHASE_OW:
		break;
	case PHASE_OGW:
		return xiPhase_OGW(Pin, Ziin);
		break;
	default:
		break;
	}
}

double BOMixture::rhoPhase(double Pin, double T, double* Ziin)
{
	switch (Mode)
	{
	case PHASE_W:
		break;
	case PHASE_OW:
		break;
	case PHASE_OGW:
		return rhoPhase_OGW(Pin, Ziin);
		break;
	default:
		break;
	}
}


double BOMixture::gammaPhaseO(double Pin, double Pbbin)
{
	switch (Mode)
	{
	case PHASE_OW:
		return gammaPhaseO_OW(Pin);
	case PHASE_OGW:
		return gammaPhaseO_OGW(Pin, Pbbin);
	default:
		ERRORcheck("Wrong Mode");
		exit(0);
	}
}

double BOMixture::gammaPhaseO_OW(double Pin)
{

	return 0;
}

double BOMixture::gammaPhaseW(double Pin)
{
	std::vector<double>		data(5, 0);
	std::vector<double>		cdata(5, 0);

	PVTW.eval_all(0, Pin, data, cdata);
	double Pw0 = data[0];
	double bw0 = data[1];
	double cbw = data[2];
	double bw = (bw0 * (1 - cbw * (Pin - Pw0)));

	return Std_GammaW / bw;
}

double BOMixture::gammaPhaseG(double Pin)
{
	if (PVDG.isempty()) {
		ERRORcheck("PVDG is missing");
	}
	double bg = (CONV1 / 1000) * PVDG.eval(0, Pin, 1);

	return Std_GammaG / bg;
}
