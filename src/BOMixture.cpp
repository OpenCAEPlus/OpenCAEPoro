#include "BOMixture.hpp"

BOMixture::BOMixture(const ParamReservoir& rs_param, const USI& PVTmode, const USI& i)
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
		len = max(len, PVTW.getCol());
	}
	if (rs_param.PVCO_T.data.size() != 0) {
		PVCO.setup(rs_param.PVCO_T.data[i]);
		len = max(len, PVCO.getCol());
	}
	if (rs_param.PVDO_T.data.size() != 0) {
		PVDO.setup(rs_param.PVDO_T.data[i]);
		len = max(len, PVDO.getCol());
	}
	if (rs_param.PVDG_T.data.size() != 0) {
		PVDG.setup(rs_param.PVDG_T.data[i]);
		len = max(len, PVDG.getCol());
	}
	data.resize(len, 0);
	cdata.resize(len, 0);
}

void BOMixture::Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin, const OCP_DBL* Sjin, const OCP_DBL& Vpore, const OCP_DBL* Ziin)
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

void BOMixture::Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin)
{
#ifdef _DEBUG
	// checkNi(Niin);
#endif // _DEBUG

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


void BOMixture::BOFlash_Ni_W(const OCP_DBL& Pin, const OCP_DBL* Niin)
{

}

void BOMixture::BOFlash_Ni_OW(const OCP_DBL& Pin, const OCP_DBL* Niin)
{

}

void BOMixture::BOFlash_Sj_W(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore)
{

}

void BOMixture::BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin, const OCP_DBL& Vpore)
{

}

OCP_DBL BOMixture::xiPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin)
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

    return 0.0; // TODO: Make sure code does not reach here!
}

OCP_DBL BOMixture::rhoPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin)
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

    return 0.0; // TODO: Make sure code does not reach here!
}


OCP_DBL BOMixture::gammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin)
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

OCP_DBL BOMixture::gammaPhaseO_OW(const OCP_DBL& Pin)
{

	return 0;
}

OCP_DBL BOMixture::gammaPhaseW(const OCP_DBL& Pin)
{

	PVTW.eval_all(0, Pin, data, cdata);
	OCP_DBL Pw0 = data[0];
	OCP_DBL bw0 = data[1];
	OCP_DBL cbw = data[2];
	OCP_DBL bw = (bw0 * (1 - cbw * (Pin - Pw0)));

	return Std_GammaW / bw;
}

OCP_DBL BOMixture::gammaPhaseG(const OCP_DBL& Pin)
{
	if (PVDG.isempty()) {
		ERRORcheck("PVDG is missing");
	}
	OCP_DBL bg = (CONV1 / 1000) * PVDG.eval(0, Pin, 1);

	return Std_GammaG / bg;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/