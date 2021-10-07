#include "BOMixture.hpp"

void BOMixture::BOFlash_Sj_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL* Sjin, const OCP_DBL& Vpore)
{
	// initialize    Np = 3, Nc = 3
	PhaseExist.assign(Np, false);
	Cij.assign(Np * Nc, 0);   // 
	Ni.assign(Nc, 0);

	P = Pin;
	for (USI j = 0; j < Np; j++)
		S[j] = Sjin[j];

	// Water Property
	PVTW.eval_all(0, P, data, cdata);
	OCP_DBL Pw0 = data[0];
	OCP_DBL bw0 = data[1];
	OCP_DBL cbw = data[2];
	OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0)); 
	OCP_DBL bwp = -cbw * bw0;

	Mu[2] = data[3];
	Rho[2] = Std_RhoW / bw;
	Xi[2] = 1 / (CONV1 * bw);
	Ni[2] = Vpore * S[2] * Xi[2];
	

	USI phasecae;

	if (1 - S[1] - S[2] < TINY) {
		if (S[1] < TINY)		phasecae = PHASE_W;		// case 1 : water, no oil, no gas
		else					phasecae = PHASE_GW;	// case 2 : water, gas, no oil
	}
	else if (S[1] < TINY)		phasecae = PHASE_OW;	// case 3 : water, oil, no gas
	else						phasecae = PHASE_OGW;	// case 4 : water, oil, gas

	switch (phasecae)
	{
	case PHASE_W:
	{
		// water
		PhaseExist[2] = true;
		S[0] = 0; S[1] = 0; S[2] = 1;
		Cij[8] = 1;

		// hypothetical Oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];

		// hypothetical Gas property
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);

		// total
		V[0] = 0;	V[1] = 0;	V[2] = CONV1 * Ni[2] * bw;
		Vf = V[2];
		Vfp = CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		

		break;
	}
	case PHASE_GW:
	{
		// water, gas
		PhaseExist[1] = true;		PhaseExist[2] = true;
		Cij[4] = 1;					Cij[8] = 1;

		// hypothetical Oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];

		// gas property	
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);
		OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[2];
		Xi[1] = 1 / bg / 1000;
		Rho[1] = Std_RhoG / bg;
		Ni[1] = Vpore * S[1] * Xi[1];

		V[0] = 0;
		V[1] = 1000 * Ni[1] * bg;   // Ni[0] = 0;
		V[2] = CONV1 * Ni[2] * bw;
		// total
		Vf = V[1] + V[2];
		S[0] = 0;	S[1] = V[1] / Vf;	  S[2] = V[2] / Vf;
		Vfp = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;

		break;
	}
	case PHASE_OW:
	{
		// water, oil
		PhaseExist[0] = true;		PhaseExist[2] = true;

		// oil property
		OCP_DBL Pbb = Pbbin;
		PVCO.eval_all(0, Pbb, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bosat = data[2];
		OCP_DBL muosat = data[3];
		OCP_DBL cbosat = data[4];
		OCP_DBL cmuosat = data[5];
		OCP_DBL bo = bosat * (1 - cbosat * (P - Pbb));
		OCP_DBL bop = -bosat * cbosat;
		OCP_DBL dBo_drs = bo / bosat * cdata[2] + bosat * (cdata[4] * (Pbb - P) + cbosat * cdata[0]);
		dBo_drs /= cdata[1];

		Ni[0] = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
		Ni[1] = Ni[0] * rs;
		Xi[0] = (1 + rs) / (CONV1 * bo);
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;
		Mu[0] = muosat * (1 + cmuosat * (P - Pbb));

		Cij[0] = Ni[0] / (Ni[0] + Ni[1]);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		// total
		V[0] = CONV1 * Ni[0] * bo;
		V[2] = CONV1 * Ni[2] * bw;
		Vf = V[0] + V[2];
		S[0] = V[0] / Vf;		S[1] = 0;		S[2] = V[2] / Vf;
		Vfp = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
		Vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
		Vfi[1] = CONV1 * dBo_drs;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	case PHASE_OGW:
	{
		PhaseExist.assign(3, true);

		// oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];
		OCP_DBL crs = cdata[1];
		OCP_DBL cbosat = cdata[2];

		Mu[0] = data[3];
		Ni[0] = Vpore * (1 - S[1] - S[2]) / (CONV1 * bo);
		Xi[0] = (1 + rs) / bo / CONV1;
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;
		
		// gas property
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);
		OCP_DBL cbg = cdata[1] * (CONV1 / 1000);
		Ni[1] = Vpore * S[1] / bg / 1000 + Ni[0] * rs;
		Xi[1] = 1 / data[1] / CONV1;
		Rho[1] = Std_RhoG / bg;
		Mu[1] = data[2];

		Cij[0] = 1 / (1 + rs);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		// total
		V[0] = CONV1 * Ni[0] * bo;
		V[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		V[2] = CONV1 * Ni[2] * bw;

		Vf = V[0] + V[1] + V[2];
		S[0] = V[0] / Vf;		S[1] = V[1] / Vf;		S[2] = V[2] / Vf;
		Vfp = CONV1 * Ni[0] * cbosat + 1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	}
}


void BOMixture::BOFlash_Ni_OGW(const OCP_DBL& Pin, const OCP_DBL* Niin)
{
	// initialize    Np = 3, Nc = 3
	PhaseExist.assign(Np, false);
	Cij.assign(Np * Nc, 0);   // 

	P = Pin;
	OCP_DBL NT = 0;
	for (USI i = 0; i < Nc; i++) {
		Ni[i] = Niin[i];
		NT += Ni[i];
	}


	// Water property
	PVTW.eval_all(0, P, data, cdata);
	OCP_DBL Pw0 = data[0];
	OCP_DBL bw0 = data[1];
	OCP_DBL cbw = data[2];
	OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0));
	OCP_DBL bwp = -cbw * bw0;

	Mu[2] = data[3];
	Xi[2] = 1 / (CONV1 * bw);
	Rho[2] = Std_RhoW / bw;
	
	USI phasecase;
	OCP_DBL Rs_sat = PVCO.eval(0, P, 1);

	if (Ni[0] < NT * TINY) {
		if (Ni[1] < Ni[0] * Rs_sat)		phasecase = PHASE_W;		// water, no oil, no gas
		else							phasecase = PHASE_GW;		// water, gas, no oil
	}
	else if (Ni[1] < Ni[0] * Rs_sat)	phasecase = PHASE_OW;		// water, oil, no gas
	else                                phasecase = PHASE_OGW;		// water, oil ,gas

	switch (phasecase)
	{
	case PHASE_W:
	{
		// water
		PhaseExist[2] = true;
		S[0] = 0; S[1] = 0; S[2] = 1;
		Cij[8] = 1;

		// hypothetical Oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];

		// hypothetical Gas property
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);

		// total
		V[0] = 0;	V[1] = 0;
		V[2] = CONV1 * Ni[2] * bw;
		Vf = V[2];
		Vfp = CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	case PHASE_GW:
	{
		// water, gas
		PhaseExist[1] = true;		PhaseExist[2] = true;
		Cij[4] = 1;					Cij[8] = 1;

		// hypothetical Oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];

		// gas property
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);
		OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[2];
		Xi[1] = 1 / bg / 1000;
		Rho[1] = Std_RhoG / bg;

		// total
		V[0] = 0;
		V[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		V[2] = CONV1 * Ni[2] * bw;
		if (V[1] < 0)
			V[1] = 0;
		Vf = V[1] + V[2];
		S[0] = 0;	S[1] = V[1] / Vf;	S[2] = V[2] / Vf;
		Vfp = 1000 * Ni[1] * cbg + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	case PHASE_OW:
	{
		// water, oil
		PhaseExist[0] = true;			PhaseExist[2] = true;
		Cij[0] = Ni[0] / (Ni[0] + Ni[1]);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		// oil property
		OCP_DBL rs = Ni[1] / Ni[0];
		PVCO.eval_all(1, rs, data, cdata);
		OCP_DBL pbb = data[0];
		OCP_DBL bosat = data[2];
		OCP_DBL muosat = data[3];
		OCP_DBL cbosat = data[4];
		OCP_DBL cmuosat = data[5];
		OCP_DBL bo = bosat * (1 - cbosat * (P - pbb));
		OCP_DBL bop = -bosat * cbosat;
		OCP_DBL dBo_drs = bo / bosat * cdata[2] + bosat * (cdata[4] * (pbb - P) + cbosat * cdata[0]);

		Mu[0] = muosat * (1 + cmuosat * (P - pbb));
		Xi[0] = (1 + rs) / (CONV1 * bo);
		Rho[0] = (Std_RhoO + (1000/CONV1) * rs * Std_RhoG) / bo;

		// total
		V[0] = CONV1 * Ni[0] * bo;
		V[2] = CONV1 * Ni[2] * bw;
		Vf = V[0] + V[2];
		S[0] = V[0] / Vf; S[1] = 0; S[2] = V[2] / Vf;
		Vfp = CONV1 * (Ni[0] * bop + Ni[2] * bwp);
		Vfi[0] = CONV1 * (bo - dBo_drs * (Ni[1] / Ni[0]));
		Vfi[1] = CONV1 * dBo_drs;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	case PHASE_OGW:
	{
		PhaseExist.assign(3, true);

		// oil property
		PVCO.eval_all(0, P, data, cdata);
		OCP_DBL rs = data[1];
		OCP_DBL bo = data[2];
		OCP_DBL crs = cdata[1];
		OCP_DBL cbosat = cdata[2];

		Mu[0] = data[3];
		Xi[0] = (1 + rs) / bo / CONV1;
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;

		// gas property
		PVDG.eval_all(0, P, data, cdata);
		OCP_DBL bg = data[1] * (CONV1 / 1000);
		OCP_DBL cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[2];
		Xi[1] = 1 / data[1] / CONV1;
		Rho[1] = Std_RhoG / bg;

		// total
		Cij[0] = 1 / (1 + rs);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		V[0] = CONV1 * Ni[0] * bo;
		V[1] = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		V[2] = CONV1 * Ni[2] * bw;
		Vf = V[0] + V[1] + V[2];
		S[0] = V[0] / Vf;  S[1] = V[1] / Vf;  S[2] = V[2] / Vf;
		Vfp = CONV1 * Ni[0] * cbosat + 1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	}
}


OCP_DBL BOMixture::xiPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin)
{
	if (Ziin[1] > 1 - TINY) {
		// inj fluid is gas
		OCP_DBL bg = PVDG.eval(0, Pin, 1);
		OCP_DBL xig = (1 / CONV1) / bg;
		return xig;
	}
	else if (Ziin[2] > 1 - TINY) {
		// inj fluid is water

		PVTW.eval_all(0, Pin, data, cdata);
		OCP_DBL Pw0 = data[0];
		OCP_DBL bw0 = data[1];
		OCP_DBL cbw = data[2];
		OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0));
		OCP_DBL xiw = (1 / CONV1) / bw;
		return xiw;
	}
	else {
		ERRORcheck("Wrong Zi!");
		exit(0);
	}
}

OCP_DBL BOMixture::rhoPhase_OGW(const OCP_DBL& Pin, const OCP_DBL* Ziin)
{
	if (Ziin[1] > 1 - TINY) {
		// inj fluid is gas
		OCP_DBL bg = PVDG.eval(0, Pin, 1);
		OCP_DBL rhog = (1000 / CONV1) * Std_RhoG / bg;
		return rhog;
	}
	else if (Ziin[2] > 1 - TINY) {
		// inj fluid is water

		PVTW.eval_all(0, Pin, data, cdata);
		OCP_DBL Pw0 = data[0];
		OCP_DBL bw0 = data[1];
		OCP_DBL cbw = data[2];
		OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0));
		OCP_DBL rhow = Std_RhoW / bw;
		return rhow;
	}
	else {
		ERRORcheck("Wrong Zi!");
		exit(0);
	}
}


OCP_DBL BOMixture::gammaPhaseO_OGW(const OCP_DBL& Pin, const OCP_DBL& Pbbin)
{

	PVCO.eval_all(0, Pbbin, data, cdata);
	OCP_DBL rs = data[1];
	OCP_DBL bosat = data[2];
	OCP_DBL cbosat = data[4];
	OCP_DBL bo = bosat * (1 - cbosat * (Pin - Pbbin));
	OCP_DBL gammaO = (Std_GammaO + (1000 / CONV1) * rs * Std_GammaG) / bo;
	
	return gammaO;
}
