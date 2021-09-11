#include "BOMixture.hpp"

void BOMixture::BOFlash_Sj_OGW(const double Pin, const double Pbbin, const double* Sjin, double Vpore)
{
	// initialize    Np = 3, Nc = 3
	PhaseExist.assign(Np, false);
	Cij.assign(Np * Nc, 0);   // 
	Ni.assign(Nc, 0);

	P = Pin;
	for (int j = 0; j < Np; j++)
		S[j] = Sjin[j];

	std::vector<double>		data(6, 0);
	std::vector<double>		cdata(6, 0);
	// Water Property
	PVTW.eval_all(0, P, data, cdata);
	double Pw0 = data[0];
	double bw0 = data[1];
	double cbw = data[2];
	double bw = bw0 * (1 - cbw * (P - Pw0)); 
	double bwp = -cbw * bw0;

	Mu[2] = data[3];
	Rho[2] = Std_RhoW / bw;
	Xi[2] = 1 / (CONV1 * bw);
	Ni[2] = Vpore * S[2] * Xi[2];
	

	int phasecae;

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
		double rs = data[1];
		double bo = data[2];

		// hypothetical Gas property
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);

		// total
		double vw = CONV1 * Ni[2] * bw;
		Vf = vw;
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
		double rs = data[1];
		double bo = data[2];

		// gas property	
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);
		double cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[2];
		Xi[1] = 1 / bg / 1000;
		Rho[1] = Std_RhoG / bg;
		Ni[1] = Vpore * S[1] * Xi[1];

		double vg = 1000 * Ni[1] * bg;   // Ni[0] = 0;
		double vw = CONV1 * Ni[2] * bw;
		// total
		Vf = vg + vw;
		S[0] = 0;	S[1] = vg / Vf;	  S[2] = vw / Vf;
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
		double Pbb = Pbbin;
		PVCO.eval_all(0, Pbb, data, cdata);
		double rs = data[1];
		double bosat = data[2];
		double muosat = data[3];
		double cbosat = data[4];
		double cmuosat = data[5];
		double bo = bosat * (1 - cbosat * (P - Pbb));
		double bop = -bosat * cbosat;
		double dBo_drs = bo / bosat * cdata[2] + bosat * (cdata[4] * (Pbb - P) + cbosat * cdata[0]);
		dBo_drs /= cdata[1];

		Ni[0] = Vpore * S[0] / (CONV1 * bo);
		Ni[1] = Ni[0] * rs;
		Xi[0] = (1 + rs) / (CONV1 * bo);
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;
		Mu[0] = muosat * (1 + cmuosat * (P - Pbb));

		Cij[0] = Ni[0] / (Ni[0] + Ni[1]);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		// total
		double vo = CONV1 * Ni[0] * bo;
		double vw = CONV1 * Ni[2] * bw;
		Vf = vo + vw;
		S[0] = vo / Vf;		S[1] = 0;		S[2] = vw / Vf;
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
		double rs = data[1];
		double bo = data[2];
		double crs = cdata[1];
		double cbosat = cdata[2];

		Mu[0] = data[3];
		Ni[0] = Vpore * S[0] / (CONV1 * bo);
		Xi[0] = (1 + rs) / bo / CONV1;
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;
		
		// gas property
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);
		double cbg = cdata[1] * (CONV1 / 1000);
		Ni[1] = Vpore * S[1] / bg / 1000 + Ni[0] * rs;
		Xi[1] = 1 / data[1] / CONV1;
		Rho[1] = Std_RhoG / bg;
		Mu[1] = data[2];

		Cij[0] = 1 / (1 + rs);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		// total
		double vo = CONV1 * Ni[0] * bo;
		double vg = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		double vw = CONV1 * Ni[2] * bw;

		Vf = vo + vg + vw;
		S[0] = vo / Vf;		S[1] = vg / Vf;		S[2] = vw / Vf;
		Vfp = CONV1 * Ni[0] * cbosat + 1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	}
}



void BOMixture::BOFlash_Ni_OGW(const double Pin, const double* Niin)
{
	// initialize    Np = 3, Nc = 3
	PhaseExist.assign(Np, false);
	Cij.assign(Np * Nc, 0);   // 

	P = Pin;
	double NT = 0;
	for (int i = 0; i < Nc; i++) {
		Ni[i] = Niin[i];
		NT += Ni[i];
	}

	std::vector<double>		data;
	std::vector<double>		cdata;

	// Water property
	PVTW.eval_all(0, P, data, cdata);
	double Pw0 = data[0];
	double bw0 = data[1];
	double cbw = data[2];
	double bw = bw0 * (1 - cbw * (P - Pw0));
	double bwp = -cbw * bw0;

	Mu[2] = data[3];
	Xi[2] = 1 / (CONV1 * bw);
	Rho[2] = Std_RhoW / bw;
	
	int phasecase;
	double Rs_sat = PVCO.eval(0, P, 1);

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
		double rs = data[1];
		double bo = data[2];

		// hypothetical Gas property
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);

		// total
		double vw = CONV1 * Ni[2] * bw;
		Vf = vw;
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
		double rs = data[1];
		double bo = data[2];

		// gas property
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);
		double cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[2];
		Xi[1] = 1 / bg / 1000;
		Rho[1] = Std_RhoG / bg;

		// total
		double vg = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		double vw = CONV1 * Ni[2] * bw;
		if (vg < 0)
			vg = 0;
		Vf = vg + vw;
		S[0] = 1; S[1] = vg / Vf; S[2] = vw / Vf;
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
		double rs = Ni[1] / Ni[0];
		PVCO.eval_all(1, rs, data, cdata);
		double pbb = data[0];
		double bosat = data[2];
		double muosat = data[3];
		double cbosat = data[4];
		double cmuosat = data[5];
		double bo = bosat * (1 - cbosat * (P - pbb));
		double bop = -bosat * cbosat;
		double dBo_drs = bo / bosat * cdata[2] + bosat * (cdata[4] * (pbb - P) + cbosat * cdata[0]);

		Mu[0] = muosat * (1 + cmuosat * (P - pbb));
		Xi[0] = (1 + rs) / (CONV1 * bo);
		Rho[0] = (Std_RhoO + (1000/CONV1) * rs * Std_RhoG) / bo;

		// total
		double vo = CONV1 * Ni[0] * bo;
		double vw = CONV1 * Ni[2] * bw;
		Vf = vo + vw;
		S[0] = vo / Vf; S[1] = 0; S[2] = vw / Vf;
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
		double rs = data[1];
		double bo = data[2];
		double crs = cdata[1];
		double cbosat = cdata[2];

		Mu[0] = data[3];
		Xi[0] = (1 + rs) / bo / CONV1;
		Rho[0] = (Std_RhoO + (1000 / CONV1) * rs * Std_RhoG) / bo;

		// gas property
		PVDG.eval_all(0, P, data, cdata);
		double bg = data[1] * (CONV1 / 1000);
		double cbg = cdata[1] * (CONV1 / 1000);

		Mu[1] = data[3];
		Xi[1] = 1 / data[1] / CONV1;
		Rho[1] = Std_RhoO / bg;

		// total
		Cij[0] = 1 / (1 + rs);
		Cij[1] = 1 - Cij[0];
		Cij[4] = 1;
		Cij[8] = 1;

		double vo = CONV1 * Ni[0] * bo;
		double vg = 1000 * (Ni[1] - rs * Ni[0]) * bg;
		double vw = CONV1 * Ni[2] * bw;
		Vf = vo + vg + vw;
		S[0] = vo / Vf;  S[1] = vg / Vf;  S[2] = vw / Vf;
		Vfp = CONV1 * Ni[0] * cbosat + 1000 * (-crs * Ni[0] * bg + (Ni[1] - rs * Ni[0]) * cbg) + CONV1 * Ni[2] * bwp;
		Vfi[0] = CONV1 * bo - 1000 * rs * bg;
		Vfi[1] = 1000 * bg;
		Vfi[2] = CONV1 * bw;
		
		break;
	}
	}
}


double BOMixture::gammaPhaseO_OGW(double Pin, double Pbbin)
{
	std::vector<double>		data(6, 0);
	std::vector<double>		cdata(6, 0);

	PVCO.eval_all(0, Pbbin, data, cdata);
	double rs = data[1];
	double bosat = data[2];
	double cbosat = data[4];
	double bo = bosat * (1 - cbosat * (Pin - Pbbin));
	double gammaO = (Std_GammaO + (1000 / CONV1) * rs * Std_GammaG) / bo;
	
	return gammaO;
}
