#include "MixtureBO.hpp"


void BOMixture::BOFlash_Sj_OW(const OCP_DBL& Pin, const OCP_DBL* Sjin,
	const OCP_DBL& Vpore)
{
	S[1] = Sjin[1];
	// Water Properties
	PVTW.Eval_All(0, P, data, cdata);
	OCP_DBL Pw0 = data[0];
	OCP_DBL bw0 = data[1];
	OCP_DBL cbw = data[2];
	OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0));
	OCP_DBL bwp = -cbw * bw0;

	mu[1] = data[3];
	xi[1] = 1 / (bw * CONV1);
	rho[1] = std_RhoW / bw;

	Ni[1] = Vpore * S[1] * xi[1];

	// Oil Properties
	PVDO.Eval_All(0, P, data, cdata);
	OCP_DBL bo = data[1];
	OCP_DBL bop = cdata[1];

	mu[0] = data[2];
	xi[0] = 1 / (CONV1 * bo);
	rho[0] = std_RhoO / bo;
	Ni[0] = Vpore * (1 - S[1]) * xi[0];

	cij[0 * 2 + 0] = 1;  cij[0 * 2 + 1] = 0;
	cij[1 * 2 + 0] = 0;  cij[1 * 2 + 1] = 1;

	v[0] = CONV1 * Ni[0] * bo;
	v[1] = CONV1 * Ni[1] * bw;
	vf = v[0] + v[1];
	S[0] = v[0] / vf;
	S[1] = v[1] / vf;
	vfi[0] = CONV1 * bo;
	vfi[1] = CONV1 * bw;
	vfp = CONV1 * Ni[0] * bop + CONV1 * Ni[1] * bwp;
}


void BOMixture::BOFlash_Ni_OW_Deriv(const OCP_DBL& Pin, const OCP_DBL* Niin)
{
	phaseExist.assign(numPhase, true);
	dSec_dPri.assign(dSec_dPri.size(), 0);

	P = Pin;
	Ni[0] = Niin[0];
	Ni[1] = Niin[1];
	Nt = Ni[0] + Ni[1];

	// Water Properties
	PVTW.Eval_All(0, P, data, cdata);
	OCP_DBL Pw0 = data[0];
	OCP_DBL bw0 = data[1];
	OCP_DBL cbw = data[2];
	OCP_DBL bw = bw0 * (1 - cbw * (P - Pw0));
	OCP_DBL bwp = -cbw * bw0;

	mu[1] = data[3];
	xi[1] = 1 / (bw * CONV1);
	rho[1] = std_RhoW / bw;

	muP[1] = cdata[3];
	xiP[1] = -bwp / (bw * bw * CONV1);
	rhoP[1] = CONV1 * xiP[1] * std_RhoW;
	
	// Oil Properties
	PVDO.Eval_All(0, P, data, cdata);
	OCP_DBL bo = data[1];     
	OCP_DBL bop = cdata[1];

	mu[0] = data[2];
	xi[0] = 1 / (CONV1 * bo);
	rho[0] = std_RhoO / bo;

	muP[0] = cdata[1];
	xiP[0] = -xi[0] * bop / bo;
	rhoP[0] = -rho[0] * bop / bo;

	cij[0 * 2 + 0] = 1;  cij[0 * 2 + 1] = 0;
	cij[1 * 2 + 0] = 0;  cij[1 * 2 + 1] = 1;

	v[0] = CONV1 * Ni[0] * bo;
	v[1] = CONV1 * Ni[1] * bw;
	vf = v[0] + v[1];
	S[0] = v[0] / vf;
	S[1] = v[1] / vf;
	vfi[0] = CONV1 * bo;
	vfi[1] = CONV1 * bw;
	vfp = CONV1 * Ni[0] * bop + CONV1 * Ni[1] * bwp;

	dSec_dPri[0] = (CONV1 * Ni[0] * bop - S[0] * vfp) / vf; // dSo / dP
	dSec_dPri[1] = (CONV1 * bo - S[0] * vfi[0]) / vf; // dSo / dNo
	dSec_dPri[2] = -S[0] * vfi[1] / vf; // dSo / dNw

	dSec_dPri[3] = (CONV1 * Ni[1] * bwp - S[1] * vfp) / vf; // dSw / dP
	dSec_dPri[4] = -S[1] * vfi[0] / vf;  // dSw / dNo
	dSec_dPri[5] = (CONV1 * bw - S[1] * vfi[1]) / vf; // dSw / dNw

}