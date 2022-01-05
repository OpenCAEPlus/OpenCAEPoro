#include"MixtureComp.hpp"



COMP::COMP(const vector<string>& comp)
{
	name = comp[0];
	Pc = stod(comp[1]);
	Tc = stod(comp[2]);
	acf = stod(comp[3]);
	MW = stod(comp[4]);
	VcMW = stod(comp[5]);
	OmegaA = stod(comp[6]);
	OmegaB = stod(comp[7]);
	Vshift = stod(comp[8]);
}


MixtureComp::MixtureComp(const EoSparam& param)
{

	numPhase = param.numPhase + 1;
	numCom = param.numComp + 1;
	Allocate();
	zi.resize(numCom);

	NC = param.numComp;
	NPmax = param.numPhase;

	comp.resize(NC);
	for (USI i = 0; i < NC; i++) {
		comp[i] = COMP(param.COM[i]);
	}

	USI len = NC * NC;
	USI count = 0;
	BIP.resize(len);
	for (USI i = 0; i < NC; i++) {
		for (USI j = 0; j < NC; j++) {
			BIP[count] = param.BIP[i][j];
			count++;
		}
	}

	EoSctrl.SSMsta.maxIt = stoi(param.SSMparamSTA[0]);
	EoSctrl.SSMsta.tol = stod(param.SSMparamSTA[1]);
	EoSctrl.SSMsta.tol2 = EoSctrl.SSMsta.tol * EoSctrl.SSMsta.tol;

	EoSctrl.NRsta.maxIt = stoi(param.NRparamSTA[0]);
	EoSctrl.NRsta.tol = stod(param.NRparamSTA[1]);
	EoSctrl.NRsta.tol2 = EoSctrl.NRsta.tol * EoSctrl.NRsta.tol;

	EoSctrl.SSMsp.maxIt = stoi(param.SSMparamSP[0]);
	EoSctrl.SSMsp.tol = stod(param.SSMparamSP[1]);
	EoSctrl.SSMsp.tol2 = EoSctrl.SSMsp.tol * EoSctrl.SSMsp.tol;

	EoSctrl.NRsp.maxIt = stoi(param.NRparamSP[0]);
	EoSctrl.NRsp.tol = stod(param.NRparamSP[1]);
	EoSctrl.NRsp.tol2 = EoSctrl.NRsp.tol * EoSctrl.NRsp.tol;

	EoSctrl.RR.maxIt = stoi(param.RRparam[0]);
	EoSctrl.RR.tol = stod(param.RRparam[1]);
	EoSctrl.RR.tol2 = EoSctrl.RR.tol * EoSctrl.RR.tol;

	AllocateEoS();
	AllocatePhase();
	AllocateMethod();
}


void MixtureComp::AllocateEoS()
{
	// Allocate Memoery for EoS variables
	Ai.resize(NC);
	Bi.resize(NC);
	Aj.resize(NPmax);
	Bj.resize(NPmax);
	Zj.resize(NPmax);
	Ztmp.resize(3);
}


void MixtureComp::SolEoS(OCP_DBL& ZjT, const OCP_DBL& AjT, const OCP_DBL& BjT) const
{
	const OCP_DBL aj = AjT;
	const OCP_DBL bj = BjT;

	const OCP_DBL a = (delta1 + delta2 - 1) * bj - 1;
	const OCP_DBL b = (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1));
	const OCP_DBL c = -(aj * bj + delta1 * delta2 * bj * bj * (bj + 1));

	USI flag = CubicRoot(a, b, c, true); // True with NT
	if (flag == 1) {
		ZjT = Ztmp[0];
	}
	else {
		OCP_DBL zj1 = Ztmp[0];
		OCP_DBL zj2 = Ztmp[2];
		OCP_DBL dG = (zj2 - zj1) + log((zj1 - bj) / (zj2 - bj))
			- aj / (bj * (delta2 - delta1))
			* log((zj1 + delta1 * bj) * (zj2 + delta2 * bj) / ((zj1 + delta2 * bj) * (zj2 + delta1 * bj)));
		if (dG > 0)
			ZjT = zj1;
		else
			ZjT = zj2;
	}
}


void MixtureComp::CalAiBi()
{
	// Calculate Ai, Bi
	OCP_DBL acf, mwi;
	OCP_DBL Pri, Tri;

	for (USI i = 0; i < NC; i++) {
		acf = comp[i].acf;
		// PR
		if (acf <= 0.49) {
			mwi = 0.37464 + 1.54226 * acf - 0.26992 * pow(acf, 2);
		}
		else {
			mwi = 0.379642 + 1.48503 * acf - 0.164423 * pow(acf, 2) + 0.016667 * pow(acf, 3);
		}

		Pri = P / comp[i].Pc;
		Tri = T / comp[i].Tc;
		Ai[i] = comp[i].OmegaA * Pri / pow(Tri, 2) * pow((1 + mwi * (1 - sqrt(Tri))), 2);
		Bi[i] = comp[i].OmegaB * Pri / Tri;
	}
}


void MixtureComp::CalAjBj(OCP_DBL& AjT, OCP_DBL& BjT, const vector<OCP_DBL>& xj) const
{
	AjT = 0;
	BjT = 0;

	for (USI i1 = 0; i1 < NC; i1++) {
		BjT += Bi[i1] * xj[i1];
		AjT += xj[i1] * xj[i1] * Ai[i1] * (1 - BIP[i1 * NC + i1]);

		for (USI i2 = 0; i2 < i1; i2++) {
			AjT += 2 * xj[i1] * xj[i2] * sqrt(Ai[i1] * Ai[i2]) * (1 - BIP[i1 * NC + i2]);
		}
	}
}

void MixtureComp::AllocatePhase()
{
	// Allocate Memoery for Phase variables
	nuj.resize(NPmax);
	
	x.resize(NPmax);
	phi.resize(NPmax);
	fug.resize(NPmax);
	n.resize(NPmax);
	for (USI j = 0; j < NPmax; j++) {
		x[j].resize(NC);
		phi[j].resize(NC);
		fug[j].resize(NC);
		n[j].resize(NC);
	}
	ln = n;
	xiC.resize(NPmax);
	MW.resize(NPmax);
	phaseLabel.resize(NPmax);
	
}


void MixtureComp::CalFugPhi(vector<OCP_DBL>& phiT, vector<OCP_DBL>& fugT, const vector<OCP_DBL>& xj)
{

	OCP_DBL aj, bj, zj;
	OCP_DBL tmp;
	CalAjBj(aj, bj, xj);
	SolEoS(zj, aj, bj);

	for (USI i = 0; i < NC; i++) {
		tmp = 0;
		for (int k = 0; k < NC; k++) {
			tmp += 2 * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
		}
		phiT[i] = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) - aj / (delta1 - delta2) / bj * (tmp / aj - Bi[i] / bj) * log((zj + delta1 * bj) / (zj + delta2 * bj)));
		fugT[i] = phiT[i] * xj[i] * P;
	}

	Asta = aj;
	Bsta = bj;
	Zsta = zj;
}

void MixtureComp::CalFugPhiAll()
{
	OCP_DBL tmp;

	for (USI j = 0; j < NP; j++) {
		vector<OCP_DBL>& xj = x[j];
		vector<OCP_DBL>& phiT = phi[j];
		vector<OCP_DBL>& fugT = fug[j];
		OCP_DBL& aj = Aj[j];
		OCP_DBL& bj = Bj[j];
		OCP_DBL& zj = Zj[j];

		CalAjBj(aj, bj, xj);
		SolEoS(zj, aj, bj);

		for (USI i = 0; i < NC; i++) {
			tmp = 0;
			for (int k = 0; k < NC; k++) {
				tmp += 2 * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
			}
			phiT[i] = exp(Bi[i] / bj * (zj - 1) - log(zj - bj) - aj / (delta1 - delta2) / bj * (tmp / aj - Bi[i] / bj) * log((zj + delta1 * bj) / (zj + delta2 * bj)));
			fugT[i] = phiT[i] * xj[i] * P;
		}
	}
}


void MixtureComp::CalMW()
{
	// Calculate Molecular Weight of phase
	MW.assign(NP, 0);
	for (USI j = 0; j < NP; j++) {
		for (USI i = 0; i < NC; i++) {
			MW[j] += x[j][i] * comp[i].MW;
		}
	}
}


void MixtureComp::CalXiRho()
{
	OCP_DBL tmp;
	for (USI j = 0; j < NP; j++){

		vector<OCP_DBL>& xj = x[j];
		tmp = Zj[j] * GAS_CONSTANT * T / P;
		for (USI i = 0; i < NC; i++) {
			tmp -= xj[i] * comp[i].Vshift;
		}
		xiC[j] = 1 / tmp;
		rhoC[j] = MW[j] * xiC[j];
	}
}


USI MixtureComp::FindMWmax()
{
	// find the phase id with the highest Molecular Weight
	USI tmpId = 0;
	OCP_DBL tmpMW = MW[0];
	for (USI j = 1; j < NP; j++) {
		if (tmpMW < MW[j]) {
			tmpMW = MW[j];
			tmpId = j;
		}
	}
	return tmpId;
}

void MixtureComp::x2n()
{
	// Total moles are supposed to be 1
	for (USI j = 0; j < NP; j++) {

		// nuj[j] = fabs(nuj[j]);
		for (USI i = 0; i < NC; i++) {
			n[j][i] = nuj[j] * x[j][i];
		}
	}

}

void MixtureComp::PrintX()
{
	for (USI j = 0; j < NP; j++) {
		for (USI i = 0; i < NC; i++) {
			cout << x[j][i] << "   ";
		}
		cout << endl;
	}
	cout << "----------------------" << endl;
}


void MixtureComp::AllocateMethod()
{
	Kw.resize(4);
	for (USI i = 0; i < 4; i++) {
		Kw[i].resize(NC);
	}
	Ks.resize(NPmax - 1);
	for (USI i = 0; i < NPmax - 1; i++) {
		Ks[i].resize(NC);
	}
	phiSta.resize(NC);
	fugSta.resize(NC);

	Y.resize(NC);
	di.resize(NC);
	resSTA.resize(NC);
	
	JmatSTA.resize(NC * NC);
	Ax.resize(NC);
	Bx.resize(NC);
	Zx.resize(NC);
	resRR.resize(NPmax - 1);
	resSP.resize(NC * NPmax);
	JmatSP.resize(NC * NC * NPmax * NPmax);
	fugX.resize(NPmax);
	fugN.resize(NPmax);
	for (USI j = 0; j < NPmax; j++) {
		fugX[j].resize(NC * NC);
		fugN[j].resize(NC * NC);
	}
	An.resize(NC);
	Bn.resize(NC);
	Zn.resize(NC);
	pivot.resize(NC * NPmax, 1);
}

void MixtureComp::CalKwilson()
{
	for (USI i = 0; i < NC; i++) {
		Kw[0][i] = (comp[i].Pc / P) * exp(5.373 * (1 + comp[i].acf) * (1 - comp[i].Tc / T));
		Kw[1][i] = 1 / Kw[0][i];
		Kw[2][i] = pow(Kw[0][i], 1.0 / 3);
		Kw[3][i] = pow(Kw[1][i], 1.0 / 3);
	}
}


void MixtureComp::PhaseEquilibrium()
{
	NP = 1;
	x[0] = zi;
	CalAiBi();
	CalKwilson();
	while (!PhaseStable()) {
		NP++;
		PhaseSplit();
		if (NP == NPmax)
			break;
	}
}


bool MixtureComp::PhaseStable()
{
	if (NP == 1) {
		testPId = 0;
	}
	else {
		CalMW();
		testPId = FindMWmax();
	}

	// Test if a phase is stable, if stable return true, else return false
	bool flag = StableSSM(testPId);
	return flag;
}


bool MixtureComp::StableSSM(const USI& Id)
{
	const vector<OCP_DBL>& xj = x[Id];
	CalFugPhi(phi[Id], fug[Id], xj);

	for (USI i = 0; i < NC; i++) {
		di[i] = phi[Id][i] * xj[i];
	}

	USI lenK = Kw.size();
	OCP_DBL Stol = EoSctrl.SSMsta.tol2;
	OCP_DBL maxIt = EoSctrl.SSMsta.maxIt;
	OCP_DBL Se;
	bool flag;
	USI iter;
	for (USI k = 0; k < 1; k++) {

		Yt = 0;
		for (USI i = 0; i < NC; i++) {
			Y[i] = xj[i] / Kw[k][i];
			Yt += Y[i];
		}
		for (USI i = 0; i < NC; i++) {
			Y[i] /= Yt;
		}
		CalFugPhi(phiSta, fugSta, Y);

		Se = 0;
		for (USI i = 0; i < NC; i++) {
			Se += pow(log(fugSta[i] / fug[Id][i] * Yt), 2);
		}

		flag = true;
		iter = 0;
		while (Se > Stol)
		{
			Yt = 0;
			for (USI i = 0; i < NC; i++) {
				Y[i] = di[i] / phiSta[i];
				Yt += Y[i];
			}
			Dscalar(NC, 1 / Yt, &Y[0]);
			CalFugPhi(phiSta, fugSta, Y);
			Se = 0;
			for (USI i = 0; i < NC; i++) {
				Se += pow(log(fugSta[i] / fug[Id][i] * Yt), 2);
			}

			iter++;
			if (iter > maxIt) {
				flag = false;
				break;
			}
		}

		// StableNR(Id);

		if (flag && Yt > 1 + 1E-3) {
			return false;
		}
	}

	/*if (!flag) {
		OCP_WARNING("SSM not converged in Stability Analysis");
	}*/

	return true;
}


bool MixtureComp::StableNR(const USI& Id)
{

#ifdef DEBUG
	cout << endl << "Stable NR Begins !" << endl << endl;
#endif // DEBUG

	for (USI i = 0; i < NC; i++) {
		resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
	}

#ifdef DEBUG
	PrintDX(NC, &resSTA[0]);
#endif // DEBUG

	USI maxIt = EoSctrl.NRsp.maxIt;
	OCP_DBL Stol = EoSctrl.NRsta.tol;
	OCP_DBL Se = Dnorm2(NC, &resSTA[0]);
	OCP_DBL alpha = 1;
	USI iter = 0;
	while (Se > Stol) {

		CalFugXSTA();
		AssembleJmatSTA();
		LUSolve(NC, &JmatSTA[0], &resSTA[0], &pivot[0]);

		Dscalar(NC, Yt, &Y[0]);
		Daxpy(NC, alpha, &resSTA[0], &Y[0]);
		Yt = 0;
		for (USI i = 0; i < NC; i++) {
			Yt += Y[i];
		}
		Dscalar(NC, 1 / Yt, &Y[0]);

		CalFugPhi(phiSta, fugSta, Y);
		for (USI i = 0; i < NC; i++) {
			resSTA[i] = log(fug[Id][i] / (fugSta[i] * Yt));
		}
		Se = Dnorm2(NC, &resSTA[0]);
		iter++;
		if (iter > maxIt) {
			break;
		}

#ifdef DEBUG
		PrintDX(NC, &resSTA[0]);
#endif // DEBUG
	}
}


void MixtureComp::CalFugXSTA()
{
	vector<OCP_DBL>& fugx = fugX[0];
	OCP_DBL aj = Asta;
	OCP_DBL bj = Bsta;
	OCP_DBL zj = Zsta;
	OCP_DBL tmp = 0;

	Bx = Bi;
	for (USI i = 0; i < NC; i++) {
		tmp = 0;
		for (USI k = 0; k < NC; k++) {
			tmp += Y[k] * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
		}
		Ax[i] = 2 * tmp;
		Zx[i] = ((bj - zj) * Ax[i] + ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj))
			+ ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj
			- (delta1 + delta2 - 1) * zj * zj) * Bx[i]) / (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj
				+ (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
	}

	OCP_DBL C, D, E, G;
	OCP_DBL Cxk, Dxk, Exk, Gxk;
	OCP_DBL* phiXT;
	G = (zj + delta1 * bj) / (zj + delta2 * bj);

	for (USI i = 0; i < NC; i++) {
		C = Y[i] * P / (zj - bj);
		// C = 1 / (zj - bj);
		// D = Bx[i] * (zj - 1) / bj;
		E = -aj / ((delta1 - delta2) * bj) * (Ax[i] / aj - Bx[i] / bj);

		for (USI k = 0; k < NC; k++) {
			// Cxk = -Y[i] * (Zx[k] - Bx[k]) / ((zj - bj) * (zj - bj));
			Cxk = ((zj - bj) * delta(i, k) - Y[i] * (Zx[k] - Bx[k])) * P / ((zj - bj) * (zj - bj));
			Dxk = Bx[i] / bj * (Zx[k] - Bx[k] * (zj - 1) / bj);
			Exk = (Ax[k] * bj - aj * Bx[k]) / (bj * bj) * (Ax[i] / aj - Bx[i] / bj) + aj / bj * (2 * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]) / aj -
				Ax[k] * Ax[i] / (aj * aj) + Bx[i] * Bx[k] / (bj * bj));
			Exk /= -(delta1 - delta2);
			Gxk = (delta1 - delta2) / (zj + delta2 * bj) / (zj + delta2 * bj) * (zj * Bx[k] - Zx[k] * bj);
			fugx[i * NC + k] = 1 / C * Cxk + Dxk + Exk * log(G) + E / G * Gxk;
		}
	}
	//cout << "PhiX" << endl;
	//for (USI i = 0; i < NC; i++) {
	//	for (USI k = 0; k < NC; k++) {
	//		cout << fugx[i * NC + k] << "      ";
	//	}
	//	cout << endl;
	//}
}


void MixtureComp::AssembleJmatSTA()
{
	vector<OCP_DBL>& fugx = fugX[0];
	JmatSTA.assign(JmatSTA.size(), 0);
	OCP_DBL tmp;
	for (USI i = 0; i < NC; i++) {

		tmp = 0;
		for (USI i1 = 0; i1 < NC; i1++) {
			tmp += Y[i1] * fugx[i * NC + i1];
		}

		for (USI k = 0; k < NC; k++) {
			// Symmetric
			JmatSTA[i * NC + k] = (fugx[i * NC + k] - tmp + 1) / Yt;
		}
	}
	//cout << "JmatSTA" << endl;
	//for (USI i = 0; i < NC; i++) {
	//	for (USI k = 0; k < NC; k++) {
	//		cout << JmatSTA[k * NC + i] << "      ";
	//	}
	//	cout << endl;
	//}
}



void MixtureComp::PhaseSplit()
{
	SplitSSM(false);
	bool flag;
	while (!SplitNR()) {
		flag = SplitSSM(true);
		SplitNR();
		if (flag)
			break;
	}
}


bool MixtureComp::SplitSSM(const bool& flag)
{
#ifdef DEBUG
	cout << "SSMSP Begins!" << endl;
#endif

	if (NP == 2) {
		return SplitSSM2(flag);
	}
	else {
		return SplitSSM3(flag);
	}
}


bool MixtureComp::SplitSSM2(const bool& flag)
{
	// NP = 2 in this case
	// Ks is very IMPORTANT!
	// flag = true : Restart SSM
	// flag = false : New SSM
	OCP_DBL Se = 1;
	OCP_DBL Stol = EoSctrl.SSMsp.tol2;
	USI maxIt = EoSctrl.SSMsp.maxIt;
	if (Yt < 1.1) {
		Stol *= 1E-1;
		maxIt *= 2;
	}

	if (!flag) {
		for (USI i = 0; i < NC; i++) {
			Ks[NP - 2][i] = phi[testPId][i] / phiSta[i];
			// Ks[NP - 2][i] = phiSta[i] / phi[testPId][i];
		}
		//for (USI i = 0; i < NC; i++) {
		//	Se += pow(fugSta[i] / fug[testPId][i] - 1, 2);
		//}
	}
	else {
		Stol *= 1E-1;
		maxIt *= 2;
	}


#ifdef DEBUG
	PrintX();
#endif // DEBUG

	USI iter = 0;
	while (Se > Stol) {

		RachfordRice2();
		UpdateXRR();
		CalFugPhiAll();
		Se = 0;
		for (USI i = 0; i < NC; i++) {
			Se += pow(fug[1][i] / fug[0][i] - 1, 2);
			Ks[0][i] = phi[1][i] / phi[0][i];
		}
#ifdef DEBUG
		PrintX();
#endif // DEBUG

		iter++;
		if (iter > maxIt) {
			// OCP_WARNING("SSM not converged in Phase Spliting!");
			break;
		}
	}
	if (iter <= maxIt) {
		return true;
	}
	else {
		return false;
	}
}

bool MixtureComp::SplitSSM3(const bool& flag)
{
	return true;
}


void MixtureComp::RachfordRice2()  ///< Used when NP = 2
{
	const vector<OCP_DBL>& Ktmp = Ks[0];
	OCP_DBL Kmin = Ktmp[0];
	OCP_DBL Kmax = Ktmp[0];

	for (USI i = 1; i < NC; i++) {
		if (Ktmp[i] < Kmin)
			Kmin = Ktmp[i];
		if (Ktmp[i] > Kmax)
			Kmax = Ktmp[i];
	}

	OCP_DBL numin = 1 / (1 - Kmax);
	OCP_DBL numax = 1 / (1 - Kmin);

	nuj[0] = 0.5 * (numin + numax);

	// Solve RR with NR
	vector<OCP_DBL> tmpRR(NC, 0);
	OCP_DBL rj, J, dnuj, tmp;

	USI iter = 0;
	OCP_DBL RRtol = EoSctrl.RR.tol;
	while (true) {

		rj = 0;
		J = 0;
		for (USI i = 0; i < NC; i++) {
			tmpRR[i] = 1 + nuj[0] * (Ktmp[i] - 1);
			rj += zi[i] * (Ktmp[i] - 1) / tmpRR[i];
			J -= zi[i] * (Ktmp[i] - 1) * (Ktmp[i] - 1) / (tmpRR[i] * tmpRR[i]);
		}

		if (fabs(rj) < RRtol || iter > EoSctrl.RR.maxIt)
			break;

		dnuj = -rj / J;
		tmp = nuj[0] + dnuj;
		if (tmp < numax && tmp > numin) {
			nuj[0] = tmp;
		}
		else {
			if (dnuj > 0) {
				nuj[0] = (nuj[0] + numax) / 2;
			}
			else {
				nuj[0] = (nuj[0] + numin) / 2;
			}
		}

		iter++;
	}

	//if (iter > EoSctrl.RR.maxIt) {
	//	OCP_WARNING("RR2 not converged!");
	//}

	nuj[1] = 1 - nuj[0];

#ifdef DEBUG
	cout << nuj[0] << endl;
#endif // DEBUG

}


void MixtureComp::RachfordRice3() ///< Used when NP > 2
{

}

void MixtureComp::UpdateXRR()
{
	OCP_DBL tmp = 0;
	for (USI i = 0; i < NC; i++) {
		tmp = 1;
		for (USI j = 0; j < NP - 1; j++) {
			tmp += nuj[j] * (Ks[j][i] - 1);
		}
		x[NP - 1][i] = zi[i] / tmp;
		for (USI j = 0; j < NP - 1; j++) {
			x[j][i] = Ks[j][i] * x[NP - 1][i];
		}
	}
}

bool MixtureComp::SplitNR()
{

	for (USI j = 0; j < NP; j++) {
		nuj[j] = fabs(nuj[j]);
	}


	USI len = NC * (NP - 1);
	x2n();
	CalResSP();
	OCP_DBL eNR0;
	OCP_DBL eNR = Dnorm2(len, &resSP[0]);
	OCP_DBL NRtol = EoSctrl.NRsp.tol;

	OCP_DBL alpha;

#ifdef DEBUG
	cout << "NRSP Begins!\n";
#endif

	OCP_DBL en;
	USI iter = 0;
	eNR0 = eNR;
	bool flag = true;
	while (eNR > NRtol) {

		// eNR0 = eNR;

		ln = n;
		CalFugNAll();
		AssembleJmatSP();

		LUSolve(len, &JmatSP[0], &resSP[0], &pivot[0]);
		// PrintDX(NC, &resSP[0]);

		// SYSSolve(&uplo, len, &JmatSP[0], &resSP[0], &pivot[0]);
		// PrintDX(NC, &resSP[0]);

		alpha = CalStepNRsp();

		n[NP - 1] = zi;
		for (USI j = 0; j < NP - 1; j++) {
			Daxpy(NC, alpha, &resSP[j * NC], &n[j][0]);
			Daxpy(NC, -1, &n[j][0], &n[NP - 1][0]);

			nuj[j] = Dnorm1(NC, &n[j][0]);
			for (USI i = 0; i < NC; i++) {
				x[j][i] = n[j][i] / nuj[j];
			}
#ifdef DEBUG
			PrintDX(NC, &x[j][0]);
#endif
		}

		for (USI i = 0; i < NC; i++) {
			n[NP - 1][i] = fabs(n[NP - 1][i]);
		}
		nuj[NP - 1] = Dnorm1(NC, &n[NP - 1][0]);
		for (USI i = 0; i < NC; i++) {
			x[NP - 1][i] = n[NP - 1][i] / nuj[NP - 1];
		}


#ifdef DEBUG
		PrintDX(NC, &x[NP - 1][0]);
		cout << "---------------------" << endl;
#endif

		CalFugPhiAll();
		CalResSP();
		eNR = Dnorm2(len, &resSP[0]);
		iter++;
		if (eNR > eNR0 || iter > EoSctrl.NRsp.maxIt) {
			flag = false;
			break;
		}

		en = 0;
		for (USI j = 0; j < NP; j++) {
			Daxpy(NC, -1, &n[j][0], &ln[j][0]);
			en += Dnorm2(NC, &ln[j][0]);
		}
		if (en / NP < 1E-8) {
			break;
		}

	}
	return flag;
}


void MixtureComp::CalResSP()
{
	// So it equals -res
	for (USI j = 0; j < NP - 1; j++) {
		for (USI i = 0; i < NC; i++) {
			// if zi is too small, resSP[j][i] = 0?
			resSP[j * NC + i] = log(fug[NP - 1][i] / fug[j][i]);
		}
	}
}

void MixtureComp::CalFugNAll()
{
	OCP_DBL C, D, E, G;
	OCP_DBL Cnk, Dnk, Enk, Gnk;
	OCP_DBL tmp, aik;

	for (USI j = 0; j < NP; j++) {
		// j th phase
		vector<OCP_DBL>& fugn = fugN[j];
		const OCP_DBL& aj = Aj[j];
		const OCP_DBL& bj = Bj[j];
		const OCP_DBL& zj = Zj[j];
		const vector<OCP_DBL>& xj = x[j];

		for (USI i = 0; i < NC; i++) {
			tmp = 0;
			for (USI m = 0; m < NC; m++) {
				tmp += (1 - BIP[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
			}
			An[i] = 2 / nuj[j] * (tmp - aj);
			Bn[i] = 1 / nuj[j] * (Bi[i] - bj);
			Zn[i] = ((bj - zj) * An[i] + ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj))
				+ ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj
				- (delta1 + delta2 - 1) * zj * zj) * Bn[i]) / (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj
					+ (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
		}

		G = (zj + delta1 * bj) / (zj + delta2 * bj);

		for (USI i = 0; i < NC; i++) {
			// i th fugacity
			C = xj[i] * P / (zj - bj);
			// D = Bi[i] / bj * (zj - 1);
			tmp = 0;
			for (USI k = 0; k < NC; k++) {
				tmp += (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]) * xj[k];
			}
			E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

			for (USI k = 0; k <= i; k++) {
				// k th components

				//if (xj[i] < 1E-10) {
				//	fugn[k * NC + i] = phi[j][i] / nuj[j] * (delta(i, k) - xj[i]);
				//	/*cout << "fnn[" << j << "][" << i << "][" << k << "] = " << fugn[i * NC + k];
				//	cout << endl;*/
				//	continue;
				//}

				aik = (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]);

				Cnk = P / (zj - bj) / (zj - bj) * ((zj - bj) / nuj[j] * (delta(i, k) - xj[i]) - xj[i] * (Zn[k] - Bn[k]));
				Dnk = Bi[i] / bj * (Zn[k] - (Bi[k] - bj) * (zj - 1) / (nuj[j] * bj));
				Gnk = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) * (Bn[k] * zj - Zn[k] * bj);
				Enk = 1 / ((delta1 - delta2) * bj * bj) * (An[k] * bj - Bn[k] * aj) * (Bi[i] / bj - 2 * tmp / aj)
					+ aj / ((delta1 - delta2) * bj) * (-Bi[i] / (bj * bj) * Bn[k] - 2 / (aj * aj) * (aj * (aik - tmp) / nuj[j] - An[k] * tmp));

				fugn[i * NC + k] = 1 / C * Cnk + Dnk + Enk * log(G) + E / G * Gnk;
				fugn[k * NC + i] = fugn[i * NC + k];
				/*cout << "fnn[" << j << "][" << i << "][" << k << "] = " << fugn[i * NC + k];
				cout << endl;*/
			}
		}
	}
	// PrintFugN();
}


void MixtureComp::PrintFugN()
{
	for (USI j = 0; j < NP; j++) {
		for (USI i = 0; i < NC; i++) {
			for (USI k = 0; k < NC; k++) {
				cout << fugN[j][i * NC + k] << "   ";
			}
			cout << endl;
		}
		cout << "---------------" << endl;
	}
}


void MixtureComp::AssembleJmatSP()
{
	JmatSP.assign(JmatSP.size(), 0);

	OCP_DBL* Jtmp;
	const vector<OCP_DBL>& fugNp = fugN[NP - 1];

	for (USI j = 0; j < NP - 1; j++) {
		const vector<OCP_DBL>& fugNj = fugN[j];

		Jtmp = &JmatSP[j * NC * NC * (NP - 1)];
		for (USI i = 0; i < NC; i++) {
			// ith components
			for (USI k = 0; k < NC; k++) {
				// kth fugacity
				Jtmp[j * NC + k] = fugNj[i * NC + k] + fugNp[i * NC + k];
			}
			Jtmp += NC * (NP - 1);
		}
	}
}

OCP_DBL MixtureComp::CalStepNRsp()
{
	OCP_DBL alpha = 1;
	OCP_DBL tmp;

	for (USI j = 0; j < NP - 1; j++) {

		const OCP_DBL* nj = &n[j][0];
		const OCP_DBL* r = &resSP[j * NC];

		for (USI i = 0; i < NC; i++) {
			tmp = nj[i] + alpha * r[i];
			if (tmp < 0) {
				alpha *= 0.9 * fabs(nj[i] / r[i]);
			}
		}
	}
	return alpha;
}


void MixtureComp::AllocateOthers()
{
	muC.resize(NPmax);
	muAux.resize(NPmax);
	for (USI i = 0; i < NPmax; i++) {
		muAux[i].resize(5);
	}
	fugP.resize(NPmax);
	for (USI j = 0; j < NPmax; j++) {
		fugP[j].resize(NC);
	}
	Zp.resize(NPmax);
}


void MixtureComp::IdentifyPhase()
{
	if (NP == 1) {
		// Critical Temperature Method
		OCP_DBL A = 0;
		OCP_DBL B = 0;
		for (USI i = 0; i < NC; i++) {
			A += x[0][i] * comp[i].VcMW * comp[i].MW * comp[i].Tc;
			B += x[0][i] * comp[i].VcMW * comp[i].MW;
		}
		OCP_DBL Tc = A / B;
		if (T > Tc) {
			phaseLabel[0] = GAS;
		}
		else {
			phaseLabel[0] = OIL;
		}
	}
	else {
		// Compare MW
		CalMW();
		if (MW[0] > MW[1]) {
			phaseLabel[0] = OIL;
			phaseLabel[1] = GAS;
		}
		else {
			phaseLabel[0] = GAS;
			phaseLabel[1] = OIL;
		}
	}
}


void MixtureComp::CalViscosity()
{

}


void MixtureComp::CalViscoLBC()
{
	OCP_DBL tmp;
	OCP_DBL Tri;
	OCP_DBL xijT;
	OCP_DBL xijP;
	OCP_DBL xijV;

	for (USI j = 0; j < NP; j++) {
		vector<OCP_DBL>& xj = x[j];
		vector<OCP_DBL>& muA = muAux[j];
		muA.assign(muA.size(), 0);
		xijT = 0;
		xijP = 0;
		xijV = 0;

		for (USI i = 0; i < NC; i++) {
			tmp = 5.4402 * pow(comp[i].Tc, 1.0 / 6) / sqrt(MW[j]) / pow(comp[i].Pc, 2.0 / 3);
			Tri = T / comp[i].Tc;
			if (Tri <= 1.5) {
				tmp = 34 * 1E-5 * pow(Tri, 0.94) / tmp;
			}
			else {
				tmp = 17.78 * 1E-5 * pow((4.58 * Tri - 1.67), 0.625) / tmp;
			}
			muA[0] += xj[i] * sqrt(comp[i].MW) * tmp;
			muA[1] += xj[i] * sqrt(comp[i].MW);
			xijT += xj[i] * comp[i].Tc;
			xijP += xj[i] * comp[i].Pc;
			xijV += xj[i] * comp[i].VcMW * comp[i].MW;
		}
		muA[2] = 5.4402 * pow(xijT, 1.0 / 6) / sqrt(MW[j]) / pow(xijP, 2.0 / 3);
		muA[3] = xiC[j] * xijV;

		if (muA[3] <= 0.18) {
			muC[j] = muA[0] / muA[1] + 2.05 * 1E-4 * muA[3] / muA[2];
		}
		else {
			// muA[4] = 1.023 + 0.23364 * muA[3] + 0.58533 * pow(muA[3], 2) - 0.40758 * pow(muA[3], 3) + 0.093324 * pow(muA[3], 4);
			muA[4] = muA[3] * (muA[3] * (muA[3] * (0.093324 * muA[3] - 0.40758) + 0.58533) + 0.23364) + 1.023;
			muC[j] = muA[0] / muA[1] + 1E-4 * (pow(muA[4], 4) - 1) / muA[2];
		}
	}
}


void MixtureComp::CalViscoHZYT()
{

}

void MixtureComp::CalFugXAll()
{
	for (USI j = 0; j < NP; j++) {

		vector<OCP_DBL>& fugx = fugX[j];
		vector<OCP_DBL>& xj = x[j];
		OCP_DBL aj = Aj[j];
		OCP_DBL bj = Bj[j];
		OCP_DBL zj = Zj[j];
		OCP_DBL tmp = 0;

		Bx = Bi;
		for (USI i = 0; i < NC; i++) {
			tmp = 0;
			for (USI k = 0; k < NC; k++) {
				tmp += xj[k] * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]);
			}
			Ax[i] = 2 * tmp;
			Zx[i] = ((bj - zj) * Ax[i] + ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj))
				+ ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj
				- (delta1 + delta2 - 1) * zj * zj) * Bx[i]) / (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj
					+ (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));
		}

		OCP_DBL C, D, E, G;
		OCP_DBL Cxk, Dxk, Exk, Gxk;
		OCP_DBL* phiXT;
		G = (zj + delta1 * bj) / (zj + delta2 * bj);

		for (USI i = 0; i < NC; i++) {
			C = xj[i] * P / (zj - bj);
			// C = 1 / (zj - bj);
			// D = Bx[i] * (zj - 1) / bj;
			E = -aj / ((delta1 - delta2) * bj) * (Ax[i] / aj - Bx[i] / bj);

			for (USI k = 0; k < NC; k++) {
				// Cxk = -xj[i] * (Zx[k] - Bx[k]) / ((zj - bj) * (zj - bj));
				Cxk = ((zj - bj) * delta(i, k) - xj[i] * (Zx[k] - Bx[k])) * P / ((zj - bj) * (zj - bj));
				Dxk = Bx[i] / bj * (Zx[k] - Bx[k] * (zj - 1) / bj);
				Exk = (Ax[k] * bj - aj * Bx[k]) / (bj * bj) * (Ax[i] / aj - Bx[i] / bj) + aj / bj * (2 * (1 - BIP[i * NC + k]) * sqrt(Ai[i] * Ai[k]) / aj -
					Ax[k] * Ax[i] / (aj * aj) + Bx[i] * Bx[k] / (bj * bj));
				Exk /= -(delta1 - delta2);
				Gxk = (delta1 - delta2) / (zj + delta2 * bj) / (zj + delta2 * bj) * (zj * Bx[k] - Zx[k] * bj);
				fugx[i * NC + k] = 1 / C * Cxk + Dxk + Exk * log(G) + E / G * Gxk;
			}
		}
	}
}


void MixtureComp::CalFugPAll()
{

	OCP_DBL C, D, E, G;
	OCP_DBL Cp, Dp, Ep, Gp;
	OCP_DBL tmp;

	for (USI j = 0; j < NP; j++) {

		vector<OCP_DBL>& fugp = fugP[j];
		vector<OCP_DBL>& xj = x[j];
		OCP_DBL& aj = Aj[j];
		OCP_DBL& bj = Bj[j];
		OCP_DBL& zj = Zj[j];

		OCP_DBL Ap = aj / P;
		OCP_DBL Bp = bj / P;
		Zp[j] = ((bj - zj) * Ap + ((aj + delta1 * delta2 * (3 * bj * bj + 2 * bj))
			+ ((delta1 + delta2) * (2 * bj + 1) - 2 * delta1 * delta2 * bj) * zj
			- (delta1 + delta2 - 1) * zj * zj) * Bp) / (3 * zj * zj + 2 * ((delta1 + delta2 - 1) * bj - 1) * zj
				+ (aj + delta1 * delta2 * bj * bj - (delta1 + delta2) * bj * (bj + 1)));

		G = (zj + delta1 * bj) / (zj + delta2 * bj);
		Gp = (delta1 - delta2) / ((zj + delta2 * bj) * (zj + delta2 * bj)) * (Bp * zj - Zp[j] * bj);
		for (USI i = 0; i < NC; i++) {

			C = xj[i] * P / (zj - bj);
			// D = Bi[i] / bj * (zj - 1);

			tmp = 0;
			for (USI m = 0; m < NC; m++) {
				tmp += (1 - BIP[i * NC + m]) * sqrt(Ai[i] * Ai[m]) * xj[m];
			}

			E = -aj / ((delta1 - delta2) * bj) * (2 * tmp / aj - Bi[i] / bj);

			Cp = xj[i] / ((zj - bj) * (zj - bj)) * ((zj - bj) - P * (Zp[j] - Bp));
			Dp = Bi[i] / bj * Zp[j];
			// Ep = 0;

			fugp[i] = 1 / C * Cp + Dp + E / G * Gp;
		}
	}
}


USI MixtureComp::CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c, const bool& NTflag) const
{

	OCP_DBL Q = (a * a - 3 * b) / 9;
	OCP_DBL R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;

	OCP_DBL Q3 = Q * Q * Q;
	OCP_DBL M = R * R - Q3;

	if (M <= 0) {
		// 3 real roots
		OCP_DBL theta = acos(R / sqrt(Q3));
		Ztmp[0] = -2 * sqrt(Q) * cos(theta / 3) - a / 3;
		Ztmp[1] = -2 * sqrt(Q) * cos((theta + 2 * PI) / 3) - a / 3;
		Ztmp[2] = -2 * sqrt(Q) * cos((theta - 2 * PI) / 3) - a / 3;

		if (NTflag) {
			NTcubicroot(Ztmp[0], a, b, c);
			NTcubicroot(Ztmp[1], a, b, c);
			NTcubicroot(Ztmp[2], a, b, c);
		}


		sort(Ztmp.begin(), Ztmp.end());

		//vector<OCP_DBL> e(3, 0);
		//for (USI i = 0; i < 3; i++) {
		//	e[i] = Ztmp[i] * (Ztmp[i] * (Ztmp[i] + a) + b) + c;
		//}
		//for (USI i = 0; i < 3; i++) {
		//	cout << scientific << e[i] << "\t";
		//}

		return 3;
	}
	else {
		OCP_DBL tmp1 = -R + sqrt(M);
		OCP_DBL tmp2 = R + sqrt(M);
		OCP_DBL S = signD(tmp1) * pow(fabs(tmp1), 1.0 / 3);
		OCP_DBL T = -signD(tmp2) * pow(fabs(tmp2), 1.0 / 3);
		Ztmp[0] = S + T - a / 3;

		if (NTflag) {
			NTcubicroot(Ztmp[0], a, b, c);
		}


		//vector<OCP_DBL> e(1, 0);
		//for (USI i = 0; i < 1; i++) {
		//	e[i] = Ztmp[i] * (Ztmp[i] * (Ztmp[i] + a) + b) + c;
		//}
		//for (USI i = 0; i < 1; i++) {
		//	cout << scientific << e[i] << "\t";
		//}

		return 1;
	}
}



/// Return the sign of double di
OCP_DBL signD(const OCP_DBL& d) {
	if (d > 0) {
		return 1.0;
	}
	else if (d < 0) {
		return -1.0;
	}
	else {
		return 0.0;
	}
}


OCP_DBL delta(const USI& i, const USI& j)
{
	if (i == j) {
		return 1.0;
	}
	return 0.0;
}


void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c)
{
	OCP_DBL e = root * (root * (root + a) + b) + c;
	OCP_DBL df;
	OCP_DBL iter = 0;
	OCP_DBL optroot = root;
	OCP_DBL opte = fabs(e);

	while (fabs(e) > 1E-8) {

		df = root * (3 * root + 2 * a) + b;
		root = root - e / df;
		iter++;
		if (iter > 10)
		{
			//std::cout << "WARNING: INEXACT ROOT FOR CUBIC EQUATIONS" << std::endl;
			break;
		}
		e = root * (root * (root + a) + b) + c;
		if (fabs(e) <= opte) {
			opte = fabs(e);
			optroot = root;
		}
	}
	root = optroot;
}