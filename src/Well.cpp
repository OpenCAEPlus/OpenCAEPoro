#include <cmath>
#include "Well.hpp"

WellOpt::WellOpt(WellOptParam& Optparam)
{
	if (Optparam.Type == "INJ") {
		Type = INJ;
	}
	else if (Optparam.Type == "PROD") {
		Type = PROD;
	}	
	else {
		ERRORcheck("WRONG Well Type");
		exit(0);
	}
	
	if (Type == INJ) {
		if (Optparam.FluidType == "OIL") {
			FluidType = OIL;
		}
		else if (Optparam.FluidType == "GAS") {
			FluidType = GAS;
		}
		else if (Optparam.FluidType == "WATER") {
			FluidType = WATER;
		}
		else if (Optparam.FluidType == "SOLVENT") {
			FluidType = SOLVENT;
		}
		else {
			ERRORcheck("WRONG Fluid type");
			exit(0);
		}
	}
	
	
	if (Optparam.State == "OPEN") {
		State = OPEN;
	}
	else if (Optparam.State == "CLOSE") {
		State = CLOSE;
	}
	else {
		ERRORcheck("WRONG State type");
		exit(0);
	}
		
	if (Optparam.OptMode == "RATE") {
		OptMode = RATE_MODE;
	}
	else if (Optparam.OptMode == "ORAT") {
		OptMode = ORATE_MODE;
	}
	else if (Optparam.OptMode == "GRAT") {
		OptMode = GRATE_MODE;
	}
	else if (Optparam.OptMode == "WRAT") {
		OptMode = WRATE_MODE;
	}
	else if (Optparam.OptMode == "BHP") {
		OptMode = BHP_MODE;
	}
	else {
		ERRORcheck("WRONG Well Opt Mode");
		exit(0);
	}

	MaxRate = Optparam.MaxRate;
	MaxBHP = Optparam.MaxBHP;
	MinBHP = Optparam.MinBHP;

}

void Well::setup(Grid& myGrid, Bulk& myBulk)
{
	// zi
	if (myBulk.BLACKOIL) {
		for (auto& opt : OptSet) {

			opt.Zi.resize(myBulk.Nc, 0);
			if (opt.Type == INJ) {
				// INJ
				switch (myBulk.PVTmode)
				{
				case PHASE_W:
				case PHASE_OW:
					opt.Zi.back() = 1;
					break;
				case PHASE_OGW:
					if (opt.FluidType == GAS)		opt.Zi[1] = 1;
					else							opt.Zi[2] = 1;
					break;
				default:
					ERRORcheck("WRONG Blackoil type!");
					exit(0);
				}
			}
			else {
				// PROD
				switch (myBulk.PVTmode)
				{
				case PHASE_W:
					opt.Zi.back() = 1;
					break;
				case PHASE_OW:
					if (opt.OptMode == ORATE_MODE)			opt.Zi[0] = 1;
					else									opt.Zi[1] = 1;
					break;
				case PHASE_OGW:
					if (opt.OptMode == ORATE_MODE)			opt.Zi[0] = 1;
					else if (opt.OptMode == GRATE_MODE)		opt.Zi[1] = 1;
					else									opt.Zi[2] = 1;
					break;
				default:
					ERRORcheck("WRONG Blackoil type!");
					exit(0);
				}
			}
		}
	}
	else if (myBulk.COMPS) {
		
	} 
	else {
		ERRORcheck("Wrong Mixture Type !");
		exit(0);
	}

	Qi_lbmol.resize(myBulk.Nc);
	// perf
	PerfNum = K2 - K1 + 1;
	dG.resize(PerfNum, 0);
	Perf.resize(PerfNum);
	for (int p = 0; p < PerfNum; p++) {
		Perf[p].State = OPEN;
		int Idg = (K1 + p) * myGrid.Nx * myGrid.Ny + J * myGrid.Nx + I;
		Perf[p].Location = myGrid.ActiveMap_G2B[Idg];
		Perf[p].Depth = myBulk.Depth[Perf[p].Location];
		Perf[p].Multiplier = 1;
		Perf[p].qi_lbmol.resize(myBulk.Nc);
	}
	calWI_Peaceman_Vertical(myBulk);
	cout << "Well::setup" << endl;
}

void Well::init(const Bulk& myBulk) {
	BHP = myBulk.P[Perf[0].Location];
}

void Well::calWI_Peaceman_Vertical(const Bulk& myBulk)
{
	// this fomular needs to be carefully checked !
	// especially the dz
	if (WI > 0) {
		for (int p = 0; p < PerfNum; p++) {
			Perf[p].WI = WI;
		}
	}
	else {
		for (int p = 0; p < PerfNum; p++) {
			int Idb = Perf[p].Location;
			double kxky = myBulk.Rock_Kx[Idb] * myBulk.Rock_Ky[Idb];
			double kx_ky = myBulk.Rock_Kx[Idb] / myBulk.Rock_Ky[Idb];
			assert(kx_ky > 0);


			double dx = myBulk.Dx[Idb];
			double dy = myBulk.Dy[Idb];
			double dz = myBulk.Dz[Idb] * myBulk.Ntg[Idb];

			double ro = 0.28 * pow((dx * dx * pow(1 / kx_ky, 0.5) + dy * dy * pow(kx_ky, 0.5)), 0.5);
			ro /= (pow(kx_ky, 0.25) + pow(1 / kx_ky, 0.25));
			if (Kh < 0) {
				Perf[p].WI = (2 * PI) * (dz * pow(kxky, 0.5)) / (log(ro / Radius) + SkinFactor);
			}
			else {
				Perf[p].WI = (2 * PI) * Kh / (log(ro / Radius) + SkinFactor);
			}
			
		}
	}
	
}


void Well::assembleMat_INJ(const Bulk& myBulk, Solver<double>& mySolver, double dt)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	int wId = mySolver.Dim;
	// important !
	mySolver.Dim++;

	for (int p = 0; p < PerfNum; p++) {
		int k = Perf[p].Location; 

		double trans = 0;
		for (int j = 0; j < np; j++) {
			if (myBulk.PhaseExist[k * np + j]) {
				trans += myBulk.Kr[k * np + j] / myBulk.Mu[k * np + j];
			}		
		}
		trans *= dt;
		double Vfi_zi = 0;
		for (int i = 0; i < nc; i++) {
			Vfi_zi += myBulk.Vfi[k * nc + i] * Opt.Zi[i];
		}
		
		int pvtnum = myBulk.PVTNUM[k];
		double xi = myBulk.Flashcal[pvtnum]->xiPhase(myBulk.P[k], myBulk.T, &Opt.Zi[0]);
		double valw = 1 * xi * trans * CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier;
		double bw = valw * dG[p];
		double valb = valw * Vfi_zi;
		double bb = valb * dG[p];

		// Bulk to Well
		
		// diag
		int ptr = mySolver.DiagPtr[k];
		mySolver.Val[k][ptr] += valb;
		// off diag
		mySolver.ColId[k].push_back(wId);
		mySolver.Val[k].push_back(-valb);
		// b
		mySolver.b[k] += bb;


		// Well to Bulk
		switch (Opt.OptMode)
		{
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// diag
			mySolver.DiagVal[wId] += valw;
			// off diag
			mySolver.ColId[wId].push_back(k);
			mySolver.Val[wId].push_back(-valw);
			// b
			mySolver.b[wId] -= bw;
			break;
		case BHP_MODE:
			mySolver.ColId[wId].push_back(k);
			mySolver.Val[wId].push_back(0);
			break;
		default:
			ERRORcheck("Wrong Well Opt mode");
			exit(0);
		}
	}

	// Well Self
	assert(mySolver.Val[wId].size() == PerfNum);
	// the order of perforation is not necessarily in order
	switch (Opt.OptMode)
	{
	case RATE_MODE:
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		mySolver.ColId[wId].push_back(wId);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(mySolver.DiagVal[wId]);
		// b
		mySolver.b[wId] += dt * Opt.MaxRate;
		break;
	case BHP_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(dt);
		// b
		mySolver.b[wId] += dt * Opt.MaxBHP;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode in function");
		exit(0);
	}
}


void Well::assembleMat_PROD_BLK(const Bulk& myBulk, Solver<double>& mySolver, double dt)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	int wId = mySolver.Dim;
	// important !
	mySolver.Dim++;
	cout << Name << endl;

	for (int p = 0; p < PerfNum; p++) {
		int k = Perf[p].Location;

		double valb = 0;	double bb = 0;
		double valw = 0;	double bw = 0;
		
		for (int j = 0; j < np; j++) {
			if (!myBulk.PhaseExist[k * np + j])
				continue;

			double tempb = 0;	
			double tempw = 0;
			
			for (int i = 0; i < nc; i++) {
				tempb += myBulk.Vfi[k * nc + i] * myBulk.Cij[k * np * nc + j * nc + i];
				tempw += Opt.Zi[i] * myBulk.Cij[k * np * nc + j * nc + i];
			}
			double trans = myBulk.Xi[k * np + j] * myBulk.Kr[k * np + j] / myBulk.Mu[k * np + j];
			valb += tempb * trans;
			valw += tempw * trans;

			double dP = dG[p] - myBulk.Pc[k * np + j];
			bb += tempb * trans * dP;
			bw += tempw * trans * dP;
		}
		double trans = dt * CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier;
		valb *= trans;
		valw *= trans;
		bb *= trans;
		bw *= trans;

		// Bulk to Well
		// diag
		int ptr = mySolver.DiagPtr[k];
		mySolver.Val[k][ptr] += valb;
		// off diag
		mySolver.ColId[k].push_back(wId);
		mySolver.Val[k].push_back(-valb);
		// b
		mySolver.b[k] += bb;


		// Well to Bulk
		switch (Opt.OptMode)
		{
		case RATE_MODE:
		case ORATE_MODE:
		case GRATE_MODE:
		case WRATE_MODE:
		case LRATE_MODE:
			// diag  !!! attention! sign is -
			mySolver.DiagVal[wId] -= valw;
			// off diag
			mySolver.ColId[wId].push_back(k);
			mySolver.Val[wId].push_back(valw);
			// b
			mySolver.b[wId] += bw;
			break;
		case BHP_MODE:
			// off diag
			mySolver.ColId[wId].push_back(k);
			mySolver.Val[wId].push_back(0);
			break;
		default:
			ERRORcheck("Wrong Well Opt mode");
			exit(0);
		}
	}

	// Well Self
	assert(mySolver.Val[wId].size() == PerfNum);
	// the order of perforation is not necessarily in order
	switch (Opt.OptMode)
	{
	case RATE_MODE:
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		mySolver.ColId[wId].push_back(wId);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(mySolver.DiagVal[wId]);
		// b
		mySolver.b[wId] += dt * Opt.MaxRate;
		break;
	case BHP_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(dt);
		// b
		mySolver.b[wId] += dt * Opt.MinBHP;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode");
		exit(0);
	}
}

void Well::calInjdG(const Bulk& myBulk)
{
	double maxlen = 10;
	int seg_num = 0;
	double seg_len = 0;
	double	qtacc = 0;
	double	rhoacc = 0;
	double	rhotmp = 0;
	vector<double>		dGperf(PerfNum, 0);
	
	if (Depth <= Perf.front().Depth) {
		// Well is higher
		for (int p = PerfNum - 1; p >= 0; p--) {
			if (p == 0) {
				seg_num = ceil((Perf[0].Depth - Depth) / maxlen);
				seg_len = (Perf[0].Depth - Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p].Depth - Perf[p - 1].Depth) / maxlen);
				seg_len = (Perf[p].Depth - Perf[p - 1].Depth) / seg_num;
			}
			int n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			double Pperf = Perf[p].P;
			double Ptmp = Pperf;

			int pvtnum = myBulk.PVTNUM[n];
			for (int i = 0; i < seg_num; i++) {
				Ptmp -= myBulk.Flashcal[pvtnum]->rhoPhase(Ptmp, myBulk.T, Opt.Zi.data()) * GRAVITY_FACTOR * seg_len;
			}
			dGperf[p] = Pperf - Ptmp;
		}
		dG[0] = dGperf[0];
		for (int p = 1; p < PerfNum; p++) {
			dG[p] = dG[p - 1] + dGperf[p];
		}
	}
	else if (Depth >= Perf[PerfNum - 1].Depth) {
		// Well is lower
		for (int p = 0; p < PerfNum; p++) {
			if (p == PerfNum - 1) {
				seg_num = ceil((Depth - Perf[PerfNum - 1].Depth) / maxlen);
				seg_len = (Depth - Perf[PerfNum - 1].Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p + 1].Depth - Perf[p].Depth) / maxlen);
				seg_len = (Perf[p + 1].Depth - Perf[p].Depth) / seg_num;
			}
			int n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			double Pperf = Perf[p].P;
			double Ptmp = Pperf;

			int pvtnum = myBulk.PVTNUM[n];
			for (int i = 0; i < seg_num; i++) {
				Ptmp += myBulk.Flashcal[pvtnum]->rhoPhase(Ptmp, myBulk.T, Opt.Zi.data()) * GRAVITY_FACTOR * seg_len;
			}
			dGperf[p] = Ptmp - Pperf;
		}
		dG[PerfNum - 1] = dGperf[PerfNum - 1];
		for (int p = PerfNum - 2; p >= 0; p--) {
			dG[p] = dG[p + 1] + dGperf[p];
		}
	}
}

void Well::calProddG(const Bulk& myBulk)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	double maxlen = 10;
	int seg_num = 0;
	double seg_len = 0;
	vector<double>		tmpNi(nc, 0);
	vector<double>		dGperf(PerfNum, 0);
	double	qtacc = 0;
	double	rhoacc = 0;
	double	rhotmp = 0;
	if (Depth <= Perf.front().Depth) {
		// Well is higher
		for (int p = PerfNum - 1; p >= 0; p--) {
			if (p == 0) {
				seg_num = ceil((Perf[0].Depth - Depth) / maxlen);
				seg_len = (Perf[0].Depth - Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p].Depth - Perf[p - 1].Depth) / maxlen);
				seg_len = (Perf[p].Depth - Perf[p - 1].Depth) / seg_num;
			}

			int n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			double Pperf = Perf[p].P;
			double Ptmp = Pperf;
			Perf[p].qi_lbmol.assign(nc, 0);
			// vector<double>	qi_lbmol(nc, 0);

			for (int j = 0; j < np; j++) {
				int id = n * np + j;
				if (myBulk.PhaseExist[id]) {
					double dP = myBulk.Pj[id] - Perf[p].P;
					double kr = myBulk.Kr[id];
					double mu = myBulk.Mu[id];
					double xi = myBulk.Xi[id];
					double qj = CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier * kr / mu * dP * xi;
					double xij;
					for (int i = 0; i < nc; i++) {
						xij = myBulk.Cij[id * nc + i];
						Perf[p].qi_lbmol[i] += qj *  xij;
					}
				}
			}
			int pvtnum = myBulk.PVTNUM[n];
			tmpNi.assign(nc, 0);
			for (int p1 = PerfNum - 1; p1 >= p; p1--) {
				for (int i = 0; i < nc; i++) {
					tmpNi[i] += Perf[p1].qi_lbmol[i];
				}
			}
			for (int k = 0; k < seg_num; k++) {
				myBulk.Flashcal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
				for (int j = 0; j < myBulk.Np; j++) {
					if (myBulk.Flashcal[pvtnum]->PhaseExist[j]) {
						rhotmp = myBulk.Flashcal[pvtnum]->Rho[j];
						qtacc += myBulk.Flashcal[pvtnum]->V[j];
						rhoacc += myBulk.Flashcal[pvtnum]->V[j] * rhotmp * GRAVITY_FACTOR;
					}
				}
				Ptmp -= rhoacc / qtacc * seg_len;
			}
			dGperf[p] = Pperf - Ptmp;
		}
		dG[0] = dGperf[0];
		for (int p = 1; p < PerfNum; p++) {
			dG[p] = dG[p - 1] + dGperf[p];
		}
	}
	else if (Depth >= Perf.back().Depth) {
		// Well is lower
		for (int p = 0; p < PerfNum; p++) {
			if (p == PerfNum - 1) {
				seg_num = ceil((Depth - Perf[PerfNum - 1].Depth) / maxlen);
				seg_len = (Depth - Perf[PerfNum - 1].Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p + 1].Depth - Perf[p].Depth) / maxlen);
				seg_len = (Perf[p + 1].Depth - Perf[p].Depth) / seg_num;
			}

			int n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			double Pperf = Perf[p].P;
			double Ptmp = Pperf;
			Perf[p].qi_lbmol.assign(nc, 0);
			// vector<double>	qi_lbmol(nc, 0);

			for (int j = 0; j < np; j++) {
				int id = n * np + j;
				if (myBulk.PhaseExist[id]) {
					double dP = myBulk.Pj[id] - Perf[p].P;
					double kr = myBulk.Kr[id];
					double mu = myBulk.Mu[id];
					double xi = myBulk.Xi[id];
					double qj = CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier * kr / mu * dP * xi;
					double xij;
					for (int i = 0; i < nc; i++) {
						xij = myBulk.Cij[id * nc + i];
						Perf[p].qi_lbmol[i] += qj * xij;
					}
				}
			}
			int pvtnum = myBulk.PVTNUM[n];
			tmpNi.assign(nc, 0);
			for (int p1 = PerfNum - 1; p1 >= p; p1--) {
				for (int i = 0; i < nc; i++) {
					tmpNi[i] += Perf[p1].qi_lbmol[i];
				}
			}
			for (int k = 0; k < seg_num; k++) {
				myBulk.Flashcal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
				for (int j = 0; j < np; j++) {
					if (myBulk.Flashcal[pvtnum]->PhaseExist[j]) {
						rhotmp = myBulk.Flashcal[pvtnum]->Rho[j];
						qtacc += myBulk.Flashcal[pvtnum]->V[j];
						rhoacc += myBulk.Flashcal[pvtnum]->V[j] * rhotmp * GRAVITY_FACTOR;
					}
				}
				Ptmp += rhoacc / qtacc * seg_len;
			}
			dGperf[p] = Ptmp - Pperf;
		}
		dG[PerfNum - 1] = dGperf[PerfNum - 1];
		for (int p = PerfNum - 2; p >= 0; p--) {
			dG[p] = dG[p + 1] + dGperf[p];
		}

	}
	else {
		ERRORcheck("Wrong Well position");
		exit(0);
	}
}

double Well::calInjRate_blk(const Bulk& myBulk, int flag)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	double qj = 0;

	double Pwell;
	if (flag == -1) {
		Pwell = Opt.MaxBHP;
	}
	else {
		Pwell = BHP;
	}

	for (int p = 0; p < PerfNum; p++) {

		double Pperf = Pwell + dG[p];
		int k = Perf[p].Location;

		double trans = 0;
		for (int j = 0; j < np; j++) {
			int id = k * np + j;
			if (myBulk.PhaseExist[id]) {
				trans += myBulk.Kr[id] / myBulk.Mu[id];
			}
		}
		int pvtnum = myBulk.PVTNUM[k];
		double xi = myBulk.Flashcal[pvtnum]->xiPhase(myBulk.P[k], myBulk.T, &Opt.Zi[0]);
		double dP = Pperf - myBulk.P[k];
		qj += CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier * trans * xi * dP;
	}
	return qj;
}

double Well::calProdRate_blk(const Bulk& myBulk)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	double qj = 0;

	for (int p = 0; p < PerfNum; p++) {

		double Pperf = Opt.MinBHP + dG[p];
		int k = Perf[p].Location;

		for (int j = 0; j < np; j++) {
			int id = k * np + j;
			if (myBulk.PhaseExist[id]) {	
				double temp = 0;
				for (int i = 0; i < nc; i++) {
					temp += Opt.Zi[i] * myBulk.Cij[id * nc + i];
				}
				double kr = myBulk.Kr[id];
				double mu = myBulk.Mu[id];
				double xi = myBulk.Xi[id];
				double dP = myBulk.Pj[id] - Pperf;
				qj += CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier * kr / mu * xi * dP * temp;
			}
		}
	}
	return qj;
}

void Well::calProdqi_blk(const Bulk& myBulk)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;

	Qi_lbmol.assign(nc, 0);

	for (int p = 0; p < PerfNum; p++) {

		Perf[p].qi_lbmol.assign(nc, 0);
		double Pperf = BHP + dG[p];
		int k = Perf[p].Location;

		for (int j = 0; j < np; j++) {
			int id = k * np + j;
			if (myBulk.PhaseExist[id]) {
				double kr = myBulk.Kr[id];
				double mu = myBulk.Mu[id];
				double xi = myBulk.Xi[id];
				double dP = myBulk.Pj[id] - Pperf;
				double xij;
				for (int i = 0; i < nc; i++) {
					xij = myBulk.Cij[id * nc + i];
					Perf[p].qi_lbmol[i] += CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier * kr / mu * xi * xij * dP;
				}
			}
		}
		for (int i = 0; i < nc; i++) {
			Qi_lbmol[i] += Perf[p].qi_lbmol[i];
		}
	}

	for (int i = 0; i < nc; i++) {
		if (myBulk.PhaseLabel[i] == OIL) {
			WOPR = Qi_lbmol[i];
			break;
		}
		else if (myBulk.PhaseLabel[i] == GAS) {
			WGPR = Qi_lbmol[i];
		}
		else if (myBulk.PhaseLabel[i] == WATER) {
			WWPR = Qi_lbmol[i];
		}
	}
}

void Well::checkOptMode(const Bulk& myBulk)
{
	if (Opt.Type == INJ) {
		if (calInjRate_blk(myBulk, -1) > Opt.MaxRate) {
			Opt.OptMode = RATE_MODE;
		}
		else {
			Opt.OptMode = BHP_MODE;
		}
	}
	else {
		if (calProdRate_blk(myBulk) > Opt.MaxRate) {
			Opt.OptMode = RATE_MODE;
		}
		else {
			Opt.OptMode = BHP_MODE;
		}
	}
}
