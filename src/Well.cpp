#include <cmath>
#include "Well.hpp"

WellOpt::WellOpt(const WellOptParam& Optparam)
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

void Well::setup(const Grid& myGrid, const Bulk& myBulk)
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
	ldG = dG;
	Perf.resize(PerfNum);
	for (USI p = 0; p < PerfNum; p++) {
		Perf[p].State = OPEN;
		OCP_USI Idg = (K1 + p) * myGrid.Nx * myGrid.Ny + J * myGrid.Nx + I;
		if (!myGrid.ActiveMap_G2B[Idg].getAct()) {
			ERRORcheck("Perforation is in inactive bulk !");
			exit(0);
		}
		Perf[p].Location = myGrid.ActiveMap_G2B[Idg].getId();
		Perf[p].Depth = myBulk.Depth[Perf[p].Location];
		Perf[p].Multiplier = 1;
		Perf[p].qi_lbmol.resize(myBulk.Nc);
		Perf[p].transj.resize(myBulk.Np);
	}
	if (Depth < 0)
		Depth = Perf[0].Depth;

	calWI_Peaceman_Vertical(myBulk);
	cout << "Well::setup" << endl;
}

void Well::init(const Bulk& myBulk) {
	BHP = myBulk.P[Perf[0].Location];
}

OCP_DBL Well::calCFL(const Bulk& myBulk, const OCP_DBL& dt) const
{
	OCP_DBL cfl = 0;
	OCP_DBL tmp = 0;
	for (USI p = 0; p < PerfNum; p++) {
		OCP_USI k = Perf[p].Location;
		tmp = fabs(Perf[p].qt_ft3) * dt;
		tmp /= myBulk.Rock_Vp[k];
		if (cfl < tmp)
			cfl = tmp;
	}
	return cfl;
}

void Well::calWI_Peaceman_Vertical(const Bulk& myBulk)
{
	// this fomular needs to be carefully checked !
	// especially the dz
	if (WI > 0) {
		for (USI p = 0; p < PerfNum; p++) {
			Perf[p].WI = WI;
		}
	}
	else {
		for (USI p = 0; p < PerfNum; p++) {
			OCP_USI Idb = Perf[p].Location;
			OCP_DBL kxky = myBulk.Rock_Kx[Idb] * myBulk.Rock_Ky[Idb];
			OCP_DBL kx_ky = myBulk.Rock_Kx[Idb] / myBulk.Rock_Ky[Idb];
			assert(kx_ky > 0);


			OCP_DBL dx = myBulk.Dx[Idb];
			OCP_DBL dy = myBulk.Dy[Idb];
			OCP_DBL dz = myBulk.Dz[Idb] * myBulk.Ntg[Idb];

			OCP_DBL ro = 0.28 * pow((dx * dx * pow(1 / kx_ky, 0.5) + dy * dy * pow(kx_ky, 0.5)), 0.5);
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


void Well::assembleMat_INJ_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, const OCP_DBL& dt) const
{
	USI nc = myBulk.Nc;
	OCP_USI wId = mySolver.Dim;
	// important !
	mySolver.Dim++;

	for (USI p = 0; p < PerfNum; p++) {
		OCP_USI k = Perf[p].Location; 

		OCP_DBL Vfi_zi = 0;
		for (USI i = 0; i < nc; i++) {
			Vfi_zi += myBulk.Vfi[k * nc + i] * Opt.Zi[i];
		}
		
		USI pvtnum = myBulk.PVTNUM[k];
		Perf[p].Xi = myBulk.Flashcal[pvtnum]->xiPhase(myBulk.P[k], myBulk.T, &Opt.Zi[0]);
		OCP_DBL valw = dt * Perf[p].Xi * Perf[p].transj[0];
		OCP_DBL bw = valw * dG[p];
		OCP_DBL valb = valw * Vfi_zi;
		OCP_DBL bb = valb * dG[p];

		// Bulk to Well
		
		// diag
		USI ptr = mySolver.DiagPtr[k];
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
		mySolver.ColId[wId].push_back(wId);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(dt);
		// b
		mySolver.b[wId] += dt * Opt.MaxBHP;
		// u   initial value
		mySolver.u[wId] = Opt.MaxBHP;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode in function");
		exit(0);
	}
}


void Well::assembleMat_PROD_BLK_IMPES(const Bulk& myBulk, Solver<OCP_DBL>& mySolver, const OCP_DBL& dt) const
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;
	OCP_USI wId = mySolver.Dim;
	// important !
	mySolver.Dim++;

	for (USI p = 0; p < PerfNum; p++) {
		OCP_USI k = Perf[p].Location;

		OCP_DBL valb = 0;	OCP_DBL bb = 0;
		OCP_DBL valw = 0;	OCP_DBL bw = 0;
		
		for (USI j = 0; j < np; j++) {
			if (!myBulk.PhaseExist[k * np + j])
				continue;

			OCP_DBL tempb = 0;	
			OCP_DBL tempw = 0;
			
			for (USI i = 0; i < nc; i++) {
				tempb += myBulk.Vfi[k * nc + i] * myBulk.Cij[k * np * nc + j * nc + i];
				tempw += Opt.Zi[i] * myBulk.Cij[k * np * nc + j * nc + i];
			}
			OCP_DBL trans = dt * Perf[p].transj[j] * myBulk.Xi[k * np + j];
			valb += tempb * trans;
			valw += tempw * trans;

			OCP_DBL dP = dG[p] - myBulk.Pc[k * np + j];
			bb += tempb * trans * dP;
			bw += tempw * trans * dP;
		}

		// Bulk to Well
		// diag
		USI ptr = mySolver.DiagPtr[k];
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
		mySolver.ColId[wId].push_back(wId);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(dt);
		// b
		mySolver.b[wId] += dt * Opt.MinBHP;
		// u   initial value
		mySolver.u[wId] = Opt.MinBHP;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode");
		exit(0);
	}
}

void Well::smoothdG()
{
	for (USI p = 0; p < PerfNum; p++) {
		dG[p] = (ldG[p] + dG[p]) / 2;      // seems better
		// dG[p] = ldG[p] + 0.618 * (dG[p] - ldG[p]);
	}
}

void Well::caldG(const Bulk& myBulk)
{
	if (Opt.Type == INJ)
		calInjdG(myBulk);
	else
		calProddG(myBulk);
}

void Well::calInjdG(const Bulk& myBulk)
{
	OCP_DBL maxlen = 10;
	USI seg_num = 0;
	OCP_DBL seg_len = 0;
	vector<OCP_DBL>		dGperf(PerfNum, 0);
	
	if (Depth <= Perf.front().Depth) {
		// Well is higher
		for (OCP_INT p = PerfNum - 1; p >= 0; p--) {
			if (p == 0) {
				seg_num = ceil((Perf[0].Depth - Depth) / maxlen);
				seg_len = (Perf[0].Depth - Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p].Depth - Perf[p - 1].Depth) / maxlen);
				seg_len = (Perf[p].Depth - Perf[p - 1].Depth) / seg_num;
			}
			OCP_USI n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			OCP_DBL Pperf = Perf[p].P;
			OCP_DBL Ptmp = Pperf;

			USI pvtnum = myBulk.PVTNUM[n];
			for (USI i = 0; i < seg_num; i++) {
				Ptmp -= myBulk.Flashcal[pvtnum]->rhoPhase(Ptmp, myBulk.T, Opt.Zi.data()) * GRAVITY_FACTOR * seg_len;
			}
			dGperf[p] = Pperf - Ptmp;
		}
		dG[0] = dGperf[0];
		for (USI p = 1; p < PerfNum; p++) {
			dG[p] = dG[p - 1] + dGperf[p];
		}
	}
	else if (Depth >= Perf[PerfNum - 1].Depth) {
		// Well is lower
		for (USI p = 0; p < PerfNum; p++) {
			if (p == PerfNum - 1) {
				seg_num = ceil((Depth - Perf[PerfNum - 1].Depth) / maxlen);
				seg_len = (Depth - Perf[PerfNum - 1].Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p + 1].Depth - Perf[p].Depth) / maxlen);
				seg_len = (Perf[p + 1].Depth - Perf[p].Depth) / seg_num;
			}
			OCP_USI n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			OCP_DBL Pperf = Perf[p].P;
			OCP_DBL Ptmp = Pperf;

			USI pvtnum = myBulk.PVTNUM[n];
			for (USI i = 0; i < seg_num; i++) {
				Ptmp += myBulk.Flashcal[pvtnum]->rhoPhase(Ptmp, myBulk.T, Opt.Zi.data()) * GRAVITY_FACTOR * seg_len;
			}
			dGperf[p] = Ptmp - Pperf;
		}
		dG[PerfNum - 1] = dGperf[PerfNum - 1];
		for (OCP_INT p = PerfNum - 2; p >= 0; p--) {
			dG[p] = dG[p + 1] + dGperf[p];
		}
	}
}

void Well::calProddG(const Bulk& myBulk)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;
	OCP_DBL maxlen = 10;
	USI seg_num = 0;
	OCP_DBL seg_len = 0;
	vector<OCP_DBL>		tmpNi(nc, 0);
	vector<OCP_DBL>		dGperf(PerfNum, 0);
	OCP_DBL	qtacc = 0;
	OCP_DBL	rhoacc = 0;
	OCP_DBL	rhotmp = 0;

	

	if (Depth <= Perf.front().Depth) {
		// Well is higher

		// check qi_lbmol   ----   test
		if (Perf[PerfNum - 1].State == CLOSE) {
			for (OCP_INT p = PerfNum - 2; p >= 0; p--) {
				if (Perf[p].State == OPEN) {
					for (USI i = 0; i < nc; i++) {
						Perf[PerfNum - 1].qi_lbmol[i] = Perf[p].qi_lbmol[i];
					}
					break;
				}
			}
		}

		for (OCP_INT p = PerfNum - 1; p >= 0; p--) {

			if (p == 0) {
				seg_num = ceil((Perf[0].Depth - Depth) / maxlen);
				seg_len = (Perf[0].Depth - Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p].Depth - Perf[p - 1].Depth) / maxlen);
				seg_len = (Perf[p].Depth - Perf[p - 1].Depth) / seg_num;
			}

			OCP_USI n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			OCP_DBL Pperf = Perf[p].P;
			OCP_DBL Ptmp = Pperf;

			USI pvtnum = myBulk.PVTNUM[n];
			tmpNi.assign(nc, 0);
			for (OCP_INT p1 = PerfNum - 1; p1 >= p; p1--) {
				for (USI i = 0; i < nc; i++) {
					tmpNi[i] += Perf[p1].qi_lbmol[i];
				}
			}
			for (USI k = 0; k < seg_num; k++) {
				myBulk.Flashcal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
				for (USI j = 0; j < myBulk.Np; j++) {
					if (myBulk.Flashcal[pvtnum]->PhaseExist[j]) {
						rhotmp = myBulk.Flashcal[pvtnum]->Rho[j];
						qtacc += myBulk.Flashcal[pvtnum]->V[j] / seg_num;
						rhoacc += myBulk.Flashcal[pvtnum]->V[j] * rhotmp * GRAVITY_FACTOR / seg_num;
					}
				}
				Ptmp -= rhoacc / qtacc * seg_len;
			}
			dGperf[p] = Pperf - Ptmp;
		}
		dG[0] = dGperf[0];
		for (USI p = 1; p < PerfNum; p++) {
			dG[p] = dG[p - 1] + dGperf[p];
		}
	}
	else if (Depth >= Perf.back().Depth) {
		// Well is lower

		// check qi_lbmol   ----   test
		if (Perf[0].State == CLOSE) {
			for (USI p = 1; p <= PerfNum; p++) {
				if (Perf[p].State == OPEN) {
					for (USI i = 0; i < nc; i++) {
						Perf[PerfNum - 1].qi_lbmol[i] = Perf[p].qi_lbmol[i];
					}
					break;
				}
			}
		}

		for (USI p = 0; p < PerfNum; p++) {
			if (p == PerfNum - 1) {
				seg_num = ceil((Depth - Perf[PerfNum - 1].Depth) / maxlen);
				seg_len = (Depth - Perf[PerfNum - 1].Depth) / seg_num;
			}
			else {
				seg_num = ceil((Perf[p + 1].Depth - Perf[p].Depth) / maxlen);
				seg_len = (Perf[p + 1].Depth - Perf[p].Depth) / seg_num;
			}

			OCP_USI n = Perf[p].Location;
			Perf[p].P = BHP + dG[p];
			OCP_DBL Pperf = Perf[p].P;
			OCP_DBL Ptmp = Pperf;

			USI pvtnum = myBulk.PVTNUM[n];
			tmpNi.assign(nc, 0);
			for (OCP_INT p1 = PerfNum - 1; p1 >= p; p1--) {
				for (USI i = 0; i < nc; i++) {
					tmpNi[i] += Perf[p1].qi_lbmol[i];
				}
			}
			for (USI k = 0; k < seg_num; k++) {
				myBulk.Flashcal[pvtnum]->Flash_Ni(Ptmp, myBulk.T, tmpNi.data());
				for (USI j = 0; j < np; j++) {
					if (myBulk.Flashcal[pvtnum]->PhaseExist[j]) {
						rhotmp = myBulk.Flashcal[pvtnum]->Rho[j];
						qtacc += myBulk.Flashcal[pvtnum]->V[j] / seg_num;
						rhoacc += myBulk.Flashcal[pvtnum]->V[j] * rhotmp * GRAVITY_FACTOR / seg_num;
					}
				}
				Ptmp += rhoacc / qtacc * seg_len;
			}
			dGperf[p] = Ptmp - Pperf;
		}
		dG[PerfNum - 1] = dGperf[PerfNum - 1];
		for (OCP_INT p = PerfNum - 2; p >= 0; p--) {
			dG[p] = dG[p + 1] + dGperf[p];
		}

	}
	else {
		ERRORcheck("Wrong Well position");
		exit(0);
	}
}

void Well::calTrans(const Bulk& myBulk)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;

	if (Opt.Type == INJ) {
		for (USI p = 0; p < PerfNum; p++) {
			Perf[p].transj.assign(np, 0);
			OCP_USI k = Perf[p].Location;
			OCP_DBL temp = CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier;
			
			// single phase
			for (USI j = 0; j < np; j++) {
				OCP_USI id = k * np + j;
				if (myBulk.PhaseExist[id]) {
					Perf[p].transj[0] += myBulk.Kr[id] / myBulk.Mu[id];
				}		
			}
			Perf[p].transj[0] *= temp;
		}
	}
	else{
		for (USI p = 0; p < PerfNum; p++) {
			Perf[p].transj.assign(np, 0);
			OCP_USI k = Perf[p].Location;
			OCP_DBL temp = CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier;

			// multi phase
			for (USI j = 0; j < np; j++) {
				OCP_USI id = k * np + j;
				if (myBulk.PhaseExist[id]) {
					Perf[p].transj[j] = temp * myBulk.Kr[id] / myBulk.Mu[id];
				}	
			}
		}
	}
}

void Well::calFlux(const Bulk& myBulk, const bool flag)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;

	if (Opt.Type == INJ) {

		for (USI p = 0; p < PerfNum; p++) {
			Perf[p].P = BHP + dG[p];
			OCP_USI k = Perf[p].Location;
			OCP_DBL dP = Perf[p].P - myBulk.P[k];
			dP *= -1.0;

			Perf[p].qt_ft3 = Perf[p].transj[0] * dP;

			USI pvtnum = myBulk.PVTNUM[k];
			if (flag)
				Perf[p].Xi = myBulk.Flashcal[pvtnum]->xiPhase(myBulk.P[k], myBulk.T, &Opt.Zi[0]);
			OCP_DBL xi = Perf[p].Xi;
			for (USI i = 0; i < nc; i++) {
				Perf[p].qi_lbmol[i] = Perf[p].qt_ft3 * xi * Opt.Zi[i];
			}
		}
	}
	else {

		for (USI p = 0; p < PerfNum; p++) {
			Perf[p].P = BHP + dG[p];
			OCP_USI k = Perf[p].Location;
			Perf[p].qt_ft3 = 0;
			Perf[p].qi_lbmol.assign(nc, 0);

			for (USI j = 0; j < np; j++) {
				OCP_USI id = k * np + j;
				if (myBulk.PhaseExist[id]) {
					OCP_DBL dP = myBulk.Pj[id] - Perf[p].P;
					Perf[p].qt_ft3 += Perf[p].transj[j] * dP;
					//cout << p << " P[" << j << "] = " << myBulk.Pj[id] << endl;
					//cout << p << " Perf = " << Perf[p].P << endl;


					OCP_DBL xi = myBulk.Xi[id];
					OCP_DBL xij;
					for (USI i = 0; i < nc; i++) {
						xij = myBulk.Cij[id * nc + i];
						Perf[p].qi_lbmol[i] += Perf[p].transj[j] * dP * xi * xij;
					}
				}
			}
		}
	}
}

void Well::massConserve(Bulk& myBulk, const OCP_DBL& dt) const
{
	USI nc = myBulk.Nc;

	for (USI p = 0; p < PerfNum; p++) {
		OCP_USI k = Perf[p].Location;
		for (USI i = 0; i < nc; i++) {
			myBulk.Ni[k * nc + i] -= Perf[p].qi_lbmol[i] * dt;
		}
	}
}


OCP_DBL Well::calInjRate_blk(const Bulk& myBulk)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;
	OCP_DBL qj = 0;

	for (USI p = 0; p < PerfNum; p++) {

		OCP_DBL Pperf = Opt.MaxBHP + dG[p];
		OCP_USI k = Perf[p].Location;

		USI pvtnum = myBulk.PVTNUM[k];
		OCP_DBL xi = myBulk.Flashcal[pvtnum]->xiPhase(myBulk.P[k], myBulk.T, &Opt.Zi[0]);
		OCP_DBL dP = Pperf - myBulk.P[k];
		qj += Perf[p].transj[0] * xi * dP;
	}
	return qj;
}

OCP_DBL Well::calProdRate_blk(const Bulk& myBulk)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;
	OCP_DBL qj = 0;

	for (USI p = 0; p < PerfNum; p++) {

		OCP_DBL Pperf = Opt.MinBHP + dG[p];
		OCP_USI k = Perf[p].Location;

		for (USI j = 0; j < np; j++) {
			OCP_USI id = k * np + j;
			if (myBulk.PhaseExist[id]) {	
				OCP_DBL temp = 0;
				for (USI i = 0; i < nc; i++) {
					temp += Opt.Zi[i] * myBulk.Cij[id * nc + i];
				}
				OCP_DBL xi = myBulk.Xi[id];
				OCP_DBL dP = myBulk.Pj[id] - Pperf;
				qj += Perf[p].transj[j] * xi * dP * temp;
			}
		}
	}
	return qj;
}

void Well::calInjqi_blk(const Bulk& myBulk, const OCP_DBL& dt)
{
	USI nc = myBulk.Nc;
	OCP_DBL qj = 0;

	for (USI p = 0; p < PerfNum; p++) {

		Perf[p].P = BHP + dG[p];
		OCP_USI k = Perf[p].Location;

		//OCP_DBL xi = Perf[p].Xi;
		//OCP_DBL dP = Perf[p].P - myBulk.P[k];
		//qj += Perf[p].transj[0] * xi * dP;

		for (USI i = 0; i < nc; i++)
			qj += Perf[p].qi_lbmol[i];
	}
	if (Opt.FluidType == WATER) {
		WWIR = -qj;
		WWIT += WWIR * dt;
	}
	else {
		WGIR = -qj;
		WGIT += WGIR * dt;
	}
}

void Well::calProdqi_blk(const Bulk& myBulk, const OCP_DBL& dt)
{
	USI np = myBulk.Np;
	USI nc = myBulk.Nc;

	Qi_lbmol.assign(nc, 0);

	for (USI p = 0; p < PerfNum; p++) {

		//Perf[p].qi_lbmol.assign(nc, 0);
		//Perf[p].P = BHP + dG[p];
		//int k = Perf[p].Location;

		//for (int j = 0; j < np; j++) {
		//	int id = k * np + j;
		//	if (myBulk.PhaseExist[id]) {
		//		OCP_DBL xi = myBulk.Xi[id];
		//		OCP_DBL dP = myBulk.Pj[id] - Perf[p].P;
		//		OCP_DBL xij;
		//		for (int i = 0; i < nc; i++) {
		//			xij = myBulk.Cij[id * nc + i];
		//			Perf[p].qi_lbmol[i] += Perf[p].transj[j] * xi * xij * dP;
		//		}
		//	}
		//}
		for (USI i = 0; i < nc; i++) {
			Qi_lbmol[i] += Perf[p].qi_lbmol[i];
		}
	}

	for (USI i = 0; i < nc; i++) {
		if (myBulk.PhaseLabel[i] == OIL) {
			WOPR = Qi_lbmol[i];
			WOPT += WOPR * dt;
		}
		else if (myBulk.PhaseLabel[i] == GAS) {
			WGPR = Qi_lbmol[i];
			WGPT += WGPR * dt;
		}
		else if (myBulk.PhaseLabel[i] == WATER) {
			WWPR = Qi_lbmol[i];
			WWPT += WWPR * dt;
		}
	}
}

void Well::checkOptMode(const Bulk& myBulk)
{
	if (Opt.Type == INJ) {
		if (calInjRate_blk(myBulk) > Opt.MaxRate) {
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

OCP_INT Well::checkP(const Bulk& myBulk)
{
	// 0 : all correct
	// 1 : negative P    
	// 2 : outlimited P
	// 3 : crossflow happens

	if (BHP < 0) 
		return 1;
	for (USI p = 0; p < PerfNum; p++) {
		if (Perf[p].State == OPEN && Perf[p].P < 0) {
#ifdef _DEBUG
			cout << "WARNING: Well " << Name << " Perf[" << p << "].P = " << Perf[p].P << endl;
#endif // _DEBUG
			return 1;
		}
			
	}


	if (Opt.Type == INJ) {
		if (Opt.OptMode != BHP_MODE && BHP > Opt.MaxBHP) {
			Opt.OptMode = BHP_MODE;
			return 2;
		}
	}
	else {
		if (Opt.OptMode != BHP_MODE && BHP < Opt.MinBHP) {
			Opt.OptMode = BHP_MODE;
			return 2;
		}
	}
	OCP_INT flag = 0;
	flag = checkCrossFlow(myBulk);

	return flag;
}

OCP_INT Well::checkCrossFlow(const Bulk& myBulk)
{
	OCP_USI  k;
	bool flagC = true;

	if (Opt.Type == PROD) {
		USI np = myBulk.Np;
		for (USI p = 0; p < PerfNum; p++) {
			k = Perf[p].Location;
			OCP_DBL minP = myBulk.P[k];
			// THINK MORE !!!
			//for (int j = 0; j < np; j++) {
			//	if (myBulk.PhaseExist[k * np + j])
			//		minP = minP < myBulk.Pj[k * np + j] ? minP : myBulk.Pj[k * np + j];
			//}
			if (Perf[p].State == OPEN && minP < Perf[p].P) {
				Perf[p].State = CLOSE;
				Perf[p].Multiplier = 0;
				flagC = false;
			}
			else if (Perf[p].State == CLOSE && minP > Perf[p].P) {
				Perf[p].State = OPEN;
				Perf[p].Multiplier = 1;
			}
		}
	}
	else {
		for (USI p = 0; p < PerfNum; p++) {
			k = Perf[p].Location;
			if (Perf[p].State == OPEN && myBulk.P[k] > Perf[p].P) {
				Perf[p].State = CLOSE;
				Perf[p].Multiplier = 0;
				flagC = false;
			}
			else if (Perf[p].State == CLOSE && myBulk.P[k] < Perf[p].P) {
				Perf[p].State = OPEN;
				Perf[p].Multiplier = 1;
			}
		}
	}


	bool flag = false;
	// check well --  if all perf are closed, open the depthest Perf
	for (USI p = 0; p < PerfNum; p++) {
		if (Perf[p].State == OPEN) {
			flag = true;
			break;
		}
			
	}
	if (!flag) {
		return 1;
		// open the depthest Perf
		Perf.back().State = OPEN;
		Perf.back().Multiplier = 1;
	}

	if (!flagC) {
		dG = ldG;
		calTrans(myBulk);
		calFlux(myBulk);
		caldG(myBulk);
		smoothdG();
		// checkOptMode(myBulk);
		return 3;
	}

	return 0;
}

void Well::showPerfStatus() const
{
	cout << "----------------------------" << endl;
	cout << Name << ":    " << Opt.OptMode << endl; 
	for (USI p = 0; p < PerfNum; p++) {
		cout << "Perf[" << p << "].State = " << Perf[p].State << endl;
	}
}
