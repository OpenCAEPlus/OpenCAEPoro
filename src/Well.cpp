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
	Opt = OptSet[0];
	PerfNum = K2 - K1 + 1;
	dG.resize(PerfNum, 0);
	Perf.resize(PerfNum);
	for (int p = 0; p < PerfNum; p++) {
		Perf[p].State = OPEN;
		int Idg = (K1 + p) * myGrid.Nx * myGrid.Ny + J * myGrid.Nx + I;
		Perf[p].Location = myGrid.ActiveMap_G2B[Idg];
		Perf[p].Depth = myBulk.Depth[Perf[p].Location];
		Perf[p].Multiplier = 1;
	}
	calWI_Peaceman_Vertical(myBulk);
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

void Well::allocateMat(Solver& mySolver)
{
	for (int p = 0; p < PerfNum; p++) {
		mySolver.RowCapacity[Perf[p].Location]++;
	}
}

void Well::assembleMat_INJ(const Bulk& myBulk, Solver& mySolver)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	int wId = mySolver.Dim;
	// important !
	mySolver.Dim++;

	for (int p = 0; p < PerfNum; p++) {
		int k = Perf[p].Location; 
		double xi = 0;

		double trans = 0;
		for (int j = 0; j < np; j++) {
			if (myBulk.PhaseExist[k * np + j]) {
				trans += myBulk.Kr[k * np + j] / myBulk.Mu[k * np + j];
			}		
		}

		double Vfi_zi = 0;
		for (int i = 0; i < nc; i++) {
			Vfi_zi += myBulk.Vfi[k * nc + i] * Opt.Zi[i];
		}
		
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
			mySolver.b[wId] += bw;
			break;
		case BHP_MODE:
			mySolver.ColId[wId].push_back(k);
			mySolver.Val[wId].push_back(0);
			break;
		default:
			ERRORcheck("Wrong Well Opt mode in function");
			exit(0);
		}
	}

	// Well Self
	assert(mySolver.Val[wId].size() == PerfNum);
	// the order of perforation is not necessarily in order
	switch (Opt.OptMode)
	{
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(mySolver.DiagVal[wId]);
		// b
		mySolver.b[wId] += 1 * Opt.OptValue;
		break;
	case BHP_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(1);
		// b
		mySolver.b[wId] += 1 * Opt.OptValue;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode in function");
		exit(0);
	}
}


void Well::assembleMat_PROD_BLK(const Bulk& myBulk, Solver& mySolver)
{
	int np = myBulk.Np;
	int nc = myBulk.Nc;
	int wId = mySolver.Dim;
	// important !
	mySolver.Dim++;

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
			bb += valb * dP;
			bw += valw * dP;
		}
		double trans = 1 * CONV1 * CONV2 * Perf[p].WI * Perf[p].Multiplier;
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
			ERRORcheck("Wrong Well Opt mode in function");
			exit(0);
		}
	}

	// Well Self
	assert(mySolver.Val[wId].size() == PerfNum);
	// the order of perforation is not necessarily in order
	switch (Opt.OptMode)
	{
	case ORATE_MODE:
	case GRATE_MODE:
	case WRATE_MODE:
	case LRATE_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(mySolver.DiagVal[wId]);
		// b
		mySolver.b[wId] += 1 * Opt.OptValue;
		break;
	case BHP_MODE:
		// diag
		mySolver.ColId[wId].push_back(PerfNum);
		mySolver.DiagPtr[wId] = PerfNum;
		mySolver.Val[wId].push_back(1);
		// b
		mySolver.b[wId] += 1 * Opt.OptValue;
		break;
	default:
		ERRORcheck("Wrong Well Opt mode in function");
		exit(0);
	}
}



