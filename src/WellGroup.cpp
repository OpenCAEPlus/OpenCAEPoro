#include "WellGroup.hpp"

void WellGroup::inputParam(ParamWell& Well_param)
{
	WellNum = Well_param.well.size();
	WellG.resize(WellNum);
	int t = Well_param.CriticalTime.size();
	vector<int>		wellCriticalTime;
	for (int w = 0; w < WellNum; w++) {
		WellG[w].Name = Well_param.well[w].Name;
		WellG[w].Direction = Well_param.well[w].Direction;
		WellG[w].Depth = Well_param.well[w].Dref;
		WellG[w].Radius = Well_param.well[w].Diameter / 2;
		WellG[w].Kh = Well_param.well[w].Kh;
		WellG[w].SkinFactor = Well_param.well[w].SkinFactor;
		WellG[w].WI = Well_param.well[w].WI;
		WellG[w].I = Well_param.well[w].I - 1;
		WellG[w].J = Well_param.well[w].J - 1;
		WellG[w].K1 = Well_param.well[w].K1 - 1;
		WellG[w].K2 = Well_param.well[w].K2 - 1;
		
		// Opt
		WellG[w].OptSet.resize(t);
		int n = Well_param.well[w].OptParam.size();
		wellCriticalTime.clear();
		wellCriticalTime.resize(n + 1);
		for (int i = 0; i < n; i++) {
			wellCriticalTime[i] = Well_param.well[w].OptParam[i].d;
		}
		wellCriticalTime.back() = t;
		for (int i = 0; i < n; i++) {
			for (int d = wellCriticalTime[i]; d < wellCriticalTime[i + 1]; d++) {
				WellG[w].OptSet[d] = WellOpt(Well_param.well[w].OptParam[i].Opt);
			}
		}
	}

	cout << "WellGroup::inputParam" << endl;
}

void WellGroup::setup(Grid& myGrid, Bulk& myBulk)
{
	setupWell(myGrid, myBulk);
	setupMixture(myBulk);
}

void WellGroup::setupWell(Grid& myGrid, Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].setup(myGrid, myBulk);
	}
}

void WellGroup::setupMixture(Bulk& myBulk)
{

	Flashcal = myBulk.getMixture();
	cout << "WellGroup::setupMixture" << endl;
}



void WellGroup::init(const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].init(myBulk);
	}
}

void WellGroup::applyControl(int i)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].Opt = WellG[w].OptSet[i];
	}
}

OCP_DBL WellGroup::calCFL(const Bulk& myBulk, OCP_DBL dt)
{
	OCP_DBL cflw = 0;
	OCP_DBL tmp = 0;
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			tmp = WellG[w].calCFL(myBulk, dt);
			if (cflw < tmp)
				cflw = tmp;
		}
	}
	return cflw;
}

void WellGroup::calFlux(const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			WellG[w].calFlux(myBulk);
		}
	}
}

void WellGroup::massConserve(Bulk& myBulk, OCP_DBL dt)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			WellG[w].massConserve(myBulk, dt);
		}
	}
}

void WellGroup::prepareWell(const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {

			WellG[w].calTrans(myBulk);
			WellG[w].calFlux(myBulk, true);
			WellG[w].caldG(myBulk);
			// test
			WellG[w].smoothdG();
			WellG[w].checkOptMode(myBulk);
		}
	}
}

void WellGroup::assemblaMat_WB(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt) const
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {

			switch (WellG[w].WellType())
			{
			case INJ:
				WellG[w].assembleMat_INJ(myBulk, mySolver, dt);
				break;
			case PROD:
				WellG[w].assembleMat_PROD_BLK(myBulk, mySolver, dt);
				break;
			default:
				ERRORcheck("Wrong Well Type in function");
				exit(0);
			}
		}
	}
}

void WellGroup::getSol_IMPES(const vector<OCP_DBL>& u, int bid)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			WellG[w].BHP = u[bid + w];
            WellG[w].updatePerfP();
		}
	}
}

void WellGroup::calIPRT(const Bulk& myBulk, OCP_DBL dt)
{
	FGIR = 0;
	FWIR = 0;
	FOPR = 0;
	FGPR = 0;
	FWPR = 0;
	for (int w = 0; w < WellNum; w++) {
		WellG[w].WGIR = 0;
		WellG[w].WWIR = 0;
		WellG[w].WOPR = 0;
		WellG[w].WGPR = 0;
		WellG[w].WWPR = 0;

		if (WellG[w].WellState()) {
			if (WellG[w].WellType() == PROD) {
				WellG[w].calProdqi_blk(myBulk, dt);
			}
			else {
				WellG[w].calInjqi_blk(myBulk, dt);
			}
		}
		FGIR += WellG[w].WGIR;
		FWIR += WellG[w].WWIR;
		FOPR += WellG[w].WOPR;
		FGPR += WellG[w].WGPR;
		FWPR += WellG[w].WWPR;
	}
	FGIT += FGIR * dt;
	FWIT += FWIR * dt;
	FOPT += FOPR * dt;
	FGPt += FGPR * dt;
	FWPT += FWPR * dt;
}


int WellGroup::checkP(const Bulk& myBulk)
{
	// 0 : All correct
	// 1   : negative P, cut the timestep and resolve
	// 2.1 : change well mode to BHP, resolve
	// 2.2 : crossflow happens, then close corresponding Perf, resolve
	bool flag2 = false;
	bool flag3 = false;

	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {

			int flag = WellG[w].checkP(myBulk);
#ifdef _DEBUG
			WellG[w].showPerfStatus();
#endif // _DEBUG
			switch (flag)
			{
			case 1:
				return 1;
			case 2:
				flag2 = true;
				break;
			case 3:
				flag3 = true;
				break;
			default:
				break;
			}
		}
	}

	if (flag2 || flag3)
		return 2;

	return 0;
}



// return the index of Specified well name
int WellGroup::getIndex(string name)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].Name == name) {
			return w;
		}
	}
	ERRORcheck("No such well name!");
	exit(0);
}

