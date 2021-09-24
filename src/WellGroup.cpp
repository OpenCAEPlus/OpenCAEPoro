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

void WellGroup::checkOptMode(const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			WellG[w].checkOptMode(myBulk);
		}
	}
}

void WellGroup::calWelldG(const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			if (WellG[w].WellType() == PROD) {
				WellG[w].calProddG(myBulk);
			}
			else {
				WellG[w].calInjdG(myBulk);
			}
		}
	}
}

void WellGroup::prepareWell(const Bulk& myBulk)
{
	calWelldG(myBulk);
	checkOptMode(myBulk);
}

void WellGroup::assemblaMat_WB(Solver<double>& mySolver, const Bulk& myBulk, double dt)
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

void WellGroup::getP_IMPES(vector<double>& u, int bid)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].WellState()) {
			WellG[w].BHP = u[bid + w];
		}
	}
}
