#include "WellGroup.hpp"

void WellGroup::inputParam(ParamWell& Well_param)
{
	WellNum = Well_param.WellNum;
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
		WellG[w].I = Well_param.well[w].I;
		WellG[w].J = Well_param.well[w].J;
		WellG[w].K1 = Well_param.well[w].K1;
		WellG[w].K2 = Well_param.well[w].K2;
		
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


void WellGroup::allocateMat(Solver& mySolver)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}

void WellGroup::assemblaMat_WB(Solver& mySolver, const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].Opt.State) {

			switch (WellG[w].Opt.Type)
			{
			case INJ:
				WellG[w].assembleMat_INJ(myBulk, mySolver);
				break;
			case PROD:
				WellG[w].assembleMat_PROD_BLK(myBulk, mySolver);
				break;
			default:
				ERRORcheck("Wrong Well Type in function");
				exit(0);
			}
		}
	}
}
