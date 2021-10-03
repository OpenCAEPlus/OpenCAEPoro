#pragma once
#include "Well.hpp"
#include "ParamWell.hpp"

using namespace std;

class WellGroup
{
public:
	WellGroup() = default;
	
	void inputParam(ParamWell& Well_param);
	void setup(Grid& myGrid, Bulk& myBulk);
	void setupWell(Grid& myGrid, Bulk& myBulk);
	void setupMixture(Bulk& myBulk);

	OCP_DBL calCFL(const Bulk& myBulk, OCP_DBL dt);

	void calFlux(const Bulk& myBulk);
	void massConserve(Bulk& myBulk, OCP_DBL dt);

	template<typename T>
	void allocateMat(Solver<T>& mySolver);

	void init(const Bulk& myBulk);
	void applyControl(int i);

	void prepareWell(const Bulk& myBulk);

	void calIPRT(const Bulk& myBulk, OCP_DBL dt);

	void assemblaMat_WB(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt);
	int getWellNum() { return WellNum; }
	string& getWellName(int i) { return WellG[i].Name; }
	int getIndex(string name);
	// Field
	OCP_DBL getFOPR() { return FOPR; }
	OCP_DBL getFOPT() { return FOPT; }
	OCP_DBL getFGPR() { return FGPR; }
	OCP_DBL getFGPT() { return FGPt; }
	OCP_DBL getFWPR() { return FWPR; }
	OCP_DBL getFWPT() { return FWPT; }
	OCP_DBL getFGIR() { return FGIR; }
	OCP_DBL getFGIT() { return FGIT; }
	OCP_DBL getFWIR() { return FWIR; }
	OCP_DBL getFWIT() { return FWIT; }

	// Well
	OCP_DBL getWOPR(int w) { return WellG[w].WOPR; }
	OCP_DBL getWOPT(int w) { return WellG[w].WOPT; }
	OCP_DBL getWGPR(int w) { return WellG[w].WGPR; }
	OCP_DBL getWGPT(int w) { return WellG[w].WGPT; }
	OCP_DBL getWWPR(int w) { return WellG[w].WWPR; }
	OCP_DBL getWWPT(int w) { return WellG[w].WWPT; }
	OCP_DBL getWGIR(int w) { return WellG[w].WGIR; }
	OCP_DBL getWGIT(int w) { return WellG[w].WGIT; }
	OCP_DBL getWWIR(int w) { return WellG[w].WWIR; }
	OCP_DBL getWWIT(int w) { return WellG[w].WWIT; }
	// BHP
	OCP_DBL getWBHP(int w) { return WellG[w].BHP; }

	void getSol_IMPES(vector<OCP_DBL>& u, int bid);

	void setLastStep() { for (auto& w : WellG)	w.ldG = w.dG; }

	int checkP(const Bulk& myBulk);

private:
	int							WellNum;
	std::vector<Well>			WellG;
	std::vector<Mixture*>		Flashcal;
	OCP_DBL						FGIR{ 0 };
	OCP_DBL						FGIT{ 0 };
	OCP_DBL						FWIR{ 0 };
	OCP_DBL						FWIT{ 0 };
	OCP_DBL						FOPR{ 0 };
	OCP_DBL						FOPT{ 0 };
	OCP_DBL						FGPR{ 0 };
	OCP_DBL						FGPt{ 0 };
	OCP_DBL						FWPR{ 0 };
	OCP_DBL						FWPT{ 0 };
};

template<typename T>
void WellGroup::allocateMat(Solver<T>& mySolver)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}
