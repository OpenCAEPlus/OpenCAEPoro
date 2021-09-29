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

	double calCFL(const Bulk& myBulk, double dt);

	void calFlux(const Bulk& myBulk);
	void massConserve(Bulk& myBulk, double dt);

	template<typename T>
	void allocateMat(Solver<T>& mySolver);

	void init(const Bulk& myBulk);
	void applyControl(int i);

	void prepareWell(const Bulk& myBulk);

	void calIPRT(const Bulk& myBulk, double dt);

	void assemblaMat_WB(Solver<double>& mySolver, const Bulk& myBulk, double dt);
	int getWellNum() { return WellNum; }
	string& getWellName(int i) { return WellG[i].Name; }
	int getIndex(string name);
	// Field
	double getFOPR() { return FOPR; }
	double getFOPT() { return FOPT; }
	double getFGPR() { return FGPR; }
	double getFGPT() { return FGPt; }
	double getFWPR() { return FWPR; }
	double getFWPT() { return FWPT; }
	double getFGIR() { return FGIR; }
	double getFGIT() { return FGIT; }
	double getFWIR() { return FWIR; }
	double getFWIT() { return FWIT; }

	// Well
	double getWOPR(int w) { return WellG[w].WOPR; }
	double getWOPT(int w) { return WellG[w].WOPT; }
	double getWGPR(int w) { return WellG[w].WGPR; }
	double getWGPT(int w) { return WellG[w].WGPT; }
	double getWWPR(int w) { return WellG[w].WWPR; }
	double getWWPT(int w) { return WellG[w].WWPT; }
	double getWGIR(int w) { return WellG[w].WGIR; }
	double getWGIT(int w) { return WellG[w].WGIT; }
	double getWWIR(int w) { return WellG[w].WWIR; }
	double getWWIT(int w) { return WellG[w].WWIT; }
	// BHP
	double getWBHP(int w) { return WellG[w].BHP; }

	void getSol_IMPES(vector<double>& u, int bid);

	void setLastStep() { for (auto& w : WellG)	w.ldG = w.dG; }

	int checkP(const Bulk& myBulk);

private:
	int							WellNum;
	std::vector<Well>			WellG;
	std::vector<Mixture*>		Flashcal;
	double						FGIR{ 0 };
	double						FGIT{ 0 };
	double						FWIR{ 0 };
	double						FWIT{ 0 };
	double						FOPR{ 0 };
	double						FOPT{ 0 };
	double						FGPR{ 0 };
	double						FGPt{ 0 };
	double						FWPR{ 0 };
	double						FWPT{ 0 };
};

template<typename T>
void WellGroup::allocateMat(Solver<T>& mySolver)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}
