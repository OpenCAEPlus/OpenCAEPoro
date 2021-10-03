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
	void allocateMat(Solver<T>& mySolver) const;

	void init(const Bulk& myBulk);
	void applyControl(int i);

	void prepareWell(const Bulk& myBulk);

	void calIPRT(const Bulk& myBulk, OCP_DBL dt);

	void assemblaMat_WB(Solver<OCP_DBL>& mySolver, const Bulk& myBulk, OCP_DBL dt) const;
	int getWellNum() const { return WellNum; }
	string& getWellName(int i) { return WellG[i].Name; }
	int getIndex(string name);
	// Field
	OCP_DBL getFOPR() const { return FOPR; }
	OCP_DBL getFOPT() const { return FOPT; }
	OCP_DBL getFGPR() const { return FGPR; }
	OCP_DBL getFGPT() const { return FGPt; }
	OCP_DBL getFWPR() const { return FWPR; }
	OCP_DBL getFWPT() const { return FWPT; }
	OCP_DBL getFGIR() const { return FGIR; }
	OCP_DBL getFGIT() const { return FGIT; }
	OCP_DBL getFWIR() const { return FWIR; }
	OCP_DBL getFWIT() const { return FWIT; }

	// Well
	OCP_DBL getWOPR(int w) const { return WellG[w].WOPR; }
	OCP_DBL getWOPT(int w) const { return WellG[w].WOPT; }
	OCP_DBL getWGPR(int w) const { return WellG[w].WGPR; }
	OCP_DBL getWGPT(int w) const { return WellG[w].WGPT; }
	OCP_DBL getWWPR(int w) const { return WellG[w].WWPR; }
	OCP_DBL getWWPT(int w) const { return WellG[w].WWPT; }
	OCP_DBL getWGIR(int w) const { return WellG[w].WGIR; }
	OCP_DBL getWGIT(int w) const { return WellG[w].WGIT; }
	OCP_DBL getWWIR(int w) const { return WellG[w].WWIR; }
	OCP_DBL getWWIT(int w) const { return WellG[w].WWIT; }
	// BHP
	OCP_DBL getWBHP(int w) const { return WellG[w].BHP; }

	void getSol_IMPES(const vector<OCP_DBL>& u, int bid);

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
void WellGroup::allocateMat(Solver<T>& mySolver) const
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}
