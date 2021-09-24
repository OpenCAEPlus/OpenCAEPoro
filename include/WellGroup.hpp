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

	template<typename T>
	void allocateMat(Solver<T>& mySolver);

	void init(const Bulk& myBulk);
	void applyControl(int i);

	void checkOptMode(const Bulk& myBulk);
	void calWelldG(const Bulk& myBulk);
	void prepareWell(const Bulk& myBulk);

	void assemblaMat_WB(Solver<double>& mySolver, const Bulk& myBulk, double dt);
	int getWellNum() { return WellNum; }
	string getWellName(int i) { return WellG[i].Name; }

	void getP_IMPES(vector<double>& u, int bid);

private:
	int							WellNum;
	std::vector<Well>			WellG;
	std::vector<Mixture*>		Flashcal;
};

template<typename T>
void WellGroup::allocateMat(Solver<T>& mySolver)
{
	for (int w = 0; w < WellNum; w++) {
		WellG[w].allocateMat(mySolver);
	}
}
