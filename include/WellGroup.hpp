#pragma once
#include "Well.hpp"
#include "ParamWell.hpp"

using namespace std;

class WellGroup
{

public:
	WellGroup() = default;
	
	void inputParam(ParamWell& Well_param);
	void setup(Bulk& myBulk);
	void setupPerf();
	void setupMixture(Bulk& myBulk);

	void allocateMat(Solver& mySolver);
	void assemblaMat_WB(Solver& mySolver, const Bulk& myBulk);
	int getWellNum() { return WellNum; };

private:
	int							WellNum;
	std::vector<Well>			WellG;
	std::vector<Mixture*>		Flashcal;
};