#pragma once
#include "Well.hpp"

using namespace std;

class WellGroup
{

public:
	WellGroup() = default;
	WellGroup(int nw) :WellNum(nw) { WellG.reserve(WellNum); };

	void assemblaMat_WB(Solver& mySolver, const Bulk& myBulk);
	int getWellNum() { return WellNum; };

private:
	int						WellNum;
	std::vector<Well>		WellG;
	Mixture*				Flash;
};