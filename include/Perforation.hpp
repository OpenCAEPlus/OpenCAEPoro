#pragma once
#include <vector>

class Perforation
{
	friend class Well;
public:
	 Perforation() = default;
	 void setState(bool flag) { State = flag; };

private:
	int						State;
	int						Location;
	double					Depth;
	double					Trans;
	double					P;
	double					WI;
	double					Multiplier;

	std::vector<double>		qi_lbmol;
};
