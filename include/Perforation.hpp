#pragma once
#include <vector>

class Perforation
{
	friend class Well;
public:
	 Perforation() = default;
	 void setState(bool flag) { State = flag; };

private:
	bool					State;
	int						Location;
	double					Depth;
	double					Trans;
	double					P;
	double					WI;
	double					Multiplier;

	double					Xi;			// inj  single phase
	std::vector<double>		qi_lbmol;
	std::vector<double>		transj;
	double					qt_ft3;
};
