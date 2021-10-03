#pragma once
#include <vector>
#include <iostream>
#include "ParamReservoir.hpp"
#include "OpenCAEPoro_consts.hpp"

using namespace std;

class Grid
{	
	friend class Bulk;
	friend class Connection_BB;
	friend class Well;
public:
	Grid() = default;

	void setup();

	int getBulkNum() { return Num; }
	int getConnNum() { return ConnNum; }
	int getActiveBulkNum() { return ActiveBulkNum; }

	void calDepthV();
	void calActiveBulk(OCP_DBL e1, OCP_DBL e2);		// fill ActiveMap_B2G and ActiveMap_G2B

	void inputParam(ParamReservoir& rs_param);

	int getIndex(int i, int j, int k);

private:
	int					Nx;					// num of bulks along x-aixs
	int					Ny;					// num of bulks along y-aixs
	int					Nz;					// num of bulks along z-aixs
	int					Num;				// Nx * Ny * Nz
	int					ConnNum;			// num of connection


	vector<OCP_DBL>		Tops;
	vector<OCP_DBL>		Depth;				// depth: Num
	vector<OCP_DBL>		Dx;					// dx: Num
	vector<OCP_DBL>		Dy;					// dy: Num
	vector<OCP_DBL>		Dz;					// dz: Num
	vector<OCP_DBL>		V;					// volume : Num
	vector<OCP_DBL>		Ntg;				// net to gross
	vector<OCP_DBL>		Poro;				// initial porosity
	vector<OCP_DBL>		Kx;					// Absolute permeability in x direction
	vector<OCP_DBL>		Ky;					// Absolute permeability in y direction
	vector<OCP_DBL>		Kz;					// Absolute permeability in z direction

	// Region
	vector<int>			SATNUM;
	vector<int>			PVTNUM;

	int					ActiveBulkNum;
	vector<int>			ActiveMap_B2G;		// size: Active Num
	vector<int>			ActiveMap_G2B;		// size: Num
	
};
