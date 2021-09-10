#pragma once
#include <vector>
#include <iostream>

using namespace std;

class Grid
{	
	friend class Bulk;
	friend class Connection_BB;
	friend class Well;
public:
	Grid() = default;
	Grid(int nx, int ny, int nz) :Nx(nx), Ny(ny), Nz(nz) {
		Num = nx*ny*nz;
		ConnNum = 3*Nx*Ny*Nz - Nx*Ny - Ny*Nz - Nz*Nx;
	}
	int getBulkNum() { return Num; }
	int getConnNum() { return ConnNum; }
	int getActiveBulkNum() { return ActiveBulkNum; }


	void calActiveBulk(double e1, double e2);		// fill ActiveMap_B2G and ActiveMap_G2B


private:
	int					Nx;					// num of bulks along x-aixs
	int					Ny;					// num of bulks along y-aixs
	int					Nz;					// num of bulks along z-aixs
	int					Num;				// Nx * Ny * Nz
	int					ConnNum;			// num of connection


	vector<double>		Depth;				// depth: Num
	vector<double>		Dx;					// dx: Num
	vector<double>		Dy;					// dy: Num
	vector<double>		Dz;					// dz: Num
	vector<double>		V;					// volume : Num
	vector<double>		Ntg;				// net to gross
	vector<double>		Poro;				// initial porosity
	vector<double>		Kx;					// Absolute permeability in x direction
	vector<double>		Ky;					// Absolute permeability in y direction
	vector<double>		Kz;					// Absolute permeability in z direction

	int					ActiveBulkNum;
	vector<int>			ActiveMap_B2G;		// size: Active Num
	vector<int>			ActiveMap_G2B;		// size: Num
	
};
