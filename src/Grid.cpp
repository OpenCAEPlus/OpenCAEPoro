#include "Grid.hpp"

void Grid::setup()
{
	calDepthV();
	calActiveBulk(1E-6, 1E-6);
}

void Grid::calDepthV()
{
	Depth.resize(Num, 0);
	int nxny = Nx * Ny;
	// 0th layer
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			int id = j * Nx + i;
			Depth[id] = Tops[id] + Dz[id] / 2;
		}
	}
	// 1th - (Nz-1)th layer
	for (int k = 1; k < Nz; k++) {
		int knxny = k * nxny;
		for (int j = 0; j < Ny; j++) {
			for (int i = 0; i < Nx; i++) {
				int id = knxny + j * Nx + i;
				Depth[id] = Depth[id - nxny] + Dz[id - nxny] / 2 + Dz[id] / 2;
			}
		}
	}

	V.resize(Num);
	for (int i = 0; i < Num; i++)
		V[i] = Dx[i] * Dy[i] * Dz[i];
	cout << "Grid::calDepthV" << endl;
}

void Grid::calActiveBulk(double e1, double e2)
{
	ActiveMap_B2G.reserve(Num);
	ActiveMap_G2B.resize(Num);
	int count = 0;
	for (int n = 0; n < Num; n++) {
		if (Poro[n] * Ntg[n] < e1 || Dx[n] * Dy[n] * Dz[n] < e2) {
			ActiveMap_G2B[n] = -1;
			continue;
		}
		ActiveMap_B2G.push_back(n);
		ActiveMap_G2B[n] = count;
		count++;
	}
	ActiveBulkNum = count;
	cout << (Num - ActiveBulkNum) / Num << "%  grids is inactive" << endl;
	cout << "Grid::calActiveBulk" << endl;
}


void Grid::inputParam(ParamReservoir& rs_param)
{
	Nx = rs_param.Dimens.Nx;
	Ny = rs_param.Dimens.Ny;
	Nz = rs_param.Dimens.Nz;
	Num = rs_param.Num;
	ConnNum = 3 * Nx * Ny * Nz - Nx * Ny - Ny * Nz - Nz * Nx;

	Tops = rs_param.Tops;
	Dx = rs_param.Dx;
	Dy = rs_param.Dy;
	Dz = rs_param.Dz;

	Ntg = rs_param.Ntg;
	Poro = rs_param.Poro;
	Kx = rs_param.PermX;
	Ky = rs_param.PermY;
	Kz = rs_param.PermZ;

	SATNUM.resize(Num, 0);
	if (rs_param.SATNUM.activity) {
		for (int i = 0; i < Num; i++) {
			SATNUM[i] = (int)(rs_param.SATNUM.data[i]) - 1;
		}
	}
	PVTNUM.resize(Num, 0);
	if (rs_param.PVTNUM.activity) {
		for (int i = 0; i < Num; i++) {
			PVTNUM[i] = (int)(rs_param.PVTNUM.data[i]) - 1;
		}
	}
	cout << "Grid::inputParam" << endl;
}
