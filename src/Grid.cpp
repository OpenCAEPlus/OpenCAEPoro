#include "Grid.hpp"

void Grid::setup()
{
	calDepthV();
	calActiveBulk(1E-6, 1E-6);
}

void Grid::calDepthV()
{
	Depth.resize(Num, 0);
	OCP_USI nxny = Nx * Ny;
	// 0th layer
	for (USI j = 0; j < Ny; j++) {
		for (USI i = 0; i < Nx; i++) {
			OCP_USI id = j * Nx + i;
			Depth[id] = Tops[id] + Dz[id] / 2;
		}
	}
	// 1th - (Nz-1)th layer
	for (USI k = 1; k < Nz; k++) {
		OCP_USI knxny = k * nxny;
		for (USI j = 0; j < Ny; j++) {
			for (USI i = 0; i < Nx; i++) {
				OCP_USI id = knxny + j * Nx + i;
				Depth[id] = Depth[id - nxny] + Dz[id - nxny] / 2 + Dz[id] / 2;
			}
		}
	}

	V.resize(Num);
	for (OCP_USI i = 0; i < Num; i++)
		V[i] = Dx[i] * Dy[i] * Dz[i];
	cout << "Grid::calDepthV" << endl;
}

void Grid::calActiveBulk(OCP_DBL e1, OCP_DBL e2)
{
	ActiveMap_B2G.reserve(Num);
	ActiveMap_G2B.resize(Num);
	OCP_USI count = 0;
	for (OCP_USI n = 0; n < Num; n++) {
		if (Poro[n] * Ntg[n] < e1 || Dx[n] * Dy[n] * Dz[n] < e2) {
			ActiveMap_G2B[n] = GB_Pair(false, 0);
			continue;
		}
		ActiveMap_B2G.push_back(n);
		ActiveMap_G2B[n] = GB_Pair(true, count);
		count++;
	}
	ActiveBulkNum = count;
	cout << (Num - ActiveBulkNum) / Num << "%  grids is inactive" << endl;
	cout << "Grid::calActiveBulk" << endl;
}


void Grid::inputParam(const ParamReservoir& rs_param)
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
		for (OCP_USI i = 0; i < Num; i++) {
			SATNUM[i] = (USI)(rs_param.SATNUM.data[i]) - 1;
		}
	}
	PVTNUM.resize(Num, 0);
	if (rs_param.PVTNUM.activity) {
		for (OCP_USI i = 0; i < Num; i++) {
			PVTNUM[i] = (USI)(rs_param.PVTNUM.data[i]) - 1;
		}
	}
	cout << "Grid::inputParam" << endl;
}


OCP_USI Grid::getIndex(USI i, USI j, USI k)
{
	OCP_USI id = k * Nx * Ny + j * Nx + i;
	bool activity = ActiveMap_G2B[id].getAct();
	if (!activity) {
		ERRORcheck("(" + to_string(i) + "," + to_string(j) + "," + to_string(k) + ") is inactive");
		exit(0);
	}
	id = ActiveMap_G2B[id].getId();
	return id;
}
