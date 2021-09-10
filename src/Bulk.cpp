#include "Bulk.hpp"
#include <ctime>
#include <algorithm>


void Bulk::init(const Grid& myGrid)
{
	Num = myGrid.ActiveBulkNum;
	Dx.resize(Num, 0);
	Dy.resize(Num, 0);
	Dz.resize(Num, 0); 
	Depth.resize(Num, 0);
	Ntg.resize(Num, 0);
	Rock_V.resize(Num, 0);
	Rock_PoroInit.resize(Num, 0);
	Rock_KxInit.resize(Num, 0);
	Rock_KyInit.resize(Num, 0);
	Rock_KzInit.resize(Num, 0);

	for (int bIdb = 0; bIdb < Num; bIdb++) {
		int bIdg = myGrid.ActiveMap_B2G[bIdb];

		Dx[bIdb] = myGrid.Dx[bIdg];
		Dy[bIdb] = myGrid.Dy[bIdg];
		Dz[bIdb] = myGrid.Dz[bIdg];
		Depth[bIdb] = myGrid.Depth[bIdg];
		Ntg[bIdb] = myGrid.Ntg[bIdg];

		Rock_V[bIdb] = myGrid.V[bIdg] * Ntg[bIdb];
		Rock_PoroInit[bIdb] = myGrid.Poro[bIdg];
		Rock_KxInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KyInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KzInit[bIdb] = myGrid.Kx[bIdg];
	}
}


void Bulk::initSjPc()
{

}

// Flash
void Bulk::flash_Sj()
{
	for (int n = 0; n < Num; n++) {
		Flashcal->Flash_Sj(P[n], Pbb[n], T[n], &S[n * Np], Rock_V[n] * Rock_Poro[n], &Ni[n * Nc]);
		for (int i = 0; i < Nc; i++) {
			Ni[n * Nc + i] = Flashcal->Ni[i];
		}
		passValue(n);
	}
	
}

void Bulk::flash_Ni()
{
	for (int n = 0; n < Num; n++) {
		Flashcal->Flash_Ni(P[n], T[n], &Ni[n * Nc]);
		passValue(n);
	}
}

void Bulk::passValue(int n)
{
	int bId = n * Np;
	for (int j = 0; j < Np; j++) {
		PhaseExist[bId + j] = Flashcal->PhaseExist[j];
		if (PhaseExist[j]) {
			for (int i = 0; i < Nc; i++) {
				Cij[bId + j * Nc + i] = Flashcal->Cij[j * Nc + i];
			}
			S[bId + j]		= Flashcal->S[j];
			Xi[bId + j]		= Flashcal->Xi[j];
			Rho[bId + j]	= Flashcal->Rho[j];
			Mu[bId + j]		= Flashcal->Mu[j];
		}
	}
	Vf[n] = Flashcal->Vf;
	Vfp[n] = Flashcal->Vfp;
	bId = n * Nc;
	for (int i = 0; i < Nc; i++) {
		Vfi[bId + i] = Flashcal->Vfi[i];
	}
}

// relative permeability and capillary pressure
void Bulk::calKrPc()
{
	for (int n = 0; n < Num; n++) {
		int bId = n * Np;
		Flow->calKrPc(&S[bId], &Kr[bId], &Pc[bId]);
	}
}

// Rock
void Bulk::calporo()
{
	for (int n = 0; n < Num; n++) {
		double dP = P[n] - Rock_Pref;
		Rock_Poro[n] = Rock_PoroInit[n] * (1 + Rock_C1 * dP + Rock_C2 * dP * dP);
	}
}
