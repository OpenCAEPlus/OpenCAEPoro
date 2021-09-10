#include "Grid.h"


void Grid::calActiveBulk(double e1, double e2)
{
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
}


