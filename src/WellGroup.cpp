#include "WellGroup.hpp"

void WellGroup::assemblaMat_WB(Solver& mySolver, const Bulk& myBulk)
{
	for (int w = 0; w < WellNum; w++) {
		if (WellG[w].State) {

			switch (WellG[w].Type)
			{
			case INJ:
				WellG[w].assembleMat_INJ(myBulk, mySolver);
				break;
			case PROD:
				WellG[w].assembleMat_PROD_BLK(myBulk, mySolver);
				break;
			default:
				ERRORcheck("Wrong Well Type in function");
				exit(0);
			}
		}
	}
}
