#include "BOMixture.h"



void BOMixture::Flash_Sj(const double Pin, const double Pbbin, const double Tin, const double* Sjin, double Vpore, const double* Ziin)
{
	switch (Mode)
	{
	case PHASE_W:
		BOFlash_Sj_W(Pin, Sjin, Vpore, Ziin);
		break;

	case PHASE_OW:
		BOFlash_Sj_OW(Pin, Sjin, Vpore, Ziin);
		break;

	case PHASE_OGW:
		BOFlash_Sj_OGW(Pin, Pbbin, Sjin, Vpore, Ziin);
		break;

	default:
		ERRORcheck("Wrong Mode");
		exit(0);
		break;
	}
}

void BOMixture::Flash_Ni(const double Pin, const double Tin, const double* Niin)
{
	switch (Mode)
	{
	case PHASE_W:
		BOFlash_Ni_W(Pin, Niin);
		break;

	case PHASE_OW:
		BOFlash_Ni_OW(Pin, Niin);
		break;

	case PHASE_OGW:
		BOFlash_Ni_OGW(Pin, Niin);
		break;

	default:
		ERRORcheck("Wrong Mode");
		exit(0);
	}
}


void BOMixture::BOFlash_Ni_W(const double Pin, const double* Niin)
{

}

void BOMixture::BOFlash_Ni_OW(const double Pin, const double* Niin)
{

}

void BOMixture::BOFlash_Sj_W(const double Pin, const double* Sjin, double Vpore, const double* Ziin)
{

}

void BOMixture::BOFlash_Sj_OW(const double Pin, const double* Sjin, double Vpore, const double* Ziin)
{

}


double BOMixture::calRhoO(double Pin, double Pbbin)
{
	switch (Mode)
	{
	case PHASE_OW:
		return calRhoO_OW(Pin);
	case PHASE_OGW:
		return calRhoO_OGW(Pin, Pbbin);
	default:
		ERRORcheck("Wrong Mode");
		exit(0);
	}
}

double BOMixture::calRhoO_OW(double Pin)
{

}
