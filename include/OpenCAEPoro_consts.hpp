#pragma once
#include <iostream>

#define ERRORcheck(exp) \
std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in " << __FILE__ ;


// general consts
const double TINY			= 1E-8;
const double PI				= 3.141592653;

// pysical consts
const double GRAVITY_FACTOR	= 0.00694444;	// 0.00694444 ft2 psi / lb
const double RHOW_STD		= 62.3664;		// lb / ft3
const double RHOAIR_STD		= 0.076362;		// lb / ft3
const double PRESSURE_STD	= 14.6959;		// psia   =   1 atm

// Units consts
const double CONV1			= 5.61458;		// 1 bbl = 5.61458 ft3
const double CONV2			= 1.12712;		// Darcy constant in Field * 1E3

// Mixture Type
const int BLKOIL			= 1;
const int EoS_PVTW			= 2;

// Phase
const int PHASE_W			= 1;
const int PHASE_GW			= 2;
const int PHASE_OW			= 3;
const int PHASE_OGW			= 4;
const int PHASE_OG			= 5;


// Well params
const int INJ			= 0;
const int PROD			= 1;
const int CLOSE			= 0;
const int OPEN			= 1;
const int HORIZONTAL	= 0;
const int VERTICAL		= 1;
// Well opt param
const int ORATE_MODE	= 0;
const int GRATE_MODE	= 1;
const int WRATE_MODE	= 2;
const int LRATE_MODE	= 3;
const int BHP_MODE		= 4;

