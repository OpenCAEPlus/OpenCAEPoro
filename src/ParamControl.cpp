/*! \file    ParamControl.cpp
 *  \brief   ParamControl class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamControl.hpp"

void ParamControl::Init(string& indir)
{
	dir = indir;
	InitMethod();
	InitTime();
	InitTuning();
}

void ParamControl::InitMethod()
{
	method = "IMPES";
	linearSolve = "./csr.dat";
}

void ParamControl::InitTuning()
{
	tuning.resize(3);

	// Timestepping controls
	tuning[0].resize(10);
	tuning[0][0] = 1.0;			// Maximum Init step length of next timestep
	tuning[0][1] = 365.0;		// Maximum length of timesteps after the next
	tuning[0][2] = 0.1;			// Minimum length of all timesteps
	tuning[0][3] = 0.15;		// Minimum choppable timestep
	tuning[0][4] = 3.0;			// Maximum timestep increase factor
	tuning[0][5] = 0.3;			// Minimum timestep cutback factor
	tuning[0][6] = 0.1;			// Factor by which timestep is cut after convergence failure
	tuning[0][7] = 1.25;		// Maximum increase factor after a convergence failure
	// Maximum throughput ratio
	tuning[0][8] = method == "IMPES" ? 0.2 : 1E20;			
	tuning[0][9] = -1;			// Maximum length of the next timestep following a well modification : no limit

	// Time truncation and convergence controls
	tuning[1].resize(13);
	// Target time truncation error
	tuning[1][0] = method == "IMPES" ? 1.0 : 0.1;			
	// Target non-linear convergence error
	tuning[1][1] = method == "IMPES" ? 0.5 : 1E-3;
	tuning[1][2] = 1E-7;		// Target material balance error
	// Target linear convergence error
	tuning[1][3] = method == "IMPES" ? 1E-5 : 1E-4;
	tuning[1][4] = 10.0;		// Maximum time truncation error
	// Maximum non-linear convergence error
	tuning[1][5] = method == "IMPES" ? 0.75 : 0.01;
	tuning[1][6] = 1E-6;		// Maximum material balance error
	// Maximum linear convergence error
	tuning[1][7] = method == "IMPES" ? 1E-4 : 1E-3;
	tuning[1][8] = 1E-3;		// Maximum well flow rate convergence error
	tuning[1][9] = 0.025;		// Target Fluid-in-place error for LGR runs
	tuning[1][10] = -1;			// Target surfactant change (Surfactant Model only)
	tuning[1][11] = 0.01;		// Threshold for damping in ion exchange calc. (Multi-Comp. Brine Model only)
	tuning[1][12] = 1;			// Weighting factor for active tracer updates when called from Newton Loop

	// Control of Newton and linear iterations
	tuning[2].resize(10);
	// Maximum number of Newton iterations in a timestep
	tuning[2][0] = method == "IMPES" ? 4 : 12;
	tuning[2][1] = 1;			// Minimum number of Newton iterations in a timestep
	tuning[2][2] = 25;			// Maximum number of linear iterations in a Newton iteration
	tuning[2][3] = 1;			// Minimum number of linear iterations in a Newton iteration
	tuning[2][4] = 8;			// Maximum number of iterations within well flow calculation
	tuning[2][5] = 8;			// Maximum number of iterations for BHP in THP controlled wells
	tuning[2][6] = 1E6;			// Maximum pressure change at last Newton iteration
	tuning[2][7] = 1E6;			// Maximum saturation change at last Newton iteration
	// Target maximum pressure change in a timestep
	tuning[2][8] = method == "IMPES" ? 100 : 1E6;
	tuning[2][9] = -1;			// Maximum tolerable pressure change in a timestep
}

void ParamControl::InputMETHOD(ifstream& ifs)
{
	vector<string>		vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	if (vbuf[0] == "FIM") {
		method = "FIM";
		linearSolve = "./brs.dat";
	}
		
	if (vbuf.size() > 1)
		linearSolve = vbuf[1];

#ifdef _DEBUG
	cout << "METHOD" << endl;
	cout << method << "  " << linearSolve << endl;
#endif // _DEBUG

}

void ParamControl::InputTUNING(ifstream& ifs)
{
	assert(criticalTime.size() == 1);

	TUNING tmp(tuning);
	USI d = criticalTime.size() - 1;

	USI row = 0;
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		/*if (vbuf[0] == "/")
			break;*/
		
		DealDefault(vbuf);
		OCP_INT len = vbuf.size();
		
		for (OCP_INT i = 0; i < len - 1; i++) {
			tmp[row][i] = stod(vbuf[i]);
		}
		if (vbuf[len - 1] != "/") {
			tmp[row][len - 1] = stod(vbuf[len - 1]);
		}
		else {
			row++;
		}
		if (row == 3)
			break;
	}
	tuning_T.push_back(TuningPair(d, tmp));
	DisplayTuning();
	cout << "TUNING" << endl;
}

void ParamControl::DisplayTuning() const
{
	for (auto v : tuning_T) {
		cout << v.d << endl;
		for (auto v1 : v.Tuning) {
			for (auto v2 : v1) {
				cout << v2 << "   ";
			}
			cout << "/ " << endl;
		}
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/