#include "ParamControl.hpp"

void ParamControl::init()
{
	initMethod();
	initTime();
	initTuning();
}

void ParamControl::initMethod()
{
	Method = "IMPES";
	LinearSolve = "./csr.dat";
}

void ParamControl::initTuning()
{
	Tuning.resize(3);

	// Timestepping controls
	Tuning[0].resize(10);
	Tuning[0][0] = 1.0;			// Maximum init step length of next timestep
	Tuning[0][1] = 365.0;		// Maximum length of timesteps after the next
	Tuning[0][2] = 0.1;			// Minimum length of all timesteps
	Tuning[0][3] = 0.15;		// Minimum choppable timestep
	Tuning[0][4] = 3.0;			// Maximum timestep increase factor
	Tuning[0][5] = 0.3;			// Minimum timestep cutback factor
	Tuning[0][6] = 0.1;			// Factor by which timestep is cut after convergence failure
	Tuning[0][7] = 1.25;		// Maximum increase factor after a convergence failure
	// Maximum throughput ratio
	Tuning[0][8] = Method == "IMPES" ? 0.2 : 1E20;			
	Tuning[0][9] = -1;			// Maximum length of the next timestep following a well modification : no limit

	// Time truncation and convergence controls
	Tuning[1].resize(13);
	// Target time truncation error
	Tuning[1][0] = Method == "IMPES" ? 1.0 : 0.1;			
	// Target non-linear convergence error
	Tuning[1][1] = Method == "IMPES" ? 0.5 : 1E-3;
	Tuning[1][2] = 1E-7;		// Target material balance error
	// Target linear convergence error
	Tuning[1][3] = Method == "IMPES" ? 1E-5 : 1E-4;
	Tuning[1][4] = 10.0;		// Maximum time truncation error
	// Maximum non-linear convergence error
	Tuning[1][5] = Method == "IMPES" ? 0.75 : 0.01;
	Tuning[1][6] = 1E-6;		// Maximum material balance error
	// Maximum linear convergence error
	Tuning[1][7] = Method == "IMPES" ? 1E-4 : 1E-3;
	Tuning[1][8] = 1E-3;		// Maximum well flow rate convergence error
	Tuning[1][9] = 0.025;		// Target Fluid-in-place error for LGR runs
	Tuning[1][10] = -1;			// Target surfactant change (Surfactant Model only)
	Tuning[1][11] = 0.01;		// Threshold for damping in ion exchange calc. (Multi-Comp. Brine Model only)
	Tuning[1][12] = 1;			// Weighting factor for active tracer updates when called from Newton Loop

	// Control of Newton and linear iterations
	Tuning[2].resize(10);
	// Maximum number of Newton iterations in a timestep
	Tuning[2][0] = Method == "IMPES" ? 4 : 12;
	Tuning[2][1] = 1;			// Minimum number of Newton iterations in a timestep
	Tuning[2][2] = 25;			// Maximum number of linear iterations in a Newton iteration
	Tuning[2][3] = 1;			// Minimum number of linear iterations in a Newton iteration
	Tuning[2][4] = 8;			// Maximum number of iterations within well flow calculation
	Tuning[2][5] = 8;			// Maximum number of iterations for BHP in THP controlled wells
	Tuning[2][6] = 1E6;			// Maximum pressure change at last Newton iteration
	Tuning[2][7] = 1E6;			// Maximum saturation change at last Newton iteration
	// Target maximum pressure change in a timestep
	Tuning[2][8] = Method == "IMPES" ? 100 : 1E6;
	Tuning[2][9] = -1;			// Maximum tolerable pressure change in a timestep
}

void ParamControl::inputMETHOD(ifstream& ifs)
{
	vector<string>		vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	if (vbuf[0] == "FIM") {
		Method = "FIM";
		LinearSolve = "./brs.dat";
	}
		
	if (vbuf.size() > 1)
		LinearSolve = vbuf[1];

	cout << "METHOD" << endl;
	cout << Method << "  " << LinearSolve << endl;

}

void ParamControl::inputTUNING(ifstream& ifs)
{
	TUNING tmp(Tuning);
	int d = CriticalTime.size() - 1;

	int row = 0;
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		/*if (vbuf[0] == "/")
			break;*/
		
		DealDefault(vbuf);
		int len = vbuf.size();
		
		for (int i = 0; i < len - 1; i++) {
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
	Tuning_T.push_back(TuningPair(d, tmp));
	showTuning();
	cout << "TUNING" << endl;
}

void ParamControl::showTuning()
{
	for (auto v : Tuning_T) {
		cout << v.d << endl;
		for (auto v1 : v.Tuning) {
			for (auto v2 : v1) {
				cout << v2 << "   ";
			}
			cout << "/ " << endl;
		}
	}
}
