#include "OpenCAEPoro.hpp"

void OpenCAEPoro::InputParam(ParamRead& param)
{
	reservoir.InputParam(param);
	control.InputParam(param.param_Control);
	output.InputParam(param.param_Output);
}

void OpenCAEPoro::SetupReservoir(ParamRead& param)
{
	InputParam(param);
	reservoir.Setup();
	output.Setup(reservoir, control.Dir);
	SetupSolver();
}

void OpenCAEPoro::SetupSolver()
{
	if (control.Method == IMPES) {
		impes.SetupParam(control.Dir, control.solveFile);
		impes.AllocateMat(reservoir);
		cout << "IMPES Method Applys !" << endl;
	}
	else if (control.Method == FIM) {
		cout << "FIM Method Applys !" << endl;
	}
}

void OpenCAEPoro::InitReservoir()
{
	reservoir.Init();
}

void OpenCAEPoro::RunSimulation()
{
	GetWallTime Timer;
	Timer.Start();


	switch (control.Method)
	{
	case IMPES:
		impes.run(reservoir, control, output);
	default:
		break;
	}


	control.TotalTime = Timer.Stop() / 1000;
}

void OpenCAEPoro::OutputResults() {

	cout << endl;
	cout << "Final time:          " << control.Current_time << " Days" << endl;
	cout << "Total linear steps:  " << control.LS_iter_total << endl;
	cout << "Linear solve time:   " << control.LS_time << "s" << endl;
	cout << "Total time steps:    " << control.NR_iter_total << endl;
	cout << "Simulation time:     " << control.TotalTime << "s" << endl;
	output.PrintInfo();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/