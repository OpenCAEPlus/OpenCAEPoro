#include "OpenCAEPoro.hpp"

void OpenCAEPoro::inputParam(ParamRead& param)
{
	reservoir.inputParam(param);
	control.inputParam(param.Control_param);
	output.inputParam(param.Output_param);
}

void OpenCAEPoro::setup(ParamRead& param)
{
	inputParam(param);
	reservoir.setup();
	output.setup(reservoir, control.Dir);
	setupSolver();
}

void OpenCAEPoro::setupSolver()
{
	if (control.Method == IMPES) {
		Impes.setupParam(control.Dir, control.SolveFile);
		Impes.allocateMat(reservoir);
		cout << "IMPES Method Applys !" << endl;
	}
	else if (control.Method == FIM) {
		cout << "FIM Method Applys !" << endl;
	}
}

void OpenCAEPoro::init()
{
	reservoir.init();
}

void OpenCAEPoro::run()
{
	GetWallTime Timer;
	Timer.Start();


	switch (control.Method)
	{
	case IMPES:
		Impes.run(reservoir, control, output);
	default:
		break;
	}


	control.TotalTime = Timer.Stop() / 1000;
}

void OpenCAEPoro::out() {

	cout << endl;
	cout << "Final time:          " << control.Current_time << " Days" << endl;
	cout << "Total linear steps:  " << control.LS_iter_total << endl;
	cout << "Linear solve time:   " << control.LS_time << "s" << endl;
	cout << "Total time steps:    " << control.NR_iter_total << endl;
	cout << "Simulation time:     " << control.TotalTime << "s" << endl;
	output.printInfo();
}
