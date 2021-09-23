#include "OpenCAEPoro.hpp"

void OpenCAEPoro::inputParam(ParamRead& param)
{
	reservoir.inputParam(param);
	control.inputParam(param.Control_param);
	output.inputParam(param.Output_param);
}

void OpenCAEPoro::setup()
{
	reservoir.setup();
	output.setup(reservoir);
}

void OpenCAEPoro::allocateMat()
{
	if (control.Method == IMPES) {
		solver.allocate(reservoir.bulk.getBulkNum() + reservoir.wellgroup.getWellNum());
		reservoir.conn.allocateMat(solver);
		reservoir.wellgroup.allocateMat(solver);
		solver.allocateColVal();
		cout << "OpenCAEPoro::allocateMat" << endl;
	}
}

void OpenCAEPoro::init()
{
	reservoir.init();
}

void OpenCAEPoro::SolveP(double dt)
{
	reservoir.assembleMat(solver, dt);
	solver.clearData();
}

void OpenCAEPoro::run()
{
	double  numdates = control.CriticalTime.size();
	for (int d = 0; d < numdates - 1; d++) {
		reservoir.wellgroup.applyControl(d);
		control.ApplyControl(d);
		control.initTime(d);
		while (control.Current_time < control.CriticalTime[d+1]) {
			
			runIMPES(control.Current_dt);

		}

	}
}

void OpenCAEPoro::runIMPES(double& dt)
{
	while (true) {
		reservoir.wellgroup.prepareWell(reservoir.bulk);
		reservoir.calCFL(dt);
		SolveP(dt);
		cout << "stop" << endl;
	}
	
}
