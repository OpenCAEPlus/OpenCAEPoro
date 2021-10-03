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
	output.setup(reservoir, control.Dir);
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
	solver.initSolver(control.Dir, control.SolveFile);
}

void OpenCAEPoro::SolveP(OCP_DBL dt)
{
	reservoir.assembleMat(solver, dt);

#ifdef _DEBUG
	solver.checkVal();
#endif // _DEBUG

#ifdef __SOLVER_FASP__
	
	solver.assemble_Fasp();
	

	GetWallTime Timer;
	Timer.Start();
	int status = solver.faspsolve();
	control.LS_time += Timer.Stop() / 1000;

#ifdef _DEBUG
	// solver.showMat_CSR("testA.dat", "testb.dat");
	// solver.showSolution("testx.dat");
#endif // _DEBUG

	
	solver.free_Fasp();

	control.LS_iter = status;
	control.LS_iter_total += status;
	
#endif // __SOLVER_FASP__

	reservoir.getSol_IMPES(solver.getSol());
	solver.clearData();
}

void OpenCAEPoro::run()
{
	GetWallTime Timer;
	Timer.Start();

	OCP_DBL  numdates = control.CriticalTime.size();
	for (int d = 0; d < numdates - 1; d++) {
		reservoir.wellgroup.applyControl(d);
		control.ApplyControl(d);
		control.initTime(d);
		while (control.CriticalTime[d+1] - control.Current_time > TINY) {
			
			runIMPES(control.Current_dt);

		}

	}

	OCP_DBL endtime = Timer.Stop();

	cout << endl;
	cout << "Final time:          " << control.Current_time << " Days" << endl;
    cout << "Total linear steps:  " << control.LS_iter_total << endl;
	cout << "Linear solve time:   " << control.LS_time << "s" << endl;
    cout << "Total time steps:    " << control.NR_iter_total << endl;
	cout << "Simulation time:     " << endtime / 1000 << "s" << endl;
	output.printInfo();
}

void OpenCAEPoro::runIMPES(OCP_DBL& dt)
{
	OCP_DBL ve = 0.01;
	OCP_DBL cfl = 1;
	int	   flagCheck = 0;

	reservoir.wellgroup.prepareWell(reservoir.bulk);

	cfl = reservoir.calCFL(dt);
	if (cfl > 1)
		dt /= (cfl + 1);

	while (true)
	{
		SolveP(dt);

		// first check : Pressure check
		flagCheck = reservoir.checkP();
		if (flagCheck == 1) {
			dt /= 2;
			continue;
		}
		else if (flagCheck == 2) {
			continue;
		}

		reservoir.conn.calFlux(reservoir.bulk);
		reservoir.wellgroup.calFlux(reservoir.bulk);

		// second check : cfl check
		cfl = reservoir.calCFL(dt);
		if (cfl > 1) {
			dt /= 2;
			reservoir.resetVal01();
			continue;
		}

		reservoir.conn.massConserve(reservoir.bulk, dt);
		reservoir.wellgroup.massConserve(reservoir.bulk, dt);

		// third check: Ni check
		if (!reservoir.checkNi()) {
			dt /= 2;
			reservoir.resetVal01();
			continue;
		}

		reservoir.bulk.flash_Ni();
		reservoir.bulk.calVporo();

		if (!reservoir.checkVe(ve)) {
			dt /= 2;
			reservoir.resetVal02();
			continue;
		}

		reservoir.bulk.calKrPc();
		reservoir.conn.calFlux(reservoir.bulk);

		break;
	}
	


	reservoir.wellgroup.calIPRT(reservoir.bulk, dt);
	control.Tstep += 1;
	control.NR_iter = 1;
	control.NR_iter_total += 1;

	reservoir.bulk.calMaxChange();
	control.setNextTstep(reservoir);
	reservoir.bulk.setLastStep();
	reservoir.wellgroup.setLastStep();

	output.setVal(reservoir, control);

#ifdef _DEBUG
	cout << fixed << setprecision(3) << control.Current_time << "Days \n";
#endif // _DEBUG

}
