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

void OpenCAEPoro::SolveP(double dt)
{
	reservoir.assembleMat(solver, dt);
#ifdef __SOLVER_FASP__
	
	solver.assemble_Fasp();
	// solver.showMat_CSR("testA.dat", "testb.dat");
	int status = solver.faspsolve();
	// solver.showSolution("testx.dat");
	solver.free_Fasp();

	control.LS_iter = status;
	control.LS_iter_total += status;
	
#endif // __SOLVER_FASP__
	reservoir.getSol_IMPES(solver.getSol());
	solver.clearData();
}

void OpenCAEPoro::run()
{
	double  numdates = control.CriticalTime.size();
	for (int d = 0; d < numdates - 1; d++) {
		reservoir.wellgroup.applyControl(d);
		control.ApplyControl(d);
		control.initTime(d);
		while (control.CriticalTime[d+1] - control.Current_time > TINY) {
			
			runIMPES(control.Current_dt);

		}

	}
	output.printInfo();
}

void OpenCAEPoro::runIMPES(double& dt)
{
	double cfl = 1;

	reservoir.wellgroup.prepareWell(reservoir.bulk);

	cfl = reservoir.calCFL(dt);
	if (cfl > 1)
		dt /= (cfl + 1);

	SolveP(dt);
	reservoir.conn.calFlux(reservoir.bulk);
	reservoir.wellgroup.calFlux(reservoir.bulk);
	reservoir.conn.massConserve(reservoir.bulk, dt);
	reservoir.wellgroup.massConserve(reservoir.bulk, dt);

	reservoir.bulk.flash_Ni();
	reservoir.bulk.calKrPc();
	reservoir.bulk.calVporo();
	reservoir.conn.calFlux(reservoir.bulk);


	reservoir.wellgroup.calIPRT(reservoir.bulk, dt);
	control.Tstep += 1;
	control.NR_iter = 1;
	control.NR_iter_total += 1;

	reservoir.bulk.calMaxChange();
	control.setNextTstep(reservoir);
	reservoir.bulk.setLastStep();

	cout << fixed << setprecision(3) << control.Current_time << "Days \n";

	output.setVal(reservoir, control);

}
