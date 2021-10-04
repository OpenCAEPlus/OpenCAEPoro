#include "Method.hpp"


void OCP_IMPES::setupParam(const string& dir, const string& file)
{
	solver.setupParam(dir, file);
}

void OCP_IMPES::allocateMat(const Reservoir& rs)
{
	solver.allocate(rs.bulk.getBulkNum() + rs.wellgroup.getWellNum());
	rs.conn.allocateMat(solver);
	rs.wellgroup.allocateMat(solver);
	solver.allocateColVal();
}

void OCP_IMPES::run(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output)
{

	unsigned int numdates = ctrl.getNumDates();
	for (unsigned int d = 0; d < numdates - 1; d++) {
		rs.wellgroup.applyControl(d);
		ctrl.ApplyControl(d);
		ctrl.initTime(d);
		while (ctrl.CriticalTime[d + 1] - ctrl.Current_time > TINY) {

			goOneStep(rs, ctrl);
			output.setVal(rs, ctrl);

		}
	}
}

void OCP_IMPES::goOneStep(Reservoir& rs, OCP_Control& ctrl)
{
	OCP_DBL ve = 0.01;
	OCP_DBL cfl = 1;
	int	   flagCheck = 0;
	double& dt = ctrl.Current_dt;

	rs.wellgroup.prepareWell(rs.bulk);

	cfl = rs.calCFL(dt);
	if (cfl > 1)
		dt /= (cfl + 1);

	while (true)
	{
		SolveP(rs, ctrl, dt);

		// first check : Pressure check
		flagCheck = rs.checkP();
		if (flagCheck == 1) {
			dt /= 2;
			continue;
		}
		else if (flagCheck == 2) {
			continue;
		}

		rs.conn.calFlux(rs.bulk);
		rs.wellgroup.calFlux(rs.bulk);

		// second check : cfl check
		cfl = rs.calCFL(dt);
		if (cfl > 1) {
			dt /= 2;
			rs.resetVal01();
			continue;
		}

		rs.conn.massConserve(rs.bulk, dt);
		rs.wellgroup.massConserve(rs.bulk, dt);

		// third check: Ni check
		if (!rs.checkNi()) {
			dt /= 2;
			rs.resetVal01();
			continue;
		}

		rs.bulk.flash_Ni();
		rs.bulk.calVporo();

		// fouth check: Volume error check
		if (!rs.checkVe(ve)) {
			dt /= 2;
			rs.resetVal02();
			continue;
		}

		rs.bulk.calKrPc();
		rs.conn.calFlux(rs.bulk);

		break;
	}


	rs.wellgroup.calIPRT(rs.bulk, dt);
	ctrl.Tstep += 1;
	ctrl.NR_iter = 1;
	ctrl.NR_iter_total += 1;

	rs.bulk.calMaxChange();
	ctrl.setNextTstep(rs);
	rs.bulk.setLastStep();
	rs.wellgroup.setLastStep();
}


void OCP_IMPES::SolveP(Reservoir& rs, OCP_Control& ctrl, const OCP_DBL& dt)
{
	rs.assembleMat(solver, dt);

#ifdef _DEBUG
	solver.checkVal();
#endif // _DEBUG

#ifdef __SOLVER_FASP__

	solver.assemble_Fasp();
	GetWallTime Timer;
	Timer.Start();
	int status = solver.faspsolve();
	ctrl.LS_time += Timer.Stop() / 1000;

#ifdef _DEBUG
	// solver.showMat_CSR("testA.dat", "testb.dat");
	// solver.showSolution("testx.dat");
#endif // _DEBUG


	solver.free_Fasp();

	ctrl.LS_iter = status;
	ctrl.LS_iter_total += status;

#endif // __SOLVER_FASP__

	rs.getSol_IMPES(solver.getSol());
	solver.clearData();
}
