#include "Solver.hpp"

void Solver::RunSimulation(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output)
{
	GetWallTime timer;
	timer.Start();

	USI numdates = ctrl.GetNumDates();
	output.PrintInfoSched(rs, ctrl, timer.Stop());
	for (USI d = 0; d < numdates - 1; d++) {
		rs.ApplyControl(d);
		ctrl.ApplyControl(d);
		ctrl.InitTime(d);
		while (!ctrl.IfCriticalTime(d+1)) {
			GoOneStep(rs, ctrl);
			output.SetVal(rs, ctrl);
		}
		output.PrintInfoSched(rs, ctrl, timer.Stop());
	}

	ctrl.RecordTotalTime(timer.Stop() / 1000);
}

void Solver::GoOneStep(Reservoir& rs, OCP_Control& ctrl)
{
	// cout << setprecision(3) << ctrl.GetCurTime() << "days\n";

	OCP_DBL& dt = ctrl.GetCurDt();
	Prepare(rs, dt);
	
	while (true)
	{
		if (dt < MIN_TIME_STEP) OCP_ABORT("Time stepsize is too small!");

		AssembleSolve(rs, ctrl, dt);
		if (!UpdateProperty(rs, dt)) {
			continue;
		}
		FinishStep(rs, ctrl);
		break;
	}
	
}



void Solver::Prepare(Reservoir& rs, OCP_DBL& dt)
{
	FSolver.Prepare(rs, dt);
}

void Solver::AssembleSolve(Reservoir& rs, OCP_Control& ctrl, const OCP_DBL& dt)
{
	// Assemble Mat.
	FSolver.AssembleMat(rs, dt);

	// Assemble Mat from above.

	// Solve linear system.
	if (true) {
		FSolver.SolveLinearSystem(rs, ctrl);
	}

}

bool Solver::UpdateProperty(Reservoir& rs, OCP_DBL& dt)
{
	bool flag;
	flag = FSolver.UpdateProperty(rs, dt);
	if (!flag) {
		return false;
	}

	return true;
}

void Solver::FinishStep(Reservoir& rs, OCP_Control& ctrl)
{
	FSolver.FinishStep(rs, ctrl);
}


void Solver::SetupParamLS(const string& dir, const string& file)
{
	// LSolver.SetupParam(dir, file);
	FSolver.SetupParamLS(dir, file);
}


void Solver::AllocateMat(const Reservoir& rs)
{
	FSolver.AllocateMat(rs);
}