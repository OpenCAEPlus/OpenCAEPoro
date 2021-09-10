#include "Solver.hpp"

void Solver::assembleSolverA()
{

}


void Solver::clearData()
{
	// ColId.clear();
	for (int i = 0; i < Dim; i++) {
		Val[i].clear();
	}
	// DiagPtr.clear();
	DiagVal.assign(Dim, 0);
	b.assign(Dim, 0);
}
