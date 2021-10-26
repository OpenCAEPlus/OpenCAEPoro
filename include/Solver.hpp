/*! \file    Solver.hpp
 *  \brief   Solver class declaration
 *  \author  Shizhe Li
 *  \date    Oct/21/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "FluidSolver.hpp"
#include "LinearSolver.hpp"


#ifndef __SOLVER_HEADER__
#define __SOLVER_HEADER__

class Solver
{
public:

	/// Simulation begins
	void RunSimulation(Reservoir& rs, OCP_Control& ctrl, OCP_Output& output);
	/// Run one time step
	void GoOneStep(Reservoir& rs, OCP_Control& ctrl);

	/// Before solve: prepare for assembling matrix.
	void Prepare(Reservoir& rs, OCP_DBL& dt);
	/// Assemble and Solve: assemble mat from different parts and solve.
	void AssembleSolve(Reservoir& rs, OCP_Control& ctrl, const OCP_DBL& dt);
	/// Update properties after solving.
	bool UpdateProperty(Reservoir& rs, OCP_DBL& dt);
	/// Finish current timestep.
	void FinishStep(Reservoir& rs, OCP_Control& ctrl);

	/// Setup linear solver params.
	void SetupParamLS(const string& dir, const string& file);
	/// Allocate mat mem
	void AllocateMat(const Reservoir& rs);
	/// Initialize the reservoir
	void InitReservoir(Reservoir& rs) const;


private:

	FluidSolver     FSolver;
	LinearSolver	LSolver;

};




#endif /* end if __SOLVER_HEADER__ */
