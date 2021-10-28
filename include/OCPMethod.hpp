/*! \file    OCPMethod.hpp
 *  \brief   Solution methods in OpenCAEPoro
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPMETHOD_HEADER__
#define __OCPMETHOD_HEADER__

 // OpenCAEPoro header files
#include "Reservoir.hpp"
#include "OCPControl.hpp"
#include "LinearSolver.hpp"

class OCPMethod
{
public:

	// General
	void CalLsColCapacity(LinearSolver& ls, const Reservoir& rs);
	void InitMatSparsity(LinearSolver& ls, const BulkConn& conn);

	// IMPEC
	void SetupIMPEC(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl);
	void AllocateRsIMPEC(Reservoir& rs);
	void SetupLsIMPEC(LinearSolver& ls, const Reservoir& rs, const OCPControl& ctrl);
	void AllocateLsIMPEC(LinearSolver& ls, const Reservoir& rs);
	void AssembleLsIMPEC(LinearSolver& ls, const Reservoir& rs, const OCP_DBL& dt);
	void AssembleLsConnIMPEC(LinearSolver& ls, const BulkConn& conn, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleLsWellIMPEC(LinearSolver& ls, const WellGroup& wellG, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleLsINJWellIMPEC(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleLsPRODWellIMPEC(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt);

	//// FIM
	void SetupFIM(Reservoir& rs, LinearSolver& ls, const OCPControl& ctrl);
	void AllocateRsFIM(Reservoir& rs);
	void SetupLsFIM(LinearSolver& ls, const Reservoir& rs, const OCPControl& ctrl);
	void AllocateLsFIM(LinearSolver& ls, const Reservoir& rs);
	void AssembleMatFIM(LinearSolver& ls, const Reservoir& rs, const OCP_DBL& dt);
	void AssembleMatConnFIM(LinearSolver& ls, const BulkConn& conn, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleMatWellFIM(LinearSolver& ls, const WellGroup& wellG, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleMatINJWellFIM(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt);
	void AssembleMatPRODWellFIM(LinearSolver& ls, const Well& well, const Bulk& bulk, const OCP_DBL& dt);

	virtual ~OCPMethod() = 0;
};




#endif /* end if __OCP_METHOD_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/