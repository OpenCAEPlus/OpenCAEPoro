#include "Reservoir.hpp"

void Reservoir::inputParam(ParamRead& param)
{
	grid.inputParam(param.Rs_param);
	bulk.inputParam(param.Rs_param);
	wellgroup.inputParam(param.Well_param);
}

void Reservoir::setup()
{
	grid.setup();
	bulk.setup(grid);
	conn.setup(grid, bulk);
	wellgroup.setup(grid, bulk);
}


void Reservoir::init()
{
	if (bulk.mixMode() == BLKOIL)
		bulk.initSjPc_blk(50);
	else if (bulk.mixMode() == EoS_PVTW)
		bulk.initSjPc_comp(50);

	bulk.calVporo();
	bulk.flash_Sj();
	bulk.calKrPc();
	bulk.setLastStep();
	conn.calFlux(bulk);
	wellgroup.init(bulk);

}


OCP_DBL Reservoir::calCFL(const OCP_DBL& dt) const
{
	OCP_DBL cflB = conn.calCFL(bulk, dt);
	OCP_DBL cflW = wellgroup.calCFL(bulk, dt);
	OCP_DBL cfl = max(cflB, cflW);

	return cfl;
}

// assemble mat
void Reservoir::assembleMat(Solver<OCP_DBL>& mysolver, const OCP_DBL& dt) const
{
	conn.initAssembleMat(mysolver);
	conn.assembleMat_IMPES(mysolver, bulk, dt);
	wellgroup.assemblaMat_WB_IMPES(mysolver, bulk, dt);
}

void Reservoir::getSol_IMPES(const vector<OCP_DBL>& u)
{
	bulk.getSol_IMPES(u);
	wellgroup.getSol_IMPES(u, bulk.getBulkNum());
}

OCP_INT Reservoir::checkP()
{
	if (!bulk.checkP())
		return 1;

	OCP_INT flag = 0;
	flag = wellgroup.checkP(bulk);
	return flag;
}


void Reservoir::resetVal01()
{
	bulk.resetP();
	bulk.resetPj();
	conn.calFlux(bulk);
}

void Reservoir::resetVal02()
{
	bulk.resetP();
	bulk.resetPj();
	conn.calFlux(bulk);

	bulk.resetNi();
	bulk.flash_Ni();

	bulk.resetVp();
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/