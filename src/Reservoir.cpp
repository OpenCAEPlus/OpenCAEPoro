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


double Reservoir::calCFL(double dt)
{
	double cflB = conn.calCFL(bulk, dt);
	double cflW = wellgroup.calCFL(bulk, dt);
	double cfl = max(cflB, cflW);

	return cfl;
}

// assemble mat
void Reservoir::assembleMat(Solver<double>& mysolver, double dt)
{
	conn.initAssembleMat(mysolver);
	conn.assembleMat(mysolver, bulk, dt);
	wellgroup.assemblaMat_WB(mysolver, bulk, dt);
}

void Reservoir::getSol_IMPES(vector<double>& u)
{
	bulk.getSol_IMPES(u);
	wellgroup.getSol_IMPES(u, bulk.getBulkNum());
}

int Reservoir::checkP()
{
	if (!bulk.checkP())
		return 1;
	return wellgroup.checkP(bulk);
}


void Reservoir::resetVal()
{
	bulk.resetVal();
	conn.calFlux(bulk);
}
