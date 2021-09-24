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
	conn.calFlux(bulk);

	wellgroup.init(bulk);
}


double Reservoir::calCFL(double dt)
{

}

// assemble mat
void Reservoir::assembleMat(Solver<double>& mysolver, double dt)
{
	conn.initAssembleMat(mysolver);
	conn.assembleMat(mysolver, bulk, dt);
	wellgroup.assemblaMat_WB(mysolver, bulk, dt);
}

void Reservoir::getP_IMPES(vector<double>& u)
{
	bulk.getP_IMPES(u);
	wellgroup.getP_IMPES(u, bulk.getBulkNum());
}
