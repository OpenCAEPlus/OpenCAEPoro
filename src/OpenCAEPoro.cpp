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
	output.setup(reservoir);
}
