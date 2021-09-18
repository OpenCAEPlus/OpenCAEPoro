#include "OpenCAEPoro.hpp"

void OpenCAEPoro::inputParam(ParamRead& param)
{
	reservoir.inputParam(param);
}

void OpenCAEPoro::setup()
{
	reservoir.setup();
}
