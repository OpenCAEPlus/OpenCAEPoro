#include "CAEOutput.hpp"

void Summary::inputParam(OutputSummary& summary_param)
{
	FPR = summary_param.FPR;
	FOPR = summary_param.FOPR;
	FOPT = summary_param.FOPT;
	FGPR = summary_param.FGPR;
	FGPT = summary_param.FGPT;
	FWPR = summary_param.FWPR;
	FWPT = summary_param.FWPT;
	FGIR = summary_param.FGIR;
	FGIT = summary_param.FGIT;
	FWIR = summary_param.FWIR;
	FWIT = summary_param.FWIT;


	WOPR = summary_param.WOPR;
	WOPT = summary_param.WOPT;
	WGPR = summary_param.WGPR;
	WGPT = summary_param.WGPT;
	WWPR = summary_param.WWPR;
	WWPT = summary_param.WWPT;
	WGIR = summary_param.WGIR;
	WGIT = summary_param.WGIT;
	WWIR = summary_param.WWIR;
	WWIT = summary_param.WWIT;

	BPR = summary_param.BPR;

	cout << "Summary::inputParam" << endl;
}


void CAEOutput::inputParam(ParamOutput& Output_param)
{
	summary.inputParam(Output_param.Summary);
}

void CAEOutput::setup()
{

}
