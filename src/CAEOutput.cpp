#include "CAEOutput.hpp"

void Summary::inputParam(OutputSummary& summary_param)
{
	FPR = summary_param.FPR;
	FOPR = summary_param.FOPR;
	FOPT = summary_param.FOPT;
	FGPR = summary_param.FGPR;
	FGPt = summary_param.FGPt;
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

	WBHP = summary_param.WBHP;

	BPR = summary_param.BPR;

	cout << "Summary::inputParam" << endl;
}

void Summary::setup(Reservoir& reservoir)
{
	Sumdata.push_back(SumPair("TIME", "  ", "DAY"));
	Sumdata.push_back(SumPair("NRiter", "  ", "  "));
	Sumdata.push_back(SumPair("LSiter", "  ", "  "));
	if (FPR)
		Sumdata.push_back(SumPair("FPR", "  ", "PSIA"));
	if (FOPR)
		Sumdata.push_back(SumPair("FOPR", "  ", "STB/DAY"));
	if (FOPT)
		Sumdata.push_back(SumPair("FOPT", "  ", "STB"));
	if (FGPR)
		Sumdata.push_back(SumPair("FGPR", "  ", "MSCF/DAY"));
	if (FGPt)
		Sumdata.push_back(SumPair("FGPT", "  ", "MSCF"));
	if (FWPR)
		Sumdata.push_back(SumPair("FWPR", "  ", "STB/DAY"));
	if (FWPT)
		Sumdata.push_back(SumPair("FWPT", "  ", "STB"));
	if (FGIR)
		Sumdata.push_back(SumPair("FGIR", "  ", "MSCF/DAY"));
	if (FGIT)
		Sumdata.push_back(SumPair("FGIT", "  ", "MSCF"));
	if (FWIR)
		Sumdata.push_back(SumPair("FWIR", "  ", "STB/DAY"));
	if (FWIT)
		Sumdata.push_back(SumPair("FWIT", "  ", "STB"));

	if (WOPR.activity) {
		if (WOPR.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WOPR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
			}
		}
		else {
			int num = WOPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WOPR", WOPR.obj[w], "STB/DAY"));
			}
		}
	}
	
	if (WOPT.activity) {
		if (WOPT.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WOPT", reservoir.wellgroup.getWellName(w), "STB"));
			}
		}
		else {
			int num = WOPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WOPT", WOPT.obj[w], "STB"));
			}
		}
	}

	if (WGPR.activity) {
		if (WGPR.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WGPR", reservoir.wellgroup.getWellName(w), "MSCF/DAY"));
			}
		}
		else {
			int num = WGPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGPR", WGPR.obj[w], "MSCF/DAY"));
			}
		}
	}

	if (WGPT.activity) {
		if (WGPT.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WGPT", reservoir.wellgroup.getWellName(w), "MSCF"));
			}
		}
		else {
			int num = WGPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGPT", WGPT.obj[w], "MSCF"));
			}
		}
	}

	if (WWPR.activity) {
		if (WWPR.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WWPR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
			}
		}
		else {
			int num = WWPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWPR", WWPR.obj[w], "STB/DAY"));
			}
		}
	}

	if (WWPT.activity) {
		if (WWPT.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WWPT", reservoir.wellgroup.getWellName(w), "STB"));
			}
		}
		else {
			int num = WWPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWPT", WWPT.obj[w], "STB"));
			}
		}
	}

	if (WGIR.activity) {
		if (WGIR.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WGIR", reservoir.wellgroup.getWellName(w), "MSCF/DAY"));
			}
		}
		else {
			int num = WGIR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGIR", WGIR.obj[w], "MSCF/DAY"));
			}
		}
	}

	if (WGIT.activity) {
		if (WGIT.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WGIT", reservoir.wellgroup.getWellName(w), "MSCF"));
			}
		}
		else {
			int num = WGIT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGIT", WGIT.obj[w], "MSCF"));
			}
		}
	}

	if (WWIR.activity) {
		if (WWIR.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WWIR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
			}
		}
		else {
			int num = WWIR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWIR", WWIR.obj[w], "STB/DAY"));
			}
		}
	}

	if (WWIT.activity) {
		if (WWIT.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WWIT", reservoir.wellgroup.getWellName(w), "STB"));
			}
		}
		else {
			int num = WWIT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWIT", WWIT.obj[w], "STB"));
			}
		}
	}

	if (WBHP.activity) {
		if (WBHP.obj[0] == "All") {
			for (int w = 0; w < reservoir.wellgroup.getWellNum(); w++) {
				Sumdata.push_back(SumPair("WBHP", reservoir.wellgroup.getWellName(w), "PSIA"));
			}
		}
		else {
			int num = WBHP.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WBHP", WBHP.obj[w], "PSIA"));
			}
		}
	}

	if (BPR.activity) {
		int n = BPR.obj.size();
		for (int i = 0; i < n; i++) {
			string temp = "( " + to_string(BPR.obj[i].I) + ", "
				+ to_string(BPR.obj[i].J) + ", "
				+ to_string(BPR.obj[i].K) + " )";
			Sumdata.push_back(SumPair("BPR", temp, "PSIA"));
		}
	}

	cout << "Summary::setup" << endl;

}

void CAEOutput::inputParam(ParamOutput& Output_param)
{
	summary.inputParam(Output_param.Summary);
}

void CAEOutput::setup(Reservoir& reservoir)
{
	summary.setup(reservoir);
}
