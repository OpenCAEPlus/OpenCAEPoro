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

	int wellnum = reservoir.wellgroup.getWellNum();

	if (WOPR.activity) {
		if (WOPR.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WOPR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
				WOPR.index.push_back(w);
			}
		}
		else {
			int num = WOPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WOPR", WOPR.obj[w], "STB/DAY"));
				WOPR.index.push_back(reservoir.wellgroup.getIndex(WOPR.obj[w]));
			}
		}
	}
	
	if (WOPT.activity) {
		if (WOPT.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WOPT", reservoir.wellgroup.getWellName(w), "STB"));
				WOPT.index.push_back(w);
			}
		}
		else {
			int num = WOPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WOPT", WOPT.obj[w], "STB"));
				WOPT.index.push_back(reservoir.wellgroup.getIndex(WOPT.obj[w]));
			}
		}
	}

	if (WGPR.activity) {
		if (WGPR.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WGPR", reservoir.wellgroup.getWellName(w), "MSCF/DAY"));
				WGPR.index.push_back(w);
			}
		}
		else {
			int num = WGPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGPR", WGPR.obj[w], "MSCF/DAY"));
				WGPR.index.push_back(reservoir.wellgroup.getIndex(WGPR.obj[w]));
			}
		}
	}

	if (WGPT.activity) {
		if (WGPT.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WGPT", reservoir.wellgroup.getWellName(w), "MSCF"));
				WGPT.index.push_back(w);
			}
		}
		else {
			int num = WGPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGPT", WGPT.obj[w], "MSCF"));
				WGPT.index.push_back(reservoir.wellgroup.getIndex(WGPT.obj[w]));
			}
		}
	}

	if (WWPR.activity) {
		if (WWPR.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WWPR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
				WWPR.index.push_back(w);
			}
		}
		else {
			int num = WWPR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWPR", WWPR.obj[w], "STB/DAY"));
				WWPR.index.push_back(reservoir.wellgroup.getIndex(WWPR.obj[w]));
			}
		}
	}

	if (WWPT.activity) {
		if (WWPT.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WWPT", reservoir.wellgroup.getWellName(w), "STB"));
				WWPT.index.push_back(w);
			}
		}
		else {
			int num = WWPT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWPT", WWPT.obj[w], "STB"));
				WWPT.index.push_back(reservoir.wellgroup.getIndex(WWPT.obj[w]));
			}
		}
	}

	if (WGIR.activity) {
		if (WGIR.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WGIR", reservoir.wellgroup.getWellName(w), "MSCF/DAY"));
				WGIR.index.push_back(w);
			}
		}
		else {
			int num = WGIR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGIR", WGIR.obj[w], "MSCF/DAY"));
				WGIR.index.push_back(reservoir.wellgroup.getIndex(WGIR.obj[w]));
			}
		}
	}

	if (WGIT.activity) {
		if (WGIT.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WGIT", reservoir.wellgroup.getWellName(w), "MSCF"));
				WGIT.index.push_back(w);
			}
		}
		else {
			int num = WGIT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WGIT", WGIT.obj[w], "MSCF"));
				WGIT.index.push_back(reservoir.wellgroup.getIndex(WGIT.obj[w]));
			}
		}
	}

	if (WWIR.activity) {
		if (WWIR.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WWIR", reservoir.wellgroup.getWellName(w), "STB/DAY"));
				WWIR.index.push_back(w);
			}
		}
		else {
			int num = WWIR.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWIR", WWIR.obj[w], "STB/DAY"));
				WWIR.index.push_back(reservoir.wellgroup.getIndex(WWIR.obj[w]));
			}
		}
	}

	if (WWIT.activity) {
		if (WWIT.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WWIT", reservoir.wellgroup.getWellName(w), "STB"));
				WWIT.index.push_back(w);
			}
		}
		else {
			int num = WWIT.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WWIT", WWIT.obj[w], "STB"));
				WWIT.index.push_back(reservoir.wellgroup.getIndex(WWIT.obj[w]));
			}
		}
	}

	if (WBHP.activity) {
		if (WBHP.obj[0] == "All") {
			for (int w = 0; w < wellnum; w++) {
				Sumdata.push_back(SumPair("WBHP", reservoir.wellgroup.getWellName(w), "PSIA"));
				WBHP.index.push_back(w);
			}
		}
		else {
			int num = WBHP.obj.size();
			for (int w = 0; w < num; w++) {
				Sumdata.push_back(SumPair("WBHP", WBHP.obj[w], "PSIA"));
				WBHP.index.push_back(reservoir.wellgroup.getIndex(WBHP.obj[w]));
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
			int I = BPR.obj[i].I - 1;
			int J = BPR.obj[i].J - 1;
			int K = BPR.obj[i].K - 1;
			BPR.index.push_back(reservoir.grid.getIndex(I, J, K));
		}
	}

	cout << "Summary::setup" << endl;

}

void Summary::setVal(Reservoir& reservoir, CAEControl& ctrl)
{
	int n = 0;

	// TIME
	Sumdata[n++].val.push_back(ctrl.getCurTime());
	// NRiter
	Sumdata[n++].val.push_back(ctrl.getNRiter());
	// LSiter
	Sumdata[n++].val.push_back(ctrl.getLSiter());

	// FPR
	if (FPR)
		Sumdata[n++].val.push_back(reservoir.bulk.calFPR());
	if (FOPR)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFOPR());
	if (FOPT)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFOPT());
	if (FGPR)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFGPR());
	if (FGPt)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFGPT());
	if (FWPR)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFWPR());
	if (FWPT)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFWPT());
	if (FGIR)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFGIR());
	if (FGIT)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFGIT());
	if (FWIR)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFWIR());
	if (FWIT)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getFWIT());

	int len = 0;
	// WOPR
	len = WOPR.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWOPR(WOPR.index[w]));
	// WOPT
	len = WOPT.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWOPT(WOPT.index[w]));
	// WGPR
	len = WGPR.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWGPR(WGPR.index[w]));
	// WGPT
	len = WGPT.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWGPT(WGPT.index[w]));
	// WWPR
	len = WWPR.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWWPR(WWPR.index[w]));
	// WWPT
	len = WWPT.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWWPT(WWPT.index[w]));
	// WGIR
	len = WGIR.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWGIR(WGIR.index[w]));
	// WGIT
	len = WGIT.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWGIT(WGIT.index[w]));
	// WWIR
	len = WWIR.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWWIR(WWIR.index[w]));
	// WWIT
	len = WWIT.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWWIT(WWIT.index[w]));
	// WBHP
	len = WBHP.index.size();
	for (int w = 0; w < len; w++)
		Sumdata[n++].val.push_back(reservoir.wellgroup.getWBHP(WBHP.index[w]));

	// BPR
	len = BPR.index.size();
	for (int i = 0; i < len; i++)
		Sumdata[n++].val.push_back(reservoir.bulk.getP(BPR.index[i]));

	cout << "Summary::setVal" << endl;
}

void CAEOutput::inputParam(ParamOutput& Output_param)
{
	summary.inputParam(Output_param.Summary);
}

void CAEOutput::setup(Reservoir& reservoir)
{
	summary.setup(reservoir);
}

void CAEOutput::setVal(Reservoir& reservoir, CAEControl& ctrl)
{
	summary.setVal(reservoir, ctrl);
}

