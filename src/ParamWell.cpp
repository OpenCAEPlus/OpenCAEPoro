#include "ParamWell.hpp"
#include "OpenCAEPoro_consts.hpp"

WellOptParam::WellOptParam(string intype, vector<string>& vbuf)
{
	type = intype;
	if (type == "INJ") {
		fluidType = vbuf[1];
		state = vbuf[2];
		optMode = vbuf[3];
		maxRate = stod(vbuf[4]);
		maxBHP = stod(vbuf[5]);
	}
	else if (type == "PROD"){
		state = vbuf[1];
		optMode = vbuf[2];
		maxRate = stod(vbuf[3]);
		minBHP = stod(vbuf[4]);
	}
	else {
		Paramcheck("WRONG Well Type");
		exit(0);
	}
}

WellParam::WellParam(vector<string>& info)
{
	name = info[0];
	if (info[1] != "DEFAULT")
		group = info[1];
	I = stoi(info[2]);
	J = stoi(info[3]);
	if (info[4] != "DEFAULT")
		depth = stod(info[4]);
}

void ParamWell::InputWELSPECS(ifstream& ifs)
{
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		well.push_back(WellParam(vbuf));
	}
	cout << "WELSPECS" << endl;
}

void ParamWell::InputCOMPDAT(ifstream& ifs)
{
	USI num = well.size();
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		string src = vbuf[0];
		bool match = (src.find("*") != string::npos);
		if (match) {
			src.pop_back();
		}
		bool tmp = false;

		for (USI w = 0; w < num; w++) {
			if (match)
				tmp = (well[w].name.find(src) != string::npos);
			else
				tmp = (well[w].name == src);

			if (tmp) {
				if (vbuf[1] == "DEFAULT" || vbuf[2] == "DEFAULT") {
					well[w].I_perf = well[w].I;
					well[w].J_perf = well[w].J;
				}
				else {
					well[w].I_perf = stoi(vbuf[1]);
					well[w].J_perf = stoi(vbuf[2]);
				}
				well[w].K1 = stoi(vbuf[3]);
				well[w].K2 = stoi(vbuf[4]);
				if (vbuf[5] != "DEFAULT")
					well[w].WI = stod(vbuf[5]);
				if (vbuf[6] != "DEFAULT")
					well[w].diameter = stod(vbuf[6]);
				if (vbuf[7] != "DEFAULT")
					well[w].kh = stod(vbuf[7]);
				if (vbuf[8] != "DEFAULT")
					well[w].skinFactor = stod(vbuf[8]);
			}
		}
	}
	cout << "COMPDAT" << endl;
}

void ParamWell::InputWCONINJE(ifstream& ifs)
{
	assert(criticalTime.size() > 0);

	USI d = criticalTime.size() - 1;
	USI num = well.size();
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		string src = vbuf[0];
		bool match = (src.find("*") != string::npos);
		if (match) {
			src.pop_back();
		}

		if (match) {
			for (USI w = 0; w < num; w++) {
				if (well[w].name.find(src) != string::npos) {
					well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
				}
			}
		}
		else {
			for (USI w = 0; w < num; w++) {
				if (well[w].name == src) {
					well[w].optParam.push_back(WellOptPair(d, "INJ", vbuf));
				}
			}
		}
		
	}
	cout << "WCONINJE" << endl;
}

void ParamWell::InputWCONPROD(ifstream& ifs)
{
	assert(criticalTime.size() > 0);

	USI d = criticalTime.size() - 1;
	USI num = well.size();
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		string src = vbuf[0];
		bool match = (src.find("*") != string::npos);
		if (match) {
			src.pop_back();
		}

		if (match) {
			for (USI w = 0; w < num; w++)
				if (well[w].name.find(src) != string::npos)
					well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
		}
		else {
			for (USI w = 0; w < num; w++)
				if (well[w].name == src)
					well[w].optParam.push_back(WellOptPair(d, "PROD", vbuf));
		}

	}
	cout << "WCONPROD" << endl;
}

void ParamWell::InputTSTEP(ifstream& ifs)
{
	assert(criticalTime.size() > 0);

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		OCP_INT len = vbuf.size();
		for (OCP_INT i = 0; i < len - 1; i++) {
			OCP_DBL t = criticalTime.back() + stod(vbuf[i]);
			criticalTime.push_back(t);
		}
		if (vbuf.back() != "/") {
			OCP_DBL t = criticalTime.back() + stod(vbuf.back());
			criticalTime.push_back(t);
		}
	}
	cout << "TSTEP" << endl;
}

void ParamWell::InputWELTARG(ifstream& ifs)
{
	assert(criticalTime.size() > 0);

	USI d = criticalTime.size() - 1;
	USI num = well.size();
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		string src = vbuf[0];
		bool match = (src.find("*") != string::npos);
		if (match) {
			src.pop_back();
		}
		bool tmp = false;

		for (USI w = 0; w < num; w++) {
			if (match)
				tmp = (well[w].name.find(src) != string::npos);
			else
				tmp = (well[w].name == src);

			if (tmp) {

				WellOptPair tar = well[w].optParam.back();
				tar.d = d;
				tar.opt.optMode = vbuf[1];
				OCP_DBL val = stod(vbuf[2]);
				if (vbuf[1] == "BHP") {
					if (tar.opt.type == "INJ")
						tar.opt.maxBHP = val;
					else
						tar.opt.minBHP = val;
				}
				else {
					tar.opt.maxRate = val;
				}
				well[w].optParam.push_back(tar);
			}
		}
	}
	cout << "WELTARG" << endl;
}

// check
void ParamWell::CheckParam() const
{
	CheckPerf();
}
void ParamWell::CheckPerf() const
{
	USI wellnum = well.size();
	for (USI w = 0; w < wellnum; w++) {
		if ((well[w].I != well[w].I_perf) || (well[w].J != well[w].J_perf)) {
			Paramcheck("This situation have not been supported!");
			exit(0);
		}
	}
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/