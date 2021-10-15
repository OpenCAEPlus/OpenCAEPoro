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

				USI k1 = stoi(vbuf[3]);
				USI k2 = stoi(vbuf[4]);

				for (USI k = k1; k <= k2; k++) {
					if (vbuf[1] == "DEFAULT" || vbuf[2] == "DEFAULT") {
						well[w].I_perf.push_back(well[w].I);
						well[w].J_perf.push_back(well[w].J);
					}
					else {
						well[w].I_perf.push_back(stoi(vbuf[1]));
						well[w].J_perf.push_back(stoi(vbuf[2]));
					}
					well[w].K_perf.push_back(k);

					if (vbuf[5] != "DEFAULT")
						well[w].WI.push_back(stod(vbuf[5]));
					else
						well[w].WI.push_back(-1.0);

					if (vbuf[6] != "DEFAULT")
						well[w].diameter.push_back(stod(vbuf[6]));
					else
						well[w].diameter.push_back(1.0);

					if (vbuf[7] != "DEFAULT")
						well[w].kh.push_back(stod(vbuf[7]));
					else
						well[w].kh.push_back(-1.0);

					if (vbuf[8] != "DEFAULT")
						well[w].skinFactor.push_back(stod(vbuf[8]));
					else
						well[w].skinFactor.push_back(0.0);

					if (vbuf[9] != "DEFAULT")
						well[w].direction.push_back(vbuf[9]);
					else
						well[w].direction.push_back("z");
				}	
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
	USI perfnum;
	for (USI w = 0; w < wellnum; w++) {
		perfnum = well[w].I_perf.size();
		if (well[w].J_perf.size() != perfnum) {
			Paramcheck("Wrong Perforations J_perf!");
			exit(0);
		}
		if (well[w].K_perf.size() != perfnum) {
			Paramcheck("Wrong Perforations K_perf!");
			exit(0);
		}
		if (well[w].diameter.size() != perfnum) {
			Paramcheck("Wrong Perforations diameter!");
			exit(0);
		}
		if (well[w].WI.size() != perfnum) {
			Paramcheck("Wrong Perforations WI!");
			exit(0);
		}
		if (well[w].kh.size() != perfnum) {
			Paramcheck("Wrong Perforations kh!");
			exit(0);
		}
		if (well[w].skinFactor.size() != perfnum) {
			Paramcheck("Wrong Perforations skinFactor!");
			exit(0);
		}
		if (well[w].direction.size() != perfnum) {
			Paramcheck("Wrong Perforations direction!");
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