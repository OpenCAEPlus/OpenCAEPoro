#include "ParamWell.hpp"
#include "OpenCAEPoro_consts.hpp"

WellOptParam::WellOptParam(string type, vector<string>& vbuf)
{
	Type = type;
	if (Type == "INJ") {
		FluidType = vbuf[1];
		State = vbuf[2];
		OptMode = vbuf[3];
		MaxRate = stod(vbuf[4]);
		MaxBHP = stod(vbuf[5]);
	}
	else if (Type == "PROD"){
		State = vbuf[1];
		OptMode = vbuf[2];
		MaxRate = stod(vbuf[3]);
		MinBHP = stod(vbuf[4]);
	}
	else {
		Paramcheck("WRONG Well Type");
		exit(0);
	}
}

WellParam::WellParam(vector<string>& info)
{
	Name = info[0];
	if (info[1] != "DEFAULT")
		Group = info[1];
	I = stoi(info[2]);
	J = stoi(info[3]);
	if (info[4] != "DEFAULT")
		Dref = stod(info[4]);
}

void ParamWell::inputWELSPECS(ifstream& ifs)
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

void ParamWell::inputCOMPDAT(ifstream& ifs)
{
	int num = well.size();
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

		for (int w = 0; w < num; w++) {
			if (match)
				tmp = (well[w].Name.find(src) != string::npos);
			else
				tmp = (well[w].Name == src);

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
					well[w].Diameter = stod(vbuf[6]);
				if (vbuf[7] != "DEFAULT")
					well[w].Kh = stod(vbuf[7]);
				if (vbuf[8] != "DEFAULT")
					well[w].SkinFactor = stod(vbuf[8]);
			}
		}
	}
	cout << "COMPDAT" << endl;
}

void ParamWell::inputWCONINJE(ifstream& ifs)
{

	int d = CriticalTime.size() - 1;
	int num = well.size();
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
			for (int w = 0; w < num; w++) {
				if (well[w].Name.find(src) != string::npos) {
					well[w].OptParam.push_back(WellOptPair(d, "INJ", vbuf));
				}
			}
		}
		else {
			for (int w = 0; w < num; w++) {
				if (well[w].Name == src) {
					well[w].OptParam.push_back(WellOptPair(d, "INJ", vbuf));
				}
			}
		}
		
	}
	cout << "WCONINJE" << endl;
}

void ParamWell::inputWCONPROD(ifstream& ifs)
{

	int d = CriticalTime.size() - 1;
	int num = well.size();
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
			for (int w = 0; w < num; w++)
				if (well[w].Name.find(src) != string::npos)
					well[w].OptParam.push_back(WellOptPair(d, "PROD", vbuf));
		}
		else {
			for (int w = 0; w < num; w++)
				if (well[w].Name == src)
					well[w].OptParam.push_back(WellOptPair(d, "PROD", vbuf));
		}

	}
	cout << "WCONPROD" << endl;
}

void ParamWell::inputTSTEP(ifstream& ifs)
{
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealDefault(vbuf);
		int len = vbuf.size();
		for (int i = 0; i < len - 1; i++) {
			OCP_DBL t = CriticalTime.back() + stod(vbuf[i]);
			CriticalTime.push_back(t);
		}
		if (vbuf.back() != "/") {
			OCP_DBL t = CriticalTime.back() + stod(vbuf.back());
			CriticalTime.push_back(t);
		}
	}
	cout << "TSTEP" << endl;
}

void ParamWell::inputWELTARG(ifstream& ifs)
{
	int d = CriticalTime.size() - 1;
	int num = well.size();
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

		for (int w = 0; w < num; w++) {
			if (match)
				tmp = (well[w].Name.find(src) != string::npos);
			else
				tmp = (well[w].Name == src);

			if (tmp) {

				WellOptPair tar = well[w].OptParam.back();
				tar.d = d;
				tar.Opt.OptMode = vbuf[1];
				OCP_DBL val = stod(vbuf[2]);
				if (vbuf[1] == "BHP") {
					if (tar.Opt.Type == "INJ")
						tar.Opt.MaxBHP = val;
					else
						tar.Opt.MinBHP = val;
				}
				else {
					tar.Opt.MaxRate = val;
				}
				well[w].OptParam.push_back(tar);
			}
		}
	}
	cout << "WELTARG" << endl;
}

// check
void ParamWell::checkParam()
{
	checkPerf();
}
void ParamWell::checkPerf()
{
	int wellnum = well.size();
	for (int w = 0; w < wellnum; w++) {
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
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/