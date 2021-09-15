#include "ParamWell.hpp"

WellOptParam::WellOptParam(string type, vector<string>& vbuf)
{
	Type = type;
	if (Type == "INJ") {
		FluidType = vbuf[1];
		State = vbuf[2];
		OptMode = vbuf[3];
		MaxRate = atof(vbuf[4].c_str());
		MaxBHP = atof(vbuf[5].c_str());
	}
	else {
		State = vbuf[1];
		OptMode = vbuf[2];
		MaxRate = atof(vbuf[3].c_str());
		MinBHP = atof(vbuf[4].c_str());
	}	
}

WellParam::WellParam(vector<string>& info)
{
	Name = info[0];
	if (info[1] != "DEFAULT")
		Group = info[1];
	I = atoi(info[2].c_str());
	J = atoi(info[3].c_str());
	if (info[4] != "DEFAULT")
		Dref = atof(info[4].c_str());
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
	cout << "hello" << endl;
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
					well[w].I_perf = atoi(vbuf[1].c_str());
					well[w].J_perf = atoi(vbuf[2].c_str());
				}
				well[w].K1 = atoi(vbuf[3].c_str());
				well[w].K2 = atoi(vbuf[4].c_str());
				if (vbuf[5] != "DEFAULT")
					well[w].Trans = atof(vbuf[5].c_str());
				if (vbuf[6] != "DEFAULT")
					well[w].Diameter = atof(vbuf[6].c_str());
				if (vbuf[7] != "DEFAULT")
					well[w].Kh = atof(vbuf[7].c_str());
				if (vbuf[8] != "DEFAULT")
					well[w].SkinFactor = atof(vbuf[8].c_str());
			}
		}
	}
	cout << "hello" << endl;
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
	cout << "hello" << endl;
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
	cout << "hello" << endl;
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
			double t = CriticalTime.back() + atof(vbuf[i].c_str());
			CriticalTime.push_back(t);
		}
		if (vbuf.back() != "/") {
			double t = CriticalTime.back() + atof(vbuf.back().c_str());
			CriticalTime.push_back(t);
		}
	}
	cout << "hello" << endl;
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
				double val = atof(vbuf[2].c_str());
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
	cout << "hello" << endl;
}
