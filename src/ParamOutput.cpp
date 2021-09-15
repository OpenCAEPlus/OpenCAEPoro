#include "ParamOutput.hpp"

void ParamOutput::inputSUMMARY(ifstream& ifs) 
{
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;
		string keyword = vbuf[0];

		switch (Map_str2int(&keyword[0], keyword.size()))
		{
		case Map_str2int("FPR", 3):
			FPR = true;
			break;

		// Field
		case Map_str2int("FOPR", 4):
			FOPR = true;
			break;

		case Map_str2int("FOPT", 4):
			FOPT = true;
			break;

		case Map_str2int("FGPR", 4):
			FGPR = true;
			break;

		case Map_str2int("FGPT", 4):
			FGPT = true;
			break;

		case Map_str2int("FWPR", 4):
			FWPR = true;
			break;

		case Map_str2int("FWPT", 4):
			FWPT = true;
			break;

		case Map_str2int("FGIR", 4):
			FGIR = true;
			break;

		case Map_str2int("FGIT", 4):
			FGIT = true;
			break;

		case Map_str2int("FWIR", 4):
			FWIR = true;
			break;

		case Map_str2int("FWIT", 4):
			FWIT = true;
			break;

		// Well
		case Map_str2int("WOPR", 4):
			inputType_A(ifs, WOPR);
			break;

		case Map_str2int("WOPT", 4):
			inputType_A(ifs, WOPT);
			break;

		case Map_str2int("WGPR", 4):
			inputType_A(ifs, WGPR);
			break;

		case Map_str2int("WGPT", 4):
			inputType_A(ifs, WGPT);
			break;

		case Map_str2int("WWPR", 4):
			inputType_A(ifs, WWPR);
			break;

		case Map_str2int("WWPT", 4):
			inputType_A(ifs, WWPT);
			break;

		case Map_str2int("WGIR", 4):
			inputType_A(ifs, WGIR);
			break;

		case Map_str2int("WGIT", 4):
			inputType_A(ifs, WGIT);
			break;

		case Map_str2int("WWIR", 4):
			inputType_A(ifs, WWIR);
			break;

		case Map_str2int("WWIT", 4):
			inputType_A(ifs, WWIT);
			break;

		case Map_str2int("WBHP", 4):
			inputType_A(ifs, WBHP);
			break;

		case Map_str2int("BPR", 3):
			inputType_B(ifs, BPR);
			break;
		}
	}
}


void ParamOutput::inputType_A(ifstream& ifs, Type_A& obj)
{
	obj.status = true;
	vector<string>		vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/") {
		obj.obj.push_back("All");
	}
	else {
		int len = vbuf.size();
		for (int i = 0; i < len - 1; i++) {
			obj.obj.push_back(vbuf[i]);
		}
		if (vbuf.back() != "/")
			obj.obj.push_back(vbuf.back());

		while (ReadLine(ifs, vbuf)) {
			if (vbuf[0] == "/")
				break;

			int len = vbuf.size();
			for (int i = 0; i < len - 1; i++) {
				obj.obj.push_back(vbuf[i]);
			}
			if (vbuf.back() != "/")
				obj.obj.push_back(vbuf.back());
		}
	}
	cout << "hello" << endl;
}

void ParamOutput::inputType_B(ifstream& ifs, Type_B& obj)
{

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf)) {
		if (vbuf[0] == "/")
			break;

		obj.status = true;
		DealDefault(vbuf);
		int i = atoi(vbuf[0].c_str());
		int j = atoi(vbuf[1].c_str());
		int k = atoi(vbuf[2].c_str());

		obj.obj.push_back(COOIJK(i, j, k));
	}
	cout << "hello" << endl;
}
