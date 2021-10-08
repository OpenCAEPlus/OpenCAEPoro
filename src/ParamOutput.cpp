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
			Summary.FPR = true;
			break;

		// Field
		case Map_str2int("FOPR", 4):
			Summary.FOPR = true;
			break;

		case Map_str2int("FOPT", 4):
			Summary.FOPT = true;
			break;

		case Map_str2int("FGPR", 4):
			Summary.FGPR = true;
			break;

		case Map_str2int("FGPT", 4):
			Summary.FGPt = true;
			break;

		case Map_str2int("FWPR", 4):
			Summary.FWPR = true;
			break;

		case Map_str2int("FWPT", 4):
			Summary.FWPT = true;
			break;

		case Map_str2int("FGIR", 4):
			Summary.FGIR = true;
			break;

		case Map_str2int("FGIT", 4):
			Summary.FGIT = true;
			break;

		case Map_str2int("FWIR", 4):
			Summary.FWIR = true;
			break;

		case Map_str2int("FWIT", 4):
			Summary.FWIT = true;
			break;

		// Well
		case Map_str2int("WOPR", 4):
			inputType_A(ifs, Summary.WOPR);
			break;

		case Map_str2int("WOPT", 4):
			inputType_A(ifs, Summary.WOPT);
			break;

		case Map_str2int("WGPR", 4):
			inputType_A(ifs, Summary.WGPR);
			break;

		case Map_str2int("WGPT", 4):
			inputType_A(ifs, Summary.WGPT);
			break;

		case Map_str2int("WWPR", 4):
			inputType_A(ifs, Summary.WWPR);
			break;

		case Map_str2int("WWPT", 4):
			inputType_A(ifs, Summary.WWPT);
			break;

		case Map_str2int("WGIR", 4):
			inputType_A(ifs, Summary.WGIR);
			break;

		case Map_str2int("WGIT", 4):
			inputType_A(ifs, Summary.WGIT);
			break;

		case Map_str2int("WWIR", 4):
			inputType_A(ifs, Summary.WWIR);
			break;

		case Map_str2int("WWIT", 4):
			inputType_A(ifs, Summary.WWIT);
			break;

		case Map_str2int("WBHP", 4):
			inputType_A(ifs, Summary.WBHP);
			break;

		case Map_str2int("BPR", 3):
			inputType_B(ifs, Summary.BPR);
			break;
		}
	}
	cout << "SUMMARY" << endl;
}


void ParamOutput::inputType_A(ifstream& ifs, Type_A_o& obj)
{
	obj.activity = true;
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
	cout << "Type_A" << endl;
}

void ParamOutput::inputType_B(ifstream& ifs, Type_B_o& obj)
{

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf)) {
		if (vbuf[0] == "/")
			break;

		obj.activity = true;
		DealDefault(vbuf);
		int i = stoi(vbuf[0]);
		int j = stoi(vbuf[1]);
		int k = stoi(vbuf[2]);

		obj.obj.push_back(COOIJK(i, j, k));
	}
	cout << "Type_B" << endl;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/08/2021      Create file                          */
/*----------------------------------------------------------------------------*/