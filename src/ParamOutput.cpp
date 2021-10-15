/*! \file    ParamOutput.cpp
 *  \brief   ParamOutput class definition
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "ParamOutput.hpp"

void ParamOutput::InputSUMMARY(ifstream& ifs) 
{
	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;
		string keyword = vbuf[0];

		switch (Map_Str2Int(&keyword[0], keyword.size()))
		{
		case Map_Str2Int("FPR", 3):
			summary.FPR = true;
			break;

		// Field
		case Map_Str2Int("FOPR", 4):
			summary.FOPR = true;
			break;

		case Map_Str2Int("FOPT", 4):
			summary.FOPT = true;
			break;

		case Map_Str2Int("FGPR", 4):
			summary.FGPR = true;
			break;

		case Map_Str2Int("FGPT", 4):
			summary.FGPt = true;
			break;

		case Map_Str2Int("FWPR", 4):
			summary.FWPR = true;
			break;

		case Map_Str2Int("FWPT", 4):
			summary.FWPT = true;
			break;

		case Map_Str2Int("FGIR", 4):
			summary.FGIR = true;
			break;

		case Map_Str2Int("FGIT", 4):
			summary.FGIT = true;
			break;

		case Map_Str2Int("FWIR", 4):
			summary.FWIR = true;
			break;

		case Map_Str2Int("FWIT", 4):
			summary.FWIT = true;
			break;

		// Well
		case Map_Str2Int("WOPR", 4):
			InputType_A(ifs, summary.WOPR);
			break;

		case Map_Str2Int("WOPT", 4):
			InputType_A(ifs, summary.WOPT);
			break;

		case Map_Str2Int("WGPR", 4):
			InputType_A(ifs, summary.WGPR);
			break;

		case Map_Str2Int("WGPT", 4):
			InputType_A(ifs, summary.WGPT);
			break;

		case Map_Str2Int("WWPR", 4):
			InputType_A(ifs, summary.WWPR);
			break;

		case Map_Str2Int("WWPT", 4):
			InputType_A(ifs, summary.WWPT);
			break;

		case Map_Str2Int("WGIR", 4):
			InputType_A(ifs, summary.WGIR);
			break;

		case Map_Str2Int("WGIT", 4):
			InputType_A(ifs, summary.WGIT);
			break;

		case Map_Str2Int("WWIR", 4):
			InputType_A(ifs, summary.WWIR);
			break;

		case Map_Str2Int("WWIT", 4):
			InputType_A(ifs, summary.WWIT);
			break;

		case Map_Str2Int("WBHP", 4):
			InputType_A(ifs, summary.WBHP);
			break;

		case Map_Str2Int("BPR", 3):
			InputType_B(ifs, summary.BPR);
			break;
		}
	}
	cout << "SUMMARY" << endl;
}


void ParamOutput::InputType_A(ifstream& ifs, Type_A_o& obj)
{
	obj.activity = true;
	vector<string>		vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/") {
		obj.obj.push_back("All");
	}
	else {
		OCP_INT len = vbuf.size();
		for (OCP_INT i = 0; i < len - 1; i++) {
			obj.obj.push_back(vbuf[i]);
		}
		if (vbuf.back() != "/")
			obj.obj.push_back(vbuf.back());

		while (ReadLine(ifs, vbuf)) {
			if (vbuf[0] == "/")
				break;

			OCP_INT len = vbuf.size();
			for (OCP_INT i = 0; i < len - 1; i++) {
				obj.obj.push_back(vbuf[i]);
			}
			if (vbuf.back() != "/")
				obj.obj.push_back(vbuf.back());
		}
	}
	cout << "Type_A" << endl;
}

void ParamOutput::InputType_B(ifstream& ifs, Type_B_o& obj)
{

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf)) {
		if (vbuf[0] == "/")
			break;

		obj.activity = true;
		DealDefault(vbuf);
		USI i = stoi(vbuf[0]);
		USI j = stoi(vbuf[1]);
		USI k = stoi(vbuf[2]);

		obj.obj.push_back(COOIJK(i, j, k));
	}
	cout << "Type_B" << endl;
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*----------------------------------------------------------------------------*/