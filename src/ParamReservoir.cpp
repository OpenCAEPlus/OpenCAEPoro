#include "ParamReservoir.hpp"

vector<double>* ParamReservoir::FindPtr(string& varName)
{
	vector<double>*		myPtr = nullptr;

	switch (Map_str2int(&varName[0], varName.size()))
	{
	case Map_str2int("'DX'", 4):
		myPtr = &Dx;
		break;

	case Map_str2int("'DY'", 4):
		myPtr = &Dy;
		break;

	case Map_str2int("'DZ'", 4):
		myPtr = &Dz;
		break;

	case Map_str2int("'PORO'", 6):
		myPtr = &Poro;
		break;

	case Map_str2int("'PERMX'", 7):
		myPtr = &PermX;
		break;

	case Map_str2int("'PERMY'", 7):
		myPtr = &PermY;
		break;

	case Map_str2int("'PERMZ'", 7):
		myPtr = &PermZ;
		break;

	case Map_str2int("'NTG'", 5):
		myPtr = &Ntg;
		break;

	case Map_str2int("'TOPS'", 6):
		myPtr = &Tops;
		break;
	}
	return myPtr;
}

TableSet* ParamReservoir::FindPtr_T(string& varName)
{
	TableSet*	myPtr = nullptr;

	switch (Map_str2int(&varName[0], varName.size()))
	{
	case Map_str2int("SWOF", 4):
		myPtr = &SWOF_T;
		break;

	case Map_str2int("SGOF", 4):
		myPtr = &SGOF_T;
		break;

	case Map_str2int("PBVD", 4):
		myPtr = &PBVD_T;
		break;

	case Map_str2int("PVDG", 4):
		myPtr = &PVDG_T;
		break;

	case Map_str2int("PVCO", 4):
		myPtr = &PVCO_T;
		break;

	case Map_str2int("PVTW", 4):
		myPtr = &PVTW_T;
		break;
	}
	return myPtr;
}

void ParamReservoir::initVar()
{
	Tops.resize(Dimens.Nx * Dimens.Ny);
	Dx.resize(Num, 0);
	Dy.resize(Num, 0);
	Dz.resize(Num, 0);
	Ntg.resize(Num, 1);
	Poro.resize(Num, 0);
	PermX.resize(Num, 0);
	PermY.resize(Num, 0);
	PermZ.resize(Num, 0);
}

void ParamReservoir::initTab()
{
	SWOF_T.Name = "SWOF";   SWOF_T.colNum = 4;
	SGOF_T.Name = "SGOF";	SGOF_T.colNum = 4;
	PBVD_T.Name = "PBVD";	PBVD_T.colNum = 2;
	PVCO_T.Name = "PVCO";	PVCO_T.colNum = 6;
	PVDG_T.Name = "PVDG";	PVDG_T.colNum = 3;
	PVTW_T.Name = "PVTW";	PVTW_T.colNum = 5;
}

void ParamReservoir::setVal(vector<double>& obj, double val, vector<int>& index)
{
	int Nx = Dimens.Nx;
	int Ny = Dimens.Ny;
	int NxNy = Nx * Ny;
	int id = 0;

	for (int k = index[4]; k <= index[5]; k++) {
		for (int j = index[2]; j <= index[3]; j++) {
			for (int i = index[0]; i <= index[1]; i++) {
				id = k * NxNy + j * Nx + i;
				obj[id] = val;		
			}
		}
	}
}

void ParamReservoir::copyVal(vector<double>& obj, vector<double>& src, vector<int>& index)
{
	int Nx = Dimens.Nx;
	int Ny = Dimens.Ny;
	int NxNy = Nx * Ny;
	int id = 0;

	for (int k = index[4]; k <= index[5]; k++) {
		for (int j = index[2]; j <= index[3]; j++) {
			for (int i = index[0]; i <= index[1]; i++) {
				id = k * NxNy + j * Nx + i;
				obj[id] = src[id];
			}
		}
	}
}

void ParamReservoir::multiplyVal(vector<double>& obj, double val, vector<int>& index)
{
	int Nx = Dimens.Nx;
	int Ny = Dimens.Ny;
	int NxNy = Nx * Ny;
	int id = 0;

	for (int k = index[4]; k <= index[5]; k++) {
		for (int j = index[2]; j <= index[3]; j++) {
			for (int i = index[0]; i <= index[1]; i++) {
				id = k * NxNy + j * Nx + i;
				obj[id] *= val;
			}
		}
	}
}

void ParamReservoir::inputDIMENS(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	Dimens.Nx = atoi(vbuf[0].c_str());
	Dimens.Ny = atoi(vbuf[1].c_str());
	Dimens.Nz = atoi(vbuf[2].c_str());
	Num = Dimens.Nx * Dimens.Ny * Dimens.Nz;

	initVar();
}

void ParamReservoir::outputDIMENS()
{
	cout << Dimens.Nx << "  " << Dimens.Ny << "  " << Dimens.Nz << endl;
}

void ParamReservoir::inputEQUALS(ifstream& ifs)
{
	vector<int>		index(6, 0);
	vector<string>		vbuf;

	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		index[0] = 0, index[1] = Dimens.Nx - 1;
		index[2] = 0, index[3] = Dimens.Ny - 1;
		index[4] = 0, index[5] = Dimens.Nz - 1;

		string objName = vbuf[0];
		double val = atof(vbuf[1].c_str());

		DealDefault(vbuf);

		for (int n = 2; n < 8; n++) {
			if (vbuf[n] != DEFAULT)
				index[n - 2] = atoi(vbuf[n].c_str()) - 1;
		}

		vector<double>* objPtr = FindPtr(objName);

		if (objPtr != nullptr) {
			if (objName == "'TOPS'") {
				index[4] = index[5] = 0;
			}
			setVal(*objPtr, val, index);
		}
		else {
			Paramcheck("Wrong Var Name: " + objName);
			exit(0);
		}
	}	
	//cout << Dx[0] << endl;
	//cout << Dy[0] << endl;
	//cout << Dz[0] << endl;
	//cout << Poro[0] << endl;
	//cout << Ntg[0] << endl;
	//cout << PermX[0] << endl;
	//cout << PermY[0] << endl;
	//cout << PermZ[0] << endl;
	//cout << Tops[0] << endl;

}

void ParamReservoir::inputCOPY(ifstream& ifs)
{
	vector<string>		vbuf;
	vector<int>			index(6, 0);

	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		index[0] = 0, index[1] = Dimens.Nx - 1;
		index[2] = 0, index[3] = Dimens.Ny - 1;
		index[4] = 0, index[5] = Dimens.Nz - 1;

		string srcName = vbuf[0];
		string objName = vbuf[1];
		DealDefault(vbuf);
		for (int n = 2; n < 8; n++) {
			if (vbuf[n] != DEFAULT)
				index[n - 2] = atoi(vbuf[n].c_str()) - 1;
		}

		vector<double>* srcPtr = FindPtr(srcName);
		vector<double>* objPtr = FindPtr(objName);
		if (srcPtr != nullptr && objPtr != nullptr) {
			copyVal(*objPtr, *srcPtr, index);
		}
		else{
			Paramcheck("Wrong Var Name: " + srcName + " " + objName);
			exit(0);
		}
	}
	//cout << PermY[0] << endl;
}

void ParamReservoir::inputMULTIPLY(ifstream& ifs)
{
	vector<string>		vbuf;
	vector<int>			index(6, 0);

	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		index[0] = 0, index[1] = Dimens.Nx - 1;
		index[2] = 0, index[3] = Dimens.Ny - 1;
		index[4] = 0, index[5] = Dimens.Nz - 1;

		string objName = vbuf[0];
		double val = atof(vbuf[1].c_str());

		DealDefault(vbuf);
		for (int n = 2; n < 8; n++) {
			if (vbuf[n] != DEFAULT)
				index[n - 2] = atoi(vbuf[n].c_str()) - 1;
		}

		vector<double>* objPtr = FindPtr(objName);
		if (objPtr != nullptr) {
			if (objName == "'TOPS'") {
				index[4] = index[5] = 0;
			}
			multiplyVal(*objPtr, val, index);
		}
		else {
			Paramcheck("Wrong Var Name: " + objName);
			exit(0);
		}
	}
}

void ParamReservoir::inputTABLE(ifstream& ifs, string tabName)
{
	TableSet*	obj;
	obj = FindPtr_T(tabName);
	if (obj == nullptr) {
		Paramcheck("Wrong Tab Name :" + tabName);
		exit(0);
	}
	
	int col = obj->colNum;
	vector<vector<double>> tmpTab(col);

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		for (int i = 0; i < col; i++) {
			tmpTab[i].push_back(atof(vbuf[i].c_str()));
		}

		if (vbuf.back() == "/") {
			obj->data.push_back(tmpTab);
			for (int j = 0; j < col; j++) {
				tmpTab[j].clear();
			}		
		}	
	}
	if (!tmpTab[0].empty())
		obj->data.push_back(tmpTab);

	obj->showTab();
	cout << "hello" << endl;
}
