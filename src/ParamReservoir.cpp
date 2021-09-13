#include "ParamReservoir.hpp"

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

void ParamReservoir::setVal(double* obj, double val, vector<int>& index)
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

void ParamReservoir::copyVal(double* obj, double* src, vector<int>& index)
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

		double* objPtr = FindPtr(objName);

		if (objPtr != nullptr) {
			if (objName == "'TOPS'") {
				index[4] = index[5] = 0;
			}
			setVal(objPtr, val, index);
		}
		else {
			Paramcheck("Wrong Var Name: " + objName);
			exit(0);
		}
	}	
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

		double* srcPtr = FindPtr(srcName);
		double* objPtr = FindPtr(objName);
		if (srcPtr != nullptr && objPtr != nullptr) {
			copyVal(objPtr, srcPtr, index);
		}
		else{
			Paramcheck("Wrong Var Name: " + srcName + " " + objName);
			exit(0);
		}
	}
	cout << PermX[0] << endl;
	cout << PermY[0] << endl;
}


double* ParamReservoir::FindPtr(string& varName)
{
	double* myPtr = nullptr;

	switch (Map_str2int(&varName[0], varName.size()))
	{
	case Map_str2int("'DX'", 4):
		myPtr = &Dx[0];
		break;

	case Map_str2int("'DY'", 4):
		myPtr = &Dy[0];
		break;

	case Map_str2int("'DZ'", 4):
		myPtr = &Dz[0];
		break;

	case Map_str2int("'PORO'", 6):
		myPtr = &Poro[0];
		break;

	case Map_str2int("'PERMX'", 7):
		myPtr = &PermX[0];
		break;

	case Map_str2int("'PERMY'", 7):
		myPtr = &PermY[0];
		break;

	case Map_str2int("'PERMZ'", 7):
		myPtr = &PermZ[0];
		break;

	case Map_str2int("'NTG'", 5):
		myPtr = &Ntg[0];
		break;

	case Map_str2int("'TOPS'", 6):
		myPtr = &Tops[0];
		break;
	}
	return myPtr;
}