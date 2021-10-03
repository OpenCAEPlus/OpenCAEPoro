#include "ParamReservoir.hpp"

vector<OCP_DBL>* ParamReservoir::FindPtr(string& varName)
{
	vector<OCP_DBL>*		myPtr = nullptr;

	switch (Map_str2int(&varName[0], varName.size()))
	{
	case Map_str2int("DX", 2):
		myPtr = &Dx;
		break;

	case Map_str2int("DY", 2):
		myPtr = &Dy;
		break;

	case Map_str2int("DZ", 2):
		myPtr = &Dz;
		break;

	case Map_str2int("PORO", 4):
		myPtr = &Poro;
		break;

	case Map_str2int("NTG", 3):
		myPtr = &Ntg;
		break;

	case Map_str2int("PERMX", 5):
		myPtr = &PermX;
		break;

	case Map_str2int("PERMY", 5):
		myPtr = &PermY;
		break;

	case Map_str2int("PERMZ", 5):
		myPtr = &PermZ;
		break;

	case Map_str2int("TOPS", 4):
		myPtr = &Tops;
		break;

	case Map_str2int("PRESSURE", 8):
		myPtr = &Pressure;
		break;

	case Map_str2int("Ni", 2):
		myPtr = &Ni;
		break;

	case Map_str2int("SATNUM", 6):
		SATNUM.activity = true;
		myPtr = &SATNUM.data;
		break;

	case Map_str2int("PVTNUM", 6):
		PVTNUM.activity = true;
		myPtr = &PVTNUM.data;
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

	case Map_str2int("PVCO", 4):
		myPtr = &PVCO_T;
		break;

	case Map_str2int("PVDO", 4):
		myPtr = &PVDO_T;
		break;

	case Map_str2int("PVDG", 4):
		myPtr = &PVDG_T;
		break;
	
	case Map_str2int("PVTW", 4):
		myPtr = &PVTW_T;
		break;
	}
	return myPtr;
}

void ParamReservoir::init()
{
	initTab();

	RTEMP = 60;

	Gravity.data.resize(3);
	Gravity.data[0] = 45.5;			// oil
	Gravity.data[1] = 1.0;			// pure water
	Gravity.data[2] = 0.7773;		// air

	Density.data.resize(3);
	Density.data[0] = 37.457;		// oil
	Density.data[1] = 62.366416;	// pure water
	Density.data[2] = 0.062428;		// air

	Rock.Pref	= 14.7;
	Rock.Cr		= 3.406E-6;
}


void ParamReservoir::initTab()
{
	SWOF_T.Name = "SWOF";   SWOF_T.colNum = 4;
	SGOF_T.Name = "SGOF";	SGOF_T.colNum = 4;
	PBVD_T.Name = "PBVD";	PBVD_T.colNum = 2;
	PVCO_T.Name = "PVCO";	PVCO_T.colNum = 6;
	PVDO_T.Name = "PVDO";	PVDO_T.colNum = 3;
	PVDG_T.Name = "PVDG";	PVDG_T.colNum = 3;
	PVTW_T.Name = "PVTW";	PVTW_T.colNum = 5;
}

template<typename T>
void ParamReservoir::setVal(vector<T>& obj, T val, vector<int>& index)
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

template<typename T>
void ParamReservoir::copyVal(vector<T>& obj, vector<T>& src, vector<int>& index)
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

void ParamReservoir::multiplyVal(vector<OCP_DBL>& obj, OCP_DBL val, vector<int>& index)
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

void ParamReservoir::inputCOMPS(ifstream& ifs)
{
	COMPS = true;
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	Nc = stoi(vbuf[0]);
}


void ParamReservoir::inputDIMENS(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	Dimens.Nx = stoi(vbuf[0]);
	Dimens.Ny = stoi(vbuf[1]);
	Dimens.Nz = stoi(vbuf[2]);
	Num = Dimens.Nx * Dimens.Ny * Dimens.Nz;
}

void ParamReservoir::outputDIMENS()
{
	std::cout << "DIMENS" << endl;
	std::cout << Dimens.Nx << "  " << Dimens.Ny << "  " << Dimens.Nz << endl;
}

void ParamReservoir::inputRTEMP(ifstream& ifs)
{
	vector<string>		vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	RTEMP = stod(vbuf[0]);
	cout << "RTEMP" << endl;
	cout << RTEMP << endl;
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
		OCP_DBL val = stod(vbuf[1]);

		DealDefault(vbuf);

		for (int n = 2; n < 8; n++) {
			if (vbuf[n] != "DEFAULT")
				index[n - 2] = stoi(vbuf[n]) - 1;
		}

		vector<OCP_DBL>* objPtr = FindPtr(objName);

		if (objPtr != nullptr) {
			if (objName == "TOPS") {
				objPtr->resize(Dimens.Nx * Dimens.Ny);
				index[4] = index[5] = 0;
			}
			else {
				objPtr->resize(Num);
			}
			setVal(*objPtr, val, index);
		}
		else {
			Paramcheck("Wrong Var Name: " + objName);
			exit(0);
		}
	}	
	std::cout << &Dx << endl;
	std::cout << &Dy << endl;
	std::cout << &Dz << endl;
	std::cout << &Poro << endl;
	std::cout << &Ntg << endl;
	std::cout << &PermX << endl;
	std::cout << &PermY << endl;
	std::cout << &PermZ << endl;
	std::cout << &Tops << endl;
	std::cout << SATNUM.activity << endl;
	std::cout << PVTNUM.activity << endl;
	cout << "EQUALS" << endl;
}

void ParamReservoir::inputGRID(ifstream& ifs, string& keyword)
{
	vector<OCP_DBL>* objPtr = nullptr;
	objPtr = FindPtr(keyword);
	if (objPtr == nullptr) {
		Paramcheck("unmatched keyword!");
		exit(0);
	}
	else {
		if (keyword == "TOPS")
			objPtr->reserve(Dimens.Nx * Dimens.Ny);
		else
			objPtr->reserve(Num);
	}

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf)) {
		if (vbuf[0] == "/")
			break;
		
		for (auto str : vbuf) {
			objPtr->push_back(stod(str));
		}
	}
	std::cout << &PermX << endl;
	std::cout << &PermY << endl;
	std::cout << &PermZ << endl;
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
			if (vbuf[n] != "DEFAULT")
				index[n - 2] = stoi(vbuf[n]) - 1;
		}

		vector<OCP_DBL>* srcPtr = FindPtr(srcName);
		vector<OCP_DBL>* objPtr = FindPtr(objName);
		if (srcPtr != nullptr && objPtr != nullptr) {
			objPtr->resize(srcPtr->size());
			copyVal(*objPtr, *srcPtr, index);
		}
		else{
			Paramcheck("Wrong Var Name: " + srcName + " " + objName);
			exit(0);
		}
	}
	std::cout << PermX[0] << endl;
	std::cout << PermY[0] << endl;
	std::cout << PermZ[0] << endl;
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
		OCP_DBL val = stod(vbuf[1]);

		DealDefault(vbuf);
		for (int n = 2; n < 8; n++) {
			if (vbuf[n] != "DEFAULT")
				index[n - 2] = stoi(vbuf[n]) - 1;
		}

		vector<OCP_DBL>* objPtr = FindPtr(objName);
		if (objPtr != nullptr) {
			if (objName == "TOPS") {
				index[4] = index[5] = 0;
			}
			multiplyVal(*objPtr, val, index);
		}
		else {
			Paramcheck("Wrong Var Name: " + objName);
			exit(0);
		}
	}
	std::cout << PermX[0] << endl;
	std::cout << PermY[0] << endl;
	std::cout << PermZ[0] << endl;
}

void ParamReservoir::inputTABLE(ifstream& ifs, string& tabName)
{
	TableSet*	obj;
	obj = FindPtr_T(tabName);
	if (obj == nullptr) {
		Paramcheck("Wrong Tab Name :" + tabName);
		exit(0);
	}
	
	int col = obj->colNum;
	vector<vector<OCP_DBL>> tmpTab(col);

	vector<string>		vbuf;
	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		for (int i = 0; i < col; i++) {
			tmpTab[i].push_back(stod(vbuf[i]));
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
	std::cout << "TABLE" << endl;
}

void ParamReservoir::inputROCK(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	Rock.Pref	= stod(vbuf[0]);
	Rock.Cr		= stod(vbuf[1]);

	std::cout << "ROCK" << endl;
	std::cout << Rock.Pref << "  " << Rock.Cr << endl;
}

void ParamReservoir::inputGRAVITY(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	for (int i = 0; i < 3; i++) {
		if (vbuf[i] != "DEFAULT") {
			Gravity.activity = true;
			Gravity.data[i] = stod(vbuf[i]);
		}		
	}

	std::cout << "GRAVITY" << endl;
	std::cout << Gravity.data[0] << "  " <<  Gravity.data[1] << "  "  << Gravity.data[2] << endl;
}

void ParamReservoir::inputDENSITY(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	DealDefault(vbuf);
	for (int i = 0; i < 3; i++) {
		if (vbuf[i] != "DEFAULT") {
			Density.activity = true;
			Density.data[i] = stod(vbuf[i]);
		}	
	}

	std::cout << "DENSITY" << endl;
	std::cout << Density.data[0] << "  " << Density.data[1] << "  " << Density.data[2] << endl;
}

void ParamReservoir::inputEQUIL(ifstream& ifs)
{
	vector<string>	vbuf;
	ReadLine(ifs, vbuf);
	if (vbuf[0] == "/")
		return;

	EQUIL.resize(6, 0);
	DealDefault(vbuf);
	for (int i = 0; i < 6; i++) {
		if (vbuf[i] != "DEFAULT")
			EQUIL[i] = stod(vbuf[i]);
	}
	
	std::cout << "EQUIL" << endl;
	for (int i = 0; i < 6; i++)
		std::cout << EQUIL[i] << "  ";
	std::cout << endl;
}

void ParamReservoir::inputTABDIMS(ifstream& ifs)
{
	vector<string>  vbuf;
	ReadLine(ifs, vbuf);
	NTSFUN = stoi(vbuf[0]);
	NTPVT = stoi(vbuf[1]);
	cout << "TABDIMS" << endl;
}

void ParamReservoir::inputRegion(ifstream& ifs, string keyword)
{
	Type_A_r<OCP_DBL>* ptr = nullptr;
	int	lim = 0;
	if (keyword == "SATNUM") {
		ptr = &SATNUM; 
		lim = NTSFUN;
	}
	else {
		ptr = &PVTNUM;
		lim = NTPVT;
	}
		
	ptr->activity = true;
	ptr->data.reserve(Num);
	vector<string>	vbuf;
	vector<int>		obj;
	vector<int>		region;

	while (ReadLine(ifs, vbuf))
	{
		if (vbuf[0] == "/")
			break;

		DealData(vbuf, obj, region);

		// check region
		for (auto r : region) {
			if (r > lim) {
				Paramcheck("Region is out of Range!");
				exit(0);
			}
		}

		int len = obj.size();
		for (int i = 0; i < len; i++) {
			for (int j = 0; j < obj[i]; j++) {
				ptr->data.push_back(region[i]);
			}
		}
	}
	cout << &SATNUM << endl;
	cout << &PVTNUM << endl;
	cout << "Region" << endl;
}


// check
void ParamReservoir::checkParam()
{
	checkGrid();
	checkEQUIL();
	checkDenGra();
	checkPhase();
	checkPhaseTab();
	checkRegion();
}

void ParamReservoir::checkGrid()
{
	if (Tops.size() != Dimens.Nx * Dimens.Ny) {
		Paramcheck("Mistakes in Tops !");
		exit(0);
	}
	if (Dx.size() != Num) {
		Paramcheck("Mistakes in Dx !");
		exit(0);
	}
	if (Dy.size() != Num) {
		Paramcheck("Mistakes in Dy !");
		exit(0);
	}
	if (Dz.size() != Num) {
		Paramcheck("Mistakes in Dz !");
		exit(0);
	}
	if (Ntg.size() != Num) {
		Paramcheck("WARNING: Mistakes in Ntg, set to 1 .");
		Ntg.resize(Num, 1);
	}
	if (Poro.size() != Num) {
		Paramcheck("Mistakes in Poro !");
		exit(0);
	}
	if (PermX.size() != Num) {
		Paramcheck("Mistakes in PermX !");
		exit(0);
	}
	if (PermY.size() != Num) {
		Paramcheck("Mistakes in PermY !");
		exit(0);
	}
	if (PermZ.size() != Num) {
		Paramcheck("Mistakes in PermZ !");
		exit(0);
	}
}

void ParamReservoir::checkEQUIL()
{
	if (EQUIL.empty()) {
		Paramcheck("EQUIL is missing !");
		exit(0);
	}
}

void ParamReservoir::checkDenGra()
{
	if (Density.activity && Gravity.activity) {
		Paramcheck("Density and Gravity have been conflict, just one of them is needed !");
		exit(0);
	}	
}

void ParamReservoir::checkPhase()
{
	if (DISGAS && (!GAS && !OIL)) {
		Paramcheck("DISGAS can only be used only if OIL and GAS are both present");
		exit(0);
	}
}

void ParamReservoir::checkPhaseTab()
{
	if (!BLACKOIL && !COMPS) {
		Paramcheck("WRONG MODEl: choose BLACKOIL or COMPS !");
		exit(0);
	}

	if (WATER && OIL && SWOF_T.data.empty()) {
		Paramcheck("SWOF is missing !");
		exit(0);
	}

	if (GAS && OIL && SGOF_T.data.empty()) {
		Paramcheck("SGOF is missing !");
		exit(0);
	}

	if (WATER && PVTW_T.data.empty()) {
		Paramcheck("PVTW is missing !");
		exit(0);
	}

	if (BLACKOIL) {
		if (OIL && DISGAS && PVCO_T.data.empty()) {
			Paramcheck("PVCO is missing !");
			exit(0);
		}

		if (OIL && (!DISGAS) && PVDO_T.data.empty()) {
			Paramcheck("PVDO is missing !");
			exit(0);
		}

		if (GAS && PVDG_T.data.empty()) {
			Paramcheck("PVDG is missing !");
			exit(0);
		}
	}
}

void ParamReservoir::checkRegion()
{
	if (SATNUM.activity && SATNUM.data.size() != Num) {
		Paramcheck("missing data in SATNUM");
		exit(0);
	}
	if (PVTNUM.activity && PVTNUM.data.size() != Num) {
		Paramcheck("missing data in PVTNUM");
		exit(0);
	}
}
void ParamReservoir::checkEqlRegion()
{
	if (PBVD_T.data.size() > 1) {
		Paramcheck("Only one equilibration region is supported now !");
		exit(0);
	}
}
