#include "Bulk.hpp"
#include <ctime>
#include <algorithm>

void Bulk::inputParam(ParamReservoir& rs_param)
{
	Rock_Pref = rs_param.Rock.Pref;
	Rock_C1 = rs_param.Rock.Cr;
	Rock_C2 = Rock_C1;

	T = rs_param.RTEMP;
	BLACKOIL = rs_param.BLACKOIL;
	COMPS = rs_param.COMPS;
	OIL = rs_param.OIL;
	GAS = rs_param.GAS;
	WATER = rs_param.WATER;
	DISGAS = rs_param.DISGAS;

	EQUIL.Dref = rs_param.EQUIL[0];
	EQUIL.Pref = rs_param.EQUIL[1];
	EQUIL.PBVD.setup(rs_param.PBVD_T.data[0]);
	
	if (BLACKOIL) {
		if (WATER && !OIL && !GAS) {
			// water
			Np = 1;	Nc = 1;
			SATmode = PHASE_W;
			PVTmode = PHASE_W;
		}
		else if (WATER && OIL && !GAS) {
			// water, dead oil
			Np = 2;	Nc = 2;
			EQUIL.DOWC = rs_param.EQUIL[2];
			EQUIL.PcOWC = rs_param.EQUIL[3];
			SATmode = PHASE_OW;
			PVTmode = PHASE_OW;
		}
		else if (WATER && OIL && GAS && !DISGAS) {
			// water, dead oil, dry gas
			Np = 3;	Nc = 3;
			EQUIL.DOWC = rs_param.EQUIL[2];
			EQUIL.PcOWC = rs_param.EQUIL[3];
			EQUIL.DGOC = rs_param.EQUIL[4];
			EQUIL.PcGOC = rs_param.EQUIL[5];
			SATmode = PHASE_OGW;
			PVTmode = PHASE_OGW;      // maybe it should be added later
		}
		else if (WATER && OIL && GAS && DISGAS) {
			// water, live oil, dry gas
			Np = 3;	Nc = 3;
			EQUIL.DOWC = rs_param.EQUIL[2];
			EQUIL.PcOWC = rs_param.EQUIL[3];
			EQUIL.DGOC = rs_param.EQUIL[4];
			EQUIL.PcGOC = rs_param.EQUIL[5];
			SATmode = PHASE_OGW;
			PVTmode = PHASE_OGW;
		}
		rs_param.Np = Np;
		rs_param.Nc = Nc;
		for (int i = 0; i < rs_param.NTSFUN; i++)
			Flow.push_back(new FlowUnit(rs_param, SATmode, i));
		if (OIL & GAS & WATER) {
			for (int i = 0; i < rs_param.NTSFUN; i++) {
				Flow[i]->generate_SWPCWG();
			}
		}
		for (int i = 0; i < rs_param.NTPVT; i++)
			Flashcal.push_back(new BOMixture(rs_param, PVTmode, i));
		cout << "Bulk::inputParam" << endl;
	}
	else if (COMPS) {
		InitZi = rs_param.InitZi;
	}
}


void Bulk::setup(const Grid& myGrid)
{
	Num = myGrid.ActiveBulkNum;
	Dx.resize(Num, 0);
	Dy.resize(Num, 0);
	Dz.resize(Num, 0); 
	Depth.resize(Num, 0);
	Ntg.resize(Num, 0);
	Rock_VpInit.resize(Num, 0);
	Rock_Vp.resize(Num, 0);
	Rock_KxInit.resize(Num, 0);
	Rock_KyInit.resize(Num, 0);
	Rock_KzInit.resize(Num, 0);
	SATNUM.resize(Num, 0);
	PVTNUM.resize(Num, 0);
	

	for (int bIdb = 0; bIdb < Num; bIdb++) {
		int bIdg = myGrid.ActiveMap_B2G[bIdb];

		Dx[bIdb] = myGrid.Dx[bIdg];
		Dy[bIdb] = myGrid.Dy[bIdg];
		Dz[bIdb] = myGrid.Dz[bIdg];
		Depth[bIdb] = myGrid.Depth[bIdg];
		Ntg[bIdb] = myGrid.Ntg[bIdg];

		Rock_VpInit[bIdb] = myGrid.V[bIdg] * myGrid.Ntg[bIdg] * myGrid.Poro[bIdg];
		// Rock_PoroInit[bIdb] = myGrid.Poro[bIdg];
		Rock_KxInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KyInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KzInit[bIdb] = myGrid.Kx[bIdg];

		SATNUM[bIdb] = myGrid.SATNUM[bIdg];
		PVTNUM[bIdb] = myGrid.PVTNUM[bIdg];
	}
	Rock_Vp = Rock_VpInit;
	// Rock_Poro = Rock_PoroInit;
	Rock_Kx = Rock_KxInit;
	Rock_Ky = Rock_KyInit;
	Rock_Kz = Rock_KzInit;

	// physical variable
	P.resize(Num);
	Pj.resize(Num * Np);
	Pc.resize(Num * Np);
	PhaseExist.resize(Num * Np);
	Ni.resize(Num * Nc);
	S.resize(Num * Np);
	Xi.resize(Num * Np);
	Cij.resize(Num * Np * Nc);
	Rho.resize(Num * Np);
	Mu.resize(Num * Np);
	Kr.resize(Num * Np);

	Vf.resize(Num, 0);
	Vfi.resize(Num * Nc, 0);
	Vfp.resize(Num, 0);

	lP = P;
	lNi = Ni;
}

void Bulk::initSjPc_blk(int tabrow)
{
	Pbub.resize(Num);

	double Dref = EQUIL.Dref;		double Pref = EQUIL.Pref;
	double DOWC = EQUIL.DOWC;		double PcOWC = EQUIL.PcOWC;
	double DOGC = EQUIL.DGOC;		double PcGOC = EQUIL.PcGOC;

	double Zmin = 1E8;
	double Zmax = 0;
	for (int n = 0; n < Num; n++) {
		double temp1 = Depth[n] - Dz[n] / 2;
		double temp2 = Depth[n] + Dz[n] / 2;
		Zmin = Zmin < temp1 ? Zmin : temp1;
		Zmax = Zmax > temp2 ? Zmax : temp2;
	}
	double tabdz = (Zmax - Zmin) / (tabrow - 1);

	// creater table
	ReservoirTable<double>	DepthP(tabrow, 4);
	vector<double>&		Ztmp = DepthP.getCol(0);
	vector<double>&		Potmp = DepthP.getCol(1);
	vector<double>&		Pgtmp = DepthP.getCol(2);
	vector<double>&		Pwtmp = DepthP.getCol(3);

	// cal Tab_Ztmp
	Ztmp[0] = Zmin;
	for (int i = 1; i < tabrow; i++) {
		Ztmp[i] = Ztmp[i - 1] + tabdz;
	}

	int beginId = 0;
	// find the RefId
	if (Dref <= Ztmp[0]) {
		beginId = 0;
	}
	else if (Dref >= Ztmp[tabrow - 1]) {
		beginId = tabrow - 1;
	}
	else {
		beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(), [s = Dref](auto t) {return t > s; }));
		beginId--;
	}


	// begin calculating oil pressure:
	double Pbb = Pref;
	double gammaOtmp, gammaWtmp, gammaGtmp;
	double Ptmp;
	int mynum = 10;  double mydz = 0;
	double Poref, Pgref, Pwref;
	double Pbegin = 0;

	if (Dref < DOGC) {
		//reference pressure is gas pressure
		if (Flow[0]->empty_SGOF()) {
			ERRORcheck("SGOF is missing !");
			exit(0);
		}

		Pgref = Pref;
		gammaGtmp = Flashcal[0]->gammaPhaseG(Pgref);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;

		//find the gas pressure
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the oil pressure in Dref by Pgref
		Poref = 0; Ptmp = Pgref;
		mydz = (DOGC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal[0]->gammaPhaseG(Ptmp);
			Ptmp += gammaGtmp * mydz;
		}
		Ptmp -= PcGOC;
		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, DOGC - i * mydz, 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		if (!EQUIL.PBVD.isempty()) {
			Pbb = EQUIL.PBVD.eval(0, Dref, 1);
		}
		gammaOtmp = Flashcal[0]->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else if (Dref > DOWC)
	{
		// reference pressure is water pressure
		Pwref = Pref;
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		//find the water pressure
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the oil pressure in Dref by Pwref
		Poref = 0; Ptmp = Pwref;
		mydz = (DOWC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp += gammaWtmp * mydz;
		}
		Ptmp += PcOWC;

		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, DOWC - i * mydz, 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		if (!EQUIL.PBVD.isempty()) {
			Pbb = EQUIL.PBVD.eval(0, Dref, 1);
		}
		gammaOtmp = Flashcal[0]->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		if (!Flow[0]->empty_SGOF()) {
			// find the gas pressure in Dref by Poref
			Pgref = 0; Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;

			for (int i = 0; i < mynum; i++) {
				if (!EQUIL.PBVD.isempty()) {
					Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
				}
				gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
				Ptmp += gammaOtmp * mydz;
			}
			Ptmp += PcGOC;
			for (int i = 0; i < mynum; i++) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Ptmp);
				Ptmp -= gammaGtmp * mydz;
			}
			Pgref = Ptmp;

			//find the gas pressure in tab
			gammaGtmp = Flashcal[0]->gammaPhaseG(Pgref);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (int id = beginId; id > 0; id--) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}
			for (int id = beginId; id < tabrow - 1; id++) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
			}
		}
		
	}
	else
	{
		// reference pressure is oil pressure
		Poref = Pref;
		if (!EQUIL.PBVD.isempty()) {
			Pbb = EQUIL.PBVD.eval(0, Dref, 1);
		}
		gammaOtmp = Flashcal[0]->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		//find the oil pressure in tab
		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		if (!Flow[0]->empty_SGOF()) {
			// find the gas pressure in Dref by Poref
			Pgref = 0; Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;

			for (int i = 0; i < mynum; i++) {
				if (!EQUIL.PBVD.isempty()) {
					Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
				}
				gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
				Ptmp += gammaOtmp * mydz;
			}
			Ptmp += PcGOC;
			for (int i = 0; i < mynum; i++) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Ptmp);
				Ptmp -= gammaGtmp * mydz;
			}
			Pgref = Ptmp;

			//find the gas pressure in tab
			gammaGtmp = Flashcal[0]->gammaPhaseG(Pgref);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (int id = beginId; id > 0; id--) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}

			for (int id = beginId; id < tabrow - 1; id++) {
				gammaGtmp = Flashcal[0]->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
			}
		}
		


		// find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
			}
			gammaOtmp = Flashcal[0]->gammaPhaseO(Ptmp, Pbb);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}

	DepthP.display();

	// calculate Pc from DepthP to calculate Sj
	std::vector<double>		data(4, 0);
	std::vector<double>		cdata(4, 0);
	for (int n = 0; n < Num; n++) {
		DepthP.eval_all(0, Depth[n], data, cdata);
		double Po = data[1];
		double Pg = data[2];
		double Pw = data[3];
		double Pcgo = Pg - Po;
		double Pcow = Po - Pw;
		double Sw = Flow[0]->evalinv_SWOF(3, Pcow, 0);
		double Sg = 0;
		if (!Flow[0]->empty_SGOF()) {
			Sg = Flow[0]->eval_SGOF(3, Pcgo, 0);
		}
		if (Sw + Sg > 1) {
			// should me modified
			double Pcgw = Pcow + Pcgo;
			Sw = Flow[0]->evalinv_SWPCWG(1, Pcgw, 0);
			Sg = 1 - Sw;
		}

		if (1 - Sw < TINY) {
			// all water
			Po = Pw + Flow[0]->eval_SWOF(0, 1.0, 3);
			// Pg = Po + Flow->eval_SGOF(0, 0.0, 3);
		}
		else if (1 - Sg < TINY) {
			// all gas
			Po = Pg - Flow[0]->eval_SGOF(0, 1.0, 3);
			// Pw = Po - Flow->eval_SWOF(0, 0.0, 3);
		}
		else if (1 - Sw - Sg < TINY) {
			// water and gas
			Po = Pg - Flow[0]->eval_SGOF(0, Sg, 3);
			// Pw = Po - Flow->eval_SWOF(0, Sw, 3);
		}
		P[n] = Po;
		lP[n] = Po;
		Pbub[n] = Pbb;
		S[n * Np + Np - 1] = Sw;
		if (!Flow[0]->empty_SGOF()) {
			S[n * Np + Np - 2] = Sg;
		}

	}
}

void Bulk::initSjPc_comp(int tabrow)
{
	InitZi.resize(Num * Nc);


	double Dref = EQUIL.Dref;		double Pref  = EQUIL.Pref;
	double DOWC = EQUIL.DOWC;		double PcOWC = EQUIL.PcOWC;
	double DOGC = EQUIL.DGOC;		double PcGOC = EQUIL.PcGOC;
	
	double Zmin = 1E8;
	double Zmax = 0;
	for (int n = 0; n < Num; n++) {
		double temp1 = Depth[n] - Dz[n] / 2;
		double temp2 = Depth[n] + Dz[n] / 2;
		Zmin = Zmin < temp1 ? Zmin : temp1;
		Zmax = Zmax > temp2 ? Zmax : temp2;
	}
	double tabdz = (Zmax - Zmin) / (tabrow - 1);

	// creater table
	ReservoirTable<double>	DepthP(tabrow, 4);
	vector<double>& Ztmp = DepthP.getCol(0);
	vector<double>& Potmp = DepthP.getCol(1);
	vector<double>& Pgtmp = DepthP.getCol(2);
	vector<double>& Pwtmp = DepthP.getCol(3);

	// cal Tab_Ztmp
	Ztmp[0] = Zmin;
	for (int i = 1; i < tabrow; i++) {
		Ztmp[i] = Ztmp[i - 1] + tabdz;
	}

	int beginId = 0;
	// find the RefId
	if (Dref <= Ztmp[0]) {
		beginId = 0;
	}
	else if (Dref >= Ztmp[tabrow - 1]) {
		beginId = tabrow - 1;
	}
	else {
		beginId = distance(Ztmp.begin(), find_if(Ztmp.begin(), Ztmp.end(), [s = Dref](auto t) {return t > s; }));
		beginId--;
	}
	

	// begin calculating oil pressure:
	double mytemp = T;
	double gammaOtmp, gammaWtmp, gammaGtmp;
	double Ptmp;
	int mynum = 10;  double mydz = 0;
	double Poref, Pgref, Pwref;
	double Pbegin = 0;

	if (Dref < DOGC) {
		//reference pressure is gas pressure
		Pgref = Pref;
		gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		//find the gas pressure
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the oil pressure in Dref by Pgref
		Poref = 0; Ptmp = Pgref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaGtmp * mydz;
		}
		Ptmp -= PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		gammaOtmp = Flashcal[0]->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else if (Dref > DOWC) 
	{
		// reference pressure is water pressure
		Pwref = Pref;
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		//find the water pressure
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the oil pressure in Dref by Pwref
		Poref = 0; Ptmp = Pwref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp += gammaWtmp * mydz;
		}
		Ptmp += PcOWC;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		gammaOtmp = Flashcal[0]->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the gas pressure in Dref by Poref
		Pgref = 0; Ptmp = Poref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp += PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaGtmp * mydz;
		}
		Pgref = Ptmp;

		//find the gas pressure in tab
		gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else
	{
		// reference pressure is oil pressure
		Poref = Pref;
		gammaOtmp = Flashcal[0]->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		//find the oil pressure
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}	
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the gas pressure in Dref by Poref
		Pgref = 0; Ptmp = Poref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp += PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaGtmp * mydz;
		}
		Pgref = Ptmp;

		//find the gas pressure in tab
		gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal[0]->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal[0]->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal[0]->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal[0]->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}

	DepthP.display();

	// calculate Pc from DepthP to calculate Sj
	std::vector<double>		data(4, 0);
	std::vector<double>		cdata(4, 0);
	for (int n = 0; n < Num; n++) {
		DepthP.eval_all(0, Depth[n], data, cdata);
		double Po = data[1];
		double Pg = data[2];
		double Pw = data[3];
		double Pcgo = Pg - Po;
		double Pcow = Po - Pw;
		double Sw = Flow[0]->evalinv_SWOF(3, Pcow, 0);
		double Sg = 0;
		if (!Flow[0]->empty_SGOF()) {
			Sg = Flow[0]->eval_SGOF(3, Pcgo, 0);
		}
		if (Sw + Sg > 1) {
			// should me modified
			double Pcgw = Pcow + Pcgo;
			Sw = Flow[0]->evalinv_SWPCWG(1, Pcgw, 0);
			Sg = 1 - Sw;
		}

		if (1 - Sw < TINY) {
			// all water
			Po = Pw + Flow[0]->eval_SWOF(0, 1.0, 3);
			// Pg = Po + Flow->eval_SGOF(0, 0.0, 3);
		}
		else if (1 - Sg < TINY) {
			// all gas
			Po = Pg - Flow[0]->eval_SGOF(0, 1.0, 3);
			// Pw = Po - Flow->eval_SWOF(0, 0.0, 3);
		}
		else if (1 - Sw - Sg < TINY) {
			// water and gas
			Po = Pg - Flow[0]->eval_SGOF(0, Sg, 3);
			// Pw = Po - Flow->eval_SWOF(0, Sw, 3);
		}
		P[n] = Po;
		lP[n] = Po;
		S[n * Np + Np - 1] = Sw;
		if (!Flow[0]->empty_SGOF()) {
			S[n * Np + Np - 2] = Sg;
		}	
	}
}


// Flash
void Bulk::flash_Sj()
{
	for (int n = 0; n < Num; n++) {
		Flashcal[PVTNUM[n]]->Flash_Sj(P[n], Pbub[n], T, &S[n * Np], Rock_Vp[n], InitZi.data());
		for (int i = 0; i < Nc; i++) {
			Ni[n * Nc + i] = Flashcal[PVTNUM[n]]->Ni[i];
		}
		passFlashValue(n);
	}
}

void Bulk::flash_Ni()
{
	for (int n = 0; n < Num; n++) {
		Flashcal[PVTNUM[n]]->Flash_Ni(P[n], T, &Ni[n * Nc]);
		passFlashValue(n);
	}
}

void Bulk::passFlashValue(int n)
{
	int bId = n * Np;
	int pvtnum = PVTNUM[n];
	for (int j = 0; j < Np; j++) {
		PhaseExist[bId + j] = Flashcal[pvtnum]->PhaseExist[j];
		if (PhaseExist[j]) {
			for (int i = 0; i < Nc; i++) {
				Cij[bId * Nc + j * Nc + i] = Flashcal[pvtnum]->Cij[j * Nc + i];
			}
			S[bId + j]		= Flashcal[pvtnum]->S[j];
			Xi[bId + j]		= Flashcal[pvtnum]->Xi[j];
			Rho[bId + j]	= Flashcal[pvtnum]->Rho[j];
			Mu[bId + j]		= Flashcal[pvtnum]->Mu[j];
		}
	}
	Vf[n] = Flashcal[pvtnum]->Vf;
	Vfp[n] = Flashcal[pvtnum]->Vfp;
	bId = n * Nc;
	for (int i = 0; i < Nc; i++) {
		Vfi[bId + i] = Flashcal[pvtnum]->Vfi[i];
	}
}

// relative permeability and capillary pressure
void Bulk::calKrPc()
{
	for (int n = 0; n < Num; n++) {
		int bId = n * Np;
		Flow[SATNUM[n]]->calKrPc(&S[bId], &Kr[bId], &Pc[bId]);
		for (int j = 0; j < Np; j++)
			Pj[n * Np + j] = P[n] + Pc[n * Np + j];
	}
}

// Rock
void Bulk::calVporo()
{
	for (int n = 0; n < Num; n++) {
		double dP = Rock_C1 * (P[n] - Rock_Pref);
		Rock_Vp[n] = Rock_VpInit[n] * (1 + dP + dP * dP / 2);
	}
}

int Bulk::mixMode()
{
	if (BLACKOIL)
		return BLKOIL;
	if (COMPS)
		return EoS_PVTW;
}

void Bulk::getP_IMPES(vector<double>& u)
{
	for (int n = 0; n < Num; n++)
		P[n] = u[n];
}
