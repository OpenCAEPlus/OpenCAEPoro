#include "Bulk.hpp"
#include <ctime>
#include <algorithm>


void Bulk::init(const Grid& myGrid)
{
	Num = myGrid.ActiveBulkNum;
	Dx.resize(Num, 0);
	Dy.resize(Num, 0);
	Dz.resize(Num, 0); 
	Depth.resize(Num, 0);
	Ntg.resize(Num, 0);
	Rock_V.resize(Num, 0);
	Rock_PoroInit.resize(Num, 0);
	Rock_KxInit.resize(Num, 0);
	Rock_KyInit.resize(Num, 0);
	Rock_KzInit.resize(Num, 0);

	for (int bIdb = 0; bIdb < Num; bIdb++) {
		int bIdg = myGrid.ActiveMap_B2G[bIdb];

		Dx[bIdb] = myGrid.Dx[bIdg];
		Dy[bIdb] = myGrid.Dy[bIdg];
		Dz[bIdb] = myGrid.Dz[bIdg];
		Depth[bIdb] = myGrid.Depth[bIdg];
		Ntg[bIdb] = myGrid.Ntg[bIdg];

		Rock_V[bIdb] = myGrid.V[bIdg] * Ntg[bIdb];
		Rock_PoroInit[bIdb] = myGrid.Poro[bIdg];
		Rock_KxInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KyInit[bIdb] = myGrid.Kx[bIdg];
		Rock_KzInit[bIdb] = myGrid.Kx[bIdg];
	}
}

void Bulk::initSjPc_blk(int tabrow)
{
	double Dref = EQUIL.Dref;		double Pref = EQUIL.Pref;
	double DOWC = EQUIL.DOWC;		double PcOWC = EQUIL.PcOWC;
	double DOGC = EQUIL.DGOC;		double PcGOC = EQUIL.PcGOC;

	double myT = T[0];

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
	vector<double>  Ztmp(tabrow, 0);
	vector<double>  Potmp(tabrow, 0);
	vector<double>	Pgtmp(tabrow, 0);
	vector<double>  Pwtmp(tabrow, 0);

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
	double Ptmp, Pcowtmp, Pcgotmp, Pcowg;
	int mynum = 10;  double mydz = 0;
	double Poref, Pgref, Pwref;
	double Pbegin = 0;

	if (Dref < DOGC) {
		//reference pressure is gas pressure
		if (Flow->empty_SGOF()) {
			ERRORcheck("SGOF is missing !");
			exit(0);
		}

		Pgref = Pref;
		gammaGtmp = Flashcal->gammaPhaseG(Pgref);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;

		//find the gas pressure
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the oil pressure in Dref by Pgref
		Poref = 0; Ptmp = Pgref;
		mydz = (DOGC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal->gammaPhaseG(Ptmp);
			Ptmp += gammaGtmp * mydz;
		}
		Ptmp -= PcGOC;
		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, DOGC - i * mydz, 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		if (!EQUIL.PBVD.isempty()) {
			Pbb = EQUIL.PBVD.eval(0, Dref, 1);
		}
		gammaOtmp = Flashcal->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else if (Dref > DOWC)
	{
		// reference pressure is water pressure
		Pwref = Pref;
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		//find the water pressure
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the oil pressure in Dref by Pwref
		Poref = 0; Ptmp = Pwref;
		mydz = (DOWC - Dref) / mynum;

		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp += gammaWtmp * mydz;
		}
		Ptmp += PcOWC;

		for (int i = 0; i < mynum; i++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, DOWC - i * mydz, 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		if (!EQUIL.PBVD.isempty()) {
			Pbb = EQUIL.PBVD.eval(0, Dref, 1);
		}
		gammaOtmp = Flashcal->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		if (!Flow->empty_SGOF()) {
			// find the gas pressure in Dref by Poref
			Pgref = 0; Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;

			for (int i = 0; i < mynum; i++) {
				if (!EQUIL.PBVD.isempty()) {
					Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
				}
				gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
				Ptmp += gammaOtmp * mydz;
			}
			Ptmp += PcGOC;
			for (int i = 0; i < mynum; i++) {
				gammaGtmp = Flashcal->gammaPhaseG(Ptmp);
				Ptmp -= gammaGtmp * mydz;
			}
			Pgref = Ptmp;

			//find the gas pressure in tab
			gammaGtmp = Flashcal->gammaPhaseG(Pgref);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (int id = beginId; id > 0; id--) {
				gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}
			for (int id = beginId; id < tabrow - 1; id++) {
				gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
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
		gammaOtmp = Flashcal->gammaPhaseO(Poref, Pbb);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;

		//find the oil pressure
		for (int id = beginId; id > 0; id--) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			if (!EQUIL.PBVD.isempty()) {
				Pbb = EQUIL.PBVD.eval(0, Ztmp[id], 1);
			}
			gammaOtmp = Flashcal->gammaPhaseO(Potmp[id], Pbb);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		if (!Flow->empty_SGOF()) {
			// find the gas pressure in Dref by Poref
			Pgref = 0; Ptmp = Poref;
			mydz = (DOGC - Dref) / mynum;

			for (int i = 0; i < mynum; i++) {
				if (!EQUIL.PBVD.isempty()) {
					Pbb = EQUIL.PBVD.eval(0, Dref + i * mydz, 1);
				}
				gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
				Ptmp += gammaOtmp * mydz;
			}
			Ptmp += PcGOC;
			for (int i = 0; i < mynum; i++) {
				gammaOtmp = Flashcal->gammaPhaseG(Ptmp);
				Ptmp -= gammaGtmp * mydz;
			}
			Pgref = Ptmp;

			//find the gas pressure in tab
			gammaGtmp = Flashcal->gammaPhaseG(Pgref);
			Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
			Pgtmp[beginId] = Pbegin;

			for (int id = beginId; id > 0; id--) {
				gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
				Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
			}

			for (int id = beginId; id < tabrow - 1; id++) {
				gammaGtmp = Flashcal->gammaPhaseG(Pgtmp[id]);
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
			gammaOtmp = Flashcal->gammaPhaseO(Ptmp, Pbb);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;

		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}

		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}

	ReservoirTable<double>	DepthP;
	DepthP.pushCol(Ztmp);
	DepthP.pushCol(Potmp);
	DepthP.pushCol(Pgtmp);
	DepthP.pushCol(Pwtmp);
	DepthP.display();

	// calculate Pc from DepthP to calculate Sj
	std::vector<double>		data;
	std::vector<double>		cdata;
	for (int n = 0; n < Num; n++) {
		DepthP.eval_all(0, Depth[n], data, cdata);
		double Po = data[1];
		double Pg = data[2];
		double Pw = data[3];
		double Pcgo = Pg - Po;
		double Pcow = Po - Pw;
		double Sw = Flow->evalinv_SWOF(3, Pcow, 0);
		double Sg = 0;
		if (!Flow->empty_SGOF()) {
			Sg = Flow->eval_SGOF(3, Pcgo, 0);
		}
		if (Sw + Sg > 1) {
			// should me modified
			double Pcgw = Pcow + Pcgo;
			Sw = Flow->evalinv_SWPCWG(1, Pcgw, 0);
			Sg = 1 - Sw;
		}

		if (1 - Sw < TINY) {
			// all water
			Po = Pw + Flow->eval_SWOF(0, 1.0, 3);
			// Pg = Po + Flow->eval_SGOF(0, 0.0, 3);
		}
		else if (1 - Sg < TINY) {
			// all gas
			Po = Pg - Flow->eval_SGOF(0, 1.0, 3);
			// Pw = Po - Flow->eval_SWOF(0, 0.0, 3);
		}
		else if (1 - Sw - Sg < TINY) {
			// water and gas
			Po = Pg - Flow->eval_SGOF(0, Sg, 3);
			// Pw = Po - Flow->eval_SWOF(0, Sw, 3);
		}
		P[n] = Po;
		lP[n] = Po;
		S[n * Np + Np - 1] = Sw;
		if (!Flow->empty_SGOF()) {
			S[n * Np + Np - 2] = Sg;
		}

	}
}

void Bulk::initSjPc_comp(int tabrow)
{

	double Dref = EQUIL.Dref;		double Pref  = EQUIL.Pref;
	double DOWC = EQUIL.DOWC;		double PcOWC = EQUIL.PcOWC;
	double DOGC = EQUIL.DGOC;		double PcGOC = EQUIL.PcGOC;

	double myT = T[0];
	
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
	vector<double>  Ztmp(tabrow, 0);
	vector<double>  Potmp(tabrow, 0);
	vector<double>	Pgtmp(tabrow, 0);
	vector<double>  Pwtmp(tabrow, 0);

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
	double mytemp = T[0];
	double gammaOtmp, gammaWtmp, gammaGtmp;
	double Ptmp, Pcowtmp, Pcgotmp, Pcowg;
	int mynum = 10;  double mydz = 0;
	double Poref, Pgref, Pwref;
	double Pbegin = 0;

	if (Dref < DOGC) {
		//reference pressure is gas pressure
		Pgref = Pref;
		gammaGtmp = Flashcal->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		//find the gas pressure
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the oil pressure in Dref by Pgref
		Poref = 0; Ptmp = Pgref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaGtmp * mydz;
		}
		Ptmp -= PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		gammaOtmp = Flashcal->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		//find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else if (Dref > DOWC) 
	{
		// reference pressure is water pressure
		Pwref = Pref;
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		//find the water pressure
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the oil pressure in Dref by Pwref
		Poref = 0; Ptmp = Pwref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp += gammaWtmp * mydz;
		}
		Ptmp += PcOWC;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaOtmp * mydz;
		}
		Poref = Ptmp;

		//find the oil pressure in tab
		gammaOtmp = Flashcal->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the gas pressure in Dref by Poref
		Pgref = 0; Ptmp = Poref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp += PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaGtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaGtmp * mydz;
		}
		Pgref = Ptmp;

		//find the gas pressure in tab
		gammaGtmp = Flashcal->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}
	else
	{
		// reference pressure is oil pressure
		Poref = Pref;
		gammaOtmp = Flashcal->gammaPhaseOG(Poref, mytemp, &InitZi[0]);
		Pbegin = Poref + gammaOtmp * (Ztmp[beginId] - Dref);
		Potmp[beginId] = Pbegin;
		
		//find the oil pressure
		for (int id = beginId; id > 0; id--) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id - 1] = Potmp[id] - gammaOtmp * (Ztmp[id] - Ztmp[id - 1]);
		}	
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Potmp[id], mytemp, &InitZi[0]);
			Potmp[id + 1] = Potmp[id] + gammaOtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the gas pressure in Dref by Poref
		Pgref = 0; Ptmp = Poref;
		mydz = (DOGC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp += PcGOC;
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp -= gammaGtmp * mydz;
		}
		Pgref = Ptmp;

		//find the gas pressure in tab
		gammaGtmp = Flashcal->gammaPhaseOG(Pgref, mytemp, &InitZi[0]);
		Pbegin = Pgref + gammaGtmp * (Ztmp[beginId] - Dref);
		Pgtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id - 1] = Pgtmp[id] - gammaGtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaGtmp = Flashcal->gammaPhaseOG(Pgtmp[id], mytemp, &InitZi[0]);
			Pgtmp[id + 1] = Pgtmp[id] + gammaGtmp * (Ztmp[id + 1] - Ztmp[id]);
		}


		// find the water pressure in Dref by Poref
		Pwref = 0; Ptmp = Poref;
		mydz = (DOWC - Dref) / mynum;
		
		for (int i = 0; i < mynum; i++) {
			gammaOtmp = Flashcal->gammaPhaseOG(Ptmp, mytemp, &InitZi[0]);
			Ptmp += gammaOtmp * mydz;
		}
		Ptmp -= PcOWC;
		for (int i = 0; i < mynum; i++) {
			gammaWtmp = Flashcal->gammaPhaseW(Ptmp);
			Ptmp -= gammaWtmp * mydz;
		}
		Pwref = Ptmp;

		//find the water pressure in tab
		gammaWtmp = Flashcal->gammaPhaseW(Pwref);
		Pbegin = Pwref + gammaWtmp * (Ztmp[beginId] - Dref);
		Pwtmp[beginId] = Pbegin;
		
		for (int id = beginId; id > 0; id--) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id - 1] = Pwtmp[id] - gammaWtmp * (Ztmp[id] - Ztmp[id - 1]);
		}
		
		for (int id = beginId; id < tabrow - 1; id++) {
			gammaWtmp = Flashcal->gammaPhaseW(Pwtmp[id]);
			Pwtmp[id + 1] = Pwtmp[id] + gammaWtmp * (Ztmp[id + 1] - Ztmp[id]);
		}
	}

	ReservoirTable<double>	DepthP;
	DepthP.pushCol(Ztmp);
	DepthP.pushCol(Potmp);
	DepthP.pushCol(Pgtmp);
	DepthP.pushCol(Pwtmp);
	DepthP.display();

	// calculate Pc from DepthP to calculate Sj
	std::vector<double>		data;
	std::vector<double>		cdata;
	for (int n = 0; n < Num; n++) {
		DepthP.eval_all(0, Depth[n], data, cdata);
		double Po = data[1];
		double Pg = data[2];
		double Pw = data[3];
		double Pcgo = Pg - Po;
		double Pcow = Po - Pw;
		double Sw = Flow->evalinv_SWOF(3, Pcow, 0);
		double Sg = 0;
		if (!Flow->empty_SGOF()) {
			Sg = Flow->eval_SGOF(3, Pcgo, 0);
		}
		if (Sw + Sg > 1) {
			// should me modified
			double Pcgw = Pcow + Pcgo;
			Sw = Flow->evalinv_SWPCWG(1, Pcgw, 0);
			Sg = 1 - Sw;
		}

		if (1 - Sw < TINY) {
			// all water
			Po = Pw + Flow->eval_SWOF(0, 1.0, 3);
			// Pg = Po + Flow->eval_SGOF(0, 0.0, 3);
		}
		else if (1 - Sg < TINY) {
			// all gas
			Po = Pg - Flow->eval_SGOF(0, 1.0, 3);
			// Pw = Po - Flow->eval_SWOF(0, 0.0, 3);
		}
		else if (1 - Sw - Sg < TINY) {
			// water and gas
			Po = Pg - Flow->eval_SGOF(0, Sg, 3);
			// Pw = Po - Flow->eval_SWOF(0, Sw, 3);
		}
		P[n] = Po;
		lP[n] = Po;
		S[n * Np + Np - 1] = Sw;
		if (!Flow->empty_SGOF()) {
			S[n * Np + Np - 2] = Sg;
		}	
	}
}


// Flash
void Bulk::flash_Sj()
{
	for (int n = 0; n < Num; n++) {
		Flashcal->Flash_Sj(P[n], Pbb[n], T[n], &S[n * Np], Rock_V[n] * Rock_Poro[n], &InitZi[0]);
		for (int i = 0; i < Nc; i++) {
			Ni[n * Nc + i] = Flashcal->Ni[i];
		}
		passFlashValue(n);
	}
}

void Bulk::flash_Ni()
{
	for (int n = 0; n < Num; n++) {
		Flashcal->Flash_Ni(P[n], T[n], &Ni[n * Nc]);
		passFlashValue(n);
	}
}

void Bulk::passFlashValue(int n)
{
	int bId = n * Np;
	for (int j = 0; j < Np; j++) {
		PhaseExist[bId + j] = Flashcal->PhaseExist[j];
		if (PhaseExist[j]) {
			for (int i = 0; i < Nc; i++) {
				Cij[bId + j * Nc + i] = Flashcal->Cij[j * Nc + i];
			}
			S[bId + j]		= Flashcal->S[j];
			Xi[bId + j]		= Flashcal->Xi[j];
			Rho[bId + j]	= Flashcal->Rho[j];
			Mu[bId + j]		= Flashcal->Mu[j];
		}
	}
	Vf[n] = Flashcal->Vf;
	Vfp[n] = Flashcal->Vfp;
	bId = n * Nc;
	for (int i = 0; i < Nc; i++) {
		Vfi[bId + i] = Flashcal->Vfi[i];
	}
}

// relative permeability and capillary pressure
void Bulk::calKrPc()
{
	for (int n = 0; n < Num; n++) {
		int bId = n * Np;
		Flow->calKrPc(&S[bId], &Kr[bId], &Pc[bId]);
	}
}

// Rock
void Bulk::calporo()
{
	for (int n = 0; n < Num; n++) {
		double dP = P[n] - Rock_Pref;
		Rock_Poro[n] = Rock_PoroInit[n] * (1 + Rock_C1 * dP + Rock_C2 * dP * dP);
	}
}
