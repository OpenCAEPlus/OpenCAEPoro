#pragma once
#include "BOMixture.hpp"
#include "FlowUnit.hpp"
#include "Grid.hpp"
#include "Mixture.hpp"
#include "OpenCAEPoro_consts.hpp"
#include "ParamReservoir.hpp"
#include "Solver.hxx"
#include <iostream>
#include <vector>

using namespace std;
// Bulk contains the infomation of each bulk
// variables are ordered according to the time consuming of program

class ParamEQUIL
{
    friend class Bulk;

private:
    //! EQUIL 6 items
    OCP_DBL Dref, Pref, DOWC, PcOWC, DGOC, PcGOC;
    //! PBVD
    ReservoirTable<OCP_DBL> PBVD;
};

class Bulk
{
    friend class Connection_BB;
    friend class Well;

public:
    Bulk() = default;

    int getBulkNum() { return Num; }

    void inputParam(ParamReservoir& rs_param);
    void setup(const Grid& myGrid);

    void initSjPc_blk(int depnum);
    void initSjPc_comp(int depnum);
    void setLastStep()
    {
        lP  = P;
        lPj = Pj;
        lNi = Ni;
        lS  = S;
        Rock_lVp = Rock_Vp;
    }
    void calMaxChange();

    // Flash
    void flash_Sj();
    void flash_Ni();
    void passFlashValue(int n);

    // relative permeability and capillary pressure
    void calKrPc();
    // Rock
    void calVporo();

    // get flash -- pass it to well
    std::vector<Mixture*>& getMixture() { return Flashcal; }

    int mixMode();

    // solver
    void getSol_IMPES(vector<OCP_DBL>& u);

    // calculate FPR
    OCP_DBL calFPR();
    OCP_DBL getP(int n) { return P[n]; }
    bool   checkP();
    bool   checkNi();
    bool   checkVe(const OCP_DBL Vlim);
    void   resetP() { P = lP; }
    void   resetPj(){ Pj = lPj; }
    void   resetNi() { Ni = lNi; }
    void   resetVp() { Rock_Vp = Rock_lVp; }


    OCP_DBL getdPmax() { return dPmax; }
    OCP_DBL getdNmax() { return dNmax; }
    OCP_DBL getdSmax() { return dSmax; }
    OCP_DBL getdVmax() { return dVmax; }

private:
    int Num; // num of active bulk
    // Physical infomation
    int Np; // num of phase
    int Nc; // num of component

    OCP_DBL              T;          // temperature : Num
    std::vector<OCP_DBL> Pbub;       // buble point pressere: Num
    std::vector<OCP_DBL> P;          // pressure: Num
    std::vector<OCP_DBL> Pj;         // phase pressure: Np*Num
    std::vector<OCP_DBL> Pc;         // capillary pressure of phase: Np*Num
    std::vector<bool>   PhaseExist; // existence of phase
    std::vector<OCP_DBL> S;          // saturation of phase j
    std::vector<OCP_DBL> Rho;        // mass density of phase: Np*Num
    std::vector<OCP_DBL> Xi;         // molar density of phase: Np*Num
    std::vector<OCP_DBL> Cij;        // Nij / Nj : Np*Nc*Num
    std::vector<OCP_DBL> Ni;         // molar of ith component in bulk: NC*Num
    std::vector<OCP_DBL> Mu;         // viscosity of phase: Np*Num
    std::vector<OCP_DBL> Kr;         // relative permeability of phase: Np*Num

    std::vector<OCP_DBL> Vj;
    std::vector<OCP_DBL> Vf;  // total fluid volume
    std::vector<OCP_DBL> Vfi; // dVt / dNi
    std::vector<OCP_DBL> Vfp; // dVt / dP

    vector<int>            PhaseLabel;
    std::vector<OCP_DBL>    InitZi; // initial component for EoS : Nc - 1
    int                    PVTmode;
    std::vector<int>       PVTNUM;
    std::vector<Mixture*>  Flashcal;
    int                    SATmode;
    std::vector<int>       SATNUM;
    std::vector<FlowUnit*> Flow;

    // last STEP
    std::vector<OCP_DBL> lP;
    std::vector<OCP_DBL> lPj;
    std::vector<OCP_DBL> lNi;
    std::vector<OCP_DBL> lS;
    std::vector<OCP_DBL> Rock_lVp;

    // max change
    OCP_DBL dPmax;
    OCP_DBL dNmax;
    OCP_DBL dVmax;
    OCP_DBL dSmax;

    // Bulk rock infomation
    std::vector<OCP_DBL> Dx;          // dx
    std::vector<OCP_DBL> Dy;          // dy
    std::vector<OCP_DBL> Dz;          // dz
    std::vector<OCP_DBL> Depth;       // depth: Num
    std::vector<OCP_DBL> Ntg;         // Ntg: Num
    std::vector<OCP_DBL> Rock_VpInit; // Vgrid * ntg * poro_init
    std::vector<OCP_DBL> Rock_Vp;     // Vgrid * ntg * poro

    // std::vector<OCP_DBL>		Rock_Poro;			// current porosity
    // std::vector<OCP_DBL>		Rock_PoroInit;		// initial porosity
    OCP_DBL              Rock_Pref;
    OCP_DBL              Rock_C1;
    OCP_DBL              Rock_C2;
    std::vector<OCP_DBL> Rock_KxInit;
    std::vector<OCP_DBL> Rock_Kx;
    std::vector<OCP_DBL> Rock_KyInit;
    std::vector<OCP_DBL> Rock_Ky;
    std::vector<OCP_DBL> Rock_KzInit;
    std::vector<OCP_DBL> Rock_Kz;

    // Reservoir information
    ParamEQUIL EQUIL;
    bool       BLACKOIL;
    bool       COMPS;
    bool       Oil;
    bool       Gas;
    bool       Water;
    bool       DisGas;
};
