/*! \file    MixtureComp.hpp
 *  \brief   MixtureComp class declaration
 *  \author  Shizhe Li
 *  \date    Nov/30/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURECOMP_HEADER__
#define __MIXTURECOMP_HEADER__

// Standard header files
#include <algorithm>
#include <math.h>
#include <vector>

// OpenCAEPoro header files
#include "DenseMat.hpp"
#include "Mixture.hpp"
#include "OCPTable.hpp"

using namespace std;

/// Params for SSM in Phase Stability Analysis
class SSMparamSTA
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL Ktol{ 1E-4 };   ///< tolerace^2 for K
	OCP_DBL dYtol{ 1E-6 };
	OCP_DBL eYt{ 1E-8 };    ///< if Yt > 1 + eYt, than single phase is unstable
	OCP_DBL tol2;   ///< tol*tol
	// test
    USI curIt; ///< current Iterations
    OCP_DBL curSk;
};

/// Params for NR in Phase Stability Analysis
class NRparamSTA
{
public:
    USI     maxIt; ///< Max Iteration
    OCP_DBL tol;   ///< Tolerance
    OCP_DBL tol2;  ///< tol*tol
    // test
    USI curIt;     ///< current Iters
};


/// Params for SSM in Phase Split
class SSMparamSP
{
public:
    USI     maxIt;   ///< Max Iteration
    OCP_DBL tol;     ///< Tolerance
    OCP_DBL tol2;    ///< tol*tol
    OCP_DBL realTol; ///< Real tol
    bool    conflag; ///< convergence flag, if converges, conflag = true
    // test
    USI curIt;     ///< current Iters
};


/// Params for NR in Phase Split
class NRparamSP
{
public:
    USI     maxIt;   ///< Max Iteration
    OCP_DBL tol;     ///< Tolerance
    OCP_DBL tol2;    ///< tol*tol
    OCP_DBL realTol; ///< Real tol
    bool    conflag; ///< convergence flag, if converges, conflag = true
    // test
    USI curIt;     ///< current Iters
};

/// Param for Solving Rachford-Rice Equations
class RRparam
{
public:
    USI     maxIt; ///< Max Iteration
    OCP_DBL tol;   ///< Tolerance
    OCP_DBL tol2;  ///< tol*tol
    // test
    USI curIt;     ///< current Iters
};

class EoScontrol
{
    friend class MixtureComp;

private:
    SSMparamSTA SSMsta;
    NRparamSTA  NRsta;
    SSMparamSP  SSMsp;
    NRparamSP   NRsp;
    RRparam     RR;
};

class COMP
{
public:
    COMP() = default;
    COMP(const vector<string>& comp);

public:
    string  name;   ///< Name of components
    OCP_DBL Pc;     ///< Critical Pressure
    OCP_DBL Tc;     ///< Critical Temperature
    OCP_DBL acf;    ///< Acentric Factor
    OCP_DBL MW;     ///< Molecular Weight
    OCP_DBL VcMW;   ///< Critical Volume / MW
    OCP_DBL Vc;     ///< VcMW * MW
    OCP_DBL OmegaA; ///< Param A of Components
    OCP_DBL OmegaB; ///< Param B of Components
    OCP_DBL Vshift; ///< shift volume
};

class MixtureComp : public Mixture
{

public:

    OCP_ULL GetSSMSTAiters() override { return SSMSTAiters; }
    OCP_ULL GetNRSTAiters() override { return NRSTAiters; }
    OCP_ULL GetSSMSPiters() override { return SSMSPiters; }
    OCP_ULL GetNRSPiters() override { return NRSPiters; }

private:
	// for dubug
	OCP_ULL SSMSTAiters{ 0 };
	OCP_ULL NRSTAiters{ 0 };
	OCP_ULL SSMSPiters{ 0 };
	OCP_ULL NRSPiters{ 0 };

public:
	MixtureComp() = default;
	MixtureComp(const ParamReservoir& rs_param, const USI& i) :MixtureComp(rs_param.EoSp, i) {
		mixtureType = EOS_PVTW;
		if (rs_param.PVTW_T.data.size() != 0) {
			PVTW.Setup(rs_param.PVTW_T.data[i]);
			if (rs_param.gravity.activity)
				std_RhoW = RHOW_STD * rs_param.gravity.data[1];
			if (rs_param.density.activity)
				std_RhoW = RHOW_STD;
			std_GammaW = GRAVITY_FACTOR * std_RhoW;
			data.resize(5);
			cdata.resize(5);
		}
	};
	MixtureComp(const EoSparam& param, const USI& i);
	void InitFlash(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
		const OCP_DBL* Sjin, const OCP_DBL& Vpore,
		const OCP_DBL* Ziin) override;
	// ftype = 0, flash from single phase
	// ftype = 1, skip phase stablity analysis and num of phase = 1
	// ftype = 1, skip phase stablity analysis and num of phase = 2
	void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
	void CalFlash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin);
    void FlashDeriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
        const OCP_DBL* Niin, const USI& ftype, const USI& lastNP,
        const OCP_DBL* lastKs) override;
	OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) override;
	OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& Tin,
		const OCP_DBL* Ziin) override;
	OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override { OCP_ABORT("Should not be used in Compositional mode!"); };
	OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override { OCP_ABORT("Should not be used in Compositional mode!"); };
	OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override;
	OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
		const OCP_DBL* Ziin) override ;
	void setPT(const OCP_DBL& p, const OCP_DBL& t) { P = p; T = t; }
	void setZi(const OCP_DBL* Ziin) { Dcopy(NC, &zi[0], Ziin); }
	void setZi() { for (USI i = 0; i < NC; i++) zi[i] = Ni[i] / Nh; }
	void setNi(const OCP_DBL* Niin) { Dcopy(numCom, &Ni[0], Niin); }
	void CallId();
    USI GetFtype() override { return ftype; }
    void CalSurfaceTension();
    OCP_DBL GetSurTen() override { return surTen; }

private:

	// Basic Components Informations
	vector<string> Cname; ///< Name of hydrocarbon components
	vector<OCP_DBL> Tc; ///< Critical temperature of hydrocarbon components
	vector<OCP_DBL> Pc; ///< Critical pressure of hydrocarbon components
	vector<OCP_DBL> Vc; ///< Critical volume of hydrocarbon components
	vector<OCP_DBL> MWC; ///< Molecular Weight of hydrocarbon components
	vector<OCP_DBL> Acf; ///< Acentric factor of hydrocarbon components
	vector<OCP_DBL> OmegaA; ///< OMEGA_A of hydrocarbon components
	vector<OCP_DBL> OmegaB; ///< OMEGA_B of hydrocarbon components
	vector<OCP_DBL> Vshift; ///< Volume shift of hydrocarbon components
	vector<OCP_DBL> Parachor; ///< PARACHOR of hydrocarbon components
	vector<OCP_DBL> Zc; ///< Critical Z-factor of hydrocarbon components
	bool ParachorAct;
	// for viscosity calculation
	vector<OCP_DBL> Vcvis; ///< Critical volume used for viscosity calculations only.
	vector<OCP_DBL> Zcvis; ///< Critical Z-factor used for viscosity calculations only.
	vector<OCP_DBL> LBCcoef; ///< LBC coefficients for viscosity calculation
    vector<OCP_DBL> BIC; ///< Binary interaction between hydrocarbon components

    // Model information
    bool miscible; ///< Miscible treatment of hydrocarbons, used in compositional Model.
    OCP_DBL surTen; ///< Surface tension

	// Initial properties
	USI	NC; ///< num of hydrocarbon components
	USI NPmax; ///< num of hydrocarbon phase
	OCP_DBL P; ///< Current Pressure
	OCP_DBL T; ///< Current Temperature
	vector<OCP_DBL>   zi; ///< mole fraction of hydrocarbon components
	vector<COMP>  comp; ///< properties of hydrocarbon components
	USI lId; ///< index of lightest components
	EoScontrol EoSctrl; ///< method params for solving phase equilibrium
	USI ftype{ 0 };
    USI tmpFtype;

	vector<OCP_DBL> Plist;
	vector<OCP_DBL> Tlist;
	vector<OCP_DBL> Ytlist;

private:
    OCPTable        PVTW;       ///< PVT table for water.
    OCP_DBL         std_RhoW;   ///< mass density of water phase in standard condition.
    OCP_DBL         std_GammaW; ///< std_RhoW * gravity factor.
    vector<OCP_DBL> data;       ///< container used to store the results of values of
                                ///< interpolation of PVT tables.
    vector<OCP_DBL> cdata;      ///< container used to store the results of slopes of
                                ///< interpolation of PVT tables.

public:
    // EoS Function
    // Allocate memory for EoS
    void AllocateEoS();
    // Setup and Solve EoS for specified A,B
    void SolEoS(OCP_DBL& ZjT, const OCP_DBL& AjT, const OCP_DBL& BjT) const;
    // Calculate Ai and Bi
    void CalAiBi();
    // Calculate Aj and Bj with specified xj
    void CalAjBj(OCP_DBL& AjT, OCP_DBL& BjT, const vector<OCP_DBL>& xj) const;
    void CalAjBj(OCP_DBL& AjT, OCP_DBL& BjT, const OCP_DBL* xj) const;
    /// Result is stored in Ztmp.
    USI CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c,
                  const bool& NTflag = false) const;
    /// test
    void PrintZtmp()
    {
        for (auto z : Ztmp) cout << z << "\t";
    }

private:
    // EoS Variables
    vector<OCP_DBL>         Ai;
    vector<OCP_DBL>         Bi;
    vector<OCP_DBL>         Aj;
    vector<OCP_DBL>         Bj;
    vector<OCP_DBL>         Zj;
    mutable vector<OCP_DBL> Ztmp; ///< Cubic root space,size: 3

    // PR default
    OCP_DBL delta1 = 2.41421356237;
    OCP_DBL delta2 = -0.41421356237;
    OCP_DBL delta1P2;
    OCP_DBL delta1M2;
    OCP_DBL delta1T2;

public:
    // Phase Function
    // Allocate memoery for phase variables
    void AllocatePhase();
    void CalFugPhi(vector<OCP_DBL>& phiT, vector<OCP_DBL>& fugT,
                   const vector<OCP_DBL>& xj);
    void CalFugPhi(OCP_DBL* phiT, OCP_DBL* fugT, const OCP_DBL* xj);
	void CalFugPhi(OCP_DBL* fugT, const OCP_DBL* xj);
	void CalFugPhiAll();
	void CalMW();
	void CalVfXiRho();
	void CalSaturation();
	USI FindMWmax();
	void x2n();  ///< x[j][i] -> n[j][i]
	void PrintX();

private:
	// Phase Variables
	USI lNP{ 0 };  ///< last num of hydrocarbon phase
	USI NP;   ///< current num of hydrocarbon phase
    OCP_DBL Nh; ///< total moles of components exclude water
	vector<OCP_DBL> vC; ///< vC represents the volume of phase
	vector<OCP_DBL> nu; ///< nu[j] represents the mole fraction of jth phase in flash calculation
	vector<vector<OCP_DBL>> x;   ///< x[j][i] represents the mole fraction of ith comp in jth phase
	vector<vector<OCP_DBL>> phi; ///< phi[j][i] represents the fugacity coefficient of ith comp in jth phase
	vector<vector<OCP_DBL>> fug; ///< fug[j][i] represents the fugacity of ith comp in jth phase
	vector<vector<OCP_DBL>> n; ///< n[j][i] represents the moles of ith comp in jth phase
	vector<vector<OCP_DBL>> ln; ///< last n in NR iterations.
	vector<OCP_DBL> xiC; ///< Molar density of phase
	vector<OCP_DBL> rhoC; ///< Mass density of phase;
	vector<OCP_DBL> MW; ///< Molecular Weight
	vector<USI> phaseLabel; ///< Label of phase

public:
	// Method Function
	// Allocate memoery for Method variables
	void AllocateMethod();
	void PhaseEquilibrium();
	void CalKwilson();
	bool PhaseStable();
	bool StableSSM(const USI& Id);
	bool StableSSM01(const USI& Id);
	bool StableNR(const USI& Id);
	void CalFugXSTA(); ///< Calculate d ln(Fug) / dx for Y
	void AssembleJmatSTA();
	bool CheckSplit();
	void PhaseSplit();
	void SplitSSM(const bool& flag);
	void SplitSSM2(const bool& flag);
	void SplitSSM3(const bool& flag);
	void RachfordRice2();  ///< Used when NP = 2
	void RachfordRice3(); ///< Used when NP > 2
	void UpdateXRR(); ///< Update X according to RR
	void SplitNR();
	void CalResSP();
	void CalFugNAll(const bool& Znflag = true);
	void PrintFugN();
	void AssembleJmatSP();
    /// Calculate d ln phi[i][j] / d n[k][j]
    void    CalPhiNSTA();
    void    AssembleSkipMatSTA();
	OCP_DBL CalStepNRsp();

    
    OCP_SIN GetMinEigenSkip() override { return eigenSkip[0]; }
    bool GetFlagSkip() override { return flagSkip; }

private:
    // Method Variables
    USI testPId;                ///< Index of the testing phase in stability analysis
    vector<vector<OCP_DBL>> Kw; ///< Equlibrium Constant of Whilson
    vector<vector<OCP_DBL>> Ks; ///< Approximation of Equilibrium Constant in SSM
    vector<OCP_DBL> lKs; ///< last Ks
    OCP_DBL                 Asta, Bsta, Zsta;
    vector<OCP_DBL> phiSta; ///< Fugacity coefficient used in phase stability analysis
    vector<OCP_DBL> fugSta; ///< Fugacity used in phase stability analysis
    // SSM in Stability Analysis
    vector<OCP_DBL> Y;  ///< x[i] / Yt
    OCP_DBL         Yt; ///< Sum Y
    vector<OCP_DBL> di; ///< phi(id) * x(id), id is the index of testing phase
    // NR in Stability Analysis
    vector<OCP_DBL>         resSTA;
    vector<OCP_DBL>         JmatSTA; ///< d g / d Y
    vector<vector<OCP_DBL>> fugX;    ///< d ln f / d X
    vector<OCP_DBL>         Ax;      ///< d Aj / d xkj, j is fixed
    vector<OCP_DBL>         Bx;      ///< d Bj / d xkj, j is fixed
    vector<OCP_DBL>         Zx;      ///< d Zj / d xkj, j is fixed
    // Skip Stability Analysis
    bool flagSkip; ///< if check skipping Stability Analysis
    vector<OCP_DBL> phiN;    ///< d ln phi[i][j] / d n[k][j]
    vector<OCP_SIN> skipMatSTA; ///< matrix for skipping Stability Analysis
    vector<OCP_SIN> eigenSkip; ///< eigen values of matrix for skipping Skip Stability Analysis
    vector<OCP_SIN> eigenWork; ///< work space for computing eigenvalues with ssyevd_
    OCP_INT         leigenWork; ///< length of eigenwork

    // SSM in Phase Split
    vector<OCP_DBL> tmpRR; ///< temp variables for solving Rachford-Rice equations.
    vector<OCP_DBL> resRR; ///< Error in Rachford-Rice equations.
    // NR in Phase Split
    vector<OCP_DBL> resSP;  ///< d G / d nij, G is Gibbs free energy: ln fij - ln fi,np
    vector<OCP_DBL> JmatSP; ///< Jacobian Matrix of (ln fij - ln fi,np) wrt. nij
    vector<vector<OCP_DBL>>
                    fugN;       ///< d ln fij / d nkj, in each subvector, ordered by k.
    vector<OCP_DBL> An;         ///< d Aj / d nkj, j is fixed
    vector<OCP_DBL> Bn;         ///< d Bj / d nkj, j is fixed
    vector<vector<OCP_DBL>> Zn; ///< d Zj / d nkj
    // for linearsolve with lapack 
    vector<OCP_INT> pivot; ///< used in dgesv_ in lapack
    vector<OCP_DBL> JmatWork; ///< work space for Jmat in STA and SP
    OCP_INT         lJmatWork; ///< length of JmatWork
    char        uplo{'U'};

public:
    // After Phase Equilibrium Calculation finishs, properties and some auxiliary
    // variables will be calculated.
    void AllocateOthers();
    void IdentifyPhase();
    /// Copy the basic properties from MixtureComp to Mixture
    void CopyPhase();
    void CalViscosity();
    void CalViscoLBC();
    void CalViscoHZYT();
    void CalFugXAll();
    void CalFugPAll(const bool& Zpflag = true);


    void CalXiPNX_partial();
    void CalRhoPX_partial();
    void CalMuPX_partial();
    void CalMuPXLBC_partial();
    void CalXiRhoMuPN_pfullx();

    void CalRhoPNX_full();

    void CalXiPNX_full01();
    void CalRhoPNX_full01();  
    void CalMuPX_full01();
    void CalMuPXLBC_full01();
    void CalVfiVfp_full01();
    void AssembleMatVfiVfp_full01();
    void AssembleRhsVfiVfp_full01();
    void CaldXsdXp01();
    void CaldXsdXpAPI01();

    void CalXiPNX_full02();
    void CalVfiVfp_full02();
    void AssembleMatVfiVfp_full02();
    void AssembleRhsVfiVfp_full02();
    void CaldXsdXpAPI02();

    void CalVjpVfpVfn_partial();
    void CalXiPn_partial();
    void CalRhoPn_partial();
    void CalMuPn_partial();
    void CalMuPnLBC_partial();
    void CalXiRhoMuPN_pfulln();

    void CaldXsdXpAPI03();
    void CaldXsdXp03();
    void CalVfiVfp_full03();

    void CalKeyDerx();
    void CalKeyDern();

private:
    // Phase properties and auxiliary variables
    vector<OCP_DBL> muC; ///< Viscosity of phase
    vector<vector<OCP_DBL>>
        muAux; ///< Auxiliary variables for Viscosity, used to calculate Derivative
    vector<OCP_DBL>
        muAux1I; ///< Auxiliary variables for Viscosity, used to calculate Derivative
    vector<OCP_DBL>         sqrtMWi;
    vector<vector<OCP_DBL>> fugP; ///< d ln fij / d P
    vector<OCP_DBL>         Zp;   ///< d Z / d P

    vector<OCP_DBL> JmatTmp; ///< Temp Mat for transpose of a matrix
    vector<OCP_DBL> JmatDer; ///< Used to store Jacobian Mat for calculating derivates
    /// rhs or d nij / d Nk, d nij / dP in calVtiVtp
    /// rhs or dXs / dXp in Cal dXsdXp
    vector<OCP_DBL> rhsDer;
    vector<OCP_DBL> xixC; ///< d xi / d xij
    vector<OCP_DBL> xiPC; ///< d xi / d P
    vector<OCP_DBL> xiNC; ///< d xi / d Nk
};

/// Return the sign of double di
OCP_DBL signD(const OCP_DBL& d);

OCP_DBL delta(const USI& i, const USI& j);

void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c);

#endif
