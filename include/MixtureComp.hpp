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
#include<vector>
#include<algorithm>


// OpenCAEPoro header files
#include "Mixture.hpp"
#include "DenseMat.hpp"
#include "OCPTable.hpp"

using namespace std;


/// Params for SSM in Phase Stability Analysis
class SSMparamSTA
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL tol2;   ///< tol*tol
};


/// Params for SSM in Phase Split
class SSMparamSP
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL tol2;   ///< tol*tol
};


/// Params for NR in Phase Stability Analysis
class NRparamSTA
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL tol2;   ///< tol*tol
};


/// Params for NR in Phase Split
class NRparamSP
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL tol2;   ///< tol*tol
};


/// Param for Solving Rachford-Rice Equations
class RRparam
{
public:
	USI maxIt;		///< Max Iteration
	OCP_DBL tol;	///< Tolerance
	OCP_DBL tol2;   ///< tol*tol
};


class EoScontrol
{
	friend class MixtureComp;

private:

	SSMparamSTA SSMsta;
	NRparamSTA  NRsta;
	SSMparamSP  SSMsp;
	NRparamSP   NRsp;
	RRparam RR;

};


class COMP
{
public:
	COMP() = default;
	COMP(const vector<string>& comp);
public:
	string name;    ///< Name of components
	OCP_DBL Pc;		///< Critical Pressure
	OCP_DBL Tc;		///< Critical Temperature
	OCP_DBL acf;	///< Acentric Factor
	OCP_DBL MW;		///< Molecular Weight
	OCP_DBL VcMW;	///< Critical Volume / MW
	OCP_DBL OmegaA; ///< Param A of Components 
	OCP_DBL OmegaB; ///< Param B of Components
	OCP_DBL Vshift; ///< shift volume
};


class MixtureComp : public Mixture
{
public:
	MixtureComp() = default;
	MixtureComp(const ParamReservoir& rs_param, const USI& i) :MixtureComp(rs_param.EoSp) {
		if (rs_param.PVTW_T.data.size() != 0) {
			PVTW.Setup(rs_param.PVTW_T.data[i]);
			data.resize(5);
			cdata.resize(5);
		}
	};
	MixtureComp(const EoSparam& param);
	void Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
		const OCP_DBL* Sjin, const OCP_DBL& Vpore,
		const OCP_DBL* Ziin) override {};
	void Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) override {};
	void Flash_Ni_Deriv(const OCP_DBL& Pin, const OCP_DBL& Tin,
		const OCP_DBL* Niin) override {};
	OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin) override {};
	OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& T,
		const OCP_DBL* Ziin) override {};
	OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) override {};
	OCP_DBL GammaPhaseG(const OCP_DBL& Pin) override {};
	OCP_DBL GammaPhaseW(const OCP_DBL& Pin) override {};
	OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
		const OCP_DBL* Ziin) override {};

private:
	// Initial properties
	USI	NC;
	USI NPmax;
	vector<OCP_DBL>   zi;
	vector<COMP>  comp;
	vector<OCP_DBL> BIP;

	EoScontrol EoSctrl;

	vector<OCP_DBL> Plist;
	vector<OCP_DBL> Tlist;
	vector<OCP_DBL> Ytlist;

private:
	OCPTable PVTW; ///< PVT table for water.
	vector<OCP_DBL> data;   ///< container used to store the results of values of
							///< interpolation of PVT tables.
	vector<OCP_DBL> cdata;  ///< container used to store the results of slopes of
							///< interpolation of PVT tables.


public:
	// Environment function
	void setPT(const OCP_DBL& p, const OCP_DBL& t) { P = p; T = t; }

private:
	// Environment variables
	OCP_DBL P; ///< Current Pressure
	OCP_DBL T; ///< Current Temperature


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
	/// Result is stored in Ztmp.
	USI CubicRoot(const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c, const bool& NTflag = false) const;
	/// test
	void PrintZtmp() { for (auto z : Ztmp) cout << z << "\t"; }

private:
	// EoS Variables
	vector<OCP_DBL> Ai;
	vector<OCP_DBL> Bi;
	vector<OCP_DBL> Aj;
	vector<OCP_DBL> Bj;
	vector<OCP_DBL> Zj;
	mutable vector<OCP_DBL> Ztmp; ///< Cubic root space,size: 3

	// PR default
	OCP_DBL delta1 = 2.41421356237;
	OCP_DBL delta2 = -0.41421356237;


public:
	// Phase Function
	// Allocate memoery for phase variables
	void AllocatePhase();
	void CalFugPhi(vector<OCP_DBL>& phiT, vector<OCP_DBL>& fugT, const vector<OCP_DBL>& xj);
	void CalFugPhiAll();
	void CalMW();
	void CalXiRho();
	USI FindMWmax();
	void x2n();  ///< x[j][i] -> n[j][i]
	void PrintX();

private:
	// Phase Variables
	USI NP;
	vector<OCP_DBL> nuj; ///< nuj[j] represents the mole fraction of jth phase
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
	bool StableNR(const USI& Id);
	void CalFugXSTA();
	void AssembleJmatSTA();
	void PhaseSplit();
	bool SplitSSM(const bool& flag);
	bool SplitSSM2(const bool& flag);
	bool SplitSSM3(const bool& flag);
	void RachfordRice2();  ///< Used when NP = 2
	void RachfordRice3(); ///< Used when NP > 2
	void UpdateXRR(); ///< Update X according to RR
	bool SplitNR();
	void CalResSP();
	void CalFugNAll();
	void PrintFugN();
	void AssembleJmatSP();
	OCP_DBL CalStepNRsp();

private:
	// Method Variables
	USI testPId; ///< Index of the testing phase in stability analysis
	vector<vector<OCP_DBL>> Kw; ///< Equlibrium Constant of Whilson
	vector<vector<OCP_DBL>> Ks;  ///< Approximation of Equilibrium Constant in SSM
	OCP_DBL Asta, Bsta, Zsta;
	vector<OCP_DBL> phiSta; ///< Fugacity coefficient used in phase stability analysis
	vector<OCP_DBL> fugSta; ///< Fugacity used in phase stability analysis
	// SSM in Stability Analysis
	vector<OCP_DBL> di;  ///< phi(id) * x(id), id is the index of testing phase
	vector<OCP_DBL> Y;   ///< x[i] / Yt
	OCP_DBL Yt;  ///< Sum Y
	// NR in Stability Analysis
	vector<OCP_DBL> resSTA;
	vector<OCP_DBL> JmatSTA; ///< d g / d Y
	vector<vector<OCP_DBL>> fugX; ///< d ln f / d X
	vector<OCP_DBL> Ax; ///< d Aj / d xkj, j is fixed
	vector<OCP_DBL> Bx; ///< d Bj / d xkj, j is fixed
	vector<OCP_DBL> Zx; ///< d Zj / d xkj, j is fixed

	// SSM in Phase Split
	vector<OCP_DBL> resRR; ///< Error in Rachford-Rice equations.
	// NR in Phase Split
	vector<OCP_DBL> resSP; ///< d G / d nij, G is Gibbs free energy: ln fij - ln fi,np
	vector<OCP_DBL> JmatSP; ///< Jacobian Matrix of (ln fij - ln fi,np) wrt. nij
	vector<vector<OCP_DBL>> fugN; ///< d ln fij / d nkj, in each subvector, ordered by k.
	vector<OCP_DBL> An; ///< d Aj / d nkj, j is fixed
	vector<OCP_DBL> Bn; ///< d Bj / d nkj, j is fixed
	vector<OCP_DBL> Zn; ///< d Zj / d nkj, j is fixed
	// for linearsolve with lapack
	vector<int> pivot; ///< used in dgesv_ in lapack
	char uplo{ 'U' };

public:
		// After Phase Equilibrium Calculation finishs, properties and some auxiliary variables
		// will be calculated.
		void AllocateOthers();
		void IdentifyPhase();
		void CalViscosity();
		void CalViscoLBC();
		void CalViscoHZYT();
		void CalFugXAll();
		void CalFugPAll();

private:
	// Phase properties and auxiliary variables
	vector<OCP_DBL> muC;    ///< Viscosity of phase
	vector<vector<OCP_DBL>> muAux; ///< Auxiliary variables for Viscosity, used to calculate Derivative

	vector<vector<OCP_DBL>> fugP; ///< d ln fij / d P
	vector<OCP_DBL> Zp; ///< d Z / d P
};


/// Return the sign of double di
OCP_DBL signD(const OCP_DBL& d);

OCP_DBL delta(const USI& i, const USI& j);

void NTcubicroot(OCP_DBL& root, const OCP_DBL& a, const OCP_DBL& b, const OCP_DBL& c);



#endif
