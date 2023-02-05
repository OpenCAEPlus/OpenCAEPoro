/*! \file    Mixture.hpp
 *  \brief   Mixture class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURE_HEADER__
#define __MIXTURE_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoro header files
#include "OCPConst.hpp"
#include "OptionalFeatures.hpp"
#include "ParamReservoir.hpp"
#include "WellOpt.hpp"

using namespace std;

/// Mixture is an abstract class, who contains all information used for flash
/// calculation including variables, functions. any properties of phases such as mass
/// density can calculated by it. it has the same data structure as the ones in bulks.
class Mixture
{
    friend class Bulk;
    friend class Well;
    friend class AllWells;

public:
    Mixture() = default;
    virtual ~Mixture(){};
    /// Allocate memory for common variables for basic class
    void Allocate()
    {
        Ni.resize(numCom);
        phaseExist.resize(numPhase);
        vj.resize(numPhase);
        S.resize(numPhase);
        xi.resize(numPhase);
        nj.resize(numPhase);
        xij.resize(numPhase * numCom);
        rho.resize(numPhase);
        mu.resize(numPhase);
        // Derivatives
        vfi.resize(numCom);
        rhoP.resize(numPhase);
        xiP.resize(numPhase);
        muP.resize(numPhase);
        rhox.resize(numPhase * numCom);
        xix.resize(numPhase * numCom);
        mux.resize(numPhase * numCom);
        dXsdXp.resize((numCom + 1) * (numPhase + numPhase * numCom));
        // Auxiliary variable
        pSderExist.resize(numPhase);
        pVnumCom.resize(numPhase);
        // Thermal model

        // used in FIM_n
        res.resize(numPhase + numPhase * numCom +
                   1); // a precomputed value stored in last position
        // water not in hydrocarbon, hydrocarbon not in water
        // keyDer.resize((numCom + 1) * ((numPhase - 1) * (numCom - 1) + 1));
    };
    virtual void SetupOptionalFeatures(OptionalFeatures& optFeatures,
                                       const OCP_USI&    numBulk) = 0;
    /// return type of mixture.
    USI GetMixtureType() const { return mixtureType; }
    /// flash calculation with saturation of phases.
    virtual void Flash(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) = 0;
    virtual void InitFlashIMPEC(const OCP_DBL& Pin,
                                const OCP_DBL& Pbbin,
                                const OCP_DBL& Tin,
                                const OCP_DBL* Sjin,
                                const OCP_DBL& Vpore,
                                const OCP_DBL* Ziin,
                                const OCP_USI& bId)                                 = 0;
    virtual void InitFlashFIM(const OCP_DBL& Pin,
                              const OCP_DBL& Pbbin,
                              const OCP_DBL& Tin,
                              const OCP_DBL* Sjin,
                              const OCP_DBL& Vpore,
                              const OCP_DBL* Ziin,
                              const OCP_USI& bId)                                   = 0;
    virtual void InitFlashFIMn(const OCP_DBL& Pin,
                               const OCP_DBL& Pbbin,
                               const OCP_DBL& Tin,
                               const OCP_DBL* Sjin,
                               const OCP_DBL& Vpore,
                               const OCP_DBL* Ziin,
                               const OCP_USI& bId)                                  = 0;
    /// Flash calculation with moles of components.
    virtual void FlashIMPEC(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Niin,
                            const USI&     lastNP,
                            const OCP_DBL* xijin,
                            const OCP_USI& bId) = 0;
    /// Flash calculation with moles of components and Calculate the derivative
    virtual void FlashFIM(const OCP_DBL& Pin,
                          const OCP_DBL& Tin,
                          const OCP_DBL* Niin,
                          const OCP_DBL* Sjin,
                          const USI&     lastNP,
                          const OCP_DBL* xijin,
                          const OCP_USI& bId)  = 0;
    virtual void FlashFIMn(const OCP_DBL& Pin,
                           const OCP_DBL& Tin,
                           const OCP_DBL* Niin,
                           const OCP_DBL* Sjin,
                           const OCP_DBL* xijin,
                           const OCP_DBL* njin,
                           const USI*     phaseExistin,
                           const USI&     lastNP,
                           const OCP_USI& bId) = 0;

    /// return mass density of phase
    // for blackoil model: if tarPhase is gas and water, Pin and tar phase is needed
    // for compositional model: if tar phase is hydrocarbon phase, Pin, Tin, Ziin is
    // needed. if tar phase is water, only Pin is needed.
    virtual OCP_DBL XiPhase(const OCP_DBL& Pin,
                            const OCP_DBL& Tin,
                            const OCP_DBL* Ziin,
                            const USI&     tarPhase) = 0;

    /// return mass density of phase
    // for blackoil model: if tarPhase is gas and water, Pin and tar phase is needed, if
    // tarPhase is oil,then Pbb is needed, too for compositional model: if tar phase is
    // hydrocarbon phase, Pin, Tin, Ziin is needed. if tar phase is water, only Pin is
    // needed.
    virtual OCP_DBL RhoPhase(const OCP_DBL& Pin,
                             const OCP_DBL& Pbb,
                             const OCP_DBL& Tin,
                             const OCP_DBL* Ziin,
                             const USI&     tarPhase) = 0;

    // for well
    /// Setup injZi, injProdPhase and factorINJ for INJ well
    /// Setup prodPhaseWeight for PROD well
    virtual void SetupWellOpt(WellOpt&                  opt,
                              const vector<SolventINJ>& sols,
                              const OCP_DBL&            Psurf,
                              const OCP_DBL&            Tsurf) = 0;
    /// Calculate ProdWeight for PROD well
    virtual void CalProdWeight(const OCP_DBL&         Pin,
                               const OCP_DBL&         Tin,
                               const OCP_DBL*         Niin,
                               const vector<OCP_DBL>& prodPhase,
                               vector<OCP_DBL>&       prodWeight) = 0;
    /// Calculate Production rate for PROD well
    virtual void    CalProdRate(const OCP_DBL&   Pin,
                                const OCP_DBL&   Tin,
                                const OCP_DBL*   Niin,
                                vector<OCP_DBL>& prodRate)                      = 0;
    virtual OCP_DBL CalInjWellEnthalpy(const OCP_DBL& Tin, const OCP_DBL* Ziin) = 0;

    /// check if Ni input from param is negative, it's used in debug mode to check
    /// Hidden trouble. actually, very small error in very short time may not make
    /// trouble.
    void CheckNi(const OCP_DBL* Ni)
    {
        OCP_BOOL flag = OCP_TRUE;
        for (USI i = 0; i < numCom; i++) {
            if (Ni[i] < 0) {
                cout << "Ni[" << i << "] = " << Ni[i] << endl;
                flag = OCP_FALSE;
                break; // skip the rest checks
            }
        }
        if (!flag) OCP_ABORT("Ni is negative!");
    }

    virtual OCP_DBL GetErrorPEC()           = 0;
    virtual void    OutMixtureIters() const = 0;

public:
    const OCP_DBL&  GetNt() const { return Nt; }
    const OCP_DBL&  GetNi(const USI& i) const { return Ni[i]; }
    const OCP_DBL&  GetVf() const { return vf; }
    const OCP_BOOL& GetPhaseExist(const USI& j) const { return phaseExist[j]; }
    const OCP_DBL&  GetS(const USI& j) const { return S[j]; }
    const OCP_DBL&  GetVj(const USI& j) const { return vj[j]; }
    const OCP_DBL&  GetNj(const USI& j) const { return nj[j]; }
    const OCP_DBL&  GetXij(const USI& j, const USI& i) const
    {
        return xij[j * numCom + i];
    }
    const OCP_DBL& GetRho(const USI& j) const { return rho[j]; }
    const OCP_DBL& GetXi(const USI& j) const { return xi[j]; }
    const OCP_DBL& GetMu(const USI& j) const { return mu[j]; }
    const OCP_DBL& GetVfP() const { return vfP; }
    const OCP_DBL& GetVfT() const { return vfT; }
    const OCP_DBL& GetVfi(const USI& i) const { return vfi[i]; }
    const OCP_DBL& GetRhoP(const USI& j) const { return rhoP[j]; }
    const OCP_DBL& GetRhoT(const USI& j) const { return rhoT[j]; }
    const OCP_DBL& GetXiP(const USI& j) const { return xiP[j]; }
    const OCP_DBL& GetXiT(const USI& j) const { return xiT[j]; }
    const OCP_DBL& GetMuP(const USI& j) const { return muP[j]; }
    const OCP_DBL& GetMuT(const USI& j) const { return muT[j]; }
    const OCP_DBL& GetRhoX(const USI& j, const USI& i) const
    {
        return rhox[j * numCom + i];
    }
    const OCP_DBL& GetXiX(const USI& j, const USI& i) const
    {
        return xix[j * numCom + i];
    }
    const OCP_DBL& GetMuX(const USI& j, const USI& i) const
    {
        return mux[j * numCom + i];
    }
    const OCP_BOOL&        GetPSderExist(const USI& j) const { return pSderExist[j]; }
    const USI&             GetPVnumCom(const USI& j) const { return pVnumCom[j]; }
    const vector<OCP_DBL>& GetDXsDXp() const { return dXsdXp; }
    const vector<OCP_DBL>& GetRes() const { return res; }
    const OCP_DBL          GetResPc() const { return resPc; }
    const OCP_DBL          GetUf() const { return Uf; }
    const OCP_DBL          GetUfP() const { return UfP; }
    const OCP_DBL          GetUfT() const { return UfT; }
    const OCP_DBL          GetUfi(const USI& i) const { return Ufi[i]; }
    const OCP_DBL          GetH(const USI& j) const { return H[j]; }
    const OCP_DBL          GetHT(const USI& j) const { return HT[j]; }
    const OCP_DBL&         GetHx(const USI& j, const USI& i) const
    {
        return Hx[j * numCom + i];
    }

protected:
    void SetBulkId(const OCP_USI& n) { bulkId = n; }

protected:
    USI mixtureType;  ///< indicates the type of mixture, black oil or compositional or
                      ///< others.
    OCP_USI bulkId;   ///< index of current bulk
    USI     numPhase; ///< num of phases.
    USI     numCom;   ///< num of components.
    OCP_DBL P;        ///< pressure when flash calculation.
    OCP_DBL T;        ///< temperature when flash calculation.

    OCP_DBL          vf;         ///< volume of total fluids.
    OCP_DBL          Nt;         ///< Total moles of Components.
    vector<OCP_DBL>  Ni;         ///< moles of component: numCom
    vector<OCP_BOOL> phaseExist; ///< existence of phase: numPhase
    vector<OCP_DBL>  S;          ///< saturation of phase: numPhase
    vector<OCP_DBL>  vj;         ///< volume of phase: numPhase;
    vector<OCP_DBL>  nj;         ///< mole number of phase j
    vector<OCP_DBL>  xij;        ///< Nij / nj: numPhase*numCom
    vector<OCP_DBL>  rho;        ///< mass density of phase: numPhase
    vector<OCP_DBL>  xi;         ///< molar density of phase: numPhase
    vector<OCP_DBL>  mu;         ///< viscosity of phase: numPhase

    // Derivatives
    OCP_DBL vfP; ///< dVf / dP, the derivative of volume of total fluids with respect to
                 ///< pressure.
    OCP_DBL         vfT; ///< d vf  / dT
    vector<OCP_DBL> vfi; ///< dVf / dNi: numCom  the derivative of volume of total
                         ///< fluids with respect to moles of components.

    vector<OCP_DBL> rhoP; ///< d rho / dP: numphase
    vector<OCP_DBL> rhoT; ///< d rho j / dT: numPhase
    vector<OCP_DBL> rhox; ///< d rho[j] / d x[i][j]: numphase * numCom
    vector<OCP_DBL> xiP;  ///< d xi / dP: numphase
    vector<OCP_DBL> xiT;  ///< d xi j / dT: numPhase
    vector<OCP_DBL> xix;  ///< d xi[j] / d x[i][j]: numphase * numCom
    vector<OCP_DBL> muP;  ///< d mu / dP: numPhase
    vector<OCP_DBL> muT;  ///< d mu j  / dT: numPhase
    vector<OCP_DBL> mux;  ///< d mu[j] / d x[i][j]: numphase * numCom

    vector<OCP_DBL> dXsdXp; ///< derivatives of second variables wrt. primary variables

    OCP_DBL         Uf;  ///< Internal energy of fluid
    OCP_DBL         UfP; ///< dUf / dP
    OCP_DBL         UfT; ///< dUf / dT
    vector<OCP_DBL> Ufi; ///< dUf / dNi
    vector<OCP_DBL> H;   ///< Enthalpy
    vector<OCP_DBL> HT;  ///< d Hj / d T
    vector<OCP_DBL> Hx;  ///< d Hj / d xij

    // Auxiliary variable for dSec_dPr
    vector<OCP_BOOL> pSderExist; ///< Existence of  derivative of phase saturation
    vector<USI>      pVnumCom;   ///< num of variable components in the phase

    vector<OCP_DBL> res;    ///< residual of a set of equations
    OCP_DBL         resPc;  ///< a precalculated value
    vector<OCP_DBL> keyDer; ///< d (xij*xi/mu) / dP or dNk
};

#endif /* end if __MIXTURE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/