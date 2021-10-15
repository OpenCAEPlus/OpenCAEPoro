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

using namespace std;

/// Mixture is an abstract class, who contains all information used for flash
/// calculation including variables, functions. any properties of phases such as mass
/// density can calculated by it. it has the same data structure as the ones in bulks.
class Mixture
{
    friend class Bulk;
    friend class Well;

public:
    Mixture() = default;
    virtual ~Mixture(){};

    /// return type of mixture.
    USI GetType() const { return mixtureType; }
    /// judge if table PVDG is empty, it will only be used in black oil model.
    virtual bool IsEmpty_PVDG() const = 0;
    /// flash calculation with saturation of phases.
    virtual void Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
                          const OCP_DBL* Sjin, const OCP_DBL& Vpore,
                          const OCP_DBL* Ziin) = 0;
    /// flash calculation with moles of components.
    virtual void Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin,
                          const OCP_DBL* Niin) = 0;

    /// return molar density of phase, it's used to calculate the molar density of
    /// injection fluids in injection wells.
    virtual OCP_DBL XiPhase(const OCP_DBL& Pin, const OCP_DBL& T,
                            const OCP_DBL* Ziin) = 0;

    /// return mass density of phase.
    virtual OCP_DBL RhoPhase(const OCP_DBL& Pin, const OCP_DBL& T,
                             const OCP_DBL* Ziin) = 0;

    /// return gamma of oil phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL GammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin) = 0;
    /// return gamma of water phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL GammaPhaseW(const OCP_DBL& Pin) = 0;
    /// return gamma of gas phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL GammaPhaseG(const OCP_DBL& Pin) = 0;
    /// return gamma of hydrocarbon mixture, gamma equals to mass density times gravity
    /// factor.
    virtual OCP_DBL GammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin,
                                 const OCP_DBL* Ziin) = 0;

    /// check if Ni input from param is negative, it's used in debug mode to check
    /// Hidden trouble. actually, very small error in very short time may not make
    /// trouble.
    void CheckNi(const OCP_DBL* Ni)
    {
        bool flag = false;
        for (USI i = 0; i < numCom; i++) {
            if (Ni[i] < 0) {
                cout << "###WARNING:  ";
                ERRORcheck("Ni < 0 ");
            }

            if (Ni[i] > 0) flag = true;
        }
        if (!flag) {
            cout << "###ERROR:  ";
            ERRORcheck("All Ni <= 0 ");
            exit(0);
        }
    }

protected:
    USI mixtureType; ///< indicates the type of mixture, black oil or compositional or
                     ///< others.

    USI     numPhase; ///< num of phases.
    USI     numCom;   ///< num of components.
    OCP_DBL P;        ///< pressure when flash calculation.
    OCP_DBL T;        ///< temperature when flash calculation.

    vector<OCP_DBL> Ni;         ///< moles of component: numCom
    vector<bool>    phaseExist; ///< existence of phase: numPhase
    vector<OCP_DBL> S;          ///< saturation of phase: numPhase
    vector<OCP_DBL> rho;        ///< mass density of phase: numPhase
    vector<OCP_DBL> xi;         ///< molar density of phase: numPhase
    vector<OCP_DBL> cij; ///< Nij / Nj: numPhase*numCom, Nij is the moles of component i
                         ///< in phase j, Nj is the moles of phase j.
    vector<OCP_DBL> mu;  ///< viscosity of phase: numPhase
    vector<OCP_DBL> v;   ///< volume of phase: numPhase;

    OCP_DBL vf;  ///< volume of total fluids.
    OCP_DBL vfp; ///< dVf / dP, the derivative of volume of total fluids with respect to
                 ///< pressure.
    vector<OCP_DBL> vfi; ///< dVf / dNi: numCom  the derivative of volume of total
                         ///< fluids with respect to moles of components.
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