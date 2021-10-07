/*! \file    Mixture.hpp
 *  \brief   Mixture class declaration
 *  \author  Shizhe Li
 *  \date    Oct/07/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __MIXTURE_HEADER__
#define __MIXTURE_HEADER__



#include "OpenCAEPoro_consts.hpp"
#include <iostream>
#include <vector>

using namespace std;

/// Mixture is an abstract class, who contains all information used for flash calculation including variables, functions.
/// it has the same data structure as the ones in bulks.
class Mixture
{
    friend class Bulk;
    friend class Well;

public:
    Mixture() = default;
    virtual ~Mixture(){};

    /// return type of mixture.
    USI getType() const { return MixtureType; }
    /// judge if table PVDG is empty, it will only be used in black oil model.
    virtual bool empty_PVDG() const = 0;
    /// flash calculation with saturation of phases.
    virtual void Flash_Sj(const OCP_DBL& Pin, const OCP_DBL& Pbbin, const OCP_DBL& Tin,
                          const OCP_DBL* Sjin, const OCP_DBL& Vpore, const OCP_DBL* Ziin)   = 0;
    /// flash calculation with moles of components.
    virtual void Flash_Ni(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Niin) = 0;

    /// return molar density of phase, it's used to calculate the molar density of injection fluids in injection wells.
    virtual OCP_DBL xiPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin) = 0;

    /// return mass density of phase.
    virtual OCP_DBL rhoPhase(const OCP_DBL& Pin, const OCP_DBL& T, const OCP_DBL* Ziin) = 0;

    /// return gamma of oil phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL gammaPhaseO(const OCP_DBL& Pin, const OCP_DBL& Pbbin)              = 0;
    /// return gamma of water phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL gammaPhaseW(const OCP_DBL& Pin)                            = 0;
    /// return gamma of gas phase, gamma equals to mass density times gravity factor.
    virtual OCP_DBL gammaPhaseG(const OCP_DBL& Pin)                            = 0;
    /// return gamma of hydrocarbon mixture, gamma equals to mass density times gravity factor.
    virtual OCP_DBL gammaPhaseOG(const OCP_DBL& Pin, const OCP_DBL& Tin, const OCP_DBL* Ziin) = 0;

    /// check if Ni input from param is negative, it's used in debug mode to check Hidden trouble.
    /// actually, very small error in very short time may not make trouble.
    void checkNi(const OCP_DBL* Ni)
    {
        bool flag = false;
        for (USI i = 0; i < Nc; i++) {
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
    USI         MixtureType;        ///< indicates the type of mixture, black oil or compositional or others. 

    USI    Np;                      ///< num of phases.
    USI    Nc;                      ///< num of components.
    OCP_DBL P;                      ///< pressure when flash calculation.
    OCP_DBL T;                      ///< temperature when flash calculation.

    vector<OCP_DBL> Ni;        ///< moles of component: Nc
    vector<bool>   PhaseExist; ///< existence of phase: Np
    vector<OCP_DBL> S;         ///< saturation of phase: Np
    vector<OCP_DBL> Rho;       ///< mass density of phase: Np
    vector<OCP_DBL> Xi;        ///< molar density of phase: Np
    vector<OCP_DBL> Cij;       ///< Nij / Nj: Np*Nc, Nij is the moles of component i in phase j, Nj is the moles of phase j.
    vector<OCP_DBL> Mu;        ///< viscosity of phase: Np
    vector<OCP_DBL> V;         ///< volume of phase: Np;

    OCP_DBL              Vf;   ///< volume of total fluids.
    OCP_DBL              Vfp;  ///< dVf / dP, the derivative of volume of total fluids with respect to pressure.
    vector<OCP_DBL>      Vfi;  ///< dVf / dNi: Nc  the derivative of volume of total fluids with respect to moles of components.
};


#endif