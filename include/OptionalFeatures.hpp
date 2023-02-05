/*! \file    OptionalFeatures.hpp
 *  \brief   OptionalFeatures class declaration
 *  \author  Shizhe Li
 *  \date    Dec/25/2022
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OPTIONALFEATURES_HEADER__
#define __OPTIONALFEATURES_HEADER__

#include "AcceleratePVT.hpp"
#include "PhasePermeability.hpp"

class OptionalFeatures
{
    friend class MixtureComp;
    friend class FlowUnit_ODGW01_Miscible;

    // For Output
    friend class Out4RPT;

public:
    void InputParam(const ParamReservoir& param)
    {
        miscible.InputParam(param.miscstr);
    };
    void ResetToLastTimeStep()
    {
        skipStaAnaly.ResetToLastTimeStep();
        miscible.ResetTolastTimeStep();
    }
    void UpdateLastTimeStep()
    {
        skipStaAnaly.UpdateLastTimeStep();
        miscible.UpdateLastTimeStep();
    }

    /////////////////////////////////////////////////////////////////////
    // Accelerate PVT
    /////////////////////////////////////////////////////////////////////

protected:
    SkipStaAnaly skipStaAnaly; ///< Skip Stability Analysis term

    /////////////////////////////////////////////////////////////////////
    // Phase Permeability Curve
    /////////////////////////////////////////////////////////////////////

protected:
    Miscible  miscible;  ///< Miscible term for Compositional Model
    ScalePcow scalePcow; ///< Scale water-oil capillary pressure term
};

#endif /* end if __OptionalFeatures_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/