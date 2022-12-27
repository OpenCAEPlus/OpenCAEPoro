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

public:
    void InputParam() {};
    void ResetToLastTimeStep() {
        skipStaAnaly.ResetToLastTimeStep();
    }
    void UpdateLastTimeStep() {
        skipStaAnaly.UpdateLastTimeStep();
    }


    /////////////////////////////////////////////////////////////////////
    // Accelerate PVT
    /////////////////////////////////////////////////////////////////////

protected:
    SkipStaAnaly skipStaAnaly;


    /////////////////////////////////////////////////////////////////////
    // Phase Permeability Curve
    /////////////////////////////////////////////////////////////////////

protected:


};





#endif /* end if __OptionalFeatures_HEADER__ */


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/25/2022      Create file                          */
/*----------------------------------------------------------------------------*/