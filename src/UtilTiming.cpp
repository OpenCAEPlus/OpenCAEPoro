/** \file    UtilTiming.cpp
 *  \brief   Elapsed wall-time and CPU-cycles definition
 *  \author  Chensong Zhang
 *  \date    Feb/22/2020
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2020--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "UtilTiming.hpp"

/// Stop timer and print out the duration time in ms, s, or m
void GetWallTime::StopInfo(const std::string& info, std::ostream& out) const
{
    const double duration = Stop();
    if (duration < CLOCK_USE_SEC) {
        std::cout << info << " costs " << std::fixed << std::setprecision(3) << duration
                  << "ms" << std::endl;
    } else if (duration < CLOCK_USE_MIN) {
        std::cout << info << " costs " << std::fixed << std::setprecision(3)
                  << duration / 1000.0 << "s" << std::endl;
    } else {
        std::cout << info << " costs " << std::fixed << std::setprecision(3)
                  << duration / 60000.0 << "m" << std::endl;
    }
}

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/