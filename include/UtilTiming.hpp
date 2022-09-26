/** \file    UtilTiming.hpp
 *  \brief   Elapsed wall-time and CPU-cycles declaration
 *  \author  Chensong Zhang
 *  \date    Sep/24/2019
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2019--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __TIMING_HEADER__ /*-- allow multiple inclusions --*/
#define __TIMING_HEADER__ /**< indicate timing.hxx has been included --*/

typedef unsigned long long uint64; ///< Unsigned long long int

#include <chrono> // For high-resolution CPU time
#include <iomanip>
#include <iostream>
#include <string> // For output string

// Definition of time units
const double CLOCK_USE_SEC = 5000;   ///< Show clock time in seconds
const double CLOCK_USE_MIN = 200000; ///< Show clock time in minutes

/*! \class GetWallTime
 *  \brief Get elapsed wall-time in millisecond
 *
 *  Read the current wall-time and return duration from start() to stop().
 */
class GetWallTime
{

private:
    std::chrono::steady_clock::time_point timeStamp; ///< Current CPU time

public:
#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)
    // for Window file system
    /// Start the timer
    inline void Start() { timeStamp = std::chrono::steady_clock::now(); }

    /// Stop the timer and return duration from start() in ms
    inline double Stop() const
    {
        auto elapsedTime = std::chrono::steady_clock::now() - timeStamp;
        return std::chrono::duration<double, std::milli>(elapsedTime).count();
    }
#else
    /// Start the timer
    __inline__ void Start() { timeStamp = std::chrono::steady_clock::now(); }

    /// Stop the timer and return duration from start() in ms
    __inline__ double Stop() const
    {
        auto elapsedTime = std::chrono::steady_clock::now() - timeStamp;
        return std::chrono::duration<double, std::milli>(elapsedTime).count();
    }
#endif
    /// Stop the timer and print out duration time
    void StopInfo(const std::string& info, std::ostream& out = std::cout) const;
};

#endif /*-- end if for __TIMING_HEADER__ --*/

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Sep/26/2022      Do not use RTC for x86               */
/*----------------------------------------------------------------------------*/