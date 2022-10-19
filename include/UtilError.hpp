/** \file    UtilError.hpp
 *  \brief   Logging error and warning messages
 *  \author  Ronghong Fan
 *  \date    Nov/01/2019
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2019--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __ERRORLOG_HXX__ /*-- allow multiple inclusions --*/
#define __ERRORLOG_HXX__ /**< indicate ErrorLog.hxx has been included before */

// Standard header files
#include <iomanip>
#include <iostream>
#include <sstream>

/// Print out location at (file, line) and function name
#define OCP_LOCATION                                                                   \
    "\n    --> function: " << __PRETTY_FUNCTION__                                      \
                           << "\n    --> file:     " << __FILE__ << "::" << __LINE__

/// Log error messages
//  msg: user-defined error message
//  We use do-while to allow the macro to be ended with ";"
#define OCP_MESSAGE(msg)                                                               \
    do {                                                                               \
        std::ostringstream info;                                                       \
        info << std::setprecision(16);                                                 \
        info << msg << OCP_LOCATION << '\n';                                           \
        std::cerr << info.str().c_str();                                               \
    } while (false)

/// Log warning messages
//  msg: user-defined warning message
//  We use do-while to allow the macro to be ended with ";"
#define OCP_WARNING(msg)                                                               \
    do {                                                                               \
        OCP_MESSAGE("### WARNING: " << (msg));                                         \
    } while (false)

/// Abort if critical error happens
//  msg: user-defined abort message
//  We use do-while to allow the macro to be ended with ";"
#define OCP_ABORT(msg)                                                                 \
    do {                                                                               \
        OCP_MESSAGE("### ABORT: " << (msg));                                           \
        std::abort();                                                                  \
    } while (false)

/// Assert condition and log user messages in DEBUG mode
//  cond: check condition
//  msg: user-defined error message
//  We use do-while to allow the macro to be ended with ";"
#ifndef DEBUG
#define OCP_ASSERT(cond, msg)                                                          \
    do {                                                                               \
    } while (false)
#else
#define OCP_ASSERT(cond, msg)                                                          \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            OCP_MESSAGE("### ASSERT: " << (msg) << " (" << #cond << ")");              \
            std::abort();                                                              \
        }                                                                              \
    } while (false)
#endif

/// Print Function Name
#ifndef OCPFUNCNAME
#define OCP_FUNCNAME                                                                   \
    do {                                                                               \
    } while (false)
#else
#define OCP_FUNCNAME                                                                   \
    do {                                                                               \
        std::cout << __FUNCTION__ << std::endl;                                        \
    } while (false)
#endif

#endif /* end if for __ERRORLOG_HXX__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*  Chensong Zhang      Nov/15/2021      Test DEBUG mode                      */
/*----------------------------------------------------------------------------*/