/*! \file    Doxygen.hpp
 *  \brief   Main page for Doxygen documentation
 *  \author  Chensong Zhang
 *  \date    Oct/15/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoro team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __DOXYGEN_HXX__ /*-- allow multiple inclusions --*/
#define __DOXYGEN_HXX__ /**< indicate Doxygen.hxx has been included before */

/** \mainpage Introduction
 *
 * OpenCAEPoro is part of the OpenCAEPlus project written in C++. It focuses on
 * multicomponent multiphase flow simulation in porous media. Check out our
 * <a href="https://github.com/OpenCAEPlus/OpenCAEPoro">Github repository</a>
 * for the source codes.
 *
 * For design goals and user manual, please refer to
 * <a href="https://porous.opencaeplus.com/">the OCP Website</a>.
 *
 * > This software distributed in the hope that it will be useful, but WITHOUT ANY
 * > WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * > PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details:
 * > https://www.gnu.org/licenses/lgpl-3.0.html
 *
 */

/**
 * \page download How to obtain OpenCAEPoro
 *
 * The code is freely available on GitHub https://github.com/OpenCAEPlus/OpenCAEPoro.
 *
 */

/**
 * \page build Building and Installation
 *
 * This is a simple instruction on building and testing. There is a top level cmake
 * for configuration and building of the FASPxx shared library and the test programs
 * suite. You can use a cmake-style way to compile the package; see https://cmake.org
 * on how to use cmake for your own operating system. To compile, you alos need a C++
 * compiler.
 *
 * ```bash
 *  $ mkdir Build; cd Build; cmake ..
 *  $ make
 * ```
 *
 * You may config with different cmake options; for example:
 *
 * (1) Build in Debug configuration:
 *
 * ```bash
 *  $ cmake -DCMAKE_BUILD_TYPE=Debug ..
 * ```
 *
 * (2) Build with verbose on (with building messages):
 * 
 * ```bash
 *  $ cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..
 * ```
 *
 * (3) Build with UMFPACK support (requires setting SUITESPARSE_DIR variable):
 * 
 * ```bash
 *  $ cmake -DUSE_UMFPACK=ON ..
 * ```
 * 
 * You can also use alternative direct solver packages like MUMPS, SUPERLU, or MKL 
 * PARDISO instead of UMFPACK. They can be compiled in a similar way. 
 */

/**
 * \page developers Developers
 *
 * ## Developers (in alphabetic order):
 *
 * - Chen, Xiaoxing
 *
 * - Li, Shizhe
 *
 * - Qiao, Changhe
 *
 * - Zhang, Chensong
 *
 * - Zhao, Li
 *
 * ## Project coordinator:
 *
 * - Zhang, Chensong
 * 
 */

/**
 * \page doxygen_comment Doxygen
 *
 * We use Doxygen as our automatically documentation generator which will make our
 * future maintainance minimized. You can obtain the software (Windows, Linux and
 * OS X) as well as its manual on the official website
 *
 * https://doxygen.nl/
 *
 * For an ordinary user, Doxygen is completely trivial to use. We only need to use
 * some special marker in the usual comment as we put in c-files.
 *
 */

#endif /* end if for __DOXYGEN_HXX__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Chensong Zhang      Oct/15/2021      Create file                          */
/*  Chensong Zhang      Sep/22/2022      Add new website address              */
/*----------------------------------------------------------------------------*/