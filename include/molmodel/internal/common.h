#ifndef SimTK_MOLMODEL_COMMON_H_
#define SimTK_MOLMODEL_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * Every Molmodel header and source file should include this header before
 * any other Molmodel header.
 */

#include "SimTKcommon.h"

#include <cassert>
#include <vector>
#include <limits>

// Shared libraries are messy in Visual Studio. We have to distinguish three
// cases:
//   (1) this header is being used to build the molmodel shared library (dllexport)
//   (2) this header is being used by a *client* of the molmodel shared
//       library (dllimport)
//   (3) we are building the molmodel static library, or the client is
//       being compiled with the expectation of linking with the
//       molmodel static library (nothing special needed)
// In the CMake script for building this library, we define one of the symbols
//     SimTK_MOLMODEL_BUILDING_{SHARED|STATIC}_LIBRARY
// Client code normally has no special symbol defined, in which case we'll
// assume it wants to use the shared library. However, if the client defines
// the symbol SimTK_USE_STATIC_LIBRARIES we'll suppress the dllimport so
// that the client code can be linked with static libraries. Note that
// the client symbol is not library dependent, while the library symbols
// affect only the molmodel library, meaning that other libraries can
// be clients of this one. However, we are assuming all-static or all-shared.

#ifdef _WIN32

    // avoid warning about use of non-standard "extern template" constructs
    // #ifdef _MSC_VER
    // #pragma warning(disable:4231)
    // #endif

    #if defined(SimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)
        #define SimTK_MOLMODEL_EXPORT __declspec(dllexport)

        // Keep MS VC++ quiet when it tries to instantiate incomplete template classes in a DLL.
        #ifdef _MSC_VER
        #pragma warning(disable:4661)

        // and lack of dll export of private members
        #pragma warning(disable:4251)
        #endif

    #elif defined(SimTK_MOLMODEL_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_MOLMODEL_EXPORT
    #else
        #define SimTK_MOLMODEL_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define SimTK_MOLMODEL_EXPORT // Linux, Mac
#endif

// Every SimTK Core library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_MOLMODEL_EXPORT void SimTK_version_molmodel(int* major, int* minor, int* build);
    SimTK_MOLMODEL_EXPORT void SimTK_about_molmodel(const char* key, int maxlen, char* value);
}

#endif // SimTK_MOLMODEL_COMMON_H_
