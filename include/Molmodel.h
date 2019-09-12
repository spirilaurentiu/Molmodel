#ifndef SimTK_MOLMODEL_MOLMODEL_H_
#define SimTK_MOLMODEL_MOLMODEL_H_

/* -------------------------------------------------------------------------- *
 *                             SimTK Molmodel(tm)                             *
 * -------------------------------------------------------------------------- *
 * Molmodel is part of the SimTK biosimulation toolkit originating from       *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/molmodel. *
 *                                                                            *
 * Portions copyright (c) 2007-11 Stanford University and the Authors.        *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * Contributors: Peter Eastman, Randy Radmer, Mark Friedrichs, Samuel Flores  *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * This sentence, the above copyright and permission notices, and the         *
 * following disclaimer shall be included in all copies or substantial        *
 * portions of the Software.                                                  *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- *
 */

/** @file
This is the header file that user code should include to pick up all Molmodel 
capabilities. Note that all symbols defined here will be in the SimTK 
namespace, or (where a namespace can't be used) prefixed by "SimTK_". 
Molmodel depends on SimTK's Simbody code having already been installed; see 
https://simtk.org/home/simbody for more information. **/

#include "Simbody.h"
#include "SimTKmolmodel.h"

#endif // SimTK_MOLMODEL_MOLMODEL_H_
