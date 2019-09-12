#ifndef SimTK_MOLMODEL_UNITS_H_
#define SimTK_MOLMODEL_UNITS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "molmodel/internal/common.h"
#include "SimTKcommon.h"

namespace SimTK {

// These typedefs are a documentation tool for method signatures

// Angles are not specific to small scales, so are not in the mdunits namespace

    /**
 * \brief angle in radians
 */
typedef Real Angle; // in radians

typedef Real LineAngle; // angle in range (0,Pi)
typedef Real CircleAngle; // angle in range (-Pi,Pi) or in range (0,2Pi)
typedef Real WindingAngle; // angle in range (-infinity, +infinity)

namespace mdunits {
/**
 * \brief distance in nanometers
 */
typedef Real Length; // in nanometers

/**
 * \brief mass in atomic mass units
 */
typedef Real Mass; // in atomic mass units

typedef Real Energy; // in kilojoules per mole

typedef Real Charge; // in elementary charge units

} // namespace mdunits


static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
static const Real Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

}

#endif // SimTK_MOLMODEL_UNITS_H_
