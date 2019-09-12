#ifndef SimTK_MOLMODEL_BONDGEOMETRY_H_
#define SimTK_MOLMODEL_BONDGEOMETRY_H_

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
#include "molmodel/internal/units.h"

namespace SimTK {

// Dihedral angles are define in terms of four atomic positions
// 
//   1
//    \
//     \
//      2------3
//              \
//               \
//                4
//
// calcDihedralAngle returns a dihedral angle in radians, 
// in the range (-Pi, Pi]
// given three unit vectors
// pointing in the direction of the 1->2 axis, the 2->3 axis, and the 3->4 axis, respectively.
// The answer is unchanged if the order of atoms is reversed from (1,2,3,4) to (4,3,2,1)
// This method can be used as a helper method for dihedral angles express using either 4 atomic
// location, or two bond center orientations.

/// Dihedral angle in radians in the range (-Pi, Pi]
Angle SimTK_MOLMODEL_EXPORT calcDihedralAngle(
        const Vec3& atomPos1,
        const Vec3& atomPos2, 
        const Vec3& atomPos3, 
        const Vec3& atomPos4);

/// Dihedral angle in radians in the range (-Pi, Pi]
Angle SimTK_MOLMODEL_EXPORT calcDihedralAngle(
        const UnitVec3& bond12, 
        const UnitVec3& bond23, 
        const UnitVec3& bond34);


} // namespace SimTK

#endif
