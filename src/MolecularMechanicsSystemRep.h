#ifndef SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_REP_H_
#define SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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
 * Define the private implementation of the MolecularMechanicsSystem
 * class (a kind of MultibodySystem).
 */

#include "SimTKcommon.h"

#include "molmodel/internal/common.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"

// #include "simbody/private/MultibodySystemRep.h"

namespace SimTK {

    /**
 * This class is a kind of MultibodySystem which is required to have exactly
 * one DuMMForceFieldSubsystem.
 */
//class MolecularMechanicsSystemRep : public MultibodySystemRep {
//public:
//    MolecularMechanicsSystemRep() : MultibodySystemRep()
//    {
//    }
//    ~MolecularMechanicsSystemRep() {
//    }
//
//    int setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
//        assert(!molecularMechanicsSub.isValid());
//        molecularMechanicsSub = addForceSubsystem(mm);
//        return molecularMechanicsSub;
//    }
//
//    const DuMMForceFieldSubsystem& getMolecularMechanicsForceSubsystem() const {
//        assert(molecularMechanicsSub.isValid());
//        return DuMMForceFieldSubsystem::downcast(getSubsystem(molecularMechanicsSub));
//    }
//    DuMMForceFieldSubsystem& updMolecularMechanicsForceSubsystem() {
//        assert(molecularMechanicsSub.isValid());
//        return DuMMForceFieldSubsystem::updDowncast(updSubsystem(molecularMechanicsSub));
//    }
//
//    SimTK_DOWNCAST(MolecularMechanicsSystemRep, System::Guts);
//private:
//    SubsystemIndex molecularMechanicsSub;
//};
//
//} // namespace SimTK

#endif // SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_REP_H_
