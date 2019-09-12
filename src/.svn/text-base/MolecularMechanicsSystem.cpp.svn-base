/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-9 Stanford University and the Authors.         *
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

/**@file
 *
 * Implementation of MolecularMechanicsSystem, a kind of MultibodySystem.
 */

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/MolecularMechanicsSystem.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"

#include <vector>

namespace SimTK {


    ////////////////////////////////
    // MOLECULAR MECHANICS SYSTEM //
    ////////////////////////////////

class DuMMForceFieldSubsystem;

 /*static*/ bool 
MolecularMechanicsSystem::isInstanceOf(const System& s) {
    // return MolecularMechanicsSystemRep::isA(s.getSystemGuts());
    return MultibodySystem::isInstanceOf(s);
}

/*static*/ const MolecularMechanicsSystem&
MolecularMechanicsSystem::downcast(const System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MolecularMechanicsSystem&>(s);
}
/*static*/ MolecularMechanicsSystem&
MolecularMechanicsSystem::updDowncast(System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MolecularMechanicsSystem&>(s);
}

MolecularMechanicsSystem::MolecularMechanicsSystem() 
{   // Give a hint to the Visualizer that it shouldn't show a ground plane.
    setUseUniformBackground(true);
}

MolecularMechanicsSystem::MolecularMechanicsSystem
   (SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem& mm)
{   // Give a hint to the Visualizer that it shouldn't show a ground plane.
    setUseUniformBackground(true);
    setMatterSubsystem(matter);
    setMolecularMechanicsForceSubsystem(mm);
}

int MolecularMechanicsSystem::setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
    assert(!molecularMechanicsSub.isValid());
    molecularMechanicsSub = SubsystemIndex( addForceSubsystem(mm) );
    return molecularMechanicsSub;
}

const DuMMForceFieldSubsystem& MolecularMechanicsSystem::getMolecularMechanicsForceSubsystem() const {
    assert(molecularMechanicsSub.isValid());
    return DuMMForceFieldSubsystem::downcast(getSubsystem(molecularMechanicsSub));
}
DuMMForceFieldSubsystem& MolecularMechanicsSystem::updMolecularMechanicsForceSubsystem() {
    assert(molecularMechanicsSub.isValid());
    return DuMMForceFieldSubsystem::updDowncast(updSubsystem(molecularMechanicsSub));
}

SpatialVec MolecularMechanicsSystem::
calcSystemRigidBodyMomentum(const State& state) const {
    return getMatterSubsystem().calcSystemCentralMomentum(state);
}

Vec3 MolecularMechanicsSystem::
calcSystemMassCenterLocation(const State& state) const {
    return getMatterSubsystem().calcSystemMassCenterLocationInGround(state);
}

void MolecularMechanicsSystem::
removeSystemRigidBodyMomentum(State& state, bool linearOnly) const {
    const SimbodyMatterSubsystem& matter = this->getMatterSubsystem();

    this->realize(state, Stage::Position);
    const Real       sysMass = matter.calcSystemMass(state);
    
    if (sysMass==0) 
        return; // nothing to do

    const Vec3       sysCOM  = matter.calcSystemMassCenterLocationInGround(state);
    const Inertia    sysICM  = matter.calcSystemCentralInertiaInGround(state);

    this->realize(state, Stage::Velocity);
    const SpatialVec sysMom  = matter.calcSystemCentralMomentum(state);

    // The angular velocity change we need is ICM^-1 * AM where ICM is the
    // system central inertia matrix and AM is the central angular momentum.
    const Vec3 dw_G = linearOnly ? Vec3(0) 
                                 : sysICM.toMat33().invert()*sysMom[0];

    // This is the change in linear velocity that would be needed for base
    // body mobilizers coincident with the ground origin.
    const Vec3 dv_G = sysMom[1]/sysMass + sysCOM % dw_G;

    std::vector< std::pair<const MobilizedBody*, SpatialVec> > baseAdjustments;
    for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
        const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
        if (mobod.getLevelInMultibodyTree() != 1)
            continue;

        // This body is connected directly to ground; i.e., it is a 
        // base body.
        const Transform&  X_GF = mobod.getDefaultInboardFrame(); // TODO: get from state
        const Transform&  X_FM = mobod.getMobilizerTransform(state);
        const SpatialVec& V_FM = mobod.getMobilizerVelocity(state); // current velocity

        const Transform X_GM = X_GF*X_FM;

        const Vec3 dw_F = ~X_GF.R() * dw_G; // re-express in F frame

        // Remove the extra linear velocity will be produced due to the base
        // body being away from the ground origin; convert to F frame.
        const Vec3 dv_F = ~X_GF.R() *(dv_G - X_GM.p() % dw_G);

        const SpatialVec desiredV_FM = V_FM - SpatialVec(dw_F, dv_F);
        baseAdjustments.push_back(std::make_pair(&mobod, desiredV_FM));
    }

    // Now make all the adjustments at once.
    for (unsigned i=0; i < baseAdjustments.size(); ++i) {
        const MobilizedBody& mobod       = *baseAdjustments[i].first;
        const SpatialVec&    desiredV_FM = baseAdjustments[i].second;
        mobod.setUToFitVelocity(state, desiredV_FM);
    }
    this->realize(state, Stage::Velocity);
}

void MolecularMechanicsSystem::
moveSystemMassCenter(State& state, const Vec3& newCOMLocation) const {
    const SimbodyMatterSubsystem& matter = this->getMatterSubsystem();

    this->realize(state, Stage::Position);
    const Vec3 sysCOM  = matter.calcSystemMassCenterLocationInGround(state);
    const Vec3 COMErr = sysCOM - newCOMLocation;

    std::vector< std::pair<const MobilizedBody*, Vec3> > baseAdjustments;
    for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
        const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
        if (mobod.getLevelInMultibodyTree() != 1)
            continue;

        // This body is connected directly to ground; i.e., it is a 
        // base body.
        const Transform&  X_GF = mobod.getDefaultInboardFrame(); // TODO: get from state
        const Transform&  X_FM = mobod.getMobilizerTransform(state);

        const Vec3 desiredP_FM = X_FM.p() - ~X_GF.R()*COMErr;
        baseAdjustments.push_back(std::make_pair(&mobod, desiredP_FM));
    }

    // Now make all the adjustments at once.
    for (unsigned i=0; i < baseAdjustments.size(); ++i) {
        const MobilizedBody& mobod       = *baseAdjustments[i].first;
        const Vec3&          desiredP_FM =  baseAdjustments[i].second;
        mobod.setQToFitTranslation(state, desiredP_FM);
    }
    this->realize(state, Stage::Position);
}


} // namespace SimTK

