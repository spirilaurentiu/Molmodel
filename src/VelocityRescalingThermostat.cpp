/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-9 Stanford University and the Authors.         *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "molmodel/internal/VelocityRescalingThermostat.h"
#include <cmath>
using namespace SimTK;

/**
 * This class is the internal implementation for VelocityRescalingThermostat.
 */

class VelocityRescalingThermostat::VelocityRescalingThermostatImpl {
public:
    // Bath temp and num excluded dofs need to be initialized after this but
    // before any calculations are performed.
    VelocityRescalingThermostatImpl(const MolecularMechanicsSystem&  system) 
    :   system(system), Tb(NaN), numExcludedDofs(-1) 
    {
    }

	// This is the number of dofs. TODO: we're ignoring constraint redundancy
	// but we shouldn't be! That could result in negative dofs, so we'll 
	// make sure that doesn't happen. But don't expect meaningful results
	// in that case. Note that it is the acceleration-level constraints that
	// matter; they remove dofs regardless of whether there is a corresponding
	// velocity constraint.
    int calcNumThermalDofs(const State& state) const {
        assert(numExcludedDofs >= 0);
        const int ndofs = state.getNU() - state.getNUDotErr() - numExcludedDofs;
        return std::max(ndofs, 0);
    }

    Real calcCurrentTemp(const State& state) {
        const int nThermalDofs = calcNumThermalDofs(state);
        if (nThermalDofs==0) return 0;
        const Real ke = system.calcKineticEnergy(state);
        return (2*ke) / (nThermalDofs*SimTK_BOLTZMANN_CONSTANT_MD);
    }

    void rescale(State& state) {
        assert(!(isNaN(Tb) || numExcludedDofs < 0));
        system.realize(state, Stage::Velocity);
        Real T = calcCurrentTemp(state);
        if (T == 0) {
            if (Tb==0 || calcNumThermalDofs(state)==0)
                return; // can't do anything

            Vector impulse(state.getNU()), du(state.getNU());
            Random::Uniform random;
            for (int i=0; i < state.getNU(); ++i)
                impulse[i] = random.getValue();
            // TODO: this should use an operator that correctly deals with
            // constraints and prescribed motion, and it probably ought
            // to remove system momentum also.
            system.realize(state, Stage::Dynamics); // TODO: shouldn't be needed
            system.getMatterSubsystem().multiplyByMInv(state, impulse, du);
            state.updU() = du;
            system.realize(state, Stage::Velocity);
            T = calcCurrentTemp(state);
        }
        if (T > 0) {
            const Real scale = std::sqrt(Tb/T); // sqrt(desiredT/currentT)
            state.updU() *= scale;
            system.realize(state, Stage::Velocity);
        }
    }


private:
    const MolecularMechanicsSystem& system;
    Real                            Tb;
    int                             numExcludedDofs;

friend class VelocityRescalingThermostat;
};


VelocityRescalingThermostat::VelocityRescalingThermostat
   (const MolecularMechanicsSystem& system, Real bathTemperature, 
    Real rescalingInterval, int numExcludedDofs) 
:   PeriodicEventHandler(rescalingInterval) {
    impl = new VelocityRescalingThermostatImpl(system);
    setBathTemperature(bathTemperature);    // these check for errors
    setNumExcludedDofs(numExcludedDofs);
}

VelocityRescalingThermostat::~VelocityRescalingThermostat() {
    delete impl;
}

Real VelocityRescalingThermostat::
getBathTemperature() const {return impl->Tb;}

VelocityRescalingThermostat& VelocityRescalingThermostat::
setBathTemperature(Real bathTemperature) {
    SimTK_ERRCHK1_ALWAYS(bathTemperature >= 0,
        "VelocityRescalingThermostat::setBathTemperature()",
        "The requested bath temperature %g was illegal.", bathTemperature);
    impl->Tb = bathTemperature; 
    return *this;
}


int VelocityRescalingThermostat::
getNumExcludedDofs() const {return impl->numExcludedDofs;}

VelocityRescalingThermostat& VelocityRescalingThermostat::
setNumExcludedDofs(int numExcludedDofs) {
    SimTK_ERRCHK1_ALWAYS(0 <= numExcludedDofs && numExcludedDofs <= 6,
        "VelocityRescalingThermostat::setNumExcludedDofs()",
        "The number of excluded rigid body degrees of freedom requested"
        " was %d but must be between 0 and 6.", numExcludedDofs);

    impl->numExcludedDofs = numExcludedDofs; 
    return *this;
}

Real VelocityRescalingThermostat::
calcCurrentTemperature(const State& state) const 
{   return impl->calcCurrentTemp(state); }

int VelocityRescalingThermostat::
calcNumThermalDofs(const State& state) const
{   return impl->calcNumThermalDofs(state); }

void VelocityRescalingThermostat::
rescale(State& state) const 
{   impl->rescale(state); }

void VelocityRescalingThermostat::
handleEvent(State& state, Real accuracy, bool& shouldTerminate) const 
{
    impl->rescale(state);
}
