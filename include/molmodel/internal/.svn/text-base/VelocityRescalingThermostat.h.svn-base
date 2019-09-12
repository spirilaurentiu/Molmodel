#ifndef SimTK_MOLMODEL_VELOCITY_RESCALING_THERMOSTAT_H_
#define SimTK_MOLMODEL_VELOCITY_RESCALING_THERMOSTAT_H_
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

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/MolecularMechanicsSystem.h"

namespace SimTK {

/**
 * This is an event handler that acts as a thermostat for controlling the 
 * temperature of a simulation efficiently, and in a qualitatively reasonable 
 * manner although not with the correct statistical temperature distribution.
 * For thermodynamically correct temperature control, use the Nose'-Hoover
 * thermostat instead. All values are assumed to be in MD units.
 *
 * The system being simulated is assumed to be embedded in a bath of 
 * infinite heat capacity that is at a constant temperature Tb. This thermostat
 * will periodically adjust the system temperature T so that T=Tb. The system
 * temperature T is defined here as follows:
 * <pre>    T = 2*KE/(N*kB)     </pre>
 * where KE is the total system kinetic energy, kB is Boltzmann's constant,
 * and N is the number of \e thermal degrees of freedom. By "thermal" 
 * degrees of freedom we mean the total number of mobilities, minus 
 * non-redundant constraints, minus the overall rigid body degrees of freedom 
 * for the system as a whole. By default, we assume there are 6 rigid body 
 * degrees of freedom for the system as a whole that are conserved and hence
 * not part of the degrees of freedom available for temperature control. This
 * will be true when no part of the system is tethered to ground, and all 
 * forces are conservative, in which case the system linear momentum and
 * angular momentum are conserved and should remain zero at all times. This
 * thermostat does not attempt to manage the system momentum; it is scaled
 * along with everything else so it is up to you to ensure that it stays zero
 * if it should be. The number of conserved rigid body degrees of freedom is
 * almost always 6; set it to 3 if your system will conserve only linear
 * momentum and to zero if there is no conservation expected, such as when
 * a system is tethered to ground.
 *
 * This thermostat works by velocity rescaling. At regular intervals, it 
 * calculates the total kinetic energy KE of the system, then rescales all 
 * the velocities so the total kinetic energy will exactly equal
 * N*kB*T/2 ensuring an average energy of kB*T/2 per thermal degree of 
 * freedom.  
 */

class SimTK_MOLMODEL_EXPORT VelocityRescalingThermostat 
:   public PeriodicEventHandler {
public:
    /**
     * Create a VelocityRescalingThermostat implemented as a periodic
     * event handler.
     * @param[in] system 
     *      The MolecularMechanicsSystem containing the bodies to be affected. 
     *      Typically this will be a Molmodel CompoundSystem (which is a
     *      kind of MolecularMechanicsSystem).
     * @param[in] bathTemperature 
     *      The temperature Tb to which the system temperature T is to 
     *      be periodically adjusted. The default value is 293.15 Kelvin 
     *      (20 C). Note that the definition of system temperature T is 
     *      influenced by the \p numExcludedDofs parameter.
     * @param[in] rescalingInterval 
     *      The time interval (in ps) at which to rescale velocities. 
     *      The default value is 1ps. Velocity will be rescaled to match
     *      the bath temperature whenever this much simulated time
     *      has elapsed; it is free to drift uncontrolled in between.
     *      This parameter is used to construct the PeriodicEventHandler
     *      that is the base class for the VelocityRescalingThermostat.
     * @param[in] numExcludedDofs
     *      This is the number of rigid body degrees of freedom that
     *      will be conserved during a simulation, in the range 0-6.
     *      The default is 6, meaning that the system as a whole is
     *      going to conserve angular and linear momentum. The degrees of 
     *      freedom here are not included when calculating temperature.
     *
     * @see class MassCenterMotionRemover
     */
    explicit VelocityRescalingThermostat
       (const MolecularMechanicsSystem& system, 
        Real bathTemperature    = 293.15, 
        Real rescalingInterval  = 1,
        int  numExcludedDofs    = 6);

    /// Get the temperature this thermostat is set to maintain.
    /// @return The currently set bath temperture Tb, in Kelvins.
    Real getBathTemperature() const;

    /// Set the bath temperature Tb that this thermostat is set to maintain.
    /// This may be changed during a simulation to effect an annealing
    /// scheduled.
    /// @param[in] bathTemperature 
    ///     The bath temperature Tb to which the system temperature T is
    ///     periodically adjusted.
    /// @return 
    ///     A writable reference to the current VelocityRescalingThermostat
    ///     object that was just modified.
    VelocityRescalingThermostat& setBathTemperature(Real bathTemperature);


    /// Obtain the current setting for the number of rigid body degrees
    /// of freedom to be excluded when determining the number of thermal
    /// degrees of freedom used in calculating system temperature.
    /// @return The currently set number of excluded dofs, 0-6.
    int getNumExcludedDofs() const;

    /// Set the number of rigid body degrees of freedom to be excluded
    /// from the count of thermal degrees of freedom used in temperature 
    /// calculation.
    /// @param[in] numExcludedDofs 
    ///     The number of rigid body dofs to be excluded at the next
    ///     rescaling request.
    /// @return 
    ///     A writable reference to the current VelocityRescalingThermostat
    ///     object that was just modified.
    VelocityRescalingThermostat& setNumExcludedDofs(int numExcludedDofs);


    /// Calculate the current system temperature in the same manner as
    /// this event handler will use when adjusting the temperature,
    /// taking account of the number of excluded degrees of freedom.
    /// @param[in] state 
    ///     The State from which the system's current velocities are to
    ///     be obtained for use in calculating the temperature.
    /// @return
    ///     The instantaneous system temperature T calculated using the
    ///     equation T=2*KE/(N*kB) where KE is the total kinetic energy,
    ///     kB is Boltzmann's constant, and N is the number of available
    ///     degrees of freedom as defined above.
    /// @see calcNumThermalDofs()
    Real calcCurrentTemperature(const State& state) const;

    /// Calculate the number of \e thermal degrees of freedom currently possessed
    /// by the system, as the total number of mobilities minus the number
    /// of non-redundant constraint equations, minus the number of excluded
    /// rigid body degrees of freedom. This is the number N that is used in
    /// calculation of the current system temperature.
    /// @param[in] state 
    ///     The State from which the count of mobilities and constraint
    ///     equations is obtained.
    /// @return The net number of thermal degrees of freedom available.
    /// @see calcCurrentTemperature()
    int calcNumThermalDofs(const State& state) const;

    /// Rescale the system velocities to adjust the system temperature T
    /// to be equal to the specified bath temperature Tb. This provides
    /// the same functionality as the event handler but is intended to
    /// be invoked directly rather than from a TimeStepper.
    /// @note
    /// If the current temperature is \e exactly zero (meaning all the
    /// generalized speeds u in the system are zero) then rescaling
    /// would have no effect. In that case, and assuming the desired
    /// temperature is greater than zero and there are some thermal
    /// degrees of freedom available, this method will generate a
    /// set of random velocities with the appropriate temperature.
    /// @param[in,out] state
    ///     The State from which the current velocities are obtained to
    ///     calculate the current temperature, and to which the revised
    ///     velocities are written.
    /// @see calcCurrentTemperature()
    void rescale(State& state) const;

    // This is the concrete implementation of the handleEvent() virtual
    // method required by every EventHandler; don't call this directly,
    // use rescale() instead.
    /// @cond -- don't show this in Doxygen
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const;
    /// @endcond

    ~VelocityRescalingThermostat();

    // Don't show these in Doxygen.
    /// @cond
    // OBSOLETE, use getBathTemperature() instead.
    Real getTemperature() const {return getBathTemperature();}
    // OBSOLETE, use setBathTemperature() instead.
    void setTemperature(Real t) {setBathTemperature(t);}
    /// @endcond
private:
    class VelocityRescalingThermostatImpl;
    VelocityRescalingThermostatImpl* impl;
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_VELOCITY_RESCALING_THERMOSTAT_H_
