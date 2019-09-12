#ifndef SimTK_MASS_CENTER_MOTION_REMOVER_H_
#define SimTK_MASS_CENTER_MOTION_REMOVER_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/MolecularMechanicsSystem.h"

namespace SimTK {

/**
 * This is an event handler that can remove from a molecular system any unwanted
 * "rigid body" motion of the sytem as a whole. This is useful systems for which
 * you expect conservation of momentum to prevent overall motion, but where
 * numerical error or non-conservative disturbances (like the 
 * VelocityRescalingThermostat) may lead to unwanted motion. The handler will
 * remove overall linear and (optionally) angular momentum, and can also 
 * reposition the center of mass if it has drifted from a desired location.
 */

class MassCenterMotionRemover 
:   public PeriodicEventHandler {
public:
    /**
     * Create a MassCenterMotionRemover implemented as a periodic
     * event handler; the default behavior is that system linear and angular
     * momentum will be removed but the mass center will not be relocated.
     * You can change the default behavior with the methods mentioned below.

     * @param[in] system 
     *      The MolecularMechanicsSystem containing the bodies to be affected. 
     *      Typically this will be a Molmodel CompoundSystem (which is a
     *      kind of MolecularMechanicsSystem).
     * @param[in] adjustmentInterval 
     *      The time interval (in ps) at which to adjust the momentum and
     *      possibly relocate the mass center. The default value is 1ps.
     *
     * @see disableAngularMomentumRemoval(), enableMassCenterCorrection()
     */
    explicit MassCenterMotionRemover(const MolecularMechanicsSystem& system, 
                            Real adjustmentInterval = 1)   // ps
    :   PeriodicEventHandler(adjustmentInterval),
        system(system), linearMomOnly(false), 
        enableCOMCorrection(false), comTol(1), comPos(0)
    {
    }

    /**
     * Enable or disable removal of system angular momentum (enabled by 
     * default). If you set this flag to false then only linear momentum
     * will be removed. There is no way to disable linear momentum removal.
     *
     * @param[in] disable
     *      Set to true if you want the handler only to remove linear momentum.
     *      Otherwise both linear and angular momentum are removed, which
     *      is the default behavior.
     * @return A reference to this object (that was just modified).
     */
    MassCenterMotionRemover& disableAngularMomentumRemoval(bool disable) 
    {   linearMomOnly = disable; return *this; }
    /**
     * Find out whether this handler is currently set to remove system
     * angular momentum; if this returns true then only linear momentum
     * is being removed.
     */
    bool isAngularMomentumRemovalDisabled() const {return linearMomOnly;}

    /**
     * Enable or disable adjustment of the system mass center location 
     * (disabled by default). If you set this flag to true then the handler
     * will check if the system mass center has drifted by more than the
     * specified tolerance, and if so will attempt to move the system mass
     * center back to the specified desired location. Note that there are
     * default values of the tolerance (1 nm) and the desired location (0,0,0)
     * if you haven't set them explicitly.
     *
     * @param[in] enabled
     *      Set to true if you want the handler to attempt to keep the system
     *      mass center fixed; false if you want to leave it alone.
     *      Mass center location correction is disabled by default.
     * @return A reference to this object (that was just modified).
     */
    MassCenterMotionRemover& enableMassCenterCorrection(bool enabled) 
    {   enableCOMCorrection = enabled; return *this; }
    /**
     * Find out whether this handler is currently set to reposition the system
     * mass center if it drifts; if this returns false then the system mass
     * center will not be monitored or controlled.
     */
    bool isMassCenterCorrectionEnabled() const {return enableCOMCorrection;}

    /**
     * Set the desired location (as a Ground frame Vec3 in nm) to which the
     * system mass center is to be repositioned if it moves away by more
     * than the currently set tolerance (no effect unless mass center 
     * correction has been enabled).
     *
     * @param[in] groundLocationInNm
     *      When the \p adjustmentInterval has elapsed, the handler will
     *      compare the current system mass center location with this
     *      location and move it if it has drifted beyone the specified 
     *      tolerance.
     * @return A reference to this object (that was just modified).
     *
     * @see setMassCenterLocationTolerance()
     */
    MassCenterMotionRemover& 
    setDesiredMassCenterLocation(const Vec3& groundLocationInNm) 
    {   comPos = groundLocationInNm; return *this; }
    /**
     * Obtain the current value of the desired mass center location.
     * This is where the mass center will be relocated to if mass center
     * relocation is enabled.
     */
    const Vec3& getDesiredMassCenterLocation() const {return comPos;}

    /**
     * Set the tolerance (in nm) for the allowable system mass center
     * drift before the handler will reposition it (no effect unless mass
     * center correction has been enabled). If you want to disable mass 
     * center correction it is more efficent to call 
     * enableMassCenterCorrection(false) than to set this to a large value.
     *
     * @param[in] toleranceInNm
     *      When the \p adjustmentInterval has elapsed, the handler will
     *      compare the current system mass center with the desired location
     *      and move it if it has drifted beyone the specified tolerance which
     *      is given in nm. The default is 1 nm (10 Angstroms). The default
     *      for the desired mass center location is (0,0,0) but you can
     *      change that using setDesiredMassCenterLocation().
     * @return A reference to this object (that was just modified).
     *
     * @see setDesiredMassCenterLocation(), enableMassCenterCorrection()
     */
    MassCenterMotionRemover& 
    setMassCenterLocationTolerance(Real toleranceInNm) {
        SimTK_APIARGCHECK1_ALWAYS(toleranceInNm > 0, 
		    "MassCenterMotionRemover","setMassCenterLocationTolerance", 
		    "Illegal mass center tolerance %g.", toleranceInNm);
        comTol = toleranceInNm; 
        return *this; 
    }
    /**
     * Obtain the current value of the mass center location tolerance
     * in nm. If mass center correction is enabled, this is how far the
     * mass center can drift before it gets corrected.
     */
    Real getMassCenterLocationTolerance() const {return comTol;}

    /**
     * Attempt to remove the system's central linear momentum and angular
     * momentum unless that has been disabled.
     *
     * @param[in,out] state
     *      The State from which the current system momentum is 
     *      determined and to which the correction is applied.
     */
    void removeSystemMomentum(State& state) const
    {   system.removeSystemRigidBodyMomentum(state, linearMomOnly); }

    /**
     * Attempt to move the system mass center in the given State to this
     * handler's currently set desired mass center location. Tolerance is
     * not checked.
     *
     * @param[in,out] state
     *      The State from which the current system mass center location is 
     *      determined and to which the correction is applied.
     */
    void correctMassCenterLocation(State& state) const 
    {   system.moveSystemMassCenter(state, comPos); }

    /**
     * Determine the current distance of the system mass center from the
     * currently set desired mass center location. The error is calculated
     * regardless of whether mass center correction is currently enabled.
     *
     * @param[in] state
     *      The State from which the current system mass center location is
     *      determined.
     * @return
     *      The distance of the current system mass center location from the
     *      desired location, measured in nm.
     */
    Real calcMassCenterError(const State& state) const
    {   return (system.calcSystemMassCenterLocation(state) - comPos).norm(); }

    // This is the concrete implementation of the handleEvent() virtual
    // method required by every EventHandler; don't call this directly,
    // use removeSystemMomentum() or correctMassCenterLocation() instead.
    /// @cond -- don't show this in Doxygen
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const
    {
        system.realize(state, Stage::Velocity);
        removeSystemMomentum(state); 
        if (enableCOMCorrection && calcMassCenterError(state) > comTol) {
            correctMassCenterLocation(state);
        }
    }
    /// @endcond

private:
    const MolecularMechanicsSystem& system;
    bool    linearMomOnly; // should we correct angular momentum?
    bool    enableCOMCorrection; // should we relocate the COM?
    Real    comTol;    // how far to let the COM drift (nm).
    Vec3    comPos;    // where the COM is expected to stay.
};

} // namespace SimTK

#endif // SimTK_MASS_CENTER_MOTION_REMOVER_H_
