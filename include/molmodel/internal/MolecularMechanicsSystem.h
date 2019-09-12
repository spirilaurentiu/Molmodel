#ifndef SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_H_
#define SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_H_

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

#include "SimTKcommon.h"
#include "molmodel/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

#include <vector>

namespace SimTK {

class SimbodyMatterSubsystem;
class DuMMForceFieldSubsystem;

/**
 * This is a particular kind of MultibodySystem, one intended for use in
 * molecular mechanics (MM). The defining feature is that in addition to
 * the mandatory MatterSubsystem common to all MultibodySystems, this one
 * will also have a single MolecularMechanicsForceSubsystem. Unlike the base
 * MultibodySystem, which can employ any consistent set of units, a 
 * MolecularMechanicsSystem always uses MD units with length in nm, time
 * in ps, and mass in Da (g/mole). 
 *
 * There are also some solvers with a distinctively molecular flavor.
 */
class SimTK_MOLMODEL_EXPORT MolecularMechanicsSystem : public MultibodySystem {
public:
    MolecularMechanicsSystem();
    MolecularMechanicsSystem(SimbodyMatterSubsystem&, DuMMForceFieldSubsystem&);

    // Steals ownership of the source; returns subsystem ID number.
    int setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem&);
    const DuMMForceFieldSubsystem& getMolecularMechanicsForceSubsystem() const;
    DuMMForceFieldSubsystem&       updMolecularMechanicsForceSubsystem();

    /// Calculate this system's overall "rigid body" momentum about the
    /// system's mass center, in the configuration and with the velocities
    /// as supplied in the state parameter.
    /// @param[in] state
    ///     The State from which the configuration and velocities are taken.
    /// @return 
    ///     The system central momentum as a SpatialVec P (a 2-vector of Vec3's)
    ///     such that P[0] is the system angular momentum vector about its mass 
    ///     center and P[1] is the system's linear momentum vector.These vectors
    ///     are expressed in the Ground frame.
    /// @see removeSystemRigidBodyMomentum()
    SpatialVec calcSystemRigidBodyMomentum(const State& state) const;

    /// This solver attempts to remove the system's overall "rigid body" 
    /// momentum in the given State about the system's mass center, by changing
    /// only the velocity of base body mobilizers (that is, mobilizers that 
    /// connect particular bodies directly to Ground). This will always work
    /// as long as all base bodies are "free" -- meaning that they are
    /// connected to Ground by six-dof mobilizers for full rigid bodies, 
    /// three for point masses, and five for "linear" bodies like a CO2 
    /// molecule.
    ///
    /// The method is as follows:
    ///     - calculate the system central momentum.
    ///     - viewing the whole system as a single rigid body, calculate
    ///       the linear and angular velocity change for that body that
    ///       would eliminate its momentum.
    ///     - at every point where the "rigid body" is connected to ground
    ///       by a mobilizer (i.e., at the base bodies) calculate the
    ///       corresponding rigid body velocity that would be consistent
    ///       with the needed velocity change, and apply that to the base
    ///       body mobilizer generalized speeds.
    ///
    /// It is possible that some base bodies are unable to change their 
    /// velocities appropriately, perhaps because the associated mobilizer
    /// is not free or because there is a constraint or prescribed motion 
    /// present that prevents the adjustment. In that case the solver will
    /// fail to achieve zero momentum. However, there is no error return.
    /// If you want to check, use the calcSystemRigidBodyMomentum() method.
    ///
    /// @param[in,out] state
    ///     This is the State that supplies the current configuration and
    ///     velocities. The configuration remains unchanged, but the velocities
    ///     are modified as described above.
    /// @param[in] linearOnly
    ///     Normally both the linear and angular momentum are removed. If 
    ///     you set this parameter to true, only the linear momentum will
    ///     be removed.
    /// @see calcSystemRigidBodyMomentum() to check the results
    void removeSystemRigidBodyMomentum(State& state, bool linearOnly=false) const;

    /// Determine the instantaneous location of the overall system "rigid body"
    /// mass center for the configuration supplied in the given State 
    /// parameter.
    /// @param[in] state
    ///     The source for current configuration information.
    /// @return 
    ///     The location of the system mass center measured from the Ground
    ///     origin and expressed in the Ground frame.
    /// @see moveSystemMassCenter()
    Vec3 calcSystemMassCenterLocation(const State& state) const;

    /// This solver attempts to move the system's overall "rigid body" mass
    /// center to a given location by changing only the generalized coordinates
    /// (positions) of base bodies (that is, mobilized bodies that are connected
    /// directly to Ground). This is commonly used to force the system mass 
    /// center to stay in one location if it has been drifting due to numerical
    /// errors or other non-conservative effects.
    /// 
    /// This solver will always succeed as long as all the base bodies
    /// have unrestricted translational freedom (3 translational degrees of 
    /// freedom). In that case the \e relative position of all atoms in the system
    /// will be unchanged by this solver; they will move as though attached to
    /// a single rigid body. If some of the base bodies are unable to make the
    /// required translation, then the system will be distorted and the mass
    /// center will not be at the desired location. You can use the method
    /// calcSystemMassCenterLocation() to see whether the solver was successful;
    /// it will return quietly with a partial result if it fails.
    ///
    /// @param[in,out] state
    ///     The State from which current configuration information is taken,
    ///     and to which new base body positions are written so that the
    ///     system mass center is moved to a desired location.
    /// @param[in] newCOMLocation
    ///     The desired location for the system mass center, as a vector from
    ///     the Ground origin, expressed in the Ground frame.
    /// @see calcSystemMassCenterLocation()
    void moveSystemMassCenter(State& state, const Vec3& newCOMLocation) const;

    SimTK_PIMPL_DOWNCAST(MolecularMechanicsSystem, System);
private:

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    SubsystemIndex molecularMechanicsSub;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

} // namespace SimTK

#endif // SimTK_MOLMODEL_MOLECULAR_MECHANICS_SYSTEM_H_
