#ifndef SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_H_
#define SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                             SimTK Molmoldel(tm)                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/molmodel. *
 *                                                                            *
 * Portions copyright (c) 2007-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Bruns, Peter Eastman, Randy Radmer, Samuel Flores      *
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
Define the public interface to DuMMForceFieldSubsystem, a subsystem which
provides molecular mechanics capability in a multibody framework. **/

#include "SimTKcommon.h"
#include "simbody/internal/ForceSubsystem.h"

#include "molmodel/internal/common.h"
#include "molmodel/internal/Biotype.h"

#include "OpenMMPlugin.h"


#include <cassert>


/**@defgroup MolecularMechanics     Molecular Mechanics in Molmodel
Once you have built a system of Compounds (molecules) in Molmodel, and modeled 
them with the degrees of freedom you want, you can apply molecular forces using
DuMM, the Molmodel molecular mechanics force field. **/

namespace SimTK {

class MolecularMechanicsSystem;

/**@ingroup MolecularMechanics
This namespace is used for symbols which are useful in conjunction with 
Molmodel's DuMMForceFieldSubsystem. **/
namespace DuMM {
/** @class SimTK::DuMM::AtomIndex
This is the unique integer index that DuMM assigns to atoms as they are 
added via addAtom(). **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomIndex);
/** @class SimTK::DuMM::IncludedAtomIndex
This is the unique integer index that DuMM assigns to atoms that are to be
included in force calculations.\ These represent a subset of all the
atoms and you can map from an IncludedAtomIndex to the corresponding 
AtomIndex. These are assigned during realizeTopology(). **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(IncludedAtomIndex);
/** @class SimTK::DuMM::NonbondAtomIndex
This is the unique integer index that DuMM assigns to included atoms that 
are involved in nonbonded force calculations (Coulomb, van der Waals,
and/or GBSA).\ These represent a subset of all the included atoms and you can
map from a NonbondAtomIndex to the corresponding IncludedAtomIndex. 
These are assigned during realizeTopology(). **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(NonbondAtomIndex);
/** @class SimTK::DuMM::BondIndex
This is the unique integer index that DuMM assigns to bonds as they are
added via addBond(). **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(BondIndex);
/** @class SimTK::DuMM::ClusterIndex
This is the unique integer index that DuMM assigns to clusters as they are
created via createCluster(). **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ClusterIndex);
/** @class SimTK::DuMM::AtomClassIndex
This is a unique integer associated with each "atom class", assigned by the
user when the atom class is first introduced via defineAtomClass(). This is 
really an \e id rather than an \e index since these need not be consecutively 
assigned. Typically these are defined by the force field being used. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomClassIndex);
/** @class SimTK::DuMM::ChargedAtomTypeIndex
This is a unique integer associated with each "charged atom type", assigned by
the user when the charged atom type is first introduced via 
defineChargedAtomType(). This is really an \e id rather than an \e index since 
these need not be consecutively assigned. Typically these are defined by the 
force field being used. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ChargedAtomTypeIndex);


/** @defgroup MMCustomTerms User defined molecular mechanics force terms
 *  @ingroup MolecularMechanics
 *
 * These classes are used with DuMMForceFieldSubsystem to permit user-defined 
 * functional forms for bonded terms.
 *
 * \b Basic features
 * - You can make custom bond stretch, bond bend, and bond torsion terms.
 * - Each term is represented by an abstract base class with a virtual 
 *   destructor.
 * - You write one scalar method to calculate energy and one to calculate force.
 * - You can associate multiple custom bond terms with the same atom classes; 
 *   in that case the terms are additive.
 * - You do not have to supply terms for every combination of atom classes that
 *   appears in your system; your methods will be called only for combinations 
 *   for which you have actually supplied a term.
 * - It is your responsibility to ensure that your force (or torque) method 
 *   returns the negative gradient of your energy method.
 *
 * \b Notes
 * - %DuMM may skip calling your function if the atoms all reside on the same 
 *   body; that can't affect the motion, but note that the internal energy will 
 *   not be counted in that case.
 * - %DuMM will take over ownership of the concrete object you supply to its
 *   defineCustomBond...() methods; don't delete them yourself.
 *
 * \b AmberTorsion example
 *
 * This demonstrates how the built-in Amber force might be applied as a custom 
 * torsion. This is a toy example because this particular form is already built
 * into DuMMForceFieldSubsystem.
 * 
 * \code
 * class AmberTorsion : public DuMM::CustomBondTorsion {
 * public:
 *    AmberTorsion(Real amplitudeInKJ, Real phaseInRad, int periodicity)
 *    :   amplitude(amplitudeInKJ), phase(phaseInRad), periodicity(periodicity) {}
 * 
 *    Real calcEnergy(Real dihedralAngleInRad) {
 *        return amplitude * 
 *                  ( 1 + std::cos(periodicity * dihedralAngleInRad - phase) );
 *    }
 *
 *    // Note: this is MINUS the gradient of the energy function!
 *    Real calcTorque(Real dihedralAngleInRad) {
 *        return periodicity * amplitude * 
 *                  std::sin(periodicity * dihedralAngleInRad - phase);
 *    }
 *
 * private:
 *    Real amplitude;       // in KJ/mol
 *    Real phase;           // in radians
 *    int  periodicity;
 * };
 * 
 * #include "SimTKmolmodel.h"
 * int main() {
 *     DuMMForceFieldSubsystem dumm;
 *     dumm.loadAmber99Parameters();
 *     // torsion for aliphatic H-C-C-H torsions like those in ethane
 *     dumm.defineCustomTorsion(
 *         dumm.getAtomClassIndex("HC"),
 *         dumm.getAtomClassIndex("CT"),
 *         dumm.getAtomClassIndex("CT"),
 *         dumm.getAtomClassIndex("HC"),
 *         new AmberTorsion(0.150*DuMM::Kcal2KJ, 0*DuMM::Deg2Rad, 3));
 *     ...
 * }
 * \endcode
 */
//@{

/**
 * Abstract base class for custom bond stretch terms, that is, 
 * functional forms that apply a distance-based force along
 * a 1-2 bond between a pair of atoms.
 *
 * @see DuMMForceFieldSubsystem::defineCustomBondStretch()
 */
class CustomBondStretch {
public:
    virtual ~CustomBondStretch() {}
    virtual Real calcEnergy(Real distance) const = 0;
    virtual Real calcForce(Real distance) const = 0;
};

/**
 * Abstract base class for custom bond bend functions, that is,
 * functional forms that apply an angle-based torque between
 * the two lines formed by a 1-2-3 triple of bonded atoms.
 *
 * @see DuMMForceFieldSubsystem::defineCustomBondBend()
 */
class CustomBondBend {
public:
    virtual ~CustomBondBend() {}
    virtual Real calcEnergy(Real bendAngle) const = 0;
    virtual Real calcTorque(Real bendAngle) const = 0;
};

/**
 * Abstract base class for custom torsion functions, that is,
 * functional forms that apply a dihedral-angle based torque about
 * the middle bond of a 1-2-3-4 quadruple of bonded atoms.
 *
 * \b AmberTorsion example
 *
 * This demonstrates how the built-in Amber force might be applied as a custom 
 * torsion. This is a toy example because this particular form is already built
 * into DuMMForceFieldSubsystem.
 * 
 * \code
 * class AmberTorsion : public DuMM::CustomBondTorsion {
 * public:
 *    AmberTorsion(Real amplitudeInKJ, Real phaseInRad, int periodicity)
 *    :   amplitude(amplitudeInKJ), phase(phaseInRad), periodicity(periodicity) {}
 * 
 *    Real calcEnergy(Real dihedralAngleInRad) {
 *        return amplitude * 
 *                  ( 1 + std::cos(periodicity * dihedralAngleInRad - phase) );
 *    }
 *
 *    // Note: this is MINUS the gradient of the energy function!
 *    Real calcTorque(Real dihedralAngleInRad) {
 *        return periodicity * amplitude * 
 *                  std::sin(periodicity * dihedralAngleInRad - phase);
 *    }
 *
 * private:
 *    Real amplitude;       // in KJ/mol
 *    Real phase;           // in radians
 *    int  periodicity;
 * };
 * 
 * #include "SimTKmolmodel.h"
 * int main() {
 *     DuMMForceFieldSubsystem dumm;
 *     dumm.loadAmber99Parameters();
 *     // torsion for aliphatic H-C-C-H torsions like those in ethane
 *     dumm.defineCustomTorsion(
 *         dumm.getAtomClassIndex("HC"),
 *         dumm.getAtomClassIndex("CT"),
 *         dumm.getAtomClassIndex("CT"),
 *         dumm.getAtomClassIndex("HC"),
 *         new AmberTorsion(0.150*DuMM::Kcal2KJ, 0*DuMM::Deg2Rad, 3));
 *     ...
 * }
 * \endcode
 *
 * @see DuMMForceFieldSubsystem::defineCustomBondTorsion()
 */
class CustomBondTorsion {
public:
    virtual ~CustomBondTorsion() {}
    virtual Real calcEnergy(Real dihedralAngle) const = 0;
    virtual Real calcTorque(Real dihedralAngle) const = 0;
};
//@}

/**@defgroup MMConversions          Handy conversion constants
 * @ingroup MolecularMechanics
 *
 * These are a set of multiplicative factors for use in converting
 * among the various commonly used molecular mechanics units.
 *
 * Note that these are compilation-unit statics, not members.
 * That way we can be sure they are initialized before being used.
 * To use these, multiply something in units on left of the "2" to get equivalent
 * in units on right. E.g., 180*Deg2Rad gives Pi radians.
 *
 * There are several conventions for giving van der Waals parameters.
 * Rmin is the radius at which the energy well minimum is seen
 * (actually it is 1/2 the distance between atom centers for a pair
 * of atoms of this class interacting with that minimum energy).
 * This is \e not the Lennard-Jones term (half-) Sigma, which is the radius 
 * (half distance between atom centers) at which the energy crosses zero, 
 * that is, a little closer together than when the energy well is at maximum depth. 
 * To convert for Lennard-Jones: Rmin = 2^(1/6) * Sigma. NOTE: in many
 * cases the Lennard-Jones data will be "diameter" Sigma, that is, twice
 * the value our conversion constant is expecting -- be careful!
 */
//@{
static const Real Ang2Nm  = (Real)0.1L; ///< angstroms to nanometers
static const Real Nm2Ang  = (Real)10.L; ///< nanometers to angstroms
static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN; ///< degrees to radians
static const Real Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE; ///< radians to degrees
static const Real KJ2Kcal = (Real)SimTK_KJOULE_TO_KCAL;   ///< kilojoules to kilocalories
static const Real Kcal2KJ = (Real)SimTK_KCAL_TO_KJOULE;   ///< kilocalories to kilojoules
/// half-Sigma to van der Waals radius; caution -- see discussion in module description.
static const Real Sigma2Radius = (Real)std::pow(2.L,  1.L/6.L);
/// van der Waals radius to half-Sigma; caution -- see discussion in module description.
static const Real Radius2Sigma = (Real)std::pow(2.L, -1.L/6.L);
//@}

} // namespace DuMM

/** @addtogroup MolecularMechanics */
//@{

/** This is a concrete subsystem that provides basic molecular mechanics 
functionality for coarse-grained molecules built in the SimTK framework.

UNITS: This subsystem requires that the system be modeled in "MD units"
of nanometers, daltons (g/mol), and picoseconds, yielding consistent
energy units of kJ/mol==(Da-nm^2/ps^2), and force in kJ/mol-nm. Charge
is in proton charge units e, and angles are in radians.
For convenience, we allow the force field to be defined in "KA" units,
that is, angstroms instead of nanometers, and energy in kcal rather
than kJ, and we also allow angles to be supplied in degrees. However,
these are immediately converted to the MD units described above. **/
class SimTK_MOLMODEL_EXPORT DuMMForceFieldSubsystem : public ForceSubsystem {
public:
/** These are the van der Waals mixing rules supported by DuMM. **/
enum VdwMixingRule {
    WaldmanHagler       = 1,    ///< Our default, Waldman & Hagler, J.Comp.Chem. 14(9) 1993
    HalgrenHHG          = 2,    ///< MMFF, AMOEBA
    Jorgensen           = 3,    ///< OPLS
    LorentzBerthelot    = 4,    ///< AMBER, CHARMM
    Kong                = 5     ///< Kong, J.Chem.Phys. 59(5) 1973
};

DuMMForceFieldSubsystem();
explicit DuMMForceFieldSubsystem(MolecularMechanicsSystem&);



/** @name               Define particular molecules
Methods in this group are used to define the atoms and bonds in the particular 
set of molecules being simulated. **/
/**@{**/

/// Add a new atom to the model. The AtomIndex number is returned; you don't
/// get to pick your own. Use the AtomIndex to identify this particular 
/// atom subsequently.
DuMM::AtomIndex addAtom(DuMM::ChargedAtomTypeIndex chargedAtomTypeIx);

/// Declare that there is a covalent bond between two atoms. Note that 
/// these are AtomIndex numbers, not AtomClasses or ChargedAtomTypes.
DuMM::BondIndex addBond(DuMM::AtomIndex atom1Ix, DuMM::AtomIndex atom2Ix);

/// For a given 1-2 bond, return the atoms which are connected by that bond.
/// You select one atom at a time by setting parameter \p which to 0 or 1.
/// 0 will return the lower-numbered AtomIndex, regardless of the order 
/// in which the atoms were specified to addBond().
DuMM::AtomIndex  getBondAtom(DuMM::BondIndex bond, int which) const;

/// How many atoms are currently in the model?
int getNumAtoms() const;
/// How many 1-2 bonds are currently in the model?
int getNumBonds() const;

/// Obtain the mass in Daltons (g/mol) of the atom indicated by the given 
/// AtomIndex.
Real getAtomMass(DuMM::AtomIndex atomIx) const;
/// Obtain the element (by atomic number) of the atom indicated by the 
/// given AtomIndex.
int  getAtomElement(DuMM::AtomIndex atomIx) const;
/// Obtain the van der Waals radius of the atom indicated by the given 
/// AtomIndex.
Real getAtomRadius(DuMM::AtomIndex atomIx) const;
/// Obtain the Simbody MobilizedBodyIndex of the rigid body on which a 
/// particular atom has been fixed. An exception will be thrown if this 
/// atom is not fixed to any body.
MobilizedBodyIndex getAtomBody(DuMM::AtomIndex atomIx) const;
/// Obtain the station at which a particular atom is fixed on its body. 
/// An exception will be thrown if this atom is not fixed to any body.
Vec3 getAtomStationOnBody(DuMM::AtomIndex atomIx) const;

// EU BEGIN
/// Set the station at which a particular atom is fixed on its body. 
/// An exception will be thrown if this atom is not fixed to any body.
void bsetAtomStationOnBody(DuMM::AtomIndex atomIx, Vec3 new_station_B);

// For CalcFullPotential Eliza
void bsetAllAtomStationOnBody(DuMM::AtomIndex atomIx, Vec3 new_station_B);

/// Set AtomPlacement station coressponding to DuMMAtom atomIx 
void bsetAtomPlacementStation(DuMM::AtomIndex atomIx, MobilizedBodyIndex inputMbx, Vec3 new_station);

// Stations computed every time
Vec3& updIncludedAtomStation(DuMM::AtomIndex atomIx);

// For CalcFullPotential Eliza
Vec3& updAllAtomStation(DuMM::AtomIndex atomIx);

// Get clusterIndex of a specified mobod
DuMM::ClusterIndex bgetMobodClusterIndex(MobilizedBodyIndex inputMbx) const;
// EU END

/// Obtain the station at which a particular atom is fixed within a 
/// particular Cluster (an atom can be in more than one Cluster).
/// An exception will be thrown if this atom is not fixed to the cluster.
Vec3 getAtomStationInCluster(DuMM::AtomIndex atomIx, DuMM::ClusterIndex clusterIx) const;

/// For display purposes, return the RGB value of a suggested color for an 
/// element given by atomic number. For example, if the atomicNumber is 8 
/// (Oxygen) the suggested color will be Red (1,0,0).
Vec3 getElementDefaultColor(int atomicNumber) const;
/// For display purposes, return the RGB value of a suggested color with 
/// which to display a particular atom.
Vec3 getAtomDefaultColor(DuMM::AtomIndex atomIx) const;

/**@}**/

/** @name               Define clusters and bodies
Methods in this group control the grouping of atoms into rigid clusters and
the placement of such clusters onto the rigid bodies of the underlying
multibody system.

Note: we use the term "station" to refer to a fixed location with respect
to a cluster or body frame, that is, a point "stationary" on that cluster
or body. **/
/**@{**/

/** Create an empty Cluster (rigid group of atoms). The Cluster index number is
returned; you don't get to pick your own. The name is just for display; you 
must use the index to reference the Cluster. Every Cluster has its own 
reference frame C. **/
DuMM::ClusterIndex createCluster(const char* clusterName);

/** Place an existing atom at a particular station in the local frame of a 
Cluster. It is fine for an atom to be in more than one Cluster as long as only
one of them ends up attached to a body. **/
// EU COMMENT BEGIN
void placeAtomInCluster(DuMM::AtomIndex atomIx, DuMM::ClusterIndex clusterIx,
                        const Vec3& station);
// EU BEGIN
//void placeAtomInCluster(DuMM::AtomIndex atomIx, DuMM::ClusterIndex clusterIx,
//                        Vec3 station);
// EU END

/** Place a Cluster (the child) in another Cluster (the parent). The child's
local frame C is placed at a given Transform with respect to the parent's frame
P. All the atoms in the child Cluster maintain their relative positioning. **/
void placeClusterInCluster(DuMM::ClusterIndex childClusterIndex, 
                           DuMM::ClusterIndex parentClusterIndex, 
                           const Transform& X_PC);

/** Calculate the composite mass properties of a Cluster, either in its own 
reference frame C or in reference frame B with the Cluster placed relative to 
B using the indicated Transform X_BC. **/
MassProperties calcClusterMassProperties(DuMM::ClusterIndex clusterIx, 
                                         const Transform& X_BC = Transform()) const;

/** Place a Cluster's local frame C at a particular location and orientation
with respect to a MobilizedBody's frame B. All the atoms within the Cluster 
will become fixed to the body while maintaining the same relative positions as
they had in the Cluster. **/
void attachClusterToBody(DuMM::ClusterIndex clusterIx, MobilizedBodyIndex body,
                         const Transform& X_BC = Transform());

/** Place an individual atom at a particular station on a body without an 
intervening cluster. **/
// EU COMMENT BEGIN
void attachAtomToBody(DuMM::AtomIndex atomIx, MobilizedBodyIndex body, 
                      const Vec3& station = Vec3(0));
// EU BEGIN
//void attachAtomToBody(DuMM::AtomIndex atomIx, MobilizedBodyIndex body, 
//                      Vec3 station = Vec3(0));
// EU END

/** Find the MobilizedBody on which a Cluster has been fixed in place. **/
MobilizedBodyIndex getClusterBody(DuMM::ClusterIndex clusterIx) const;

/** Find where on its body a Cluster has been placed by returning the 
transform X_BC giving the orientation and position of Cluster frame C in body
frame B. **/
Transform getClusterPlacementOnBody(DuMM::ClusterIndex clusterIx) const;

/** Find where on parent cluster P a child cluster C has been placed, by 
returning the transform X_PC. **/
Transform getClusterPlacementInCluster(DuMM::ClusterIndex childClusterIndex, 
                                       DuMM::ClusterIndex parentClusterIndex) const;
/**@}**/


/** @name       Atom/bond exclusion methods (advanced users only!)
Methods in this group give you fine control over which atoms in the defined
molecules are actually used in the calculation of nonbonded forces, and which
bonds are included in the calculation of bonded forces. By default, any atom 
or bond is included if it can contribute to forces that affect motion.

@warning
You can cause strange effects with these methods if you're not careful -- for
example excluding atoms while using GBSA will cause odd surface area 
calculations to be performed. Don't use these methods unless you know exactly
what you're doing!

During realizeTopology() %DuMM studies the defined molecules and force field
terms to determine (for example) which bonds cross bodies and which atoms
can generate forces. Then the list of included atoms and bonds is finalized 
based on the instructions you give by calling these methods. **/
/**@{**/

/** Clear the list of atoms to be included in nonbonded force 
calculations.\ This is usually called prior to
including selected groups of atoms. Any previous calls to include or 
exclude atoms from the nonbonded list will be forgotten. **/
void clearIncludedNonbondAtomList();

/** Clear the list of bonds to be included in bonded force calculations.\ This
is usually called prior to including selected bonds. Any previous calls to 
include or exclude bonds will be forgotten. **/
void clearIncludedBondList();

/** Restore the included nonbond atoms list to its %DuMM-selected default.\ This
will include all the atoms that can affect interbody nonbonded forces,
including Coulomb, Van der Waals, and solvent (GBSA) forces. Any previous 
calls to include or exclude nonbond atoms will be forgotten. **/
void resetIncludedNonbondAtomListToDefault();

/** Request that %DuMM set the included bonds list to all bonds that can
generate forces that affect motion, meaning bonds whose bonded force terms
have a non-zero amplitude and that involve atoms from two or more bodies. Any 
previous calls to include or exclude bonds will be forgotten. **/
void resetIncludedBondListToDefault();

/** Add the given atom to the included nonbond atom list if it isn't already 
there.\ This does not include bonded terms that involve this atom; you have to
request that explicitly with includeAllInterbodyBondsForOneAtom(). 
@see includeAllNonbondAtomsForOneBody() **/
void includeNonbondAtom(DuMM::AtomIndex atom);

/** Add to the included nonbond atom list all the atoms that are fixed to the
given body. This produces the same result as if includeNonbondAtom() were called 
on each of this body's atoms individually. Nothing happens if there are no 
atoms on the given body. @see includeNonbondAtom() **/
void includeAllNonbondAtomsForOneBody(MobilizedBodyIndex mobod);

/** Given an atom, include all the bonded force terms that (a) include
this atom, and (b) involve a body other than the one to which this atom
is fixed.\ This does not include the atom in nonbonded force calculations;
you have to request that explicitly with includeNonbondAtom(). **/
void includeAllInterbodyBondsForOneAtom(DuMM::AtomIndex atom);

/** Given two atoms that may appear together in some bonded force term,
include all the bonded terms that involve both of them. Note that even if
these atoms are on the same body, some of the bonded terms involving them
might still contain atoms which are on other bodies. An alternate signature
allows you to specify the two atoms by supplying a Bond instead. **/
void includeAllInterbodyBondsWithBothAtoms(DuMM::AtomIndex atom1, 
                                           DuMM::AtomIndex atom2);

/** Given a bond (that is, a 1-2 connection between atoms), extract the
two atoms from it and then make sure all interbody bond terms that involve 
both of those atoms (in any position) are included. See the other 
signature of this method for details. **/
void includeAllInterbodyBondsWithBothAtoms(DuMM::BondIndex bond);

/** Include any bonds that are needed to calculate interbody bonded forces that
involve any atoms fixed to the given mobilized body.\ This will involve a few 
atoms on neighboring bodies. This does not include the involved atoms in the
nonbonded atom list; you have to do that explicitly. 
@see includeAllInterbodyBondsBetweenTwoBodies() **/
void includeAllInterbodyBondsForOneBody(MobilizedBodyIndex mobod);

/** Include any bonds that are needed to calculate interbody bonded forces that
act between the given pair of mobilized bodies.\ Any bonded term that involves
both an atom from \a mobod1 and an atom from \a mobod2 will be included.
Nothing will happen if there is no bond that connects these two bodies. 
@see includeAllInterbodyBondsForOneBody **/
void includeAllInterbodyBondsBetweenTwoBodies
   (MobilizedBodyIndex mobod1, MobilizedBodyIndex mobod2);

/** Atoms to be included in force calculations are numbered from 
0 to getNumIncludedAtoms()-1 after realizeTopology() has been called. Use
DuMM::IncludedAtomIndex to refer to included atoms by their sequential
index numbers. **/
int getNumIncludedAtoms() const;
/** Given a DuMM::IncludedAtomIndex, return the corresponding 
DuMM::AtomIndex.\ You must already have called realizeTopology(). **/
DuMM::AtomIndex getAtomIndexOfIncludedAtom
   (DuMM::IncludedAtomIndex incAtomIx) const;

/** A subset of the included atoms are used in nonbond calculations.\ These
are numbered from 0 to getNumNonbondAtoms()-1 after realizeTopology() has been
called. Use the unique integer type DuMM::NonbondAtomIndex to refer to 
included nonbond atoms by their sequential index numbers. Note that some or
all of these nonbond atoms may also be involved in bonded force 
calculations. **/
int getNumNonbondAtoms() const;
/** Given a DuMM::NonbondAtomIndex, return the corresponding 
DuMM::IncludedAtomIndex.\ You must already have called realizeTopology().
See getAtomIndexOfNonbondAtom() if you want its DuMM::AtomIndex instead. **/
DuMM::IncludedAtomIndex getIncludedAtomIndexOfNonbondAtom
   (DuMM::NonbondAtomIndex nonbondAtomIx) const;
/** Given a DuMM::NonbondAtomIndex, return the corresponding 
DuMM::AtomIndex.\ You must already have called realizeTopology(). 
See getIncludedAtomIndexOfNonbondAtom() if you want its 
DuMM::IncludedAtomIndex instead. **/
DuMM::AtomIndex getAtomIndexOfNonbondAtom
   (DuMM::NonbondAtomIndex nonbondAtomIx) const
{   return getAtomIndexOfIncludedAtom
            (getIncludedAtomIndexOfNonbondAtom(nonbondAtomIx)); }
/**@}**/


    // DEFINE FORCE FIELD PARAMETERS
    
/** @name                 Define atom categories
An AtomClass is used to collect together a set of properties which are expected
to be shared by many individual atoms. The properties are: the element (as 
atomic number), expected valence, and van der Waals parameters. Charge is not 
included in AtomClass but in a second more detailed classification level called
ChargedAtomType. **/
/**@{**/

/** Define a new atom class for this force field, for identifying atoms of a 
particular element, number of bonds, and van der Waals parameters. You must
assign a unique index number and name, either of which can be used to 
reference this AtomClass subsequently. Typically the number and name
are defined by the force field you are using.
@param[in]      atomClassIx
    A unique integer index identifying this AtomClass.
@param[in]      atomClassName
    A unique string identifying this AtomClass.
@param[in]      atomicNumber
    The kind of element associated with this AtomClass, that is, the number of
    protons in the nucleus of atoms that are members of this class.
@param[in]      expectedValence
    The number of bonds expected for atoms in this class; this is a defining
    feature of the atom class so must be satisfied exactly. For example, this
    would be 4 for a class representing a tetrahedral carbon atom.
@param[in]      vdwRadiusInNm
    Van der Waals radius is given as Rmin, the \e radius (in nm, \e not A!) at
    which the energy well minimum is seen (actually it is 1/2 the distance 
    between atom centers for a pair of atoms of this class). This is \e not 
    Sigma, which we define as the \e radius (half distance) at which the 
    energy crosses zero, that is, a little closer together than when the 
    energy well is at maximum depth. See warning below.
@param[in]      vdwWellDepthInKJ
    This is the minimum energy value for the van der Waals force as a positive
    number in kJ/mol (\e not kcal/mol!). This is the value used in the 
    definition of the \a vdwRadiusInNm parameter.
    
@warning It is also common in force fields that van der Waals size terms are 
given by \e diameter rather than \e radius -- be sure you know what convention
was used by your source so that you can convert to radius here. To convert for
a Lennard-Jones force: Rmin = 2^(1/6)Sigma. The radius is in nm, the well 
depth in kJ/mol.

@see defineAtomClass_KA() to supply values in Angstroms and kcals. **/
void defineAtomClass(DuMM::AtomClassIndex     atomClassIx, 
                     const char*              atomClassName,
                     int  atomicNumber,  int  expectedValence,
                     Real vdwRadiusInNm, Real vdwWellDepthInKJ) 
{   defineIncompleteAtomClass(atomClassIx, atomClassName, atomicNumber, 
                              expectedValence);
    setAtomClassVdwParameters(atomClassIx, vdwRadiusInNm, vdwWellDepthInKJ); }

/** Same routine as defineAtomClass() but in Kcal/Angstrom (KA) unit system, 
that is, radius (still not sigma) is in Angstroms, and well depth in kcal/mol.
**/
void defineAtomClass_KA(DuMM::AtomClassIndex atomClassIx, 
                        const char* atomClassName,
                        int element, int valence,
                        Real vdwRadiusInAng, Real vdwWellDepthInKcal)
{
    defineAtomClass(atomClassIx, atomClassName, element, valence,
        vdwRadiusInAng*DuMM::Ang2Nm, vdwWellDepthInKcal*DuMM::Kcal2KJ);
}

/** Obsolete method -- use the other signature. Same routine as 
defineAtomClass_KA() but for backwards compatibility and compactness of 
expression in some contexts, this one accepts an integer for the class index 
rather than requiring the safer DuMM::AtomClassIndex type. **/
void defineAtomClass_KA(int atomClassIx, const char* atomClassName,
                        int element, int valence,
                        Real vdwRadiusInAng, Real vdwWellDepthInKcal)
{
    defineAtomClass_KA((DuMM::AtomClassIndex)atomClassIx, atomClassName, 
                       element, valence, vdwRadiusInAng, vdwWellDepthInKcal);
}

/** Check whether an atom class has been defined using this index. **/
bool hasAtomClass(DuMM::AtomClassIndex) const;
/** Check whether an atom class has been defined using this name. **/
bool hasAtomClass(const String& atomClassName) const;
/** Obtain the atom class index corresponding to this atom class name. **/
DuMM::AtomClassIndex getAtomClassIndex(const String& atomClassName) const;
/** Obtain an atom class index that is numerically larger than the largest
currently-defined atom class index. **/
DuMM::AtomClassIndex getNextUnusedAtomClassIndex() const;
/** Get the index number of the atom class associated with this atom. **/	
DuMM::AtomClassIndex getAtomClassIndex(DuMM::AtomIndex atomIx) const;
/** Get the van der Waals radius shared by all atoms that belong to the 
indicated atom class.\ See comments for this group for a precise definition
of what this means; there are ambiguities so don't assume you already know. **/
Real getVdwRadius(DuMM::AtomClassIndex atomClassIx) const;
/** Get the van der Waals energy well depth shared by all atoms that belong to
the indicated atom class.\ See comments for this group for a precise definition
of what this means; there are ambiguities so don't assume you already know. **/
Real getVdwWellDepth(DuMM::AtomClassIndex atomClassIx) const;

/** Define a new ChargedAtomType for this force field, for identifying atoms
of a particular AtomClass that have a particular partial charge. You must
assign a unique index number and name, either of which can be used to 
reference this ChargedAtomType subsequently. Typically the number and name
are defined by the force field you are using.
@param[in]      atomTypeIx
    A unique integer index identifying this ChargedAtomType.
@param[in]      atomTypeName
    A unique string identifying this ChargedAtomType.
@param[in]      atomClassIx
    The index number of a previously defined AtomClass.
@param[in]      partialChargeInE
    A partial charge in units of e (charge on a proton); this value is the
    same in both the MD and KA units system. **/
void defineChargedAtomType(DuMM::ChargedAtomTypeIndex atomTypeIx, 
                           const char*                atomTypeName,
                           DuMM::AtomClassIndex       atomClassIx, 
                           Real                       partialChargeInE) 
{   defineIncompleteChargedAtomType(atomTypeIx, atomTypeName, atomClassIx);
    setChargedAtomTypeCharge(atomTypeIx, partialChargeInE); }

/** This is an alternate name for defineChargedAtomType() but since partial
charge uses the same unit (e) in both the MD and KA unit systems, this is 
identical to defineChargedAtomType(). **/
void defineChargedAtomType_KA(DuMM::ChargedAtomTypeIndex    atomTypeIx, 
                              const char*                   atomTypeName,
                              DuMM::AtomClassIndex          atomClassIx, 
                              Real                          partialChargeInE)
{   defineChargedAtomType(atomTypeIx, atomTypeName, atomClassIx, 
                          partialChargeInE); } // easy!
/** Obsolete method -- use the other signature. Same routine as 
defineChargedAtomType_KA() but for backwards compatibility and compactness of 
expression in some contexts, this one accepts an integer for the charged atom
type index and the atom class index rather than requiring the safer 
DuMM::ChargedAtomTypeIndex and DuMM::AtomClassIndex types. **/
void defineChargedAtomType_KA(int atomTypeIx, const char* atomTypeName,
                              int atomClassIx, Real partialChargeInE) 
{   defineChargedAtomType_KA((DuMM::ChargedAtomTypeIndex)atomTypeIx, 
                             atomTypeName, (DuMM::AtomClassIndex)atomClassIx, 
                             partialChargeInE); }

/** Check whether a charged atom type has been defined using this index. **/
bool hasChargedAtomType(DuMM::ChargedAtomTypeIndex) const;
/** Check whether a charged atom type has been defined using this name. **/
bool hasChargedAtomType(const String& chargedTypeName) const;
/** Obtain the charged atom type index corresponding to this charged atom
type name. **/
DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex(const String& chargedTypeName) const; // TODO
/** Obtain a charged atom type index that is numerically larger than the 
largest currently-defined charged atom type index. **/
DuMM::ChargedAtomTypeIndex getNextUnusedChargedAtomTypeIndex() const; // TODO

/**@}**/


/** @name Control force field nonbonded behavior in special circumstances
These methods permit setting overall force field behavior in special
circumstances, including van der Waals mixing behavior for dissimilar
atom pairs and scaling of non-bonded terms for closely-bonded atoms. **/
/**@{**/

/** Obtain a human-readable name for one of our van der Waals mixing rules. **/
const char* getVdwMixingRuleName(VdwMixingRule) const;

/** Set the van der Waals mixing rule -- our default is Waldman-Hagler. **/
void setVdwMixingRule(VdwMixingRule);
/** Get the van der Waals mixing rule currently in effect.
@see getVdwMixingRuleName() for a human-readable version. **/
VdwMixingRule getVdwMixingRule() const;

void setVdw12ScaleFactor(Real); ///< default 0
void setVdw13ScaleFactor(Real); ///< default 0
void setVdw14ScaleFactor(Real); ///< default 1
void setVdw15ScaleFactor(Real); ///< default 1

void setCoulomb12ScaleFactor(Real); ///< default 0
void setCoulomb13ScaleFactor(Real); ///< default 0
void setCoulomb14ScaleFactor(Real); ///< default 1
void setCoulomb15ScaleFactor(Real); ///< default 1
/**@}**/


/** @name   Tinker biotypes and pre-defined force field parameter sets
DuMM understands Tinker-format parameter files that can be used to load 
in a whole force field description. This requires assigning Tinker 
Biotypes to the atoms in the molecule. **/
/**@{**/

/** Use Amber99 force field parameters. This duplicates the Tinker Amber99 
parameter set in pre-built code, so you don't need to load the parameters from
a file. **/
void loadAmber99Parameters();

/** Load force field parameters from a TINKER format force field parameter 
file.\ (Only the Amber99 force field has been tested.) **/
void populateFromTinkerParameterFile(
    std::istream& ///< input stream for TINKER format force field parameter file
    );

/** Associate a Tinker Biotype with a ChargedAtomType in this subsystem. A 
Biotype association is required for every Biotype found on atoms in the current
System. **/
void setBiotypeChargedAtomType(
    DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, ///< Prexisting charged atom type index in this subsystem
    BiotypeIndex biotypeIx ///< Preexisting BiotypeIndex defined in the Biotype class
    );

/** Get charged atom type index in this force field associated with a 
particular Biotype. **/
DuMM::ChargedAtomTypeIndex getBiotypeChargedAtomType(BiotypeIndex biotypeIx) const;

/** Generate C++ code from the current contents of this DuMM force 
field object. **/
std::ostream& generateBiotypeChargedAtomTypeSelfCode(std::ostream& os) const;
/**@}**/


/** @name                   Bond stretch terms
Bond stretch parameters (between 2 atom classes). You can use the standard,
built-in functional form (a harmonic) or define your own. **/
/**@{**/

/** Define a harmonic bond stretch term between two atom classes using the
built-in functional form. You are only allowed to specify a built-in term
once per pair of classes, regardless of the order they are listed. It is
permitted to redefine the same term as long as the parameters are the same
both times, otherwise this method will throw an exception if (class1,class2) or 
(class2,class1) has already been assigned.
Stiffness (energy per length^2) must be given in (kJ/mol)/nm^2, and nominal
bond length is in nanometers. CAUTION: Energy is kx^2 using this definition,
while force is -2kx; note factor of 2 in force -- conventions vary.

@param[in] class1, class2        The atom class pair to which this term applies; 
                                 order doesn't matter.
@param[in] stiffnessInKJperNmSq  The bond stiffness in kilojoules per nm^2.
@param[in] nominalLengthInNm     Bond length at which energy and force are zero, in nm.

@see defineBondStretch_KA() for kcal-angstrom units
@see defineCustomBondStretch() to define your own functional form **/
void defineBondStretch(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
                        Real stiffnessInKJperNmSq, Real nominalLengthInNm);

/** Same as defineBondStretch() except that for convenience this takes stiffness 
in (kcal/mol)/A^2, and nominal length is in A (angstroms). Note that these are
immediately converted to our standard MD units internally.

@see defineBondStretch() for a complete description

@param[in] class1, class2           The atom class pair to which this term applies; order doesn't matter.
@param[in] stiffnessInKcalPerAngSq  The bond stiffness in kilocalories per Angstrom^2.
@param[in] nominalLengthInAng       Bond length at which energy and force are zero, in Angstrom.
**/
void defineBondStretch_KA(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
                            Real stiffnessInKcalPerAngSq, Real nominalLengthInAng)
{
    defineBondStretch(class1, class2, 
                        stiffnessInKcalPerAngSq * DuMM::Kcal2KJ/square(DuMM::Ang2Nm),
                        nominalLengthInAng      * DuMM::Ang2Nm);
}
/** Same as defineBondStretch_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondStretch_KA(int class1, int class2,
                          Real stiffnessInKcalPerAngSq, Real nominalLengthInAng)
{   defineBondStretch_KA((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2,
                         stiffnessInKcalPerAngSq, nominalLengthInAng); }

/** Define a custom bond stretch term to be applied to the indicated pair of atom 
classes whenever they are found in a 1-2 bond. Pass in a pointer to a newly-
allocated CustomBondStretch object; the subsystem takes over ownership of the
object so it should NOT be deleted by the caller.

@param[in] class1, class2    The atom class pair to which this term applies; order doesn't matter.
@param[in] bondStretchTerm   Pointer to a heap-allocated object of a concrete type
                             derived from DuMM::CustomBondStretch.

@see DuMM::CustomBondStretch
@see defineBondStretch() to use the built-in harmonic functional form **/
void defineCustomBondStretch(DuMM::AtomClassIndex       class1, 
                                DuMM::AtomClassIndex       class2,
                                DuMM::CustomBondStretch*   bondStretchTerm);
/**@}**/



/** @name                       Bond bending terms
Bond bending parameters (for 3 atom classes bonded 1-2-3). You can use the 
standard, built-in functional form (a harmonic based on the 1-2-3 angle) or 
define your own. **/
/**@{**/

/** Define a harmonic bond bending term applying to a triple of atom classes 
using the built-in functional form. You are only allowed to specify a built-in 
term once per class triple, regardless of the order they are listed. It is
permitted to redefine the same term as long as the parameters are the same
both times, otherwise this method will throw an exception if (class1,class2,
class3) or (class3,class2,class1) has already been assigned. Stiffness (energy
per angle^2) must be given in (kJ/mol)/rad^2, and nominal bond angle is in 
<em>degrees</em>(!). Note that the nominal angle is in degrees while the 
stiffness is in radians. Odd, I know, but that seems to be how it's done!
CAUTION: Energy is ka^2 using this definition, while torque is -2ka; note 
factor of 2 in torque -- conventions vary.

@param[in] class1, class2, class3    The atom class triple to which this term applies; 
                                     1-2-3 or 3-2-1 order doesn't matter.
@param[in] stiffnessInKJPerRadSq     The bond stiffness in kilojoules per radian^2.
@param[in] nominalAngleInDeg         Bond angle at which energy and force are zero, in <em>degrees</em>(!).

@see defineBondBend_KA() for kcal-angstrom units
@see defineCustomBondBend() to define your own functional form **/
void defineBondBend(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3,
                    Real stiffnessInKJPerRadSq, Real nominalAngleInDeg);

/** Same as defineBondBend() except that for convenience this takes stiffness 
in (kcal/mol)/rad^2 (nominal angle is still in degrees). Note that the stiffness
is immediately converted to our standard MD units internally.

@see defineBondBend() for a complete description

@param[in] class1, class2, class3    The atom class triple to which this term applies; 
                                     1-2-3 or 3-2-1 order doesn't matter.
@param[in] stiffnessInKcalPerRadSq   The bond stiffness in kilocalories per radian^2.
@param[in] nominalAngleInDeg         Bond angle at which energy and force are zero, in <em>degrees</em>(!).
**/
void defineBondBend_KA(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3,
                        Real stiffnessInKcalPerRadSq, Real nominalAngleInDeg) 
{
    defineBondBend(class1,class2,class3,
                    stiffnessInKcalPerRadSq * DuMM::Kcal2KJ,
                    nominalAngleInDeg);
}
/** Same as defineBondBend_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondBend_KA(int class1, int class2, int class3,
                       Real stiffnessInKcalPerRadSq, Real nominalAngleInDeg)
{   defineBondBend_KA((DuMM::AtomClassIndex)class1,(DuMM::AtomClassIndex)class2,(DuMM::AtomClassIndex)class3,
                      stiffnessInKcalPerRadSq,nominalAngleInDeg); }

/** Define a custom bond bend term to be applied to the indicated triple of atom 
classes whenever they are found in a 1-2-3 bonded sequence. Pass in a pointer to a newly-
allocated CustomBondBend object; the subsystem takes over ownership of the
object so it should NOT be deleted by the caller.

@param[in] class1, class2, class3    The atom class triple to which this term applies; 
                                     1-2-3 or 3-2-1 order doesn't matter.
@param[in] bondBendTerm              Pointer to a heap-allocated object of a concrete type
                                     derived from DuMM::CustomBondBend.

@see DuMM::CustomBondBend
@see defineBondBend() to use the built-in harmonic functional form **/
void defineCustomBondBend(DuMM::AtomClassIndex      class1, 
                            DuMM::AtomClassIndex      class2, 
                            DuMM::AtomClassIndex      class3,
                            DuMM::CustomBondBend*     bondBendTerm);
/**@}**/



/** @name                   Bond torsion terms
Bond torsion (dihedral) parameters (for 4 atom classes bonded 1-2-3-4). You can
use the standard, built-in functional form (a combination of sinusoids) or 
define your own. Bond torsion terms produce energy and a scalar torque that 
are dependent only on the rotation angle about the 2-3 bond in the specified 
atom class sequence.

\b Theory

Label the four bonded atoms r-x-y-s. Rotation occurs about the axis v=y-x, that
is, a vector from x to y. We define a torsion angle theta using the "polymer 
convention" rather than the IUPAC one which is 180 degrees different. Ours is 
like this:
<pre>
            r                         r      s
  theta=0    \             theta=180   \    / 
              x--y                      x--y
                  \
                   s
</pre>
The sign convention is the same for IUPAC and polymer: A positive angle is 
defined by considering r-x fixed in space. Then using the right hand rule 
around v (that is, thumb points from x to y) a positive rotation rotates y->s 
in the direction of your fingers.

We use a periodic energy function like this:
<pre>
      E(theta) = sum E_n(1 + cos(n*theta - theta0_n))
</pre>
where n is the periodicity, E_n is the amplitude (kJ/mol) for term n, and 
theta0_n is the phase offset for term n. The torque term (applied about the 
v axis) is then
<pre>
      T(theta) = -[sum -n*E_n*sin(n*theta - theta0_n)]
</pre>

Note that in the functional forms above the amplitude paramter E_n could be 
considered a "half amplitude" since the total excursion of the energy and 
torque functions is 2*E_n. When entering this  parameter, be sure you 
understand whether the data you have represents the (half) amplitude as above
(most common) or the full excursion range of the function, in which case you 
should provide only 1/2 the excursion as the value for DuMM's amplitude 
parameter.

For a custom torsion bond, DuMM will provide the angle theta to the
user-written energy and torque routines, expecting a scalar value E(theta) and
T(theta) back as above. However, the functional form is specified by the user
in that case, and it is the user's responsibility to ensure that the returned
torque is the negative gradient of the energy with respect to the angle 
parameter theta.

In either case, DuMM takes care of translating the returned scalar torque
into an appropriate mobility torque when possible, or an appropriate set
of forces on each of the atoms r-x-y-s such that the desired pure torque is
realized as the resultant of the forces. **/
/**@{**/

/** Define bond torsion terms applying to a quadruple of atom classes using the
built-in functional form, a combination of sinusoids with different periods
(see full discussion in the group header for bond torsion terms). You are only 
allowed to specify a built-in term of a given period once per class quadruple, 
whether specified 1-2-3-4 or 4-3-2-1. It is permitted to redefine the same term
as long as the parameters are the same both times, otherwise this method will 
throw an exception if (class1,class2,class3,class4) or (class4,class3,class2,
class1) has already been assigned. The amplitude of the energy functions must 
be given in (kJ/mol), the period is an integer which will multiply the dihedral
angle, and the phase (offset) angle is given in degrees (not radians). 
CAUTION: various conventions exist for specifying the parameters for bond 
torsion terms. Pay particular attention to the amplitude -- in our convention
it is really a half-amplitude since the sinusoids range from -1 to 1 making the
full energy and torque "excursion" twice the amplitude.

@param[in] class1, class2, class3, class4    
                             The atom class quadruple to which this term applies; 
                             1-2-3-4 or 4-3-2-1 order doesn't matter.
@param[in] periodicity       The periodicity (dependence on input angle theta) of the sinusoidal term.
@param[in] ampInKJ           The amplitude (half-range) of the energy function, in kilojoules.
@param[in] phaseInDegrees    Angle at which periodicity*theta yields maximum energy (2*\p ampInKJ).

@see defineBondTorsion_KA()    for kcal-angstrom units
@see defineCustomBondTorsion() to define your own functional form **/
void defineBondTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity, Real ampInKJ, Real phaseInDegrees);
/** Same as defineBondTorsion() but permits two torsion terms (with different 
periods) to be specified simultaneously. **/
void defineBondTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,  
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees);
/** Same as defineBondTorsion() but permits three torsion terms (with different
periods) to be specified simultaneously. **/
void defineBondTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees);

/** Same as defineBondTorsion() but permits four torsion terms (with different
periods) to be specified simultaneously. **/
void defineBondTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees,
    int periodicity4, Real amp4InKJ, Real phase4InDegrees);

/** Same as defineBondTorsion_KA() but takes five integer class arguments for 
backwards compatibility. **/
void defineBondTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees,
    int periodicity4, Real amp4InKcal, Real phase4InDegrees,
    int periodicity5, Real amp5InKcal, Real phase5InDegrees);

/** Same as defineBondTorsion() but permits takes the amplitude in kcal/mol 
(but note that this is converted immediately to our MD unit system of 
kJ/mol). **/
void defineBondTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees)
{ 
    defineBondTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees);
}
/** Same as defineBondTorsion_KA() but permits two torsion terms (with 
different periods) to be specified simultaneously. **/
void defineBondTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,  
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees)
{
    defineBondTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees);
}
/** Same as defineBondTorsion_KA() but permits three torsion terms (with 
different periods) to be specified simultaneously. **/
void defineBondTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees)
{
    defineBondTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                        periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees);
}
/** Same as defineBondTorsion_KA() but permits three torsion terms (with 
different periods) to be specified simultaneously. **/
void defineBondTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees,
    int periodicity4, Real amp4InKcal, Real phase4InDegrees)
{
    defineBondTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                        periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees,
                        periodicity4, amp4InKcal * DuMM::Kcal2KJ, phase4InDegrees);
}

/** Same as defineBondTorsion_KA() but permits three torsion terms (with 
different periods) to be specified simultaneously. **/
void defineBondTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees,
    int periodicity4, Real amp4InKcal, Real phase4InDegrees,
    int periodicity5, Real amp5InKcal, Real phase5InDegrees)
{
    defineBondTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                        periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees,
                        periodicity4, amp4InKcal * DuMM::Kcal2KJ, phase4InDegrees,
                        periodicity5, amp5InKcal * DuMM::Kcal2KJ, phase5InDegrees);
}


/** Same as defineBondTorsion_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondTorsion_KA
    (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees)
{
    defineBondTorsion_KA
            ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
            periodicity1, amp1InKcal, phase1InDegrees);
}
/** Same as defineBondTorsion_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondTorsion_KA
    (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees)
{
    defineBondTorsion_KA
            ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
            periodicity1, amp1InKcal, phase1InDegrees,
            periodicity2, amp2InKcal, phase2InDegrees);
}
/** Same as defineBondTorsion_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondTorsion_KA
    (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees)
{
    defineBondTorsion_KA
            ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
            periodicity1, amp1InKcal, phase1InDegrees,
            periodicity2, amp2InKcal, phase2InDegrees,
            periodicity3, amp3InKcal, phase3InDegrees);
}


/** Same as defineBondTorsion_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondTorsion_KA
    (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees,
    int periodicity4, Real amp4InKcal, Real phase4InDegrees)
{
    defineBondTorsion_KA
            ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
            periodicity1, amp1InKcal, phase1InDegrees,
            periodicity2, amp2InKcal, phase2InDegrees,
            periodicity3, amp3InKcal, phase3InDegrees,
            periodicity4, amp4InKcal, phase4InDegrees);
}

/** Same as defineBondTorsion_KA() but takes integer class arguments for 
backwards compatibility. **/
void defineBondTorsion_KA
    (int class1, int class2, int class3, int class4, 
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees,
    int periodicity4, Real amp4InKcal, Real phase4InDegrees,
    int periodicity5, Real amp5InKcal, Real phase5InDegrees)
{
    defineBondTorsion_KA
            ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
            periodicity1, amp1InKcal, phase1InDegrees,
            periodicity2, amp2InKcal, phase2InDegrees,
            periodicity3, amp3InKcal, phase3InDegrees,
            periodicity4, amp4InKcal, phase4InDegrees,
            periodicity5, amp5InKcal, phase5InDegrees);
}

/** Define a custom bond torsion term to be applied to the indicated quadruple 
of atom classes whenever they are found in a 1-2-3-4 bonded sequence. Pass in a
pointer to a newly-allocated CustomBondTorsion object; the subsystem takes over
ownership of the object so it should NOT be deleted by the caller.

@param[in] class1, class2, class3, class4    
                             The atom class quadruple to which this term applies; 
                             1-2-3-4 or 4-3-2-1 order doesn't matter.
@param[in] bondTorsionTerm   Pointer to a heap-allocated object of a concrete 
                             type derived from DuMM::CustomBondTorsion.
@see DuMM::CustomBondTorsion
@see defineBondTorsion() to use the built-in sinusoidal functional form **/
void defineCustomBondTorsion(DuMM::AtomClassIndex       class1, 
                                DuMM::AtomClassIndex       class2, 
                                DuMM::AtomClassIndex       class3,
                                DuMM::AtomClassIndex       class4,
                                DuMM::CustomBondTorsion*   bondTorsionTerm);
/**@}**/



/** @name               Amber-style improper torsions
As with normal torsions, (see defineBondTorsion()), only one term may have
a given periodicity. The amplitudes are in kJ/mol. The third atom is the 
central one to which the other three are bonded; this is not the same in 
reverse order. **/
/**@{**/

/** Provide one torsion term in MD units, using kilojoules for amplitude. **/
void defineAmberImproperTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity, Real ampInKJ, Real phaseInDegrees);
/** Provide two torsion terms in MD units, using kilojoules/mole for 
amplitude. **/
void defineAmberImproperTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees);
/** Provide three torsion terms in MD units, using kilojoules/mole for 
amplitude. **/
void defineAmberImproperTorsion
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees);

/** Provide one torsion term in KA units, using kilocalories/mole for 
amplitude. **/
void defineAmberImproperTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKcal, Real phase1InDegrees)
{
    defineAmberImproperTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees);
}
/** Provide two torsion terms in KA units, using kilocalories/mole for 
amplitude. **/
void defineAmberImproperTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees)
{
    defineAmberImproperTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees);
}
/** Provide three torsion terms in KA units, using kilocalories/mole for 
amplitude. **/
void defineAmberImproperTorsion_KA
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKcal, Real phase1InDegrees,
    int periodicity2, Real amp2InKcal, Real phase2InDegrees,
    int periodicity3, Real amp3InKcal, Real phase3InDegrees)
{
    defineAmberImproperTorsion(class1,class2,class3,class4,
                        periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                        periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                        periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees);
}
/**@}**/


/** @name                   GBSA implicit solvation
These methods are used to set GBSA terms. **/
/**@{**/
void setSolventDielectric(Real);    // typically 80 for water
void setSoluteDielectric(Real);     // typically 1 or 2 for protein
Real getSolventDielectric() const;
Real getSoluteDielectric() const;

void setGbsaIncludeAceApproximation(bool);
void setGbsaIncludeAceApproximationOn()  {setGbsaIncludeAceApproximation(true );}
void setGbsaIncludeAceApproximationOff() {setGbsaIncludeAceApproximation(false);}
/**@}**/

/** @name                   Global scale factors
These <em>non-physical</em> parameters can be used to weaken or disable (or 
magnify) individual force field terms. These are always 1 for correct 
implementation of any force field; other values are primarily useful for 
testing the effects of individual terms on results or performance. Set to 0 to 
disable the corresponding term altogether. **/
/**@{**/
void setVdwGlobalScaleFactor(Real);                 ///< scale all van der Waals terms
void setCoulombGlobalScaleFactor(Real);             ///< scale all Coulomb terms
void setGbsaGlobalScaleFactor(Real);                ///< scale all GBSA terms
void setBondStretchGlobalScaleFactor(Real);         ///< scale all built-in bond stretch terms
void setBondBendGlobalScaleFactor(Real);            ///< scale all built-in bond bending terms
void setBondTorsionGlobalScaleFactor(Real);         ///< scale all built-in bond torsion terms
void setAmberImproperTorsionGlobalScaleFactor(Real);///< scale all improper torsion terms
void setCustomBondStretchGlobalScaleFactor(Real);   ///< scale all custom bond stretch terms
void setCustomBondBendGlobalScaleFactor(Real);      ///< scale all custom bond bending terms
void setCustomBondTorsionGlobalScaleFactor(Real);   ///< scale all custom bond torsion terms

Real getVdwGlobalScaleFactor() const;                   ///< get current scale factor for van der Waals terms
Real getCoulombGlobalScaleFactor() const;               ///< get current scale factor for Coulomb terms
Real getGbsaGlobalScaleFactor() const;                  ///< get current scale factor for GBSA terms
Real getBondStretchGlobalScaleFactor() const;           ///< get current scale factor for built-in bond stretch terms
Real getBondBendGlobalScaleFactor() const;              ///< get current scale factor for built-in bond bending terms
Real getBondTorsionGlobalScaleFactor() const;           ///< get current scale factor for built-in bond torsion terms
Real getAmberImproperTorsionGlobalScaleFactor() const;  ///< get current scale factor for improper torsion terms
Real getCustomBondStretchGlobalScaleFactor() const;     ///< get current scale factor for custom bond stretch terms
Real getCustomBondBendGlobalScaleFactor() const;        ///< get current scale factor for custom bond bending terms
Real getCustomBondTorsionGlobalScaleFactor() const;     ///< get current scale factor for custom bond torsion terms

/** Set all the global scale factors to the same value. This is commonly used 
to turn everything on or off, followed by selectively disabling or enabling 
individual terms. **/
void setAllGlobalScaleFactors(Real s) {
	setVdwGlobalScaleFactor(s);
	setCoulombGlobalScaleFactor(s);
	setGbsaGlobalScaleFactor(s);
	setBondStretchGlobalScaleFactor(s);
	setBondBendGlobalScaleFactor(s);
	setBondTorsionGlobalScaleFactor(s);
	setAmberImproperTorsionGlobalScaleFactor(s);
	setCustomBondStretchGlobalScaleFactor(s);
	setCustomBondBendGlobalScaleFactor(s);
	setCustomBondTorsionGlobalScaleFactor(s);

}


// Added extra functions to customize OpenMM usage (Eliza)
Real getNonbondedCutoff() const;                        ///< get current nonbonded cutoff (nm) used for LJ and Coulomb calculations
void setNonbondedCutoff(Real);                          ///< set nonbonded cutoff (nm) used for LJ and Coulomb calculations
int getNonbondedMethod() const;                         ///< get current nonbonded method used by OpenMM.
void setNonbondedMethod(int);                           ///< set nonbonded nonbonded method used by OpenMM. (0 = nocutoff; 1=cutoffnonperiodic).




/**@}**/




/** @name               Computational options
These methods control how DuMM performs its computations. **/
/**@{**/

/** For debugging, you can ask DuMM to dump some information to std::clog
about its attempts to figure out the best way to use the available
hardware. **/
void setTracing(bool);

/** Enable or disable multithreaded computation (enabled by default). **/
void setUseMultithreadedComputation(bool);
/** Is multithreaded computation enabled? **/
bool getUseMultithreadedComputation() const;

/** Request how many threads we want to use if multithreading is enabled; zero
(the default) means let DuMM choose, which will likely be one thread per 
processor. DuMM will choose to use single threaded code if there is only one 
processor, but requesting 1 thread here explicitly will cause it to use 
multithreaded code but with a single thread (useful for debugging sometimes). **/
void setNumThreadsRequested(int);
/** What was the last value passed to setNumThreadsRequested()? The default is
zero meaning DuMM chooses the number of threads. This doesn't tell you how many
threads are actually in use -- see getNumThreadsInUse() for that. **/
int getNumThreadsRequested() const;

/** Is DuMM using the multithreaded code? This could return true even
if there is just one thread, if you forced it with setNumThreadsToUse(). **/
bool isUsingMultithreadedComputation() const;

/** Find out how many threads DuMM is actually using; will be zero until after
realizeTopology(). **/
int getNumThreadsInUse() const;

/** This determines whether we use OpenMM GPU acceleration if it is available. 
By default, this is set false because OpenMM will compute only to single 
precision. Note that even if you set this flag, we won't use OpenMM unless 
(a) it is installed correctly on your machine, and (b) it can run with GPU 
acceleration. If you want to allow use of the non-accelerated Reference 
Platform provided by OpenMM, use setAllowOpenMMReference(). **/
void setUseOpenMMAcceleration(bool);
/** Return the current setting of the flag set by setUseOpenMMAcceleration(). **/
bool getUseOpenMMAcceleration() const;



/* Customize OpenMMPlugin */

/** This determines whether we use OpenMM for calculating the nonbonded
 energy & forces (van der Waals, Coulomb and GBSA ) or to calculate the
 all the energy & forces terms (included the bonded terms: bonds, bends,
 dihedral and impropers) **/
void setUseOpenMMCalcOnlyNonBonded(bool);
/** Return the current setting of the flag set by setUseOpenMMCalcOnlyNonBonded(). **/
bool getUseOpenMMCalcOnlyNonBonded() const;

// Used for OpenMM integration
    void setUseOpenMMIntegration(bool);
    bool getUseOpenMMIntegration() const;
    void setOpenMMstepsize(float);
    float getOpenMMstepsize() const;
    void setOpenMMtemperature(float);
    float getOpenMMtemperature() const;
    void setOpenMMvelocities(float);
    void OMM_updatePositions(const std::vector<SimTK::Vec3>& positions);
    const std::vector<OpenMM::Vec3>& OMM_getPositions() const;

    SimTK::Vec3 calcAtomLocationInGroundFrameThroughOMM( DuMM::AtomIndex daix ) const;
    void OMM_integrateTrajectory( int steps );
    Real OMM_calcPotentialEnergy() const;
    Real OMM_calcKineticEnergy() const;

/** Return OpennMMPluginIterface pointer**/
OpenMMPluginInterface*  getOpenMMPluginIfc() const;

/** This allows us to use OpenMM even if only the Reference platform is 
available. This is for testing/debugging; one should never use the Reference
platform in production since it will likely be slower than the CPU 
implementation. **/
void setAllowOpenMMReference(bool);
/** Return the current setting of the flag set by setAllowOpenMMReference(). **/
bool getAllowOpenMMReference() const;

/** Return true if DuMM is currently using OpenMM for its computations. **/
bool isUsingOpenMM() const;
/** Return the OpenMM Platform currently in use, or the empty string
if we're not using OpenMM. **/
std::string getOpenMMPlatformInUse() const;
/**@}**/

/** @name Bookkeeping, debugging, and internal-use-only methods
Hopefully you won't need these. **/
/**@{**/

/** How many times has the forcefield been evaluated? **/
long long getForceEvaluationCount() const;

/** Produce an ugly but comprehensive dump of the contents of DuMM's internal
data structures, sent to std::cout (stdout). **/
void dump() const;

/** Generate C++ code to reproduce forceField parameters presently in 
memory. **/
void dumpCForceFieldParameters
   (std::ostream& os, const String& methodName = "loadParameters") const;

/** Load test parameters. **/
void loadTestMoleculeParameters();

/**@}**/

/// OBSOLETE NAME -- just use setTracing().
void setTraceOpenMM(bool shouldTrace) {setTracing(shouldTrace);}
protected:
/** @cond **/ // don't let doxygen see these
SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, ForceSubsystem);
SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, Subsystem);
/** @endcond **/

void defineIncompleteAtomClass(
    DuMM::AtomClassIndex classIx, 
    const char* name, 
    int elementNumber, 
    int valence);

void defineIncompleteAtomClass_KA(
    DuMM::AtomClassIndex classIx, 
    const char* name, 
    int elementNumber, 
    int valence) 
{
    defineIncompleteAtomClass(
        classIx, 
        name, 
        elementNumber, 
        valence
        );
}

    void setAtomClassVdwParameters(DuMM::AtomClassIndex atomClassIx, Real vdwRadiusInNm, Real vdwWellDepthInKJPerMol);
    void setAtomClassVdwParameters_KA(DuMM::AtomClassIndex atomClassIx, Real radiusInAng, Real wellDepthInKcal) {
        setAtomClassVdwParameters(atomClassIx, radiusInAng*DuMM::Ang2Nm, wellDepthInKcal*DuMM::Kcal2KJ);
    }

    bool isValidAtomClass(DuMM::AtomClassIndex) const;
        
    void defineIncompleteChargedAtomType(
        DuMM::ChargedAtomTypeIndex typeIx, 
        const char* name,
        DuMM::AtomClassIndex classIx);

    void defineIncompleteChargedAtomType_KA(
        DuMM::ChargedAtomTypeIndex typeIx, 
        const char* name,
        DuMM::AtomClassIndex classIx) 
    {
        defineIncompleteChargedAtomType(typeIx, name, classIx);
    }

    void setChargedAtomTypeCharge(DuMM::ChargedAtomTypeIndex, Real charge);
    void setChargedAtomTypeCharge_KA(DuMM::ChargedAtomTypeIndex chargedAtomTypeIx, Real charge) {
        setChargedAtomTypeCharge(chargedAtomTypeIx, charge);
    }

private:
    class DuMMForceFieldSubsystemRep& updRep();
    const DuMMForceFieldSubsystemRep& getRep() const;

    friend class MolecularMechanicsSystem;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//		        GMOLMODEL - EXTRA FUNCTIONALITIES
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
public:
        Real CalcFullPotEnergyIncludingRigidBodies ( const State& state ) const;

//------------------------------------------------------------------------------



};

/** This class is just a DuMMForceFieldSubsystem for which the constructor 
pre-loads the definitions of the Amber99 force field. **/
class Amber99ForceSubsystem : public DuMMForceFieldSubsystem {
public:
    explicit Amber99ForceSubsystem(MolecularMechanicsSystem& system) 
        : DuMMForceFieldSubsystem(system)
    {
        loadAmber99Parameters();
    }
};


} // namespace SimTK


/**@}**/ // End of MolecularMechanics module

#endif // SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_H_
