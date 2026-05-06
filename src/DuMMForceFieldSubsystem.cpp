/* -------------------------------------------------------------------------- *
 *                             SimTK Molmodel(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Christopher Bruns, Randy Radmer                              *
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
 * Private implementation of DuMMForceFieldSubsystem. Units here are uniformly
 * MD units: nanometers, daltons, picoseconds, with energy in kilojoules/mole.
 * We accept angles from users in degrees, but use only radians internally.
 */

#include "molmodel/internal/DuMMForceFieldSubsystem.h"

#include <cstddef>
#include <utility>

#include "molmodel/internal/common.h"

#include "DuMMForceFieldSubsystemRep.h"
#include "OpenMM.hpp"
#include "units.h"


using namespace SimTK;

constexpr Real EPSILON = 1e-2;

auto angleIsInRange(double angle) -> bool {
    // Check if angle >= MIN_ANGLE - EPSILON AND angle <= MAX_ANGLE + EPSILON.
    // However, if the boundaries themselves have floating-point error (e.g., from calculations),
    // you might need a more complex comparison, but for fixed bounds (-180.0, 180.0),
    // the following is generally acceptable.

    const double MIN_ANGLE = -180.0;
    const double MAX_ANGLE = 180.0;

    return (angle >= MIN_ANGLE - EPSILON) && (angle <= MAX_ANGLE + EPSILON);
}

inline auto almostEqual(Real first, Real second) -> bool {
    return std::fabs(first - second) <= EPSILON * std::max({1.0, std::fabs(first), std::fabs(second)});
}

inline auto anglesAlmostEqual(Real first, Real second) -> bool {
    double diff = std::fmod(first - second + 180.0, 360.0);
    if (diff < 0) {
        diff += 360.0;
    }
    diff -= 180.0;
    return std::fabs(diff) <= EPSILON;
}

////////////////////////////////
// DUMM FORCE FIELD SUBSYSTEM //
////////////////////////////////

auto DuMMForceFieldSubsystem::isInstanceOf(const Subsystem& subsystem) -> bool {
    return DuMMForceFieldSubsystemRep::isA(subsystem.getSubsystemGuts());
}
auto DuMMForceFieldSubsystem::downcast(const Subsystem& subsystem) -> const DuMMForceFieldSubsystem& {
    assert(isInstanceOf(subsystem));
    return reinterpret_cast<const DuMMForceFieldSubsystem&>(subsystem);
}
auto DuMMForceFieldSubsystem::updDowncast(Subsystem& subsystem) -> DuMMForceFieldSubsystem& {
    assert(isInstanceOf(subsystem));
    return reinterpret_cast<DuMMForceFieldSubsystem&>(subsystem);
}
auto DuMMForceFieldSubsystem::getRep() const -> const DuMMForceFieldSubsystemRep& {
    return dynamic_cast<const DuMMForceFieldSubsystemRep&>(getSubsystemGuts());
}
auto DuMMForceFieldSubsystem::updRep() -> DuMMForceFieldSubsystemRep& {
    return dynamic_cast<DuMMForceFieldSubsystemRep&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
DuMMForceFieldSubsystem::DuMMForceFieldSubsystem() {
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
}

DuMMForceFieldSubsystem::DuMMForceFieldSubsystem(MolecularMechanicsSystem& mms) {
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
    mms.setMolecularMechanicsForceSubsystem(*this); // steal ownership
}

auto DuMMForceFieldSubsystem::getAtomClassIndex(DuMM::AtomIndex atomIx) const -> DuMM::AtomClassIndex {
    DuMM::ChargedAtomTypeIndex typeIx = getRep().atoms[atomIx].chargedAtomTypeIndex;
    return getRep().chargedAtomTypes[typeIx].atomClassIx;
}
auto DuMMForceFieldSubsystem::getVdwRadius(DuMM::AtomClassIndex atomClassIx) const -> Real {
    return getRep().atomClasses[atomClassIx].vdwRadius;
}
auto DuMMForceFieldSubsystem::getVdwWellDepth(DuMM::AtomClassIndex atomClassIx) const -> Real {
    return getRep().atomClasses[atomClassIx].vdwWellDepth;
}

/*! <!-- desk_mass_related --> */
auto DuMMForceFieldSubsystem::getAtomMass(DuMM::AtomIndex dAIx) const -> SimTK::mdunits::Mass {
    const DuMMForceFieldSubsystemRep& rep = getRep();
    return rep.getAtomMass(dAIx);
}

/*! <!-- desk_mass_related --> */
void DuMMForceFieldSubsystem::setDuMMAtomMass(SimTK::DuMM::AtomIndex dAIx, SimTK::mdunits::Mass atomicMass) {
    static const char* MethodName = "setDuMMAtomMass";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class1),
                       rep.ApiClassName,
                       MethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);

    rep.setAtomMass(dAIx, atomicMass);
}

void DuMMForceFieldSubsystem::defineIncompleteAtomClass(DuMM::AtomClassIndex atomClassIx,
                                                        const char* atomClassName,
                                                        int elementNumber,
                                                        int valence) {
    static const char* MethodName = "defineIncompleteAtomClass";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Catch nonsense arguments.
    SimTK_APIARGCHECK1(atomClassIx.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "atom class Index %d invalid: must be nonnegative",
                       (int)atomClassIx);
    SimTK_APIARGCHECK1(rep.isValidElement(elementNumber),
                       rep.ApiClassName,
                       MethodName,
                       "element %d invalid: must be a valid atomic number and have an entry here",
                       elementNumber);
    SimTK_APIARGCHECK1(valence >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "expected valence %d invalid: must be nonnegative",
                       valence);

    // Make sure there is a slot available for this atom class.
    if (atomClassIx >= (DuMM::AtomClassIndex)rep.atomClasses.size()) {
        rep.atomClasses.resize(atomClassIx + 1);
    }

    // Make sure this atom class hasn't already been defined.
    SimTK_APIARGCHECK2(!rep.atomClasses[atomClassIx].isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "atom class Index %d is already in use for '%s'",
                       (int)atomClassIx,
                       rep.atomClasses[atomClassIx].name.c_str());

    if (rep.atomClassIndicesByName.find(atomClassName) != rep.atomClassIndicesByName.end()) {
        DuMM::AtomClassIndex oldAtomClassIx = rep.atomClassIndicesByName.find(atomClassName)->second;
        if (oldAtomClassIx != atomClassIx) {
            throw(std::runtime_error(String("Duplicate atom class name: ") + atomClassName));
        }
    }

    rep.insertNewAtomClass(AtomClass(atomClassIx, atomClassName, elementNumber, valence, NaN, NaN));
}

void DuMMForceFieldSubsystem::setAtomClassVdwParameters(DuMM::AtomClassIndex atomClassIx,
                                                        Real vdwRadiusInNm,
                                                        Real vdwWellDepthInKJPerMol) {
    static const char* MethodName = "setAtomClassVdwParameters";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    SimTK_APIARGCHECK1(atomClassIx.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "atom class Index %d invalid: must be nonnegative",
                       (int)atomClassIx);
    SimTK_APIARGCHECK1(vdwRadiusInNm >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "van der Waals radius %g invalid: must be nonnegative",
                       vdwRadiusInNm);
    SimTK_APIARGCHECK1(vdwWellDepthInKJPerMol >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "van der Waals energy well depth %g invalid: must be nonnegative",
                       vdwWellDepthInKJPerMol);

    AtomClass& atomClass = rep.atomClasses[atomClassIx];
    atomClass.vdwRadius = vdwRadiusInNm;
    atomClass.vdwWellDepth = vdwWellDepthInKJPerMol;
}

auto DuMMForceFieldSubsystem::isValidAtomClass(DuMM::AtomClassIndex atomClassIx) const -> bool {
    return getRep().isValidAtomClass(atomClassIx);
}

void DuMMForceFieldSubsystem::defineIncompleteChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex,
                                                              const char* typeName,
                                                              DuMM::AtomClassIndex atomClassIx) {
    static const char* MethodName = "defineIncompleteChargedAtomType";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Check for nonsense arguments.
    SimTK_APIARGCHECK1(chargedAtomTypeIndex.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "charged atom type index %d invalid: must be nonnegative",
                       (int)chargedAtomTypeIndex);
    SimTK_APIARGCHECK1(atomClassIx.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "atom class index %d invalid: must be nonnegative",
                       (int)atomClassIx);
    // partialCharge is a signed quantity

    // Make sure the referenced atom class has already been defined.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(atomClassIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom class %d is undefined",
                       (int)atomClassIx);

    // Make sure there is a slot available for the new chargedAtomType.
    if (chargedAtomTypeIndex >= (int)rep.chargedAtomTypes.size()) {
        rep.chargedAtomTypes.resize(chargedAtomTypeIndex + 1);
    }

    // Check that this slot is not already in use.
    SimTK_APIARGCHECK2(!rep.chargedAtomTypes[chargedAtomTypeIndex].isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "charged atom type index %d is already in use for '%s'",
                       (int)chargedAtomTypeIndex,
                       rep.chargedAtomTypes[chargedAtomTypeIndex].name.c_str());

    rep.insertNewChargedAtomType(ChargedAtomType(chargedAtomTypeIndex, typeName, atomClassIx, NaN));
}

auto DuMMForceFieldSubsystem::hasAtomClass(DuMM::AtomClassIndex atomClassIndex) const -> bool {
    return getRep().hasAtomClass(atomClassIndex);
}
auto DuMMForceFieldSubsystem::hasAtomClass(const String& atomClassName) const -> bool {
    return getRep().hasAtomClass(atomClassName);
}
auto DuMMForceFieldSubsystem::getAtomClassIndex(const String& atomClassName) const -> DuMM::AtomClassIndex {
    return getRep().getAtomClassIndex(atomClassName);
}
auto DuMMForceFieldSubsystem::getNextUnusedAtomClassIndex() const -> DuMM::AtomClassIndex {
    return getRep().getNextUnusedAtomClassIndex();
}

auto DuMMForceFieldSubsystem::hasChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex) const
    -> bool {
    return getRep().hasChargedAtomType(chargedAtomTypeIndex);
}
auto DuMMForceFieldSubsystem::hasChargedAtomType(const String& chargedTypeName) const -> bool {
    return getRep().hasChargedAtomType(chargedTypeName);
}
auto DuMMForceFieldSubsystem::getChargedAtomTypeIndex(const String& chargedTypeName) const
    -> DuMM::ChargedAtomTypeIndex {
    return getRep().getChargedAtomTypeIndex(chargedTypeName);
}
auto DuMMForceFieldSubsystem::getNextUnusedChargedAtomTypeIndex() const -> DuMM::ChargedAtomTypeIndex {
    return getRep().getNextUnusedChargedAtomTypeIndex();
}

void DuMMForceFieldSubsystem::setChargedAtomTypeCharge(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex,
                                                       Real charge) {
    static const char* MethodName = "setChargedAtomTypeCharge";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Check for nonsense arguments.
    SimTK_APIARGCHECK1(chargedAtomTypeIndex.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "charged atom type index %d invalid: must be nonnegative",
                       (int)chargedAtomTypeIndex);

    auto& chargedAtomType = rep.chargedAtomTypes[chargedAtomTypeIndex];
    chargedAtomType.partialCharge = charge;
}

void DuMMForceFieldSubsystem::defineBondStretch(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                Real stiffnessInKJPerNmSq,
                                                Real nominalLengthInNm) {
    static const char* MethodName = "defineBondStretch";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class1),
                       rep.ApiClassName,
                       MethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class2),
                       rep.ApiClassName,
                       MethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(stiffnessInKJPerNmSq >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "stiffness %g is not valid: must be nonnegative",
                       stiffnessInKJPerNmSq);
    SimTK_APIARGCHECK1(nominalLengthInNm >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "nominal length %g is not valid: must be nonnegative",
                       nominalLengthInNm);

    // We canonicalize the key so that the atom class pair has the
    // lower class Index first.
    const AtomClassIndexPair key(class1, class2, true);

    // Attempt to create a new bond stretch entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = rep.bondStretch.emplace(key, key);
    auto& bondStretchEntry = ret.first->second;

    if (bondStretchEntry.hasBuiltinTerm()) {
        SimTK_APIARGCHECK2(
            bondStretchEntry.k == stiffnessInKJPerNmSq && bondStretchEntry.d0 == nominalLengthInNm,
            rep.ApiClassName,
            MethodName,
            "There was already a different built-in bond stretch term for atom class pair (%d,%d); only one "
            "is allowed."
            "\nUse a CustomBondStretch term if you need another term for the same atom class pair.",
            (int)key[0],
            (int)key[1]);
    } else {
        bondStretchEntry.setBuiltinTerm(stiffnessInKJPerNmSq, nominalLengthInNm);
    }
}

void DuMMForceFieldSubsystem::defineCustomBondStretch(DuMM::AtomClassIndex class1,
                                                      DuMM::AtomClassIndex class2,
                                                      DuMM::CustomBondStretch* customBondStretch) {
    SimTK_ASSERT(false, "DuMMForceFieldSubsystem::defineCustomBondStretch is no longer supported.");

    // static const char* MethodName = "defineCustomBondStretch";

    // invalidateSubsystemTopologyCache();

    // DuMMForceFieldSubsystemRep& rep = updRep();

    //     // Watch for nonsense arguments.
    // SimTK_APIARGCHECK1(rep.isValidAtomClass(class1), rep.ApiClassName, MethodName,
    //     "class1=%d which is not a valid atom class Index", (int) class1);
    // SimTK_APIARGCHECK1(rep.isValidAtomClass(class2), rep.ApiClassName, MethodName,
    //     "class2=%d which is not a valid atom class Index", (int) class2);
    // SimTK_APIARGCHECK(customBondStretch, rep.ApiClassName, MethodName,
    //     "CustomBondStretch pointer was null");

    //     // We canonicalize the key so that the atom class pair has the
    //     // lower class Index first.
    // const AtomClassIndexPair key(class1,class2,true);

    //     // Attempt to create a new bond stretch entry containing no valid
    //     // terms. If there was already an entry it will be returned instead
    //     // and no insertion is performed.
    // std::pair<std::map<AtomClassIndexPair,BondStretch>::iterator, bool> ret =
    //   rep.bondStretch.insert(std::pair<AtomClassIndexPair,BondStretch>
    //     (key, BondStretch(key)));

    // BondStretch& bondStretchEntry = ret.first->second;
    // bondStretchEntry.addCustomTerm(customBondStretch);
}

void DuMMForceFieldSubsystem::defineBondBend(DuMM::AtomClassIndex class1,
                                             DuMM::AtomClassIndex class2,
                                             DuMM::AtomClassIndex class3,
                                             Real stiffnessInKJPerRadSq,
                                             Real nominalAngleInDeg) {
    static const char* MethodName = "defineBondBend";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class1),
                       rep.ApiClassName,
                       MethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class2),
                       rep.ApiClassName,
                       MethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class3),
                       rep.ApiClassName,
                       MethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK1(stiffnessInKJPerRadSq >= 0,
                       rep.ApiClassName,
                       MethodName,
                       "stiffness %g is not valid: must be nonnegative",
                       stiffnessInKJPerRadSq);
    SimTK_APIARGCHECK1(0 <= nominalAngleInDeg && nominalAngleInDeg <= 180,
                       rep.ApiClassName,
                       MethodName,
                       "nominal angle %g is not valid: must be between 0 and 180 degrees, inclusive",
                       nominalAngleInDeg);

    // We canonicalize the key so that the first classIndex is no larger than the third.
    const AtomClassIndexTriple key(class1, class2, class3, true);

    // Attempt to create a new bond bend entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = rep.bondBend.insert(std::make_pair(key, BondBend(key)));
    auto& bondBendEntry = ret.first->second;

    if (bondBendEntry.hasBuiltinTerm()) {
        SimTK_APIARGCHECK3(
            bondBendEntry.k == stiffnessInKJPerRadSq
                && bondBendEntry.theta0 == nominalAngleInDeg * DuMM::Deg2Rad,
            rep.ApiClassName,
            MethodName,
            "There was already a different built-in bond bend term for atom class triple (%d,%d,%d); only "
            "one is allowed."
            "\nUse a CustomBondBend term if you need another term for the same atom class triple.",
            (int)key[0],
            (int)key[1],
            (int)key[2]);
    } else {
        bondBendEntry.setBuiltinTerm(stiffnessInKJPerRadSq, nominalAngleInDeg);
    }
}

void DuMMForceFieldSubsystem::defineCustomBondBend(DuMM::AtomClassIndex class1,
                                                   DuMM::AtomClassIndex class2,
                                                   DuMM::AtomClassIndex class3,
                                                   DuMM::CustomBondBend* customBondBend) {
    static const char* MethodName = "defineCustomBondBend";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class1),
                       rep.ApiClassName,
                       MethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class2),
                       rep.ApiClassName,
                       MethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class3),
                       rep.ApiClassName,
                       MethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK(customBondBend, rep.ApiClassName, MethodName, "CustomBondBend pointer was null");

    // We canonicalize the key so that the first classIndex is no larger than the third.
    const AtomClassIndexTriple key(class1, class2, class3, true);

    // Attempt to create a new bond bend entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = rep.bondBend.insert(std::make_pair(key, BondBend(key)));
    auto& bondBendEntry = ret.first->second;
    bondBendEntry.addCustomTerm(customBondBend);
}

//
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions, and then inserts the built in torsion terms
// if they are legitimate.
//
void DuMMForceFieldSubsystemRep::defineAnyTorsion(DuMM::AtomClassIndex class1,
                                                  DuMM::AtomClassIndex class2,
                                                  DuMM::AtomClassIndex class3,
                                                  DuMM::AtomClassIndex class4,
                                                  bool shouldCanonicalizeClassOrder,
                                                  int periodicity1,
                                                  Real amp1InKJ,
                                                  Real phase1InDegrees,
                                                  int periodicity2,
                                                  Real amp2InKJ,
                                                  Real phase2InDegrees,
                                                  int periodicity3,
                                                  Real amp3InKJ,
                                                  Real phase3InDegrees,
                                                  std::map<AtomClassIndexQuad, BondTorsion>& torsionMap,
                                                  const char* CallingMethodName) const {
    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(isValidAtomClass(class1),
                       ApiClassName,
                       CallingMethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(isValidAtomClass(class2),
                       ApiClassName,
                       CallingMethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(isValidAtomClass(class3),
                       ApiClassName,
                       CallingMethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK1(isValidAtomClass(class4),
                       ApiClassName,
                       CallingMethodName,
                       "class4=%d which is not a valid atom class Index",
                       (int)class4);
    SimTK_APIARGCHECK(periodicity1 != -1 || periodicity2 != -1 || periodicity3 != -1,
                      ApiClassName,
                      CallingMethodName,
                      "must be at least one torsion term supplied");

    if (periodicity1 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity1 && periodicity1 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity1);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);*/
        // scf changed 0 to -180 to allow NAST right handed helices

        SimTK_APIARGCHECK1(angleIsInRange(phase1InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle1(%g) is not valid: must be between -180 and 180 degrees, inclusive",
                           phase1InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            (periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity2 && periodicity2 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity2);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase2InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase2InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            periodicity3 != periodicity2,
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity3 && periodicity3 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity3);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase3InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase3InDegrees);
        // (we've already checked for any possible repeats)
    }

    // Canonicalize atom class quad by reversing order if necessary so that the
    // first class Index is numerically no larger than the fourth. Amber improper
    // torsions should not be canonicalized because order matters.
    const AtomClassIndexQuad key(class1, class2, class3, class4, shouldCanonicalizeClassOrder);

    // Attempt to create a new bond torsion entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = torsionMap.insert(std::make_pair(key, BondTorsion(key)));
    auto& bondTorsionEntry = ret.first->second;

    // A new entry or one that just had a custom term in it won't have a built in
    // term so we can load it up and we're done.
    if (!bondTorsionEntry.hasBuiltinTerm()) {
        if (periodicity1 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
        if (periodicity2 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
        if (periodicity3 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
        return;
    }

    // If we get here we have discovered that there is already a built in torsion
    // term present for this atom class quad. We can still insert new terms, and we'll
    // allow duplicates if they are identical.
    if (periodicity1 != -1) {
        const TorsionTerm& term1 = bondTorsionEntry.getTermWithPeriod(periodicity1);
        if (term1.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term1.amplitude, amp1InKJ)
                    && anglesAlmostEqual(term1.theta0 * DuMM::Rad2Deg, phase1InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity1);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
    }
    if (periodicity2 != -1) {
        const TorsionTerm& term2 = bondTorsionEntry.getTermWithPeriod(periodicity2);
        if (term2.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term2.amplitude, amp2InKJ)
                    && anglesAlmostEqual(term2.theta0 * DuMM::Rad2Deg, phase2InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity2);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
    }
    if (periodicity3 != -1) {
        const TorsionTerm& term3 = bondTorsionEntry.getTermWithPeriod(periodicity3);
        if (term3.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term3.amplitude, amp3InKJ)
                    && anglesAlmostEqual(term3.theta0 * DuMM::Rad2Deg, phase3InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity3);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
    }
}

//
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions, and then inserts the built in torsion terms
// if they are legitimate.
//
void DuMMForceFieldSubsystemRep::defineAnyTorsion(DuMM::AtomClassIndex class1,
                                                  DuMM::AtomClassIndex class2,
                                                  DuMM::AtomClassIndex class3,
                                                  DuMM::AtomClassIndex class4,
                                                  bool shouldCanonicalizeClassOrder,
                                                  int periodicity1,
                                                  Real amp1InKJ,
                                                  Real phase1InDegrees,
                                                  int periodicity2,
                                                  Real amp2InKJ,
                                                  Real phase2InDegrees,
                                                  int periodicity3,
                                                  Real amp3InKJ,
                                                  Real phase3InDegrees,
                                                  int periodicity4,
                                                  Real amp4InKJ,
                                                  Real phase4InDegrees,
                                                  std::map<AtomClassIndexQuad, BondTorsion>& torsionMap,
                                                  const char* CallingMethodName) const {
    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(isValidAtomClass(class1),
                       ApiClassName,
                       CallingMethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(isValidAtomClass(class2),
                       ApiClassName,
                       CallingMethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(isValidAtomClass(class3),
                       ApiClassName,
                       CallingMethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK1(isValidAtomClass(class4),
                       ApiClassName,
                       CallingMethodName,
                       "class4=%d which is not a valid atom class Index",
                       (int)class4);
    SimTK_APIARGCHECK(periodicity1 != -1 || periodicity2 != -1 || periodicity3 != -1 || periodicity4 != -1,
                      ApiClassName,
                      CallingMethodName,
                      "must be at least one torsion term supplied");


    if (periodicity1 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity1 && periodicity1 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity1);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);*/
        // scf changed 0 to -180 to allow NAST right handed helices

        SimTK_APIARGCHECK1(angleIsInRange(phase1InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle1(%g) is not valid: must be between -180 and 180 degrees, inclusive",
                           phase1InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            (periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity2 && periodicity2 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity2);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase2InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase2InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            periodicity3 != periodicity2,
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity3 && periodicity3 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity3);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase3InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase3InDegrees);
        // (we've already checked for any possible repeats)
    }
    if (periodicity4 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity4 && periodicity4 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity4(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity4);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase4InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle4(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase4InDegrees);
        // (we've already checked for any possible repeats)
    }

    // Canonicalize atom class quad by reversing order if necessary so that the
    // first class Index is numerically no larger than the fourth. Amber improper
    // torsions should not be canonicalized because order matters.
    const AtomClassIndexQuad key(class1, class2, class3, class4, shouldCanonicalizeClassOrder);

    // Attempt to create a new bond torsion entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = torsionMap.insert(std::make_pair(key, BondTorsion(key)));
    auto& bondTorsionEntry = ret.first->second;

    // A new entry or one that just had a custom term in it won't have a built in
    // term so we can load it up and we're done.
    if (!bondTorsionEntry.hasBuiltinTerm()) {
        if (periodicity1 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
        if (periodicity2 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
        if (periodicity3 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
        if (periodicity4 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity4, amp4InKJ, phase4InDegrees));
        }
        return;
    }

    // If we get here we have discovered that there is already a built in torsion
    // term present for this atom class quad. We can still insert new terms, and we'll
    // allow duplicates if they are identical.
    if (periodicity1 != -1) {
        const TorsionTerm& term1 = bondTorsionEntry.getTermWithPeriod(periodicity1);
        if (term1.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term1.amplitude, amp1InKJ)
                    && anglesAlmostEqual(term1.theta0 * DuMM::Rad2Deg, phase1InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity1);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
    }

    if (periodicity2 != -1) {
        const TorsionTerm& term2 = bondTorsionEntry.getTermWithPeriod(periodicity2);
        if (term2.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term2.amplitude, amp2InKJ)
                    && anglesAlmostEqual(term2.theta0 * DuMM::Rad2Deg, phase2InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity2);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
    }

    if (periodicity3 != -1) {
        const TorsionTerm& term3 = bondTorsionEntry.getTermWithPeriod(periodicity3);
        if (term3.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term3.amplitude, amp3InKJ)
                    && anglesAlmostEqual(term3.theta0 * DuMM::Rad2Deg, phase3InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity3);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
    }

    if (periodicity4 != -1) {
        const TorsionTerm& term4 = bondTorsionEntry.getTermWithPeriod(periodicity4);
        if (term4.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term4.amplitude, amp4InKJ)
                    && anglesAlmostEqual(term4.theta0 * DuMM::Rad2Deg, phase4InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity4);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity4, amp4InKJ, phase4InDegrees));
        }
    }
}

//
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions, and then inserts the built in torsion terms
// if they are legitimate. Written by S.A.T. for dihedral with 5 periodicities.
//
void DuMMForceFieldSubsystemRep::defineAnyTorsion(DuMM::AtomClassIndex class1,
                                                  DuMM::AtomClassIndex class2,
                                                  DuMM::AtomClassIndex class3,
                                                  DuMM::AtomClassIndex class4,
                                                  bool shouldCanonicalizeClassOrder,
                                                  int periodicity1,
                                                  Real amp1InKJ,
                                                  Real phase1InDegrees,
                                                  int periodicity2,
                                                  Real amp2InKJ,
                                                  Real phase2InDegrees,
                                                  int periodicity3,
                                                  Real amp3InKJ,
                                                  Real phase3InDegrees,
                                                  int periodicity4,
                                                  Real amp4InKJ,
                                                  Real phase4InDegrees,
                                                  int periodicity5,
                                                  Real amp5InKJ,
                                                  Real phase5InDegrees,
                                                  std::map<AtomClassIndexQuad, BondTorsion>& torsionMap,
                                                  const char* CallingMethodName) const {
    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(isValidAtomClass(class1),
                       ApiClassName,
                       CallingMethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(isValidAtomClass(class2),
                       ApiClassName,
                       CallingMethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(isValidAtomClass(class3),
                       ApiClassName,
                       CallingMethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK1(isValidAtomClass(class4),
                       ApiClassName,
                       CallingMethodName,
                       "class4=%d which is not a valid atom class Index",
                       (int)class4);
    SimTK_APIARGCHECK(periodicity1 != -1 || periodicity2 != -1 || periodicity3 != -1 || periodicity4 != -1
                          || periodicity5 != -1,
                      ApiClassName,
                      CallingMethodName,
                      "must be at least one torsion term supplied");

    if (periodicity1 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity1 && periodicity1 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity1);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);*/
        // scf changed 0 to -180 to allow NAST right handed helices

        SimTK_APIARGCHECK1(angleIsInRange(phase1InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle1(%g) is not valid: must be between -180 and 180 degrees, inclusive",
                           phase1InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            (periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity2 && periodicity2 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity2);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase2InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase2InDegrees);

        // No repeats.
        SimTK_APIARGCHECK1(
            periodicity3 != periodicity2,
            ApiClassName,
            CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity3 && periodicity3 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity3);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase3InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase3InDegrees);
        // (we've already checked for any possible repeats)
    }
    if (periodicity4 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity4 && periodicity4 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity4(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity4);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase4InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle4(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase4InDegrees);
        // (we've already checked for any possible repeats)
    }
    if (periodicity5 != -1) {
        // No nonsense.
        SimTK_APIARGCHECK1(1 <= periodicity5 && periodicity5 <= 6,
                           ApiClassName,
                           CallingMethodName,
                           "periodicity5(%d) is invalid: we require 1 <= periodicity <= 6",
                           periodicity4);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1(angleIsInRange(phase5InDegrees),
                           ApiClassName,
                           CallingMethodName,
                           "phaseAngle5(%g) is not valid: must be between 0 and 180 degrees, inclusive",
                           phase5InDegrees);
        // (we've already checked for any possible repeats)
    }


    // Canonicalize atom class quad by reversing order if necessary so that the
    // first class Index is numerically no larger than the fourth. Amber improper
    // torsions should not be canonicalized because order matters.
    const AtomClassIndexQuad key(class1, class2, class3, class4, shouldCanonicalizeClassOrder);

    // Attempt to create a new bond torsion entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    auto ret = torsionMap.insert(std::make_pair(key, BondTorsion(key)));
    auto& bondTorsionEntry = ret.first->second;

    // A new entry or one that just had a custom term in it won't have a built in
    // term so we can load it up and we're done.
    if (!bondTorsionEntry.hasBuiltinTerm()) {
        if (periodicity1 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
        if (periodicity2 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
        if (periodicity3 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
        if (periodicity4 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity4, amp4InKJ, phase4InDegrees));
        }
        if (periodicity5 != -1) {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity5, amp5InKJ, phase5InDegrees));
        }
        return;
    }

    // If we get here we have discovered that there is already a built in torsion
    // term present for this atom class quad. We can still insert new terms, and we'll
    // allow duplicates if they are identical.
    if (periodicity1 != -1) {
        const TorsionTerm& term1 = bondTorsionEntry.getTermWithPeriod(periodicity1);
        if (term1.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term1.amplitude, amp1InKJ)
                    && anglesAlmostEqual(term1.theta0 * SimTK::Rad2Deg, phase1InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity1);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        }
    }

    if (periodicity2 != -1) {
        const TorsionTerm& term2 = bondTorsionEntry.getTermWithPeriod(periodicity2);
        if (term2.isValid()) {
            const bool amplitudesEqual = almostEqual(term2.amplitude, amp2InKJ);
            if (!amplitudesEqual) {
                const std::string errorMsg =
                    "Atom class quad " + std::to_string((int)class1) + "," + std::to_string((int)class2) + ","
                    + std::to_string((int)class3) + "," + std::to_string((int)class4)
                    + " already had a different term with periodicity " + std::to_string(periodicity2)
                    + ": existing amplitude=" + std::to_string(term2.amplitude)
                    + " new amplitude=" + std::to_string(amp2InKJ);
                SimTK_ASSERT(amplitudesEqual, errorMsg.c_str());
            }

            const bool phasesEqual = anglesAlmostEqual(term2.theta0 * SimTK::Rad2Deg, phase2InDegrees);
            if (!phasesEqual) {
                const std::string errorMsg =
                    "Atom class quad " + std::to_string((int)class1) + "," + std::to_string((int)class2) + ","
                    + std::to_string((int)class3) + "," + std::to_string((int)class4)
                    + " already had a different term with periodicity " + std::to_string(periodicity2)
                    + ": existing phase=" + std::to_string(term2.theta0)
                    + " new phase=" + std::to_string(phase2InDegrees);
                SimTK_ASSERT(phasesEqual, errorMsg.c_str());
            }
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        }
    }

    if (periodicity3 != -1) {
        const TorsionTerm& term3 = bondTorsionEntry.getTermWithPeriod(periodicity3);
        if (term3.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term3.amplitude, amp3InKJ)
                    && anglesAlmostEqual(term3.theta0 * SimTK::Rad2Deg, phase3InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity3);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        }
    }

    if (periodicity4 != -1) {
        const TorsionTerm& term4 = bondTorsionEntry.getTermWithPeriod(periodicity4);
        if (term4.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term4.amplitude, amp4InKJ)
                    && anglesAlmostEqual(term4.theta0 * SimTK::Rad2Deg, phase4InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity4);
        } else {
            // bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
            // //Laurentiu Code
            bondTorsionEntry.addBuiltinTerm(
                TorsionTerm(periodicity4, amp4InKJ, phase4InDegrees)); // Teodor Code
        }
    }

    if (periodicity5 != -1) {
        const TorsionTerm& term5 = bondTorsionEntry.getTermWithPeriod(periodicity5);
        if (term5.isValid()) {
            SimTK_APIARGCHECK5(
                almostEqual(term5.amplitude, amp5InKJ)
                    && anglesAlmostEqual(term5.theta0 * SimTK::Rad2Deg, phase5InDegrees),
                ApiClassName,
                CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,
                (int)class2,
                (int)class3,
                (int)class4,
                periodicity5);
        } else {
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity5, amp5InKJ, phase5InDegrees));
        }
    }
}

// We allow up to 3 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                DuMM::AtomClassIndex class3,
                                                DuMM::AtomClassIndex class4,
                                                int periodicity1,
                                                Real amp1InKJ,
                                                Real phase1InDegrees,
                                                int periodicity2,
                                                Real amp2InKJ,
                                                Real phase2InDegrees,
                                                int periodicity3,
                                                Real amp3InKJ,
                                                Real phase3InDegrees) {
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();
    rep.defineAnyTorsion(class1,
                         class2,
                         class3,
                         class4,
                         true, // canonicalize
                         periodicity1,
                         amp1InKJ,
                         phase1InDegrees,
                         periodicity2,
                         amp2InKJ,
                         phase2InDegrees,
                         periodicity3,
                         amp3InKJ,
                         phase3InDegrees,
                         rep.bondTorsion,
                         MethodName);
}

//
// We allow up to 4 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                DuMM::AtomClassIndex class3,
                                                DuMM::AtomClassIndex class4,
                                                int periodicity1,
                                                Real amp1InKJ,
                                                Real phase1InDegrees,
                                                int periodicity2,
                                                Real amp2InKJ,
                                                Real phase2InDegrees,
                                                int periodicity3,
                                                Real amp3InKJ,
                                                Real phase3InDegrees,
                                                int periodicity4,
                                                Real amp4InKJ,
                                                Real phase4InDegrees) {
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();
    rep.defineAnyTorsion(class1,
                         class2,
                         class3,
                         class4,
                         true, // canonicalize
                         periodicity1,
                         amp1InKJ,
                         phase1InDegrees,
                         periodicity2,
                         amp2InKJ,
                         phase2InDegrees,
                         periodicity3,
                         amp3InKJ,
                         phase3InDegrees,
                         periodicity4,
                         amp4InKJ,
                         phase4InDegrees,
                         rep.bondTorsion,
                         MethodName);
}

//
// Torsion with up to five terms (Added by S.A.T., as it is needed
// when simulating lipids).
//
void DuMMForceFieldSubsystem::defineBondTorsion(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                DuMM::AtomClassIndex class3,
                                                DuMM::AtomClassIndex class4,
                                                int periodicity1,
                                                Real amp1InKJ,
                                                Real phase1InDegrees,
                                                int periodicity2,
                                                Real amp2InKJ,
                                                Real phase2InDegrees,
                                                int periodicity3,
                                                Real amp3InKJ,
                                                Real phase3InDegrees,
                                                int periodicity4,
                                                Real amp4InKJ,
                                                Real phase4InDegrees,
                                                int periodicity5,
                                                Real amp5InKJ,
                                                Real phase5InDegrees) {
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();
    rep.defineAnyTorsion(class1,
                         class2,
                         class3,
                         class4,
                         true, // canonicalize
                         periodicity1,
                         amp1InKJ,
                         phase1InDegrees,
                         periodicity2,
                         amp2InKJ,
                         phase2InDegrees,
                         periodicity3,
                         amp3InKJ,
                         phase3InDegrees,
                         periodicity4,
                         amp4InKJ,
                         phase4InDegrees,
                         periodicity5,
                         amp5InKJ,
                         phase5InDegrees,
                         rep.bondTorsion,
                         MethodName);
}

void DuMMForceFieldSubsystem::defineCustomBondTorsion(DuMM::AtomClassIndex class1,
                                                      DuMM::AtomClassIndex class2,
                                                      DuMM::AtomClassIndex class3,
                                                      DuMM::AtomClassIndex class4,
                                                      DuMM::CustomBondTorsion* customBondTorsion) {
    static const char* MethodName = "defineCustomBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Watch for nonsense arguments.
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class1),
                       rep.ApiClassName,
                       MethodName,
                       "class1=%d which is not a valid atom class Index",
                       (int)class1);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class2),
                       rep.ApiClassName,
                       MethodName,
                       "class2=%d which is not a valid atom class Index",
                       (int)class2);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class3),
                       rep.ApiClassName,
                       MethodName,
                       "class3=%d which is not a valid atom class Index",
                       (int)class3);
    SimTK_APIARGCHECK1(rep.isValidAtomClass(class3),
                       rep.ApiClassName,
                       MethodName,
                       "class4=%d which is not a valid atom class Index",
                       (int)class4);
    SimTK_APIARGCHECK(customBondTorsion, rep.ApiClassName, MethodName, "CustomBondTorsion pointer was null");

    // Canonicalize atom class quad by reversing order if necessary so that the
    // first class Index is numerically no larger than the fourth.
    const AtomClassIndexQuad key(class1, class2, class3, class4, true);

    // Attempt to create a new bond torsion entry containing no valid
    // terms. If there was already an entry it will be returned instead
    // and no insertion is performed.
    std::pair<std::map<AtomClassIndexQuad, BondTorsion>::iterator, bool> ret =
        rep.bondTorsion.insert(std::pair<AtomClassIndexQuad, BondTorsion>(key, BondTorsion(key)));

    BondTorsion& bondTorsionEntry = ret.first->second;
    bondTorsionEntry.addCustomTerm(customBondTorsion);
}

// Convenient signature for a bond torsion with only one term.
void DuMMForceFieldSubsystem::defineBondTorsion(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                DuMM::AtomClassIndex class3,
                                                DuMM::AtomClassIndex class4,
                                                int periodicity1,
                                                Real amp1InKJ,
                                                Real phase1InDegrees) {
    defineBondTorsion(class1,
                      class2,
                      class3,
                      class4,
                      periodicity1,
                      amp1InKJ,
                      phase1InDegrees,
                      -1,
                      0.,
                      0.,
                      -1,
                      0.,
                      0.);
}

// Convenient signature for a bond torsion with two terms.
void DuMMForceFieldSubsystem::defineBondTorsion(DuMM::AtomClassIndex class1,
                                                DuMM::AtomClassIndex class2,
                                                DuMM::AtomClassIndex class3,
                                                DuMM::AtomClassIndex class4,
                                                int periodicity1,
                                                Real amp1InKJ,
                                                Real phase1InDegrees,
                                                int periodicity2,
                                                Real amp2InKJ,
                                                Real phase2InDegrees) {
    defineBondTorsion(class1,
                      class2,
                      class3,
                      class4,
                      periodicity1,
                      amp1InKJ,
                      phase1InDegrees,
                      periodicity2,
                      amp2InKJ,
                      phase2InDegrees,
                      -1,
                      0.,
                      0.);
}

//
// This function is based on the defineBondTorsion function.
// As with the normal bond torsions, we allow up to 3 terms in a single torsion function,
// with three different periodicities. If any of these are unused, set the corresponding
// periodicity to -1.
//
void DuMMForceFieldSubsystem::defineAmberImproperTorsion(DuMM::AtomClassIndex class1,
                                                         DuMM::AtomClassIndex class2,
                                                         DuMM::AtomClassIndex class3,
                                                         DuMM::AtomClassIndex class4,
                                                         int periodicity1,
                                                         Real amp1InKJ,
                                                         Real phase1InDegrees,
                                                         int periodicity2,
                                                         Real amp2InKJ,
                                                         Real phase2InDegrees,
                                                         int periodicity3,
                                                         Real amp3InKJ,
                                                         Real phase3InDegrees) {
    static const char* MethodName = "defineAmberImproperTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();
    rep.defineAnyTorsion(class1,
                         class2,
                         class3,
                         class4,
                         false, // don't canonicalize
                         periodicity1,
                         amp1InKJ,
                         phase1InDegrees,
                         periodicity2,
                         amp2InKJ,
                         phase2InDegrees,
                         periodicity3,
                         amp3InKJ,
                         phase3InDegrees,
                         rep.amberImproperTorsion,
                         MethodName);
}

// Convenient signature for an amber improper torsion with only one term.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion(DuMM::AtomClassIndex class1,
                                                         DuMM::AtomClassIndex class2,
                                                         DuMM::AtomClassIndex class3,
                                                         DuMM::AtomClassIndex class4,
                                                         int periodicity1,
                                                         Real amp1InKJ,
                                                         Real phase1InDegrees) {
    defineAmberImproperTorsion(class1,
                               class2,
                               class3,
                               class4,
                               periodicity1,
                               amp1InKJ,
                               phase1InDegrees,
                               -1,
                               0.,
                               0.,
                               -1,
                               0.,
                               0.);
}

// Convenient signature for an amber improper torsion with two terms.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion(DuMM::AtomClassIndex class1,
                                                         DuMM::AtomClassIndex class2,
                                                         DuMM::AtomClassIndex class3,
                                                         DuMM::AtomClassIndex class4,
                                                         int periodicity1,
                                                         Real amp1InKJ,
                                                         Real phase1InDegrees,
                                                         int periodicity2,
                                                         Real amp2InKJ,
                                                         Real phase2InDegrees) {
    defineAmberImproperTorsion(class1,
                               class2,
                               class3,
                               class4,
                               periodicity1,
                               amp1InKJ,
                               phase1InDegrees,
                               periodicity2,
                               amp2InKJ,
                               phase2InDegrees,
                               -1,
                               0.,
                               0.);
}

void DuMMForceFieldSubsystem::clearIncludedNonbondAtomList() {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.includedNonbondAtoms.clear();
    inclList.includedNonbondBodies.clear();
    inclList.useDefaultNonbondList = false; // i.e., now there is nothing
}

void DuMMForceFieldSubsystem::clearIncludedBondList() {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.atomsWhoseBondsAreIncluded.clear();
    inclList.atomPairsWhoseConnectingBondsAreIncluded.clear();
    inclList.bodiesWhoseBondsAreIncluded.clear();
    inclList.bodyPairsWhoseConnectingBondsAreIncluded.clear();
    inclList.useDefaultBondList = false; // i.e., now there is nothing
}

void DuMMForceFieldSubsystem::resetIncludedNonbondAtomListToDefault() {
    clearIncludedNonbondAtomList();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = true;
}

void DuMMForceFieldSubsystem::resetIncludedBondListToDefault() {
    clearIncludedBondList();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = true;
}

void DuMMForceFieldSubsystem::includeNonbondAtom(DuMM::AtomIndex dAIx) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = false;
    inclList.includedNonbondAtoms.insert(dAIx); // ignores duplicates
}

void DuMMForceFieldSubsystem::includeAllNonbondAtomsForOneBody(MobilizedBodyIndex mobodIx) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = false;
    inclList.includedNonbondBodies.insert(mobodIx); // ignores duplicates
}

void DuMMForceFieldSubsystem::includeAllInterbodyBondsForOneAtom(DuMM::AtomIndex dAIx) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.atomsWhoseBondsAreIncluded.insert(dAIx); // ignores duplicates
}

void DuMMForceFieldSubsystem::includeAllInterbodyBondsWithBothAtoms(DuMM::AtomIndex atom1,
                                                                    DuMM::AtomIndex atom2) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.atomPairsWhoseConnectingBondsAreIncluded.insert(
        AtomIndexPair(atom1, atom2, true)); // canonicalize order; ignore dups
}

void DuMMForceFieldSubsystem::includeAllInterbodyBondsWithBothAtoms(DuMM::BondIndex bond) {
    includeAllInterbodyBondsWithBothAtoms(getBondAtom(bond, 0), getBondAtom(bond, 1));
}

void DuMMForceFieldSubsystem::includeAllInterbodyBondsForOneBody(MobilizedBodyIndex mobod) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.bodiesWhoseBondsAreIncluded.insert(mobod); // ignores dups
}

void DuMMForceFieldSubsystem::includeAllInterbodyBondsBetweenTwoBodies(MobilizedBodyIndex mobod1,
                                                                       MobilizedBodyIndex mobod2) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.bodyPairsWhoseConnectingBondsAreIncluded.insert(
        MobodIndexPair(mobod1, mobod2, true)); // canonicalize order
}

auto DuMMForceFieldSubsystem::getNumIncludedAtoms() const -> int {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
                                       "getNumIncludedAtoms",
                                       "Subsystem",
                                       "DuMMForceFieldSubsystem");
    return getRep().getNumIncludedAtoms();
}

auto DuMMForceFieldSubsystem::getAtomIndexOfIncludedAtom(DuMM::IncludedAtomIndex incAtomIx) const
    -> DuMM::AtomIndex {
    // Don't check in Release mode since this might get called a lot and
    // presumably we just checked in getNumIncludedAtoms().
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
                                       "getAtomIndexOfIncludedAtom",
                                       "Subsystem",
                                       "DuMMForceFieldSubsystem");

    return getRep().getAtomIndexOfIncludedAtom(incAtomIx);
}

auto DuMMForceFieldSubsystem::getNumNonbondAtoms() const -> int {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
                                       "getNumNonbondAtoms",
                                       "Subsystem",
                                       "DuMMForceFieldSubsystem");
    return getRep().getNumNonbondAtoms();
}

auto DuMMForceFieldSubsystem::getIncludedAtomIndexOfNonbondAtom(DuMM::NonbondAtomIndex nonbondAtomIx) const
    -> DuMM::IncludedAtomIndex {
    // Don't check in Release mode since this might get called a lot and
    // presumably we just checked in getNumNonbondAtoms().
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
                                       "getIncludedAtomIndexOfNonbondAtom",
                                       "Subsystem",
                                       "DuMMForceFieldSubsystem");

    return getRep().getIncludedAtomIndexOfNonbondAtom(nonbondAtomIx);
}

auto DuMMForceFieldSubsystem::getNonbondAtomIndex(DuMM::AtomIndex dAIx) const -> DuMM::NonbondAtomIndex {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
                                       "getIncludedAtomIndexOfNonbondAtom",
                                       "Subsystem",
                                       "DuMMForceFieldSubsystem");
    const auto& rep = getRep();
    const auto& atom = rep.getAtom(dAIx);
    return atom.getNonbondAtomIndex();
}

void DuMMForceFieldSubsystem::setTracing(bool shouldTrace) {
    updRep().tracing = shouldTrace;
}

void DuMMForceFieldSubsystem::updateOpenMMPositionsFromState(const State& state) const {
    OPENMM::get().updatePositionsCache(getRep().getNonBondedMappings(),
                                       getRep().getIncludedAtomPositionsInG(state));
}

void DuMMForceFieldSubsystem::evaluateEnergiesFromState(SimTK::Real& newPotentialEnergy,
                                                        SimTK::Real& newKineticEnergy) {
    OPENMM::get().evaluateEnergiesFromPositionCache(newPotentialEnergy, newKineticEnergy);
}

auto DuMMForceFieldSubsystem::integrateTrajectoryWithOpenMM(int steps, SimTK::Real timeStepInPicoseconds)
    -> bool {
    return OPENMM::get().integrateTrajectory(steps, timeStepInPicoseconds);
}

auto DuMMForceFieldSubsystem::createCluster(const char* clusterName) -> DuMM::ClusterIndex {
    invalidateSubsystemTopologyCache();
    // Currently there is no error checking to do. We don't insist on unique cluster names.
    return updRep().addCluster(Cluster(clusterName));
}

auto DuMMForceFieldSubsystem::addAtom(DuMM::ChargedAtomTypeIndex chargedAtomTypeIx) -> DuMM::AtomIndex {
    static const char* MethodName = "addAtom";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    SimTK_APIARGCHECK1(rep.isValidChargedAtomType(chargedAtomTypeIndex),
                       rep.ApiClassName,
                       MethodName,
                       "charged atom type %d is not valid",
                       (int)chargedAtomTypeIx);

    const auto atomIndex = (const DuMM::AtomIndex)rep.atoms.size();
    rep.atoms.push_back(DuMMAtom(chargedAtomTypeIx, atomIndex));

    return atomIndex;
}

void DuMMForceFieldSubsystem::placeAtomInCluster(DuMM::AtomIndex dAIx,
                                                 DuMM::ClusterIndex clusterIx,
                                                 const Vec3& station) {
    static const char* MethodName = "placeAtomInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure that we've seen both the atomIndex and clusterIx before.
    SimTK_APIARGCHECK1(rep.isValidAtom(atomIndex),
                       rep.ApiClassName,
                       MethodName,
                       "atom index %d is not valid",
                       (int)atomIndex);
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster index %d is not valid",
                       (int)clusterIx);

    Cluster& cluster = rep.updCluster(clusterIx);

    // Make sure that this cluster doesn't already contain this atom, either directly
    // or recursively through its subclusters.
    SimTK_APIARGCHECK3(!cluster.containsAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') already contains atom %d",
                       (int)clusterIx,
                       cluster.name.c_str(),
                       (int)dAIx);

    // Add the atom to the cluster.
    cluster.placeAtom(dAIx, station, rep);
}

void DuMMForceFieldSubsystem::placeClusterInCluster(DuMM::ClusterIndex childClusterIndex,
                                                    DuMM::ClusterIndex parentClusterIndex,
                                                    const Transform& placementInNm) {
    static const char* MethodName = "placeClusterInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1(rep.isValidCluster(childClusterIndex),
                       rep.ApiClassName,
                       MethodName,
                       "child cluster Index %d is not valid",
                       (int)childClusterIndex);
    SimTK_APIARGCHECK1(rep.isValidCluster(parentClusterIndex),
                       rep.ApiClassName,
                       MethodName,
                       "parent cluster Index %d is not valid",
                       (int)parentClusterIndex);

    Cluster& parent = rep.updCluster(parentClusterIndex);
    const Cluster& child = rep.getCluster(childClusterIndex);

    // TODO: for now, make sure the parent is a top-level cluster, meaning that it does
    // not have any parent clusters (although it can be attached to a body). This restriction
    // should be relaxed but it is tricky to get all the parents' and ancestors' content
    // lists updated correctly so I'm deferring that for now (sherm 060928).
    SimTK_APIARGCHECK2(
        parent.isTopLevelCluster(),
        rep.ApiClassName,
        MethodName,
        "parent cluster %d('%s') is not a top-level cluster so you cannot add a child cluster to it now",
        (int)parentClusterIndex,
        parent.name.c_str());

    // Child must not already be attached to a body.
    SimTK_APIARGCHECK2(
        !child.isAttachedToBody(),
        rep.ApiClassName,
        MethodName,
        "child cluster %d('%s') is already attached to a body so cannot now be placed in another cluster",
        (int)childClusterIndex,
        child.name.c_str());

    // Make sure that parent cluster doesn't already contain child cluster, either directly
    // or recursively through its subclusters.
    SimTK_APIARGCHECK4(!parent.containsCluster(childClusterIndex),
                       rep.ApiClassName,
                       MethodName,
                       "parent cluster %d('%s') already contains child cluster %d('%s')",
                       (int)parentClusterIndex,
                       parent.name.c_str(),
                       (int)childClusterIndex,
                       child.name.c_str());

    // Make sure the new child cluster doesn't contain any atoms which are already in
    // any of the trees to which the parent cluster is associated.
    // TODO: for now we need only look at the parent since we know it is top level.
    DuMM::AtomIndex atomIndex;
    SimTK_APIARGCHECK5(!parent.overlapsWithCluster(child, atomIndex),
                       rep.ApiClassName,
                       MethodName,
                       "parent cluster %d('%s') and would-be child cluster %d('%s') both contain atom %d"
                       " so they cannot have a parent/child relationship",
                       (int)parentClusterIndex,
                       parent.name.c_str(),
                       (int)childClusterIndex,
                       child.name.c_str(),
                       (int)atomIndex);

    // Add the child cluster to the parent.
    parent.placeCluster(childClusterIndex, placementInNm, rep);
}

void DuMMForceFieldSubsystem::attachClusterToBody(DuMM::ClusterIndex clusterIx,
                                                  MobilizedBodyIndex mobodIx,
                                                  const Transform& placementInNm) {
    static const char* MethodName = "attachClusterToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this cluster before, and that the body number is well formed.
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster Index %d is not valid",
                       (int)clusterIx);
    SimTK_APIARGCHECK1(mobodIx.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "body number %d is not valid: must be nonnegative",
                       (int)mobodIx);

    const Cluster& child = rep.getCluster(clusterIx);

    // Child must not already be attached to a body.
    SimTK_APIARGCHECK3(!child.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') is already attached to body %d so cannot now be"
                       " attached to a body",
                       (int)clusterIx,
                       child.name.c_str(),
                       (int)child.getMobodIndex());

    // None of the atoms in the child can be attached to any body.
    DuMM::AtomIndex tempAtomIndex;
    MobilizedBodyIndex tempBodyIndex;
    SimTK_APIARGCHECK4(!child.containsAnyAtomsAttachedToABody(tempAtomIndex, tempBodyIndex, rep),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') contains atom %d which is already attached to body %d"
                       " so the cluster cannot now be attached to another body",
                       (int)clusterIx,
                       child.name.c_str(),
                       (int)tempAtomIndex,
                       (int)tempBodyIndex);

    // Create an entry for the body if necessary, and its corresponding cluster.
    DuMMBodyIndex duMMBodyIndex = rep.ensureDuMMBodyEntryExists(mobodIx);
    Cluster& bodyCluster = rep.updCluster(rep.getDuMMBody(duMMBodyIndex).getClusterIndex());

    // Make sure that body cluster doesn't already contain child cluster, either directly
    // or recursively through its subclusters.
    SimTK_APIARGCHECK3(!bodyCluster.containsCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') is already attached (directly or indirectly) to"
                       " body %d",
                       (int)clusterIx,
                       child.name.c_str(),
                       (int)mobodIx);

    // OK, attach the cluster to the body's cluster.
    bodyCluster.placeCluster(clusterIx, placementInNm, rep);
}

// EU COMMENT BEGIN
void DuMMForceFieldSubsystem::attachAtomToBody(DuMM::AtomIndex atomIndex,
                                               MobilizedBodyIndex bodyIndex,
                                               const Vec3& stationInNm)
// EU BEGIN
// void DuMMForceFieldSubsystem::attachAtomToBody
//   (DuMM::AtomIndex atomIndex, MobilizedBodyIndex bodyIndex,
//    Vec3 stationInNm)
// EU END
{
    static const char* MethodName = "attachAtomToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this atom before, and that the body number is well formed.
    SimTK_APIARGCHECK1(rep.isValidAtom(atomIndex),
                       rep.ApiClassName,
                       MethodName,
                       "atom index %d is not valid",
                       (int)atomIndex);
    SimTK_APIARGCHECK1(bodyIndex.isValid(),
                       rep.ApiClassName,
                       MethodName,
                       "body number %d is not valid: must be nonnegative",
                       (int)bodyIndex);

    // The atom must not already be attached to a body, even this one.
    SimTK_APIARGCHECK2(!rep.getAtom(atomIndex).isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is already attached to body %d so cannot now be attached"
                       " to a body",
                       (int)atomIndex,
                       (int)rep.getAtom(atomIndex).getMobodIndex());

    // Create an entry for the body if necessary, and its corresponding cluster.
    auto duMMBodyIndex = rep.ensureDuMMBodyEntryExists(bodyIndex);
    auto& bodyCluster = rep.updCluster(rep.getDuMMBody(duMMBodyIndex).getClusterIndex());

    // Attach the atom to the body's cluster.
    bodyCluster.placeAtom(atomIndex, stationInNm, rep);
}

auto DuMMForceFieldSubsystem::calcClusterMassProperties(DuMM::ClusterIndex clusterIx,
                                                        const Transform& X_BC) const -> MassProperties {
    static const char* MethodName = "calcClusterMassProperties";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster Index %d is not valid",
                       (int)clusterIx);

    return rep.getCluster(clusterIx).calcMassProperties(X_BC, rep);
}


DuMM::BondIndex DuMMForceFieldSubsystem::addBond(DuMM::AtomIndex atom1Ix, DuMM::AtomIndex atom2Ix) {
    static const char* MethodName = "addBond";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen these atoms before.
    SimTK_APIARGCHECK1(rep.isValidAtom(atom1Ix),
                       rep.ApiClassName,
                       MethodName,
                       "atom1(%d) is not valid",
                       (int)atom1Ix);
    SimTK_APIARGCHECK1(rep.isValidAtom(atom2Ix),
                       rep.ApiClassName,
                       MethodName,
                       "atom2(%d) is not valid",
                       (int)atom2Ix);

    // An atom can't be bonded to itself.
    SimTK_APIARGCHECK1(atom1Ix != atom2Ix,
                       rep.ApiClassName,
                       MethodName,
                       "the same atom index (%d) was given for both atoms, which makes no sense",
                       (int)atom1Ix);

    // Ensure that atom1 < atom2
    if (atom1Ix > atom2Ix) {
        std::swap(atom1Ix, atom2Ix);
    }

    DuMMAtom& a1 = rep.updAtom(atom1Ix);
    DuMMAtom& a2 = rep.updAtom(atom2Ix);

    SimTK_APIARGCHECK2(!a1.isBondedTo(atom2Ix),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is already bonded to atom %d; you can only do that once",
                       (int)atom1Ix,
                       (int)atom2Ix);

    rep.bonds.push_back(Bond(atom1Ix, atom2Ix));
    a1.bond12.push_back(atom2Ix);
    a2.bond12.push_back(atom1Ix);
    return (DuMM::BondIndex)(rep.bonds.size() - 1);
}

auto DuMMForceFieldSubsystem::getNumAtoms() const -> int {
    return getRep().getNumAtoms();
}
auto DuMMForceFieldSubsystem::getNumBonds() const -> int {
    return getRep().getNumBonds();
}

// 'which' is 0 or 1 to pick one of the two atoms whose index we return.
auto DuMMForceFieldSubsystem::getBondAtom(DuMM::BondIndex bondIx, int which) const -> DuMM::AtomIndex {
    static const char* MethodName = "getBondAtom";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this bond before.
    SimTK_APIARGCHECK1(rep.isValidBond(bondIx),
                       rep.ApiClassName,
                       MethodName,
                       "bond %d is not valid",
                       (int)bondIx);

    SimTK_APIARGCHECK1(which == 0 || which == 1,
                       rep.ApiClassName,
                       MethodName,
                       "'which' was %d but must be 0 or 1 to choose one of the two atoms",
                       which);

    return rep.bonds[bondIx].atoms[which];
}

// Returns the atomic number (number of protons in nucleus).
auto DuMMForceFieldSubsystem::getAtomElement(DuMM::AtomIndex dAIx) const -> int {
    static const char* MethodName = "getAtomElement";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(atomIndex),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    return rep.getAtomElementNum(dAIx);
}

// Returned radius is in nm.
auto DuMMForceFieldSubsystem::getAtomRadius(DuMM::AtomIndex dAIx) const -> Real {
    static const char* MethodName = "getAtomRadius";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    const AtomClass& atomClass = rep.atomClasses[rep.getAtomClassIndex(dAIx)];
    return atomClass.vdwRadius;
}

// Returned station is in nm.
auto DuMMForceFieldSubsystem::getAtomStationOnBody(DuMM::AtomIndex dAIx) const -> const Vec3& {
    static const char* MethodName = "getAtomStationOnBody";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    const DuMMAtom& atom = rep.getAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(atom.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    return atom.station_B;
}

struct AtomPlacementMap {
    DuMM::ClusterIndex clusterIndex;
    std::size_t atomPlacementIx;
};

void DuMMForceFieldSubsystem::updateClustersCacheList(DuMM::AtomIndex dAIx, MobilizedBodyIndex inMbx) {
    DuMMForceFieldSubsystemRep& rep = updRep();

    const auto numAtoms = rep.getNumAtoms();
    dAIxToClusterIndex.resize(numAtoms);
    dAIxToAtomPlacementIx.resize(numAtoms);

    const auto numClusters = rep.getNumClusters();
    mbxToClusterIndex.resize(numClusters);

    for (DuMMBodyIndex dBIx(0); dBIx < rep.duMMSubsetOfBodies.size(); ++dBIx) {
        DuMMBody& body = rep.duMMSubsetOfBodies[dBIx];
        const MobodIndex mbx = body.getMobilizedBodyIndex();

        // Found the DuMM body with required mbx
        if (mbx == inMbx) {
            // Get the corresponding cluster
            const Cluster& cluster = rep.getCluster(body.clusterIndex);

            // Save cluster index
            mbxToClusterIndex[mbx] = body.clusterIndex;

            // Iterate AtomPlacementArray of this Cluster
            const auto& atomPlacements = cluster.getAllContainedAtoms();
            const auto numPlacements = atomPlacements.size();

            for (std::size_t atomPlacementIx = 0; atomPlacementIx < numPlacements; ++atomPlacementIx) {
                const auto& atomPlacement = atomPlacements[atomPlacementIx];
                assert(atomPlacement.isValid());

                if (atomPlacement.atomIndex == dAIx) {
                    dAIxToClusterIndex[dAIx] = body.clusterIndex;
                    dAIxToAtomPlacementIx[dAIx] = atomPlacementIx;
                    break;
                }
            }
        } // found the mobod
    } // every DummBody
}

/// Set the station at which a particular atom is fixed on its body.
/// An exception will be thrown if this atom is not fixed to any body.
void DuMMForceFieldSubsystem::bsetAtomStationOnBody(DuMM::AtomIndex dAIx, const Vec3& newStationB) {
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    auto& atom = rep.updAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(atom.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    atom.station_B = newStationB;
}

/// Set the station at which a particular atom is fixed on its body.
/// An exception will be thrown if this atom is not fixed to any body.
void DuMMForceFieldSubsystem::bsetAllAtomStationOnBody(DuMM::AtomIndex dAIx, const Vec3& newStationB) {
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    auto& atom = rep.updAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(atom.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    atom.station_B_All = newStationB;
}

/*!
 * <!-- Set the atom station in Atom placements of the inputMbx
 * corresponding cluster.
 * used in realizeSubsystemTopology
 * -->
 */
void DuMMForceFieldSubsystem::bsetAtomPlacementStation(DuMM::AtomIndex dAIx,
                                                       MobilizedBodyIndex inMbx,
                                                       const Vec3& newStation) {
    DuMMForceFieldSubsystemRep& rep = updRep();

    // Get cluster
    const auto& clusterIndex = dAIxToClusterIndex[dAIx];
    auto& cluster = rep.updCluster(clusterIndex);

    // Get atom placement in the cluster
    const auto& atomPlacementIx = dAIxToAtomPlacementIx[dAIx];
    auto& atomPlacements = cluster.updAllContainedAtoms();
    auto& atomPlacement = atomPlacements[atomPlacementIx];

    atomPlacement.setStation(newStation);
}

// Stations computed every time
auto DuMMForceFieldSubsystem::updIncludedAtomStation(DuMM::AtomIndex dAIx) -> Vec3& {
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    auto& atom = rep.updAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(atom.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    return rep.updIncludedAtomStation(atom.inclAtomIndex);
}

// Stations computed every time - CalcFullPotential
auto DuMMForceFieldSubsystem::updAllAtomStation(DuMM::AtomIndex dAIx) -> Vec3& {
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& rep = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not valid",
                       (int)dAIx);

    auto& atom = rep.updAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(atom.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    return rep.updAllAtomStation(atom.inclAtomIndex);
}

// Get ClusterIndex corresponding to specified Mobod
auto DuMMForceFieldSubsystem::bgetMobodClusterIndex(MobilizedBodyIndex mbx) const -> DuMM::ClusterIndex {
    return mbxToClusterIndex[mbx];
}
// EU END


// Returned placement is in nm.
auto DuMMForceFieldSubsystem::getClusterPlacementOnBody(DuMM::ClusterIndex clusterIx) const
    -> const Transform& {
    static const char* MethodName = "getClusterPlacementOnBody";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster Index %d is not valid",
                       (int)clusterIx);

    const Cluster& cluster = rep.getCluster(clusterIx);

    // Cluster must be attached to a body.
    SimTK_APIARGCHECK2(cluster.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') is not attached to a body",
                       (int)clusterIx,
                       cluster.name.c_str());

    return cluster.placement_B;
}

// Returned station is in nm.
auto DuMMForceFieldSubsystem::getAtomStationInCluster(DuMM::AtomIndex dAIx,
                                                      DuMM::ClusterIndex clusterIx) const -> const Vec3& {
    static const char* MethodName = "getAtomStationInCluster";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure that we've seen both the atomIndex and clusterIx before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom index %d is not valid",
                       (int)dAIx);
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster index %d is not valid",
                       (int)clusterIx);

    const Cluster& cluster = rep.getCluster(clusterIx);
    const AtomPlacementArray& atoms = cluster.getAllContainedAtoms();
    const auto atomPlacement = std::find(atoms.begin(), atoms.end(), AtomPlacement(dAIx, Vec3(0)));

    // We're going to be upset of this cluster doesn't contain this atom.
    SimTK_APIARGCHECK3(atomPlacement != atoms.end(),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') does not contain atom %d",
                       (int)clusterIx,
                       cluster.name.c_str(),
                       (int)dAIx);

    return atomPlacement->station;
}

// Returned placement is in nm.
auto DuMMForceFieldSubsystem::getClusterPlacementInCluster(DuMM::ClusterIndex childClusterIndex,
                                                           DuMM::ClusterIndex parentClusterIndex) const
    -> const Transform& {
    static const char* MethodName = "getClusterPlacementInCluster";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1(rep.isValidCluster(childClusterIndex),
                       rep.ApiClassName,
                       MethodName,
                       "child cluster Index %d is not valid",
                       (int)childClusterIndex);
    SimTK_APIARGCHECK1(rep.isValidCluster(parentClusterIndex),
                       rep.ApiClassName,
                       MethodName,
                       "parent cluster Index %d is not valid",
                       (int)parentClusterIndex);

    const Cluster& parent = rep.getCluster(parentClusterIndex);
    const Cluster& child = rep.getCluster(childClusterIndex);

    const auto& clusters = parent.getAllContainedClusters();
    const auto clusterPlacement =
        std::find(clusters.begin(), clusters.end(), ClusterPlacement(childClusterIndex, Transform()));

    // We're going to be upset of the parent cluster doesn't contain the child.
    SimTK_APIARGCHECK4(clusterPlacement != clusters.end(),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') does not contain cluster %d('%d')",
                       (int)parentClusterIndex,
                       parent.name.c_str(),
                       (int)childClusterIndex,
                       child.name.c_str());

    return clusterPlacement->placement;
}

auto DuMMForceFieldSubsystem::getAtomBody(DuMM::AtomIndex dAIx) const -> MobilizedBodyIndex {
    static const char* MethodName = "getAtomBody";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure that we've seen this atomIndex before.
    SimTK_APIARGCHECK1(rep.isValidAtom(dAIx),
                       rep.ApiClassName,
                       MethodName,
                       "atom index %d is not valid",
                       (int)dAIx);

    const DuMMAtom& atom = rep.getAtom(dAIx);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1(a.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "atom %d is not attached to a body",
                       (int)dAIx);

    return atom.getMobodIndex();
}

auto DuMMForceFieldSubsystem::getClusterBody(DuMM::ClusterIndex clusterIx) const -> MobilizedBodyIndex {
    static const char* MethodName = "getClusterBody";
    const DuMMForceFieldSubsystemRep& rep = getRep();

    // Make sure that we've seen this atomIndex before.
    SimTK_APIARGCHECK1(rep.isValidCluster(clusterIx),
                       rep.ApiClassName,
                       MethodName,
                       "cluster Index %d is not valid",
                       (int)clusterIx);

    const Cluster& cluster = rep.getCluster(clusterIx);

    // Cluster must be attached to a body.
    SimTK_APIARGCHECK2(cluster.isAttachedToBody(),
                       rep.ApiClassName,
                       MethodName,
                       "cluster %d('%s') is not attached to a body",
                       (int)clusterIx,
                       cluster.name.c_str());

    return cluster.getMobodIndex();
}

void DuMMForceFieldSubsystem::dump() const {
    getRep().dump();
}

// How many times has the forcefield been evaluated?
auto DuMMForceFieldSubsystem::getForceEvaluationCount() const -> long long {
    return getRep().getForceEvaluationCount();
}

auto DuMMForceFieldSubsystemRep::generateBiotypeChargedAtomTypeSelfCode(std::ostream& os) const
    -> std::ostream& {
    std::map<BiotypeIndex, DuMM::ChargedAtomTypeIndex>::const_iterator i;
    for (i = chargedAtomTypesByBiotype.begin(); i != chargedAtomTypesByBiotype.end(); ++i) {
        generateBiotypeChargedAtomTypeSelfCode(os, i->first);
    }

    return os;
}

auto DuMMForceFieldSubsystem::generateBiotypeChargedAtomTypeSelfCode(std::ostream& ostream) const
    -> std::ostream& {
    return getRep().generateBiotypeChargedAtomTypeSelfCode(ostream);
}

void DuMMForceFieldSubsystem::setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex,
                                                        BiotypeIndex biotypeIx) {
    updRep().setBiotypeChargedAtomType(chargedAtomTypeIndex, biotypeIx);
}

DuMM::ChargedAtomTypeIndex DuMMForceFieldSubsystem::getBiotypeChargedAtomType(BiotypeIndex biotypeIx) const {
    return getRep().getBiotypeChargedAtomType(biotypeIx);
}
