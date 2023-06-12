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

#include "molmodel/internal/common.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"

#include "DuMMForceFieldSubsystemRep.h"
#include "TinkerAmber99.h"

//#ifndef DEBUG
//#define DEBUG 1
//#endif

#ifdef DEBUG
#define TRACE(STR) printf("%s", STR);
#else
#define TRACE(STR)
#endif

using namespace SimTK;

    ////////////////////////////////
    // DUMM FORCE FIELD SUBSYSTEM //
    ////////////////////////////////

/*static*/ bool 
DuMMForceFieldSubsystem::isInstanceOf(const Subsystem& s) {
    return DuMMForceFieldSubsystemRep::isA(s.getSubsystemGuts());
}
/*static*/ const DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const DuMMForceFieldSubsystem&>(s);
}
/*static*/ DuMMForceFieldSubsystem&
DuMMForceFieldSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<DuMMForceFieldSubsystem&>(s);
}

const DuMMForceFieldSubsystemRep& 
DuMMForceFieldSubsystem::getRep() const {
    return dynamic_cast<const DuMMForceFieldSubsystemRep&>( getSubsystemGuts() );
}
DuMMForceFieldSubsystemRep&       
DuMMForceFieldSubsystem::updRep() {
    return dynamic_cast<DuMMForceFieldSubsystemRep&>( updSubsystemGuts() );
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
DuMMForceFieldSubsystem::DuMMForceFieldSubsystem() 
  : ForceSubsystem()
{
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
}

DuMMForceFieldSubsystem::DuMMForceFieldSubsystem(MolecularMechanicsSystem& mms) 
  : ForceSubsystem()
{
    adoptSubsystemGuts(new DuMMForceFieldSubsystemRep());
    mms.setMolecularMechanicsForceSubsystem(*this); // steal ownership
}

DuMM::AtomClassIndex DuMMForceFieldSubsystem::getAtomClassIndex(DuMM::AtomIndex atomIx) const {
	DuMM::ChargedAtomTypeIndex typeIx = getRep().atoms[atomIx].chargedAtomTypeIndex;
	return getRep().chargedAtomTypes[typeIx].atomClassIx;
}
Real DuMMForceFieldSubsystem::getVdwRadius(DuMM::AtomClassIndex atomClassIx) const {
	return getRep().atomClasses[atomClassIx].vdwRadius;
}
Real DuMMForceFieldSubsystem::getVdwWellDepth(DuMM::AtomClassIndex atomClassIx) const {
	return getRep().atomClasses[atomClassIx].vdwWellDepth;
}

void DuMMForceFieldSubsystem::dumpCForceFieldParameters(std::ostream& os, const String& methodName) const {
    const DuMMForceFieldSubsystemRep& mm = getRep();

    os << "void " << methodName << "(DuMMForceFieldSubsystem& dumm)" << std::endl;
    os << "{" << std::endl; // open method

    // 1) define atom classes
    for (DuMM::AtomClassIndex i(0); i < (int)mm.atomClasses.size(); ++i) {
        if (!mm.atomClasses[i].isValid()) continue;
        const AtomClass& atomClass = mm.atomClasses[i];

        os << "    dumm.";
        atomClass.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 2) define charged atom types
    for (DuMM::ChargedAtomTypeIndex i(0); i < (int)mm.chargedAtomTypes.size(); ++i) {
        if (!mm.chargedAtomTypes[i].isValid()) continue;

        const ChargedAtomType& chargedAtomType = mm.chargedAtomTypes[i];
        os << "    dumm.";
        chargedAtomType.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 3) bond stretch parameters
    std::map<AtomClassIndexPair, BondStretch>::const_iterator b;
    for (b = mm.bondStretch.begin(); b != mm.bondStretch.end(); ++b) {
        os << "    dumm.";
        b->second.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 4) bond bend parameters
    std::map<AtomClassIndexTriple, BondBend>::const_iterator bendI;
    for (bendI = mm.bondBend.begin(); bendI != mm.bondBend.end(); ++bendI) {
        os << "    dumm.";
        bendI->second.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 5) bond torsion parameters
    std::map<AtomClassIndexQuad, BondTorsion>::const_iterator t;
    for (t = mm.bondTorsion.begin(); t != mm.bondTorsion.end(); ++t) {
        os << "    dumm.";
        t->second.generateSelfCode(os);
        os << std::endl;
    }

    os << std::endl;

    // 6) amber-style improper torsion parameters
    for (t = mm.amberImproperTorsion.begin(); t != mm.amberImproperTorsion.end(); ++t) {
        os << "    dumm.";
        t->second.generateSelfCode(os, 2);
        os << std::endl;
    }

    os << std::endl;

    // 7) global parameters

    // van der Waals mixing rule
    os << "    dumm.setVdwMixingRule(";
    switch (getVdwMixingRule()) {
        case WaldmanHagler:
            os << "DuMMForceFieldSubsystem::WaldmanHagler";
            break;
        case HalgrenHHG:
            os << "DuMMForceFieldSubsystem::HalgrenHHG";
            break;
        case Jorgensen:
            os << "DuMMForceFieldSubsystem::Jorgensen";
            break;
        case LorentzBerthelot:
            os << "DuMMForceFieldSubsystem::LorentzBerthelot";
            break;
        case Kong:
            os << "DuMMForceFieldSubsystem::Kong";
            break;
        default:
            assert(false);
            os << "DuMMForceFieldSubsystem::WaldmanHagler";
            break;
    }
    os << ");" << std::endl;

    os << "    dumm.setVdw12ScaleFactor(" << mm.vdwScale12 << ");" << std::endl;
    os << "    dumm.setVdw13ScaleFactor(" << mm.vdwScale13 << ");" << std::endl;
    os << "    dumm.setVdw14ScaleFactor(" << mm.vdwScale14 << ");" << std::endl;
    os << "    dumm.setVdw15ScaleFactor(" << mm.vdwScale15 << ");" << std::endl;

    os << "    dumm.setCoulomb12ScaleFactor(" << mm.coulombScale12 << ");" << std::endl;
    os << "    dumm.setCoulomb13ScaleFactor(" << mm.coulombScale13 << ");" << std::endl;
    os << "    dumm.setCoulomb14ScaleFactor(" << mm.coulombScale14 << ");" << std::endl;
    os << "    dumm.setCoulomb15ScaleFactor(" << mm.coulombScale15 << ");" << std::endl;

    os << "    dumm.setVdwGlobalScaleFactor("     << mm.vdwGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setCoulombGlobalScaleFactor(" << mm.coulombGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setGbsaGlobalScaleFactor("    << mm.gbsaGlobalScaleFactor << ");" << std::endl;

    os << "    dumm.setBondStretchGlobalScaleFactor(" << mm.bondStretchGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setBondBendGlobalScaleFactor("    << mm.bondBendGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setBondTorsionGlobalScaleFactor(" << mm.bondTorsionGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setAmberImproperTorsionGlobalScaleFactor(" << mm.amberImproperTorsionGlobalScaleFactor << ");" << std::endl;

    os << "    dumm.setCustomBondStretchGlobalScaleFactor(" << mm.customBondStretchGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setCustomBondBendGlobalScaleFactor("    << mm.customBondBendGlobalScaleFactor << ");" << std::endl;
    os << "    dumm.setCustomBondTorsionGlobalScaleFactor(" << mm.customBondTorsionGlobalScaleFactor << ");" << std::endl;

    os << "    dumm.setIncludeGbsaAceApproximation(" << mm.gbsaIncludeAceApproximation << ");" << std::endl;

    os << "}" << std::endl; // end of method
}

void DuMMForceFieldSubsystem::defineIncompleteAtomClass
   (DuMM::AtomClassIndex atomClassIx, const char* atomClassName, int element, int valence)
{
    static const char* MethodName = "defineIncompleteAtomClass";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Catch nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(atomClassIx.isValid(), mm.ApiClassName, MethodName,
        "atom class Index %d invalid: must be nonnegative", (int) atomClassIx);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidElement(element), mm.ApiClassName, MethodName,
        "element %d invalid: must be a valid atomic number and have an entry here",element);
    SimTK_APIARGCHECK1_ALWAYS(valence >= 0, mm.ApiClassName, MethodName, 
        "expected valence %d invalid: must be nonnegative", valence);

        // Make sure there is a slot available for this atom class.
    if (atomClassIx >= (DuMM::AtomClassIndex)mm.atomClasses.size())
        mm.atomClasses.resize(atomClassIx+1);

        // Make sure this atom class hasn't already been defined.
    SimTK_APIARGCHECK2_ALWAYS(!mm.atomClasses[atomClassIx].isValid(), mm.ApiClassName, MethodName, 
        "atom class Index %d is already in use for '%s'", (int) atomClassIx, 
        mm.atomClasses[atomClassIx].name.c_str());

	if (mm.atomClassIndicesByName.find(atomClassName) != mm.atomClassIndicesByName.end()) {
		DuMM::AtomClassIndex oldAtomClassIx = mm.atomClassIndicesByName.find(atomClassName)->second;
		if (oldAtomClassIx != atomClassIx) {
			throw(std::runtime_error(String("Duplicate atom class name: ") + atomClassName));
		}
	}

	mm.insertNewAtomClass( AtomClass(atomClassIx, atomClassName, element, valence, 
                                            NaN, NaN) );


}

void DuMMForceFieldSubsystem::setAtomClassVdwParameters(DuMM::AtomClassIndex atomClassIx, Real vdwRadiusInNm, Real vdwWellDepthInKJPerMol) 
{
    static const char* MethodName = "setAtomClassVdwParameters";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(atomClassIx.isValid(), mm.ApiClassName, MethodName,
        "atom class Index %d invalid: must be nonnegative", (int) atomClassIx);

    SimTK_APIARGCHECK1_ALWAYS(vdwRadiusInNm >= 0, mm.ApiClassName, MethodName, 
        "van der Waals radius %g invalid: must be nonnegative", vdwRadiusInNm);
    SimTK_APIARGCHECK1_ALWAYS(vdwWellDepthInKJPerMol >= 0, mm.ApiClassName, MethodName, 
        "van der Waals energy well depth %g invalid: must be nonnegative", vdwWellDepthInKJPerMol);

    AtomClass& atomClass = mm.atomClasses[atomClassIx];
    atomClass.vdwRadius = vdwRadiusInNm;
    atomClass.vdwWellDepth = vdwWellDepthInKJPerMol;
}

bool DuMMForceFieldSubsystem::isValidAtomClass(DuMM::AtomClassIndex atomClassIx) const {
    return getRep().isValidAtomClass(atomClassIx);
}

void DuMMForceFieldSubsystem::defineIncompleteChargedAtomType
    (DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, const char* typeName, DuMM::AtomClassIndex atomClassIx)
{
    static const char* MethodName = "defineIncompleteChargedAtomType";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Check for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(chargedAtomTypeIndex.isValid(), mm.ApiClassName, MethodName,
        "charged atom type index %d invalid: must be nonnegative", (int) chargedAtomTypeIndex);
    SimTK_APIARGCHECK1_ALWAYS(atomClassIx.isValid(), mm.ApiClassName, MethodName,
        "atom class index %d invalid: must be nonnegative", (int) atomClassIx);
    // partialCharge is a signed quantity

        // Make sure the referenced atom class has already been defined.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(atomClassIx), mm.ApiClassName, MethodName,
        "atom class %d is undefined", (int) atomClassIx);

        // Make sure there is a slot available for the new chargedAtomType.
    if (chargedAtomTypeIndex >= (int)mm.chargedAtomTypes.size())
        mm.chargedAtomTypes.resize(chargedAtomTypeIndex+1);

        // Check that this slot is not already in use.
    SimTK_APIARGCHECK2_ALWAYS(!mm.chargedAtomTypes[chargedAtomTypeIndex].isValid(), mm.ApiClassName, MethodName, 
        "charged atom type index %d is already in use for '%s'", (int) chargedAtomTypeIndex, 
        mm.chargedAtomTypes[chargedAtomTypeIndex].name.c_str());

	mm.insertNewChargedAtomType(ChargedAtomType(chargedAtomTypeIndex, typeName, atomClassIx, NaN));
}

bool DuMMForceFieldSubsystem::hasAtomClass(DuMM::AtomClassIndex atomClassIndex) const {
	return getRep().hasAtomClass(atomClassIndex);
}
bool DuMMForceFieldSubsystem::hasAtomClass(const String& atomClassName) const {
	return getRep().hasAtomClass(atomClassName);
}
DuMM::AtomClassIndex DuMMForceFieldSubsystem::getAtomClassIndex(const String& atomClassName) const {
	return getRep().getAtomClassIndex(atomClassName);
}
DuMM::AtomClassIndex DuMMForceFieldSubsystem::getNextUnusedAtomClassIndex() const {
	return getRep().getNextUnusedAtomClassIndex();
}

bool DuMMForceFieldSubsystem::hasChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex) const {
	return getRep().hasChargedAtomType(chargedAtomTypeIndex);
}
bool DuMMForceFieldSubsystem::hasChargedAtomType(const String& chargedTypeName) const {
	return getRep().hasChargedAtomType(chargedTypeName);
}
DuMM::ChargedAtomTypeIndex DuMMForceFieldSubsystem::getChargedAtomTypeIndex(const String& chargedTypeName) const {
	return getRep().getChargedAtomTypeIndex(chargedTypeName);
}
DuMM::ChargedAtomTypeIndex DuMMForceFieldSubsystem::getNextUnusedChargedAtomTypeIndex() const {
	return getRep().getNextUnusedChargedAtomTypeIndex();
}

void DuMMForceFieldSubsystem::setChargedAtomTypeCharge(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, Real charge) {
    static const char* MethodName = "setChargedAtomTypeCharge";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Check for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(chargedAtomTypeIndex.isValid(), mm.ApiClassName, MethodName,
        "charged atom type index %d invalid: must be nonnegative", (int) chargedAtomTypeIndex);

    ChargedAtomType& chargedAtomType = mm.chargedAtomTypes[chargedAtomTypeIndex];
    chargedAtomType.partialCharge = charge;
}

void DuMMForceFieldSubsystem::defineBondStretch
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, Real stiffnessInKJPerNmSq, Real nominalLengthInNm)
{
    static const char* MethodName = "defineBondStretch";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerNmSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerNmSq);
    SimTK_APIARGCHECK1_ALWAYS(nominalLengthInNm >= 0, mm.ApiClassName, MethodName, 
        "nominal length %g is not valid: must be nonnegative", nominalLengthInNm);

        // We canonicalize the key so that the atom class pair has the 
        // lower class Index first.
    const AtomClassIndexPair key(class1,class2,true);

        // Attempt to create a new bond stretch entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexPair,BondStretch>::iterator, bool> ret = 
      mm.bondStretch.insert(std::pair<AtomClassIndexPair,BondStretch>
        (key, BondStretch(key)));

    BondStretch& bondStretchEntry = ret.first->second;

    if (bondStretchEntry.hasBuiltinTerm()) {
        SimTK_APIARGCHECK2_ALWAYS(
            bondStretchEntry.k==stiffnessInKJPerNmSq && bondStretchEntry.d0==nominalLengthInNm, 
            mm.ApiClassName, MethodName, 
            "There was already a different built-in bond stretch term for atom class pair (%d,%d); only one is allowed."
            "\nUse a CustomBondStretch term if you need another term for the same atom class pair.",
            (int)key[0], (int)key[1]);
    } else
        bondStretchEntry.setBuiltinTerm(stiffnessInKJPerNmSq,nominalLengthInNm);
}

void DuMMForceFieldSubsystem::defineCustomBondStretch
    (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::CustomBondStretch* customBondStretch)
{
    static const char* MethodName = "defineCustomBondStretch";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK_ALWAYS(customBondStretch, mm.ApiClassName, MethodName, 
        "CustomBondStretch pointer was null");

        // We canonicalize the key so that the atom class pair has the 
        // lower class Index first.
    const AtomClassIndexPair key(class1,class2,true);

        // Attempt to create a new bond stretch entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexPair,BondStretch>::iterator, bool> ret = 
      mm.bondStretch.insert(std::pair<AtomClassIndexPair,BondStretch>
        (key, BondStretch(key)));

    BondStretch& bondStretchEntry = ret.first->second;
    bondStretchEntry.addCustomTerm(customBondStretch);
}

void DuMMForceFieldSubsystem::defineBondBend
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, Real stiffnessInKJPerRadSq, Real nominalAngleInDeg)
{
    static const char* MethodName = "defineBondBend";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Index", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(stiffnessInKJPerRadSq >= 0, mm.ApiClassName, MethodName, 
        "stiffness %g is not valid: must be nonnegative", stiffnessInKJPerRadSq);
    SimTK_APIARGCHECK1_ALWAYS(0 <= nominalAngleInDeg && nominalAngleInDeg <= 180, 
        mm.ApiClassName, MethodName, 
        "nominal angle %g is not valid: must be between 0 and 180 degrees, inclusive", 
        nominalAngleInDeg);

        // We canonicalize the key so that the first classIndex is no larger than the third.
    const AtomClassIndexTriple key(class1,class2,class3,true);

        // Attempt to create a new bond bend entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexTriple,BondBend>::iterator, bool> ret = 
        mm.bondBend.insert(std::pair<AtomClassIndexTriple,BondBend>
            (key, BondBend(key)));

    BondBend& bondBendEntry = ret.first->second;

    if (bondBendEntry.hasBuiltinTerm()) {
        SimTK_APIARGCHECK3_ALWAYS(
               bondBendEntry.k==stiffnessInKJPerRadSq 
            && bondBendEntry.theta0==nominalAngleInDeg*DuMM::Deg2Rad, 
            mm.ApiClassName, MethodName, 
            "There was already a different built-in bond bend term for atom class triple (%d,%d,%d); only one is allowed."
            "\nUse a CustomBondBend term if you need another term for the same atom class triple.",
            (int)key[0], (int)key[1], (int)key[2]);
    } else
        bondBendEntry.setBuiltinTerm(stiffnessInKJPerRadSq,nominalAngleInDeg);
}

void DuMMForceFieldSubsystem::defineCustomBondBend
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, 
    DuMM::CustomBondBend* customBondBend)
{
    static const char* MethodName = "defineCustomBondBend";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Index", (int) class3);
    SimTK_APIARGCHECK_ALWAYS(customBondBend, mm.ApiClassName, MethodName, 
        "CustomBondBend pointer was null");

        // We canonicalize the key so that the first classIndex is no larger than the third.
    const AtomClassIndexTriple key(class1,class2,class3,true);

        // Attempt to create a new bond bend entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexTriple,BondBend>::iterator, bool> ret = 
      mm.bondBend.insert(std::pair<AtomClassIndexTriple,BondBend>
        (key, BondBend(key)));

    BondBend& bondBendEntry = ret.first->second;
    bondBendEntry.addCustomTerm(customBondBend);
}

// 
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions, and then inserts the built in torsion terms
// if they are legitimate.
//
void DuMMForceFieldSubsystemRep::defineAnyTorsion 
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    bool shouldCanonicalizeClassOrder,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees,
    std::map<AtomClassIndexQuad,BondTorsion>& torsionMap,
    const char* CallingMethodName) const
{
        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class1), ApiClassName, CallingMethodName,
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class2), ApiClassName, CallingMethodName,
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class3), ApiClassName, CallingMethodName,
        "class3=%d which is not a valid atom class Index", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class4), ApiClassName, CallingMethodName,
        "class4=%d which is not a valid atom class Index", (int) class4);
    SimTK_APIARGCHECK_ALWAYS(periodicity1!=-1 || periodicity2!=-1 || periodicity3!=-1,
        ApiClassName, CallingMethodName, "must be at least one torsion term supplied");


    if (periodicity1 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity1 && periodicity1 <= 6, ApiClassName, CallingMethodName,
            "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6", periodicity1);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);*/
        //scf changed 0 to -180 to allow NAST right handed helices

        SimTK_APIARGCHECK1_ALWAYS(-180 <= phase1InDegrees && phase1InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle1(%g) is not valid: must be between -180 and 180 degrees, inclusive", phase1InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS((periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity2 && periodicity2 <= 6, ApiClassName, CallingMethodName,
            "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6", periodicity2);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);*/

        SimTK_APIARGCHECK1_ALWAYS(0 <= phase2InDegrees && phase2InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase2InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS(periodicity3 != periodicity2, ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity3 && periodicity3 <= 6, ApiClassName, CallingMethodName,
            "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6", periodicity3);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1_ALWAYS(0 <= phase3InDegrees && phase3InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase3InDegrees);
            // (we've already checked for any possible repeats)
    }


        // Canonicalize atom class quad by reversing order if necessary so that the
        // first class Index is numerically no larger than the fourth. Amber improper
        // torsions should not be canonicalized because order matters.
    const AtomClassIndexQuad key(class1, class2, class3, class4, shouldCanonicalizeClassOrder);

        // Attempt to create a new bond torsion entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexQuad,BondTorsion>::iterator, bool> ret = 
      torsionMap.insert(std::pair<AtomClassIndexQuad,BondTorsion>
        (key, BondTorsion(key)));

    BondTorsion& bondTorsionEntry = ret.first->second;

    // A new entry or one that just had a custom term in it won't have a built in
    // term so we can load it up and we're done.
    if (!bondTorsionEntry.hasBuiltinTerm()) {
        if (periodicity1 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        if (periodicity2 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        if (periodicity3 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        return;
    }

    // If we get here we have discovered that there is already a built in torsion
    // term present for this atom class quad. We can still insert new terms, and we'll
    // allow duplicates if they are identical.

    if (periodicity1 != -1) {
        const TorsionTerm& term1 = bondTorsionEntry.getTermWithPeriod(periodicity1);
        if (term1.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term1.amplitude==amp1InKJ && term1.theta0==phase1InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity1);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
    }
    if (periodicity2 != -1) {
        const TorsionTerm& term2 = bondTorsionEntry.getTermWithPeriod(periodicity2);
        if (term2.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term2.amplitude==amp2InKJ && term2.theta0==phase2InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity2);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
    }
    if (periodicity3 != -1) {
        const TorsionTerm& term3 = bondTorsionEntry.getTermWithPeriod(periodicity3);
        if (term3.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term3.amplitude==amp3InKJ && term3.theta0==phase3InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity3);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
    }
}

// 
// This is a utility method that checks for invalid inputs to the defineBondTorsion() and
// defineAmberImproperTorsion() functions, and then inserts the built in torsion terms
// if they are legitimate.
//
void DuMMForceFieldSubsystemRep::defineAnyTorsion 
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    bool shouldCanonicalizeClassOrder,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees,
    int periodicity4, Real amp4InKJ, Real phase4InDegrees,
    std::map<AtomClassIndexQuad,BondTorsion>& torsionMap,
    const char* CallingMethodName) const
{
        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class1), ApiClassName, CallingMethodName,
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class2), ApiClassName, CallingMethodName,
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class3), ApiClassName, CallingMethodName,
        "class3=%d which is not a valid atom class Index", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(isValidAtomClass(class4), ApiClassName, CallingMethodName,
        "class4=%d which is not a valid atom class Index", (int) class4);
    SimTK_APIARGCHECK_ALWAYS(periodicity1!=-1 || periodicity2!=-1 || periodicity3!=-1,
        ApiClassName, CallingMethodName, "must be at least one torsion term supplied");


    if (periodicity1 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity1 && periodicity1 <= 6, ApiClassName, CallingMethodName,
            "periodicity1(%d) is invalid: we require 1 <= periodicity <= 6", periodicity1);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp1InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude1(%g) is not valid: must be nonnegative", amp1InKJ);*/
        //scf changed 0 to -180 to allow NAST right handed helices

        SimTK_APIARGCHECK1_ALWAYS(-180 <= phase1InDegrees && phase1InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle1(%g) is not valid: must be between -180 and 180 degrees, inclusive", phase1InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS((periodicity2 != periodicity1) && (periodicity3 != periodicity1),
            ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity1);
    }
    if (periodicity2 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity2 && periodicity2 <= 6, ApiClassName, CallingMethodName,
            "periodicity2(%d) is invalid: we require 1 <= periodicity <= 6", periodicity2);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp2InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude2(%g) is not valid: must be nonnegative", amp2InKJ);*/

        SimTK_APIARGCHECK1_ALWAYS(0 <= phase2InDegrees && phase2InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle2(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase2InDegrees);

            // No repeats.
        SimTK_APIARGCHECK1_ALWAYS(periodicity3 != periodicity2, ApiClassName, CallingMethodName,
            "only one term with a given periodicity may be specified (periodicity %d was repeated)",
            periodicity2);
    }
    if (periodicity3 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity3 && periodicity3 <= 6, ApiClassName, CallingMethodName,
            "periodicity3(%d) is invalid: we require 1 <= periodicity <= 6", periodicity3);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1_ALWAYS(0 <= phase3InDegrees && phase3InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase3InDegrees);
            // (we've already checked for any possible repeats)
    }
    if (periodicity4 != -1) {
            // No nonsense.
        SimTK_APIARGCHECK1_ALWAYS(1 <= periodicity4 && periodicity4 <= 6, ApiClassName, CallingMethodName,
            "periodicity4(%d) is invalid: we require 1 <= periodicity <= 6", periodicity4);

        // GMOL Amber allows negative dihedral energy
        /*        SimTK_APIARGCHECK1_ALWAYS(amp3InKJ >= 0, ApiClassName, CallingMethodName,
            "amplitude3(%g) is not valid: must be nonnegative", amp3InKJ);*/

        SimTK_APIARGCHECK1_ALWAYS(0 <= phase4InDegrees && phase4InDegrees <= 180, ApiClassName, CallingMethodName,
            "phaseAngle3(%g) is not valid: must be between 0 and 180 degrees, inclusive", phase4InDegrees);
            // (we've already checked for any possible repeats)
    }


        // Canonicalize atom class quad by reversing order if necessary so that the
        // first class Index is numerically no larger than the fourth. Amber improper
        // torsions should not be canonicalized because order matters.
    const AtomClassIndexQuad key(class1, class2, class3, class4, shouldCanonicalizeClassOrder);

        // Attempt to create a new bond torsion entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexQuad,BondTorsion>::iterator, bool> ret = 
      torsionMap.insert(std::pair<AtomClassIndexQuad,BondTorsion>
        (key, BondTorsion(key)));

    BondTorsion& bondTorsionEntry = ret.first->second;

    // A new entry or one that just had a custom term in it won't have a built in
    // term so we can load it up and we're done.
    if (!bondTorsionEntry.hasBuiltinTerm()) {
        if (periodicity1 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
        if (periodicity2 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
        if (periodicity3 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
        if (periodicity4 != -1)
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity4, amp4InKJ, phase4InDegrees));
        return;
    }

    // If we get here we have discovered that there is already a built in torsion
    // term present for this atom class quad. We can still insert new terms, and we'll
    // allow duplicates if they are identical.

    if (periodicity1 != -1) {
        const TorsionTerm& term1 = bondTorsionEntry.getTermWithPeriod(periodicity1);
        if (term1.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term1.amplitude==amp1InKJ && term1.theta0==phase1InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity1);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity1, amp1InKJ, phase1InDegrees));
    }
    if (periodicity2 != -1) {
        const TorsionTerm& term2 = bondTorsionEntry.getTermWithPeriod(periodicity2);
        if (term2.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term2.amplitude==amp2InKJ && term2.theta0==phase2InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity2);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity2, amp2InKJ, phase2InDegrees));
    }
    if (periodicity3 != -1) {
        const TorsionTerm& term3 = bondTorsionEntry.getTermWithPeriod(periodicity3);
        if (term3.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term3.amplitude==amp3InKJ && term3.theta0==phase3InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity3);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
    }
    if (periodicity4 != -1) {
        const TorsionTerm& term4 = bondTorsionEntry.getTermWithPeriod(periodicity4);
        if (term4.isValid()) {
            SimTK_APIARGCHECK5_ALWAYS(term4.amplitude==amp4InKJ && term4.theta0==phase4InDegrees,
                ApiClassName, CallingMethodName,
                "atom class quad (%d,%d,%d,%d) already had a different term with periodicity %d",
                (int)class1,(int)class2,(int)class3,(int)class4,periodicity4);
        } else
            bondTorsionEntry.addBuiltinTerm(TorsionTerm(periodicity3, amp3InKJ, phase3InDegrees));
    }
}


// 
// We allow up to 3 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees)
{
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.defineAnyTorsion(class1, class2, class3, class4, true, // canonicalize 
                     periodicity1, amp1InKJ, phase1InDegrees,
                     periodicity2, amp2InKJ, phase2InDegrees,
                     periodicity3, amp3InKJ, phase3InDegrees,
                     mm.bondTorsion,
                     MethodName);
}

// 
// We allow up to 3 terms in a single torsion function, with three different
// periodicities. If any of these are unused, set the corresponding periodicity
// to -1.
//
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees,
    int periodicity4, Real amp4InKJ, Real phase4InDegrees)
{
    static const char* MethodName = "defineBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.defineAnyTorsion(class1, class2, class3, class4, true, // canonicalize 
                     periodicity1, amp1InKJ, phase1InDegrees,
                     periodicity2, amp2InKJ, phase2InDegrees,
                     periodicity3, amp3InKJ, phase3InDegrees,
                     periodicity4, amp4InKJ, phase4InDegrees,
                     mm.bondTorsion,
                     MethodName);
}


void DuMMForceFieldSubsystem::defineCustomBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    DuMM::CustomBondTorsion* customBondTorsion)
{
    static const char* MethodName = "defineCustomBondTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Watch for nonsense arguments.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class1), mm.ApiClassName, MethodName, 
        "class1=%d which is not a valid atom class Index", (int) class1);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class2), mm.ApiClassName, MethodName, 
        "class2=%d which is not a valid atom class Index", (int) class2);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class3=%d which is not a valid atom class Index", (int) class3);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtomClass(class3), mm.ApiClassName, MethodName, 
        "class4=%d which is not a valid atom class Index", (int) class4);
    SimTK_APIARGCHECK_ALWAYS(customBondTorsion, mm.ApiClassName, MethodName, 
        "CustomBondTorsion pointer was null");

        // Canonicalize atom class quad by reversing order if necessary so that the
        // first class Index is numerically no larger than the fourth.
    const AtomClassIndexQuad key(class1, class2, class3, class4, true);

        // Attempt to create a new bond torsion entry containing no valid
        // terms. If there was already an entry it will be returned instead
        // and no insertion is performed.
    std::pair<std::map<AtomClassIndexQuad,BondTorsion>::iterator, bool> ret = 
      mm.bondTorsion.insert(std::pair<AtomClassIndexQuad,BondTorsion>
        (key, BondTorsion(key)));

    BondTorsion& bondTorsionEntry = ret.first->second;
    bondTorsionEntry.addCustomTerm(customBondTorsion);
}

// Convenient signature for a bond torsion with only one term.
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      -1,0.,0., -1,0.,0.);
}

// Convenient signature for a bond torsion with two terms.
void DuMMForceFieldSubsystem::defineBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees)
{
    defineBondTorsion(class1, class2, class3, class4, 
                      periodicity1,amp1InKJ,phase1InDegrees,
                      periodicity2,amp2InKJ,phase2InDegrees,
                      -1,0.,0.);
}

// 
// This function is based on the defineBondTorsion function.
// As with the normal bond torsions, we allow up to 3 terms in a single torsion function,
// with three different periodicities. If any of these are unused, set the corresponding
// periodicity to -1.
//
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees,
    int periodicity3, Real amp3InKJ, Real phase3InDegrees)
{
    static const char* MethodName = "defineAmberImproperTorsion";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.defineAnyTorsion(class1, class2, class3, class4, false, // don't canonicalize 
                 periodicity1, amp1InKJ, phase1InDegrees,
                 periodicity2, amp2InKJ, phase2InDegrees,
                 periodicity3, amp3InKJ, phase3InDegrees,
                 mm.amberImproperTorsion,
                 MethodName);
}

// Convenient signature for an amber improper torsion with only one term.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, 
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees)
{
    defineAmberImproperTorsion(class1, class2, class3, class4,
                               periodicity1,amp1InKJ,phase1InDegrees,
                               -1,0.,0., -1,0.,0.);
}

// Convenient signature for an amber improper torsion with two terms.
void DuMMForceFieldSubsystem::defineAmberImproperTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    int periodicity1, Real amp1InKJ, Real phase1InDegrees,
    int periodicity2, Real amp2InKJ, Real phase2InDegrees)
{
    defineAmberImproperTorsion(class1, class2, class3, class4,
                               periodicity1,amp1InKJ,phase1InDegrees,
                               periodicity2,amp2InKJ,phase2InDegrees,
                               -1,0.,0.);
}

void DuMMForceFieldSubsystem::setVdwMixingRule(VdwMixingRule rule) {
    //static const char* MethodName = "setVdwMixingRule";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();
    mm.vdwMixingRule = rule; 
}

DuMMForceFieldSubsystem::VdwMixingRule 
DuMMForceFieldSubsystem::getVdwMixingRule() const {
    //static const char* MethodName = "getVdwMixingRule";
    const DuMMForceFieldSubsystemRep& mm = getRep();
    return mm.vdwMixingRule; 
}

const char*
DuMMForceFieldSubsystem::getVdwMixingRuleName(VdwMixingRule rule) const {
    static const char* MethodName = "getVdwMixingRuleName";
    switch(rule) {
    case WaldmanHagler:     return "Waldman-Hagler";
    case HalgrenHHG:        return "Halgren-HHG";        
    case Jorgensen:         return "Jorgensen";        
    case LorentzBerthelot:  return "Lorentz-Berthelot"; 
    case Kong:              return "Kong";          
    default:
        SimTK_APIARGCHECK1_ALWAYS(false, "DuMMForceFieldSubsystem", MethodName,
        "Unknown van der Waals mixing rule %d", (int)rule);
    };
}

void DuMMForceFieldSubsystem::setVdw12ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw12ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale12=fac;
}
void DuMMForceFieldSubsystem::setVdw13ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw13ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale13=fac;
}
void DuMMForceFieldSubsystem::setVdw14ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw14ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale14=fac;
}
void DuMMForceFieldSubsystem::setVdw15ScaleFactor(Real fac) {
    static const char* MethodName = "setVdw15ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "van der Waals energy scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.vdwScale15=fac;
}

void DuMMForceFieldSubsystem::setCoulomb12ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb12ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-2 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale12=fac;
}

void DuMMForceFieldSubsystem::setCoulomb13ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb13ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-3 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale13=fac;
}
void DuMMForceFieldSubsystem::setCoulomb14ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb14ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-4 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale14=fac;
}
void DuMMForceFieldSubsystem::setCoulomb15ScaleFactor(Real fac) {
    static const char* MethodName = "setCoulomb15ScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac && fac <= 1, mm.ApiClassName, MethodName,
        "Coulomb scale factor (%g) for 1-5 bonded atoms was invalid: must be between 0 and 1, inclusive",
        fac);

    mm.coulombScale15=fac;
}

void DuMMForceFieldSubsystem::setVdwGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setVdwScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global van der Waals scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.vdwGlobalScaleFactor=fac;
}

void DuMMForceFieldSubsystem::setCoulombGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCoulombScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global Coulomb scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.coulombGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondStretchGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondStretchScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond stretch scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondStretchGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondBendGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondBendScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond bend scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondBendGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setBondTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setBondTorsionScaleFactor";
 
    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global bond torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.bondTorsionGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setAmberImproperTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setAmberImproperTorsionScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global amber improper torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.amberImproperTorsionGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setCustomBondStretchGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCustomBondStretchScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global custom bond stretch scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.customBondStretchGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setCustomBondBendGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCustomBondBendScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global custom bond bend scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.customBondBendGlobalScaleFactor=fac;
}
void DuMMForceFieldSubsystem::setCustomBondTorsionGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setCustomBondTorsionScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global custom bond torsion scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.customBondTorsionGlobalScaleFactor=fac;
}


void DuMMForceFieldSubsystem::setNonbondedCutoff(Real cutoff) {
    static const char* MethodName = "setNonbondedCutoff";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= cutoff, mm.ApiClassName, MethodName,
                              "Nonbonded cutoff (nm) (%g) was invalid: must be nonnegative",
                              cutoff);

    mm.nonbondedCutoff=cutoff;
}


void DuMMForceFieldSubsystem::setNonbondedMethod(int methodIndex) {
    static const char* MethodName = "setNonbondedMethod";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= methodIndex && methodIndex <= 1, mm.ApiClassName, MethodName,
                              "Nonbonded method index should be 0 (NoCutoff) or 1 (CutoffNonPeriodic): (%g) was invalid",
                              methodIndex);

    mm.nonbondedMethod=methodIndex;
}



void DuMMForceFieldSubsystem::setSolventDielectric(Real dielectric) {
    static const char* MethodName = "setSolventDielectric";
    DuMMForceFieldSubsystemRep& mm = updRep();
    SimTK_APIARGCHECK1_ALWAYS(dielectric > 0, mm.ApiClassName, MethodName,
        "Solvent dielectric (%g) was invalid: must be greater than zero", dielectric);
    invalidateSubsystemTopologyCache();
    mm.gbsaSolventDielectric = dielectric;
}
void DuMMForceFieldSubsystem::setSoluteDielectric(Real dielectric) {
    static const char* MethodName = "setSolutetDielectric";
    DuMMForceFieldSubsystemRep& mm = updRep();
    SimTK_APIARGCHECK1_ALWAYS(dielectric > 0, mm.ApiClassName, MethodName,
        "Solute dielectric (%g) was invalid: must be greater than zero", dielectric);
    invalidateSubsystemTopologyCache();
    mm.gbsaSoluteDielectric = dielectric;
}
Real DuMMForceFieldSubsystem::getSolventDielectric() const {
    return getRep().gbsaSolventDielectric;
}
Real DuMMForceFieldSubsystem::getSoluteDielectric() const {
    return getRep().gbsaSoluteDielectric;
}

void DuMMForceFieldSubsystem::setGbsaIncludeAceApproximation(bool doInclude)
{
    //static const char* MethodName = "setGbsaIncludeAceApproximation";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    mm.gbsaIncludeAceApproximation=doInclude;
}

void DuMMForceFieldSubsystem::setGbsaGlobalScaleFactor(Real fac) {
    static const char* MethodName = "setGbsaGlobalScaleFactor";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(0 <= fac, mm.ApiClassName, MethodName,
        "Global generalized Born scale factor (%g) was invalid: must be nonnegative",
        fac);

    mm.gbsaGlobalScaleFactor=fac;
}

void DuMMForceFieldSubsystem::
clearIncludedNonbondAtomList() {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.includedNonbondAtoms.clear();
    inclList.includedNonbondBodies.clear();
    inclList.useDefaultNonbondList = false; // i.e., now there is nothing
}

void DuMMForceFieldSubsystem::
clearIncludedBondList() {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.atomsWhoseBondsAreIncluded.clear();
    inclList.atomPairsWhoseConnectingBondsAreIncluded.clear();
    inclList.bodiesWhoseBondsAreIncluded.clear();
    inclList.bodyPairsWhoseConnectingBondsAreIncluded.clear();
    inclList.useDefaultBondList = false; // i.e., now there is nothing
}

void DuMMForceFieldSubsystem::
resetIncludedNonbondAtomListToDefault() {
    clearIncludedNonbondAtomList();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = true;
}

void DuMMForceFieldSubsystem::
resetIncludedBondListToDefault() {
    clearIncludedBondList();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = true;
}

void DuMMForceFieldSubsystem::
includeNonbondAtom(DuMM::AtomIndex atomIx) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = false;
    inclList.includedNonbondAtoms.insert(atomIx); // ignores duplicates
}

void DuMMForceFieldSubsystem::
includeAllNonbondAtomsForOneBody(MobilizedBodyIndex mobodIx) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultNonbondList = false;
    inclList.includedNonbondBodies.insert(mobodIx); // ignores duplicates
}

void DuMMForceFieldSubsystem::
includeAllInterbodyBondsForOneAtom(DuMM::AtomIndex ax) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.atomsWhoseBondsAreIncluded.insert(ax); // ignores duplicates
}



void DuMMForceFieldSubsystem::
includeAllInterbodyBondsWithBothAtoms(DuMM::AtomIndex ax1, DuMM::AtomIndex ax2)
{   invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.atomPairsWhoseConnectingBondsAreIncluded.insert
        (AtomIndexPair(ax1,ax2,true)); // canonicalize order; ignore dups
}

void DuMMForceFieldSubsystem::
includeAllInterbodyBondsWithBothAtoms(DuMM::BondIndex bond) {
   includeAllInterbodyBondsWithBothAtoms(getBondAtom(bond,0),
                                         getBondAtom(bond,1));
}


void DuMMForceFieldSubsystem::
includeAllInterbodyBondsForOneBody(MobilizedBodyIndex mobod) {
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.bodiesWhoseBondsAreIncluded.insert(mobod); // ignores dups
}

void DuMMForceFieldSubsystem::
includeAllInterbodyBondsBetweenTwoBodies
   (MobilizedBodyIndex mobod1, MobilizedBodyIndex mobod2) 
{
    invalidateSubsystemTopologyCache();
    InclusionListSpec& inclList = updRep().inclList;
    inclList.useDefaultBondList = false;
    inclList.bodyPairsWhoseConnectingBondsAreIncluded
        .insert(MobodIndexPair(mobod1,mobod2,true)); // canonicalize order
}

int DuMMForceFieldSubsystem::
getNumIncludedAtoms() const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(subsystemTopologyHasBeenRealized(),
        "getNumIncludedAtoms", "Subsystem", "DuMMForceFieldSubsystem");
    return getRep().getNumIncludedAtoms();
}

DuMM::AtomIndex DuMMForceFieldSubsystem::
getAtomIndexOfIncludedAtom
   (DuMM::IncludedAtomIndex incAtomIndex) const 
{
    // Don't check in Release mode since this might get called a lot and
    // presumably we just checked in getNumIncludedAtoms().
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
        "getAtomIndexOfIncludedAtom", "Subsystem", "DuMMForceFieldSubsystem");

    return getRep().getAtomIndexOfIncludedAtom(incAtomIndex);
}

int DuMMForceFieldSubsystem::
getNumNonbondAtoms() const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(subsystemTopologyHasBeenRealized(),
        "getNumNonbondAtoms", "Subsystem", "DuMMForceFieldSubsystem");
    return getRep().getNumNonbondAtoms();
}

DuMM::IncludedAtomIndex DuMMForceFieldSubsystem::
getIncludedAtomIndexOfNonbondAtom
   (DuMM::NonbondAtomIndex nbAtomIndex) const 
{
    // Don't check in Release mode since this might get called a lot and
    // presumably we just checked in getNumNonbondAtoms().
    SimTK_STAGECHECK_TOPOLOGY_REALIZED(subsystemTopologyHasBeenRealized(),
        "getIncludedAtomIndexOfNonbondAtom", "Subsystem", "DuMMForceFieldSubsystem");

    return getRep().getIncludedAtomIndexOfNonbondAtom(nbAtomIndex);
}



Real DuMMForceFieldSubsystem::getVdwGlobalScaleFactor()     const {return getRep().vdwGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getCoulombGlobalScaleFactor() const {return getRep().coulombGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getGbsaGlobalScaleFactor()    const {return getRep().gbsaGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getBondStretchGlobalScaleFactor() const {return getRep().bondStretchGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getBondBendGlobalScaleFactor()    const {return getRep().bondBendGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getBondTorsionGlobalScaleFactor() const {return getRep().bondTorsionGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getAmberImproperTorsionGlobalScaleFactor() const {return getRep().amberImproperTorsionGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getCustomBondStretchGlobalScaleFactor() const {return getRep().customBondStretchGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getCustomBondBendGlobalScaleFactor()    const {return getRep().customBondBendGlobalScaleFactor;}
Real DuMMForceFieldSubsystem::getCustomBondTorsionGlobalScaleFactor() const {return getRep().customBondTorsionGlobalScaleFactor;}

Real DuMMForceFieldSubsystem::getNonbondedCutoff()     const {return getRep().nonbondedCutoff;}
int DuMMForceFieldSubsystem::getNonbondedMethod() const {return getRep().nonbondedMethod;}







void DuMMForceFieldSubsystem::setTracing(bool shouldTrace) {
    updRep().tracing = shouldTrace;
}

bool DuMMForceFieldSubsystem::getUseMultithreadedComputation() const
{   return getRep().useMultithreadedComputation; }

void DuMMForceFieldSubsystem::setUseMultithreadedComputation(bool use)
{   invalidateSubsystemTopologyCache();
    updRep().useMultithreadedComputation = use; }

bool DuMMForceFieldSubsystem::isUsingMultithreadedComputation() const
{   return getRep().usingMultithreaded; }

int DuMMForceFieldSubsystem::getNumThreadsRequested() const
{   return getRep().numThreadsRequested; }

void DuMMForceFieldSubsystem::setNumThreadsRequested(int nThreads)
{   invalidateSubsystemTopologyCache();
    updRep().numThreadsRequested = nThreads > 0 ? nThreads : 0; }

int DuMMForceFieldSubsystem::getNumThreadsInUse() const 
{   return getRep().numThreadsInUse; }


bool DuMMForceFieldSubsystem::getUseOpenMMAcceleration() const
{   return getRep().wantOpenMMAcceleration; }


void DuMMForceFieldSubsystem::setUseOpenMMAcceleration(bool use)
{   invalidateSubsystemTopologyCache();
    updRep().wantOpenMMAcceleration = use; }



bool DuMMForceFieldSubsystem::getUseOpenMMCalcOnlyNonBonded() const
{   return getRep().wantOpenMMCalcOnlyNonBonded; }

void DuMMForceFieldSubsystem::setUseOpenMMCalcOnlyNonBonded(bool use)
{   invalidateSubsystemTopologyCache();
    updRep().wantOpenMMCalcOnlyNonBonded = use; }



bool DuMMForceFieldSubsystem::getUseOpenMMIntegration() const
{   return getRep().wantOpenMMIntegration; }

void DuMMForceFieldSubsystem::setUseOpenMMIntegration(bool use)
{   invalidateSubsystemTopologyCache();
    updRep().wantOpenMMIntegration = use; }


Real DuMMForceFieldSubsystem::OMM_calcPotentialEnergy() const
{ return getRep().openMMPlugin.calcPotentialEnergy();} 

Real DuMMForceFieldSubsystem::OMM_calcKineticEnergy() const
{ return getRep().openMMPlugin.calcKineticEnergy();} 

void DuMMForceFieldSubsystem::OMM_integrateTrajectory( int steps )
{ return updRep().openMMPlugin.integrateTrajectory(steps);}  

SimTK::Vec3 DuMMForceFieldSubsystem::calcAtomLocationInGroundFrameThroughOMM( DuMM::AtomIndex daix ) const
{
    return getRep().openMMPlugin.getAtomPosition(daix);
}

const std::vector<OpenMM::Vec3>& DuMMForceFieldSubsystem::OMM_getPositions() const
{
    return getRep().openMMPlugin.getPositions();
}

float DuMMForceFieldSubsystem::getOpenMMstepsize() const
{   return getRep().stepsize; }

void DuMMForceFieldSubsystem::setOpenMMstepsize(float value)
{   invalidateSubsystemTopologyCache();
    updRep().stepsize = value; }

float DuMMForceFieldSubsystem::getOpenMMtemperature() const
{   return getRep().temperature; }

void DuMMForceFieldSubsystem::setOpenMMtemperature(float value)
{   
    invalidateSubsystemTopologyCache();
    updRep().temperature = value;
    updRep().openMMPlugin.setVelocitiesToTemperature(value);
    
    // std::cout<<"SETTING TEMPERATURE in DUMM "<<std::endl << getRep().temperature <<std::endl<< std::flush;
}

void DuMMForceFieldSubsystem::setOpenMMvelocities(float value)
{
    updRep().openMMPlugin.setVelocitiesToTemperature(value);
}

void DuMMForceFieldSubsystem::OMM_updatePositions(const std::vector<SimTK::Vec3>& positions)
{
    updRep().openMMPlugin.updatePositions(positions);
}

bool DuMMForceFieldSubsystem::getAllowOpenMMReference() const
{   return getRep().allowOpenMMReference; }

void DuMMForceFieldSubsystem::setAllowOpenMMReference(bool allow)
{   
    invalidateSubsystemTopologyCache();
    updRep().allowOpenMMReference = allow; }

bool DuMMForceFieldSubsystem::isUsingOpenMM() const 
{   return getRep().usingOpenMM; }

std::string DuMMForceFieldSubsystem::getOpenMMPlatformInUse() const {
    return getRep().openMMPlatformInUse;
}

DuMM::ClusterIndex DuMMForceFieldSubsystem::createCluster(const char* groupName)
{
    invalidateSubsystemTopologyCache();
    // Currently there is no error checking to do. We don't insist on unique group names.
    return updRep().addCluster(Cluster(groupName));
}

DuMM::AtomIndex DuMMForceFieldSubsystem::addAtom(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex)
{
    static const char* MethodName = "addAtom";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

    SimTK_APIARGCHECK1_ALWAYS(mm.isValidChargedAtomType(chargedAtomTypeIndex), mm.ApiClassName, MethodName, 
        "charged atom type %d is not valid", (int) chargedAtomTypeIndex);

    const DuMM::AtomIndex atomIndex = (const DuMM::AtomIndex)mm.atoms.size();
    mm.atoms.push_back(DuMMAtom(chargedAtomTypeIndex, atomIndex));
    TRACE((std::string("DuMM::addAtom") + std::to_string(atomIndex) + std::string("\n")).c_str());
    return atomIndex;
}

// EU COMMENT BEGIN
void DuMMForceFieldSubsystem::placeAtomInCluster(DuMM::AtomIndex atomIndex, DuMM::ClusterIndex clusterIndex, const Vec3& stationInNm)
// EU BEGIN
//void DuMMForceFieldSubsystem::placeAtomInCluster(DuMM::AtomIndex atomIndex, DuMM::ClusterIndex clusterIndex, Vec3 stationInNm)
// EU END
{
    static const char* MethodName = "placeAtomInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both the atomIndex and clusterIndex before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom index %d is not valid", (int) atomIndex);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
        "cluster index %d is not valid", (int) clusterIndex);

    Cluster& cluster = mm.updCluster(clusterIndex);

        // Make sure that this cluster doesn't already contain this atom, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!cluster.containsAtom(atomIndex), mm.ApiClassName, MethodName,
        "cluster %d('%s') already contains atom %d", (int) clusterIndex, cluster.name.c_str(), (int) atomIndex);

        // Add the atom to the cluster.
    cluster.placeAtom(atomIndex, stationInNm, mm);
}

void DuMMForceFieldSubsystem::placeClusterInCluster
   (DuMM::ClusterIndex childClusterIndex, DuMM::ClusterIndex parentClusterIndex, const Transform& placementInNm)
{
    static const char* MethodName = "placeClusterInCluster";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterIndex), mm.ApiClassName, MethodName,
        "child cluster Index %d is not valid", (int) childClusterIndex);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterIndex), mm.ApiClassName, MethodName,
        "parent cluster Index %d is not valid", (int) parentClusterIndex);

    Cluster&       parent = mm.updCluster(parentClusterIndex);
    const Cluster& child  = mm.getCluster(childClusterIndex);

        // TODO: for now, make sure the parent is a top-level cluster, meaning that it does
        // not have any parent clusters (although it can be attached to a body). This restriction
        // should be relaxed but it is tricky to get all the parents' and ancestors' content
        // lists updated correctly so I'm deferring that for now (sherm 060928).
    SimTK_APIARGCHECK2_ALWAYS(parent.isTopLevelCluster(), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') is not a top-level cluster so you cannot add a child cluster to it now",
        (int) parentClusterIndex, parent.name.c_str());

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(!child.isAttachedToBody(), mm.ApiClassName, MethodName,
        "child cluster %d('%s') is already attached to a body so cannot now be placed in another cluster",
        (int) childClusterIndex, child.name.c_str());

        // Make sure that parent cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK4_ALWAYS(!parent.containsCluster(childClusterIndex), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') already contains child cluster %d('%s')", 
        (int) parentClusterIndex, parent.name.c_str(), (int) childClusterIndex, child.name.c_str());

        // Make sure the new child cluster doesn't contain any atoms which are already in
        // any of the trees to which the parent cluster is associated.
        // TODO: for now we need only look at the parent since we know it is top level.
    DuMM::AtomIndex atomIndex;
    SimTK_APIARGCHECK5_ALWAYS(!parent.overlapsWithCluster(child, atomIndex), mm.ApiClassName, MethodName,
        "parent cluster %d('%s') and would-be child cluster %d('%s') both contain atom %d"
        " so they cannot have a parent/child relationship",
        (int) parentClusterIndex, parent.name.c_str(), (int) childClusterIndex, child.name.c_str(), (int) atomIndex);

        // Add the child cluster to the parent.
    parent.placeCluster(childClusterIndex, placementInNm, mm);
}

void DuMMForceFieldSubsystem::attachClusterToBody
   (DuMM::ClusterIndex clusterIndex, MobilizedBodyIndex mobodIx, 
    const Transform& placementInNm) 
{
    static const char* MethodName = "attachClusterToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this cluster before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), 
        mm.ApiClassName, MethodName,
        "cluster Index %d is not valid", (int) clusterIndex);
    SimTK_APIARGCHECK1_ALWAYS(mobodIx.isValid(), mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", (int)mobodIx);

    const Cluster& child  = mm.getCluster(clusterIndex);

        // Child must not already be attached to a body.
    SimTK_APIARGCHECK3_ALWAYS(!child.isAttachedToBody(), 
        mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached to body %d so cannot now be"
        " attached to a body",
        (int)clusterIndex, child.name.c_str(), (int)child.getMobodIndex());

        // None of the atoms in the child can be attached to any body.
    DuMM::AtomIndex    tempAtomIndex;
    MobilizedBodyIndex tempBodyIndex;
    SimTK_APIARGCHECK4_ALWAYS(
        !child.containsAnyAtomsAttachedToABody(tempAtomIndex,tempBodyIndex,mm), 
        mm.ApiClassName, MethodName,
        "cluster %d('%s') contains atom %d which is already attached to body %d"
        " so the cluster cannot now be attached to another body",
        (int)clusterIndex, child.name.c_str(), (int)tempAtomIndex, 
        (int)tempBodyIndex);

    // Create an entry for the body if necessary, and its corresponding cluster.
    DuMMBodyIndex duMMBodyIndex = mm.ensureDuMMBodyEntryExists(mobodIx);
    Cluster&      bodyCluster   = 
        mm.updCluster(mm.getDuMMBody(duMMBodyIndex).getClusterIndex());

        // Make sure that body cluster doesn't already contain child cluster, either directly
        // or recursively through its subclusters.
    SimTK_APIARGCHECK3_ALWAYS(!bodyCluster.containsCluster(clusterIndex), 
        mm.ApiClassName, MethodName,
        "cluster %d('%s') is already attached (directly or indirectly) to"
        " body %d", (int)clusterIndex, child.name.c_str(), (int)mobodIx);

        // OK, attach the cluster to the body's cluster.
    bodyCluster.placeCluster(clusterIndex, placementInNm, mm);
}

// EU COMMENT BEGIN
void DuMMForceFieldSubsystem::attachAtomToBody
   (DuMM::AtomIndex atomIndex, MobilizedBodyIndex bodyIndex, 
    const Vec3& stationInNm)
// EU BEGIN
//void DuMMForceFieldSubsystem::attachAtomToBody
//   (DuMM::AtomIndex atomIndex, MobilizedBodyIndex bodyIndex, 
//    Vec3 stationInNm)
// EU END
{
    static const char* MethodName = "attachAtomToBody";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this atom before, and that the body number is well formed.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom index %d is not valid", (int) atomIndex);
    SimTK_APIARGCHECK1_ALWAYS(bodyIndex.isValid(), mm.ApiClassName, MethodName,
        "body number %d is not valid: must be nonnegative", (int)bodyIndex);

        // The atom must not already be attached to a body, even this one.
    SimTK_APIARGCHECK2_ALWAYS(!mm.getAtom(atomIndex).isAttachedToBody(), 
        mm.ApiClassName, MethodName,
        "atom %d is already attached to body %d so cannot now be attached"
        " to a body",
        (int) atomIndex, (int)mm.getAtom(atomIndex).getMobodIndex());

        // Create an entry for the body if necessary, and its corresponding cluster.
    DuMMBodyIndex duMMBodyIndex = mm.ensureDuMMBodyEntryExists(bodyIndex);
    Cluster& bodyCluster = mm.updCluster(mm.getDuMMBody(duMMBodyIndex).getClusterIndex());

        // Attach the atom to the body's cluster.
    bodyCluster.placeAtom(atomIndex, stationInNm, mm);
}

MassProperties DuMMForceFieldSubsystem::calcClusterMassProperties
   (DuMM::ClusterIndex clusterIndex, const Transform& placementInNm) const
{
    static const char* MethodName = "calcClusterMassProperties";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
        "cluster Index %d is not valid", (int) clusterIndex);

    return mm.getCluster(clusterIndex).calcMassProperties(placementInNm, mm);
}


DuMM::BondIndex DuMMForceFieldSubsystem::addBond(DuMM::AtomIndex atom1Ix, DuMM::AtomIndex atom2Ix)
{
    static const char* MethodName = "addBond";

    invalidateSubsystemTopologyCache();

    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen these atoms before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom1Ix), mm.ApiClassName, MethodName,
        "atom1(%d) is not valid", (int) atom1Ix);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atom2Ix), mm.ApiClassName, MethodName,
        "atom2(%d) is not valid", (int) atom2Ix);

        // An atom can't be bonded to itself.
    SimTK_APIARGCHECK1_ALWAYS(atom1Ix != atom2Ix, mm.ApiClassName, MethodName,
        "the same atom index (%d) was given for both atoms, which makes no sense", (int) atom1Ix);

    // Ensure that atom1 < atom2
    if (atom1Ix > atom2Ix)
        std::swap(atom1Ix,atom2Ix);

    DuMMAtom& a1 = mm.updAtom(atom1Ix);
    DuMMAtom& a2 = mm.updAtom(atom2Ix);

    SimTK_APIARGCHECK2_ALWAYS(!a1.isBondedTo(atom2Ix), mm.ApiClassName, MethodName,
        "atom %d is already bonded to atom %d; you can only do that once",
        (int) atom1Ix, (int) atom2Ix);

    mm.bonds.push_back(Bond(atom1Ix,atom2Ix));
    a1.bond12.push_back(atom2Ix);
    a2.bond12.push_back(atom1Ix);
    return (DuMM::BondIndex)(mm.bonds.size() - 1);
}

int DuMMForceFieldSubsystem::getNumAtoms() const {
    return getRep().getNumAtoms();
}
int DuMMForceFieldSubsystem::getNumBonds() const {
    return getRep().getNumBonds();
}

// 'which' is 0 or 1 to pick one of the two atoms whose index we return.
DuMM::AtomIndex DuMMForceFieldSubsystem::getBondAtom(DuMM::BondIndex bondIx, int which) const {
    static const char* MethodName = "getBondAtom";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this bond before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidBond(bondIx), mm.ApiClassName, MethodName,
        "bond %d is not valid", (int) bondIx);

    SimTK_APIARGCHECK1_ALWAYS(which==0 || which==1, mm.ApiClassName, MethodName,
        "'which' was %d but must be 0 or 1 to choose one of the two atoms", which);

    return mm.bonds[bondIx].atoms[which];
}

// Returned mass is in daltons (g/mol).
Real DuMMForceFieldSubsystem::getAtomMass(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomMass";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    const Element e = Element::getByAtomicNumber(mm.getAtomElementNum(atomIndex));
    return e.getMass();
}

// Returns the atomic number (number of protons in nucleus).
int DuMMForceFieldSubsystem::getAtomElement(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomElement";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    return mm.getAtomElementNum(atomIndex);
}

// Last vestige of Element class that was internal to DuMMForceFieldSubsystem.
// I don't want to add a "color" field to the external Element class.
// Properly, there should be an ElementColorer class or something.
Vec3 DuMMForceFieldSubsystem::getElementDefaultColor(int atomicNumber) const {
    switch (atomicNumber) {
        case 1: // hydrogen
            return Green;
        case 7: // nitrogen
            return Blue;
        case 8: // oxygen
            return Red;
        case 15: // phosphorus
            return Magenta;
        case 16: // sulfur
            return Yellow;
        case 79: // gold
            return Yellow;
        default:
            return Gray;
    }
}

Vec3 DuMMForceFieldSubsystem::getAtomDefaultColor(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomDefaultColor";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    return getElementDefaultColor(mm.getAtomElementNum(atomIndex));
}

// Returned radius is in nm.
Real DuMMForceFieldSubsystem::getAtomRadius(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomRadius";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    const AtomClass& cl = mm.atomClasses[mm.getAtomClassIndex(atomIndex)];
    return cl.vdwRadius;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationOnBody(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomStationOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    const DuMMAtom& a = mm.getAtom(atomIndex);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomIndex);

    return a.station_B;
}


// EU BEGIN
/// Set the station at which a particular atom is fixed on its body. 
/// An exception will be thrown if this atom is not fixed to any body.
void DuMMForceFieldSubsystem::bsetAtomStationOnBody(DuMM::AtomIndex atomIndex, Vec3 new_station_B){
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    DuMMAtom& a = mm.updAtom(atomIndex);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomIndex);

    a.station_B = new_station_B;
}

/// Set the station at which a particular atom is fixed on its body. 
/// An exception will be thrown if this atom is not fixed to any body.
void DuMMForceFieldSubsystem::bsetAllAtomStationOnBody(DuMM::AtomIndex atomIndex, Vec3 new_station_B){
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

        // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom %d is not valid", (int) atomIndex);

    DuMMAtom& a = mm.updAtom(atomIndex);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomIndex);

    a.station_B_All = new_station_B;
}

// Atom placements in clusters used in realizeSubsytemTopology
void DuMMForceFieldSubsystem::bsetAtomPlacementStation(DuMM::AtomIndex atomIndex, MobilizedBodyIndex inputMbx, Vec3 new_station){
  DuMMForceFieldSubsystemRep& mm = updRep();

  //Cluster& cluster = mm.clusters[DuMM::ClusterIndex(0)];
  //(cluster.allAtomPlacements.begin())->setStation(Vec3(1, 2, 3));
  //std::cout<<"cluster station after "<<(cluster.allAtomPlacements.begin())->station<<std::endl;

  //AtomPlacement ap(DuMM::AtomIndex(0), Vec3(0));
  //ap.station = Vec3(1);
  for (DuMMBodyIndex bnum(0); bnum < mm.duMMSubsetOfBodies.size(); ++bnum) {
    DuMMBody& b = mm.duMMSubsetOfBodies[bnum];
    const MobodIndex mbx = b.getMobilizedBodyIndex();
    if(mbx == inputMbx){ 
      const Cluster& cluster = mm.getCluster(b.clusterIndex);
      
      for (AtomPlacementSet::iterator app = cluster.getAllContainedAtoms().begin();
        app != cluster.getAllContainedAtoms().end();
        ++app)
      {
        assert(app->isValid());
        if(app->atomIndex == atomIndex){
          //std::cout<<"bnum "<<bnum<<" b.cIx "<<b.clusterIndex<<" b.mbx "<<b.mobilizedBodyIndex<<" app.station "<<app->station;
          app->setStation(new_station);
          //std::cout<<" app.station "<<app->station<<std::endl; fflush(stdout);
          break;
        }
      }
      
    }
  }
  
}


// Stations computed every time
Vec3& DuMMForceFieldSubsystem::updIncludedAtomStation(DuMM::AtomIndex atomIndex){
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
      "atom %d is not valid", (int) atomIndex);

    DuMMAtom& a = mm.updAtom(atomIndex);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
      "atom %d is not attached to a body", (int) atomIndex);

    return mm.updIncludedAtomStation(a.inclAtomIndex);
}

// Stations computed every time - CalcFullPotential
Vec3& DuMMForceFieldSubsystem::updAllAtomStation(DuMM::AtomIndex atomIndex){
    static const char* MethodName = "getAtomStationOnBody";
    DuMMForceFieldSubsystemRep& mm = updRep();

    // Make sure we've seen this atom before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
      "atom %d is not valid", (int) atomIndex);

    DuMMAtom& a = mm.updAtom(atomIndex);

    // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
      "atom %d is not attached to a body", (int) atomIndex);

    return mm.updAllAtomStation(a.inclAtomIndex);
}

// Get ClusterIndex corresponding to specified Mobod
DuMM::ClusterIndex DuMMForceFieldSubsystem::bgetMobodClusterIndex(MobilizedBodyIndex inputMbx) const {
  //static const char* MethodName = "bgetMobodClusterIndex";
  const DuMMForceFieldSubsystemRep& mm = getRep();

  for (DuMMBodyIndex bnum(0); bnum < mm.duMMSubsetOfBodies.size(); ++bnum) {
    const DuMMBody& b = mm.duMMSubsetOfBodies[bnum];
    const MobodIndex mbx = b.getMobilizedBodyIndex();
    if(mbx == inputMbx){ 
      // Make sure that we've seen this clusters before.
      //SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
      //  "cluster Index %d is not valid", (int) b.clusterIndex);
      return b.clusterIndex;
    }
  }

  // Should never get here, but the compiler keeps complaining.
  assert(false);
  return {};
}
// EU END



// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementOnBody(DuMM::ClusterIndex clusterIndex) const {
    static const char* MethodName = "getClusterPlacementOnBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure we've seen this cluster before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
        "cluster Index %d is not valid", (int) clusterIndex);

    const Cluster& c = mm.getCluster(clusterIndex);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", (int) clusterIndex, c.name.c_str());

    return c.placement_B;
}

// Returned station is in nm.
Vec3 DuMMForceFieldSubsystem::getAtomStationInCluster(DuMM::AtomIndex atomIndex, DuMM::ClusterIndex clusterIndex) const {
    static const char* MethodName = "getAtomStationInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both the atomIndex and clusterIndex before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom index %d is not valid", (int) atomIndex);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
        "cluster index %d is not valid", (int) clusterIndex);

    const Cluster& c = mm.getCluster(clusterIndex);
    const AtomPlacementSet& atoms = c.getAllContainedAtoms();
    const AtomPlacementSet::const_iterator ap = 
        atoms.find(AtomPlacement(atomIndex,Vec3(0)));

        // We're going to be upset of this cluster doesn't contain this atom.
    SimTK_APIARGCHECK3_ALWAYS(ap != atoms.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain atom %d", (int) clusterIndex, c.name.c_str(), (int) atomIndex);

    return ap->station;
}

// Returned placement is in nm.
Transform DuMMForceFieldSubsystem::getClusterPlacementInCluster(DuMM::ClusterIndex childClusterIndex, DuMM::ClusterIndex parentClusterIndex) const {
    static const char* MethodName = "getClusterPlacementInCluster";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen both of these clusters before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(childClusterIndex), mm.ApiClassName, MethodName,
        "child cluster Index %d is not valid", (int) childClusterIndex);
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(parentClusterIndex), mm.ApiClassName, MethodName,
        "parent cluster Index %d is not valid", (int) parentClusterIndex);

    const Cluster& parent = mm.getCluster(parentClusterIndex);
    const Cluster& child  = mm.getCluster(childClusterIndex);

    const ClusterPlacementSet& clusters = parent.getAllContainedClusters();
    const ClusterPlacementSet::const_iterator cp = 
        clusters.find(ClusterPlacement(childClusterIndex,Transform()));

        // We're going to be upset of the parent cluster doesn't contain the child.
    SimTK_APIARGCHECK4_ALWAYS(cp != clusters.end(), mm.ApiClassName, MethodName,
        "cluster %d('%s') does not contain cluster %d('%d')", 
        (int) parentClusterIndex, parent.name.c_str(), (int) childClusterIndex, child.name.c_str());

    return cp->placement;
}

MobilizedBodyIndex DuMMForceFieldSubsystem::getAtomBody(DuMM::AtomIndex atomIndex) const {
    static const char* MethodName = "getAtomBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomIndex before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidAtom(atomIndex), mm.ApiClassName, MethodName,
        "atom index %d is not valid", (int) atomIndex);

    const DuMMAtom& a = mm.getAtom(atomIndex);

        // Atom must be attached to a body.
    SimTK_APIARGCHECK1_ALWAYS(a.isAttachedToBody(), mm.ApiClassName, MethodName,
        "atom %d is not attached to a body", (int) atomIndex);

    return a.getMobodIndex();
}


MobilizedBodyIndex DuMMForceFieldSubsystem::getClusterBody(DuMM::ClusterIndex clusterIndex) const {
    static const char* MethodName = "getClusterBody";
    const DuMMForceFieldSubsystemRep& mm = getRep();

        // Make sure that we've seen this atomIndex before.
    SimTK_APIARGCHECK1_ALWAYS(mm.isValidCluster(clusterIndex), mm.ApiClassName, MethodName,
        "cluster Index %d is not valid", (int) clusterIndex);

    const Cluster& c = mm.getCluster(clusterIndex);

        // Cluster must be attached to a body.
    SimTK_APIARGCHECK2_ALWAYS(c.isAttachedToBody(), mm.ApiClassName, MethodName,
        "cluster %d('%s') is not attached to a body", (int) clusterIndex, c.name.c_str());

    return c.getMobodIndex();
}

void DuMMForceFieldSubsystem::dump() const {
    return getRep().dump();
}



// How many times has the forcefield been evaluated?
long long DuMMForceFieldSubsystem::getForceEvaluationCount() const
{
	return getRep().getForceEvaluationCount();
}

std::ostream& DuMMForceFieldSubsystemRep::generateBiotypeChargedAtomTypeSelfCode(std::ostream& os) const 
{
    std::map<BiotypeIndex, DuMM::ChargedAtomTypeIndex>::const_iterator i;
    for (i = chargedAtomTypesByBiotype.begin(); i != chargedAtomTypesByBiotype.end(); ++i)
    {
        generateBiotypeChargedAtomTypeSelfCode(os, i->first);
    }

    return os;
}

std::ostream& DuMMForceFieldSubsystem::generateBiotypeChargedAtomTypeSelfCode(std::ostream& os) const 
{
    return getRep().generateBiotypeChargedAtomTypeSelfCode(os);
}

void DuMMForceFieldSubsystem::setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, BiotypeIndex biotypeIx) 
{
    updRep().setBiotypeChargedAtomType(chargedAtomTypeIndex, biotypeIx);
}

DuMM::ChargedAtomTypeIndex DuMMForceFieldSubsystem::getBiotypeChargedAtomType(BiotypeIndex biotypeIx) const {
    return getRep().getBiotypeChargedAtomType(biotypeIx);
}

void DuMMForceFieldSubsystem::loadAmber99Parameters() 
{
    Biotype::initializePopularBiotypes();

    populateAmber99Params(*this);
}

void DuMMForceFieldSubsystem::loadTestMoleculeParameters()
{
    Biotype::initializePopularBiotypes();

    // TODO - these hard-coded chargedAtomTypeIndexs are not too cool
    // TODO - these charges are made up
    defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
    defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), "Methane H",  DuMM::AtomClassIndex(34),  -0.01);
    setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());
    setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), Biotype::MethaneH().getIndex());

    defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), "Ethane C",   DuMM::AtomClassIndex(1),  0.03);
    defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), "Ethane H",  DuMM::AtomClassIndex(34),  -0.01);
    setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), Biotype::EthaneC().getIndex());
    setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), Biotype::EthaneH().getIndex());
}


void DuMMForceFieldSubsystem::populateFromTinkerParameterFile(std::istream& tinkerStream) 
{

    //////////////////////////////////////////////////////
    // 1) Read Tinker parameter file one line at a time //
    //////////////////////////////////////////////////////

    Real radiusSizeScale = 1.0; // set to 0.5 for diameter parameter sets vs. 1.0 for radius
    Real radiusTypeScale = 1.0; // set to 2^(1/6) for sigma(R0) vs. 1.0 for Rmin

    // Bookkeeping to retrieve element and valence during biotype definition
    std::map<DuMM::AtomClassIndex, int> atomicNumberByAtomClassIndex;
    std::map<DuMM::AtomClassIndex, int> valenceByAtomClassIndex;
    std::map<DuMM::ChargedAtomTypeIndex, DuMM::AtomClassIndex> atomClassIdByChargedAtomTypeIndex;

    std::string tinkerParamFileLine;
    while( getline( tinkerStream, tinkerParamFileLine ) ) {
        std::stringstream lineStream(tinkerParamFileLine);

        // Get the first word on the line of text from the file
        // Significant words include "atom", "vdw", "bond", etc.
        String recordType;
        lineStream >> recordType;

        // forcefield name
        if (recordType == "forcefield") {
            lineStream >> updRep().forcefieldName;
        }

        // VDWTYPE [LENNARD-JONES/BUCKINGHAM/BUFFERED-14-7/MM3-HBOND/GAUSSIAN]
        else if (recordType == "vdwtype") {
            String vdwType;
            lineStream >> vdwType;

            if (vdwType == "LENNARD-JONES") ; // OK
            else { // DuMMForcefieldSubsystem doesn't know about other vdw models
                SimTK_THROW1( Exception::Cant,"Parse Exception: Can't use van der Waals model other than LENNARD-JONES" );
            }
        }

        // RADIUSRULE [ARITHMETIC/GEOMETRIC/CUBIC-MEAN]
        else if (recordType == "radiusrule") {
            String rule;
            lineStream >> rule;

            if (rule == "ARITHMETIC")
                setVdwMixingRule(LorentzBerthelot);
            else if (rule == "GEOMETRIC")
                setVdwMixingRule(Jorgensen);
            else if (rule == "CUBIC-MEAN")
                setVdwMixingRule(HalgrenHHG);
            else { // DuMMForcefieldSubsystem doesn't know about other vdw models
                SimTK_THROW1( Exception::Cant,"Parse Exception: Unrecognized radius rule" );
            }
        }

        // RADIUSTYPE [R-MIN/SIGMA]
        else if (recordType == "radiustype") {
            String radiusType;
            lineStream >> radiusType;
            
            if (radiusType == "R-MIN")
                radiusTypeScale = 1.0;
            else if (radiusType == "SIGMA") 
                radiusTypeScale = pow(2.0, (1.0/6.0));
        }

        else if (recordType == "radiussize") {
            String radiusSize;
            lineStream >> radiusSize;

            if (radiusSize == "RADIUS")
                radiusSizeScale = 1.0;
            else if (radiusSize == "DIAMETER")
                radiusSizeScale = 0.5;
            else {
                SimTK_THROW1( Exception::Cant, "Parse Error: unrecognized radius size" );
            }
        }

        // EPSILONRULE [GEOMETRIC/ARITHMETIC/HARMONIC/HHG]
        // TODO - DuMM currently lumps this with RADIUSRULE
        else if (recordType == "epsilonrule") {
            String epsilonRule; 
            lineStream >> epsilonRule;

            if (epsilonRule == "GEOMETRIC") {
                if (getVdwMixingRule() == LorentzBerthelot) continue; // already geometric
                else if (getVdwMixingRule() == Jorgensen) continue; // already geometric
                else setVdwMixingRule(Jorgensen);
            }
            else if (epsilonRule == "HHG") {
                setVdwMixingRule(HalgrenHHG);
            }
            else {
                SimTK_THROW1( Exception::Cant, "Parse Error: unrecognized epsilon rule" );
            }
        }

        else if (recordType == "vdw-14-scale") {
            Real vdw14Scale;
            lineStream >> vdw14Scale;

            setVdw14ScaleFactor(1.0 / vdw14Scale);
        }

        else if (recordType == "chg-14-scale") {
            Real chg14Scale;
            lineStream >> chg14Scale;

            setCoulomb14ScaleFactor(1.0 / chg14Scale);
        }

        else if (recordType == "dielectric") {
            Real dielectric;
            lineStream >> dielectric;

            if (dielectric != 1.0) {
                SimTK_THROW1( Exception::Cant, "Can't use dielectric other than 1.0" );
            }
        }

        // "atom" records
        // 'atom      1    14    N       "Glycine N"                 7     14.010     3'
        else if (recordType == "atom") 
        {
            int integer;
            lineStream >> integer;
            DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId(integer);

            String atomClassName;
            lineStream >> atomClassName;

            // Use getline() to get field between quotation marks
            std::string chargedAtomTypeName;
            // First try grabs spaces
            std::getline(lineStream, chargedAtomTypeName, '"');
            // Second try grabs string
            std::getline(lineStream, chargedAtomTypeName, '"');

            int elementNumber; 
            lineStream >> elementNumber;

            Real atomMass;
            lineStream >> atomMass;

            int valence;
            lineStream >> valence;

            // we don't yet know vdwRadius and vdwWellDepth
            if (!isValidAtomClass(atomClassId)) {
                defineIncompleteAtomClass_KA(
                    atomClassId, 
                    atomClassName.c_str(), 
                    elementNumber, 
                    valence);

            }
            atomicNumberByAtomClassIndex[atomClassId] = elementNumber;
            valenceByAtomClassIndex[atomClassId] = valence;

            // we don't yet know atomic partial charge
            defineIncompleteChargedAtomType_KA(
                chargedAtomTypeIndex, 
                chargedAtomTypeName.c_str(),
                atomClassId);

            atomClassIdByChargedAtomTypeIndex[chargedAtomTypeIndex] = atomClassId;
        }
       
        // "vdw" records, e.g.
        // 'vdw          1              1.9080     0.1094'
        // RecordType   AtomClass      Radius     WellDepth
        else if (recordType == "vdw") 
        {
            int integer;
            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId(integer);

            Real radius;
            lineStream >> radius;
            radius *= radiusSizeScale;
            radius *= radiusTypeScale;

            Real wellDepth;
            lineStream >> wellDepth;

            setAtomClassVdwParameters_KA(atomClassId, radius, wellDepth);
        }

        // "bond" records
        // bond stretching parameters
        // 'bond         1   22          320.0     1.4100'
        else if (recordType == "bond")
        {
            int integer;
            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId1(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId2(integer);

            Real stiffness;
            lineStream >> stiffness;

            Real length;
            lineStream >> length;

            if ( isValidAtomClass(atomClassId1) && isValidAtomClass(atomClassId2) )
                defineBondStretch_KA(atomClassId1, atomClassId2, stiffness, length);
        }

        // "angle" records
        // angle bending parameters
        // 'angle       10    1   34     50.00     109.50'
        else if (recordType == "angle") {
            int integer;
            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId1(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId2(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId3(integer);

            Real stiffness;
            lineStream >> stiffness;

            Real angle;
            lineStream >> angle;

            defineBondBend_KA(atomClassId1, atomClassId2, atomClassId3, stiffness, angle);
        }

        // improper torsions
        // imptors      1   14    2   24         10.500  180.0  2
        else if (recordType == "imptors") 
        {
            int integer;
            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId1(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId2(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId3(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId4(integer);

            Real amplitudeInKcal;
            lineStream >> amplitudeInKcal;

            Real phaseAngleInDegrees;
            lineStream >> phaseAngleInDegrees;

            int periodicity;
            lineStream >> periodicity;

            defineAmberImproperTorsion_KA
               (atomClassId1, atomClassId2, atomClassId3, atomClassId4,
                periodicity, amplitudeInKcal, phaseAngleInDegrees);
        }

        // "torsion" records
        // 'torsion      1    1    1    1     0.200 180.0 1   0.250 180.0 2   0.180   0.0 3'
        // OR
        // 'torsion     23    1    1   23          0.144    0.0  3      1.175    0.0  2'
        // OR
        // 'torsion      1    1    1    2          0.156    0.0  3'
        else if (recordType == "torsion") 
        {
            int integer;
            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId1(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId2(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId3(integer);

            lineStream >> integer;
            DuMM::AtomClassIndex atomClassId4(integer);

            // figure out how many fields are on this line
            int numberOfFields = 0;
            size_t pos = 0;
            while(pos != String::npos) {
                pos = tinkerParamFileLine.find_first_not_of(" ", pos);
                pos = tinkerParamFileLine.find_first_of(" ", pos);
                ++numberOfFields;
            }

            assert(numberOfFields >= 8);

            Real amplitude1;
            lineStream >> amplitude1;

            Real phase1;
            lineStream >> phase1;

            int periodicity1;
            lineStream >> periodicity1;

            if (numberOfFields == 8) {
                defineBondTorsion_KA(
                    atomClassId1, atomClassId2, atomClassId3, atomClassId4,
                    periodicity1, amplitude1, phase1
                    );
                continue;
            }

            assert(numberOfFields >= 11);

            Real amplitude2;
            lineStream >> amplitude2;

            Real phase2;
            lineStream >> phase2;

            int periodicity2;
            lineStream >> periodicity2;

            if (numberOfFields == 11) {
                defineBondTorsion_KA(
                    atomClassId1, atomClassId2, atomClassId3, atomClassId4,
                    periodicity1, amplitude1, phase1,
                    periodicity2, amplitude2, phase2
                    );
                continue;
            }

            assert(numberOfFields == 14);

            Real amplitude3;
            lineStream >> amplitude3;

            Real phase3;
            lineStream >> phase3;

            int periodicity3;
            lineStream >> periodicity3;

            defineBondTorsion_KA(
                atomClassId1, atomClassId2, atomClassId3, atomClassId4,
                periodicity1, amplitude1, phase1,
                periodicity2, amplitude2, phase2,
                periodicity3, amplitude3, phase3
                );
            continue;

        }

        // atom partial charge records, e.g.
        // 'charge       1       -0.4157'
        // RecordType AtomChargedType Charge
        else if ( recordType == "charge" ) {
            int integer;
            lineStream >> integer;
            DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex(integer);

            Real charge;
            lineStream >> charge;

            setChargedAtomTypeCharge(chargedAtomTypeIndex, charge);
        }

        // "biotype" records, e.g.
        // 'biotype      1    N       "Glycine"                   1'
        // RecordType Biotype AtomName ResidueName          AtomClass
        else if ( recordType == "biotype" ) 
        {
            int integer;
            lineStream >> integer;
            TinkerBiotypeIndex tinkerBiotypeIndex(integer);

            String atomName;
            lineStream >> atomName;

            String residueName;
            // Use getline() to get field between quotation marks
            // First try grabs spaces
            std::getline(lineStream, residueName, '"');
            // Second try grabs string
            std::getline(lineStream, residueName, '"');

            lineStream >> integer;
            DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex(integer);

            // Resolve element and valence
            DuMM::AtomClassIndex atomClassId = atomClassIdByChargedAtomTypeIndex[chargedAtomTypeIndex];
            int atomicNumber = atomicNumberByAtomClassIndex[atomClassId];
            int valence = valenceByAtomClassIndex[atomClassId];

            // Resolve ordinality
            // Parse residueName into residueName and context
            // TODO - include nucleotide nomenclature, if necessary
            Ordinality::Residue ordinality = Ordinality::Any;

            // Separate residue name from ordinality indicator
            // case of "Acetyl N-Terminus"
            // Handle "N-Terminal GLY"
            auto pos1 = residueName.find("N-Term");
            if (pos1 != std::string::npos) {
                ordinality = Ordinality::Initial;

                // e.g. "N-Terminal GLY"
                if (pos1 == 0) 
                    residueName.replace(pos1, 11, "");
                // e.g. "Acetyl N-Terminus"
                else 
                    residueName.replace(pos1-1, 11, "");
            }

            pos1 = residueName.find("C-Term");
            if (pos1 != std::string::npos) {
                ordinality = Ordinality::Final;

                // e.g. "C-Terminal GLY"
                if (pos1 == 0) 
                    residueName.replace(pos1, 11, "");
                // e.g. "N-MeAmide C-Terminus"
                else 
                    residueName.replace(pos1-1, 11, "");
            }
            
            // Nucleic acid biotypes, e.g.
            // "5'-Hydroxyl, RNA"
            // "5'-Phosphate OS, DNA"
            pos1 = residueName.find("5'-");
            if (pos1 != std::string::npos) {
            	ordinality = Ordinality::Initial;
            	residueName.replace(pos1, 3, "");
            	
                // get rid of that annoying wrongly placed atom name in the terminal phosphates
            	pos1 = residueName.find(" OP,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 4, ",");
            	pos1 = residueName.find(" OS,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 4, ",");
            	pos1 = residueName.find(" P,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 3, ",");
            	
            }
            pos1 = residueName.find("3'-");
            if (pos1 != std::string::npos) {
            	ordinality = Ordinality::Final;
            	residueName.replace(pos1, 3, "");
            	
                // get rid of that annoying wrongly placed atom name in the terminal phosphates
            	pos1 = residueName.find(" OP,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 4, ",");
            	pos1 = residueName.find(" OS,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 4, ",");
            	pos1 = residueName.find(" P,");
            	if (pos1 != std::string::npos) residueName.replace(pos1, 3, ",");
            }

            BiotypeIndex biotypeIx;

            if ( Biotype::exists(residueName.c_str(), atomName.c_str(), ordinality) ) {
                Biotype& biotype = Biotype::upd( residueName.c_str(), atomName.c_str(), ordinality );
                biotypeIx = biotype.getIndex();

                biotype.setTinkerBiotypeIndex(tinkerBiotypeIndex);
            }

            else biotypeIx = Biotype::defineTinkerBiotype(tinkerBiotypeIndex
                                                          , Element::getByAtomicNumber(atomicNumber)
                                                          , valence
                                                          , residueName.c_str()
                                                          , atomName.c_str()
                                                          , ordinality
            );

            setBiotypeChargedAtomType(chargedAtomTypeIndex, biotypeIx);
        }

    }

}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//		        GMOLMODEL - EXTRA FUNCTIONALITIES
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


Real DuMMForceFieldSubsystem::
CalcFullPotEnergyIncludingRigidBodies ( const State& state ) const {
    static const char* MethodName = "CalcFullPotEnergyIncludingRigidBodies";
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(subsystemTopologyHasBeenRealized(),
                                              MethodName, "Subsystem", "DuMMForceFieldSubsystem");
    return getRep().CalcFullPotEnergyIncludingRigidBodiesRep( state );
}
