#ifndef SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_REP_H_
#define SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_REP_H_

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
 * This header declares various internal classes used by the private
 * implementation of DuMMForceFieldSubsystem.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "molmodel/internal/DuMMForceFieldSubsystem.h"
#include "molmodel/internal/Element.h"
#include "molmodel/internal/Index.hpp"
#include "molmodel/internal/MolecularMechanicsSystem.h"
#include "molmodel/internal/common.h"
#include "simbody/internal/ForceSubsystemGuts.h"

#include "SimTKsimbody.h"


using namespace SimTK;


//-----------------------------------------------------------------------------
//                               ATOM CLASS
//-----------------------------------------------------------------------------
class AtomClass {
    public:
    AtomClass()
        : element(-1)
        , valence(-1)
        , vdwRadius(-1)
        , vdwWellDepth(-1) {
    }
    AtomClass(int id, const char* nm, int e, int v, Real radInNm, Real wellDepthInKJ)
        : atomClassIx(id)
        , name(nm)
        , element(e)
        , valence(v)
        , vdwRadius(radInNm)
        , vdwWellDepth(wellDepthInKJ) {
        // Permit incomplete construction, i.e. radius and depth not yet set
        assert(isValid());
    }

    bool isValid() const {
        return atomClassIx.isValid() && element > 0 && valence >= 0;
    }

    bool isComplete() const {
        return isValid() && vdwRadius >= 0 && vdwWellDepth >= 0;
    }

    void invalidateTopologicalCache() {
        vdwDij.clear();
        vdwEij.clear();
    }

    void dump() const {
        printf("   %d(%s): element=%d, valence=%d vdwRad=%g nm, vdwDepth(kJ)=%g (%g kcal)\n",
               (int)atomClassIx,
               name.c_str(),
               element,
               valence,
               vdwRadius,
               vdwWellDepth,
               vdwWellDepth * DuMM::KJ2Kcal);
        printf("    vdwDij (nm):");
        for (int i = 0; i < (int)vdwDij.size(); ++i) {
            printf(" %g", vdwDij[i]);
        }
        printf("\n    vdwEij (kJ):");
        for (int i = 0; i < (int)vdwEij.size(); ++i) {
            printf(" %g", vdwEij[i]);
        }
        printf("\n");
    }

    std::ostream& generateSelfCode(std::ostream& os) const {
        os << "defineAtomClass((DuMM::AtomClassIndex)";
        os << atomClassIx << ", ";
        os << "\"" << name << "\", ";
        os << element << ", ";
        os << valence << ", ";
        os << vdwRadius << ", ";
        os << vdwWellDepth << ");";

        return os;
    }

    // TOPOLOGICAL STATE VARIABLES
    //   Filled in during construction.

    DuMM::AtomClassIndex atomClassIx;
    std::string name;

    int element;
    int valence;       // # of direct bonds expected
    Real vdwRadius;    // ri, nm
    Real vdwWellDepth; // ei, kJ=Da-nm^2/ps^2


    // TOPOLOGICAL CACHE ENTRIES
    //   These are calculated in realizeTopology() from topological
    //   state variables (from here or others in the DuMM class).

    // After all types have been defined, we can calculate vdw
    // combining rules for dmin and well depth energy. We only fill
    // in entries for pairings of this class with itself and with
    // higher-numbered atom types, so to find the entry for class c,
    // index these arrays by c-atomClassIx where atomClassIx is the
    // class Index of the present AtomClass.
    // Note that different combining rules may be used but they
    // will always result in a pair of vdw parameters.
    Array_<Real> vdwDij; // nm
    Array_<Real> vdwEij; // kJ=Da-A^2/ps^2
};


//-----------------------------------------------------------------------------
//                         CHARGED ATOM TYPE
//-----------------------------------------------------------------------------
class ChargedAtomType {
    public:
    ChargedAtomType()
        : chargedAtomTypeIndex(DuMM::InvalidChargedAtomTypeIndex)
        , atomClassIx(DuMM::InvalidAtomIndex)
        , partialCharge(NaN) {
    }
    ChargedAtomType(DuMM::ChargedAtomTypeIndex id, const char* nm, DuMM::AtomClassIndex aclass, Real chg)
        : chargedAtomTypeIndex(id)
        , name(nm)
        , atomClassIx(aclass)
        , partialCharge(chg) {
        assert(isValid());
    }
    bool isValid() const {
        return chargedAtomTypeIndex.isValid() && atomClassIx.isValid();
    }

    void dump() const {
        printf("    %d(%s): atomClassIx=%d, chg=%g e\n",
               (int)chargedAtomTypeIndex,
               name.c_str(),
               (int)atomClassIx,
               partialCharge);
    }

    std::ostream& generateSelfCode(std::ostream& os) const {
        os << "defineChargedAtomType((DuMM::ChargedAtomTypeIndex)";
        os << chargedAtomTypeIndex << ", ";
        os << "\"" << name << "\", ";
        os << "(DuMM::AtomClassIndex)" << atomClassIx << ", ";
        os << partialCharge << ");";

        return os;
    }

    // These are all Topological state variables, filled in during construction.
    // There are no calculations to be performed.
    DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    std::string name;

    DuMM::AtomClassIndex atomClassIx;
    Real partialCharge; // qi, in e (charge on proton)
};


//-----------------------------------------------------------------------------
//                               BOND STRETCH
//-----------------------------------------------------------------------------
// This represents bond-stretch information for a pair of atom types.
// Use an AtomIndexPair as a key.
// There is always an internal bond stretch term, and there may also be one
// or more custom bond stretch terms as well. If you only have custom terms,
// set the stiffness k to 0.
class BondStretch {
    public:
    BondStretch() = default;
    explicit BondStretch(AtomClassIndexPair key)
        : classes(key)
        , k(-1)
        , d0(-1) {
        assert(classes.isValid());
    }
    BondStretch(AtomClassIndexPair key, Real stiffnessInKJperNmSq, Real lengthInNm)
        : classes(key)
        , k(stiffnessInKJperNmSq)
        , d0(lengthInNm) {
        assert(classes.isValid() && k >= 0 && d0 >= 0);
    }

    ~BondStretch() = default;

    bool hasBuiltinTerm() const {
        return k >= 0;
    }
    bool hasCustomTerm() const {
        return !customTerms.empty();
    }

    void setBuiltinTerm(Real stiffnessInKJperNmSq, Real lengthInNm) {
        assert(!hasBuiltinTerm());
        k = stiffnessInKJperNmSq;
        d0 = lengthInNm;
        assert(classes.isValid() && k >= 0 && d0 >= 0);
    }

    int addCustomTerm(std::unique_ptr<DuMM::CustomBondStretch> custom) {
        customTerms.push_back(std::move(custom));
        return static_cast<int>(customTerms.size() - 1);
    }

    bool isValid() const {
        return (classes.isValid() && hasBuiltinTerm()) || hasCustomTerm();
    } // TODO what is the correct order here?

    std::ostream& generateSelfCode(std::ostream& os) const {
        if (hasBuiltinTerm()) {
            os << "defineBondStretch((DuMM::AtomClassIndex)";
            os << (int)classes[0] << ", (DuMM::AtomClassIndex)";
            os << (int)classes[1] << ", ";
            os << k << ", ";
            os << d0 << ");";
        }
        if (hasCustomTerm()) {
            if (hasBuiltinTerm()) {
                os << std::endl;
            }
            os << "// atom class pair " << classes << " has " << customTerms.size() << " custom term(s).";
        }

        return os;
    }

    AtomClassIndexPair classes;
    Real k;  // in energy units (kJ=Da-nm^2/ps^2) per nm^2, i.e. Da/ps^2
    Real d0; // distance at which force is 0 (in nm)

    // We own the heap space here. Don't forget to clean up!
    // Array_< DuMM::CustomBondStretch* > customTerms;
    std::vector<std::unique_ptr<DuMM::CustomBondStretch>> customTerms;
};


//-----------------------------------------------------------------------------
//                                BOND BEND
//-----------------------------------------------------------------------------
class BondBend {
    public:
    // no default constructor
    explicit BondBend(AtomClassIndexTriple key)
        : classes(key)
        , k(-1)
        , theta0(-1) {
        assert(classes.isValid());
    }
    BondBend(AtomClassIndexTriple key, Real stiffnessInKJPerRadSq, Real angleInDeg)
        : classes(key)
        , k(stiffnessInKJPerRadSq)
        , theta0(angleInDeg * DuMM::Deg2Rad) {
        assert(classes.isValid() && k >= 0 && (0 <= angleInDeg && angleInDeg <= 180));
    }

    // Clean up CustomBondBend terms.
    ~BondBend() {
        for (int i = 0; i < (int)customTerms.size(); ++i) {
            delete customTerms[i];
        }
    }

    bool hasBuiltinTerm() const {
        return k >= 0;
    }
    bool hasCustomTerm() const {
        return !customTerms.empty();
    }

    void setBuiltinTerm(Real stiffnessInKJPerRadSq, Real angleInDeg) {
        assert(!hasBuiltinTerm());
        k = stiffnessInKJPerRadSq;
        theta0 = angleInDeg * DuMM::Deg2Rad;
        assert(classes.isValid() && k >= 0 && (0 <= angleInDeg && angleInDeg <= 180));
    }

    int addCustomTerm(DuMM::CustomBondBend* custom) {
        const int index = (int)customTerms.size();
        customTerms.push_back(custom);
        return index;
    }

    bool isValid() const {
        return (classes.isValid() && hasBuiltinTerm()) || hasCustomTerm();
    } // TODO what is the correct order here?

    // Given a central atom location c bonded to atoms at r and s,
    // calculate the angle between them, the potential energy,
    // and forces on each of the three atoms.
    // TODO: no need for all this if there is a mobility that
    // corresponds directly to this bond bend angle.
    void calculateAtomForces(const Vec3& cG,
                             const Vec3& rG,
                             const Vec3& sG,
                             const Real& scale,
                             const Real& customScale,
                             Real& theta,
                             Real& pe,
                             Vec3& cf,
                             Vec3& rf,
                             Vec3& sf) const;

    std::ostream& generateSelfCode(std::ostream& os) const {
        if (hasBuiltinTerm()) {
            os << "defineBondBend((DuMM::AtomClassIndex)";
            os << (int)classes[0] << ", DuMM::AtomClassIndex(";
            os << (int)classes[1] << "), DuMM::AtomClassIndex(";
            os << (int)classes[2] << "), ";
            os << k << ", ";
            os << (theta0 * DuMM::Rad2Deg) << ");";
        }
        if (hasCustomTerm()) {
            if (hasBuiltinTerm()) {
                os << std::endl;
            }
            os << "// atom class triple " << classes << " has " << customTerms.size() << " custom term(s).";
        }

        return os;
    }

    AtomClassIndexTriple classes;
    Real k;      // energy units kJ per rad^2, i.e. Da-nm^2/(ps^2-rad^2)
    Real theta0; // unstressed angle in radians

    // We own the heap space here. Don't forget to clean up!
    Array_<DuMM::CustomBondBend*> customTerms;
};


/*
//-----------------------------------------------------------------------------
//                              TORSION TERM
//-----------------------------------------------------------------------------
// Torsion term for atoms bonded r-x-y-s. Rotation occurs about
// the axis v=y-x, that is, a vector from x to y. We define a torsion
// angle theta using the "polymer convention" rather than the IUPAC
// one which is 180 degrees different. Ours is like this:
//             r                         r      s
//   theta=0    \             theta=180   \    /
//               x--y                      x--y
//                   \
//                    s
// The sign convention is the same for IUPAC and polymer:
// A positive angle is defined by considering r-x fixed in space. Then
// using the right hand rule around v (that is, thumb points from x to y)
// a positive rotation rotates y->s in the direction of your fingers.
//
// We use a periodic energy function like this:
//       E(theta) = sum E_n(1 + cos(n*theta - theta0_n))
// where n is the periodicity, E_n is the amplitude (kJ/mol) for
// term n, and theta0_n is the phase offset for term n. The torque
// term (applied about the v axis) is then
//       T(theta) = -[sum -n*E_n*sin(n*theta - theta0_n)]
// We have to translate this into forces on the four atoms.
//
*/
class TorsionTerm {
    public:
    TorsionTerm()
        : periodicity(-1)
        , amplitude(-1)
        , theta0(-1) {
    }
    TorsionTerm(int n, Real ampInKJ, Real th0InDeg)
        : periodicity(n)
        , amplitude(ampInKJ)
        , theta0(th0InDeg * DuMM::Deg2Rad) {
        // theta0 might be numerically a hair over Pi
        while (theta0 > Pi) {
            theta0 -= 2 * Pi;
        }
        while (theta0 <= -Pi) {
            theta0 += 2 * Pi;
        }


        /*        std::cout << "periodicity = " << periodicity << std::endl;
                std::cout << "amplitude = " << amplitude << std::endl;
                std::cout << "theta0 = " << theta0 << std::endl;
                std::cout << "Pi = " << Pi << std::endl;*/


        assert(isValid());
    }
    bool isValid() const {
        return periodicity > 0
               //&& amplitude >= 0 // GMOL Amber allows negative dihedral energy
               && -Pi < theta0 && theta0 <= Pi;
    }
    Real energy(Real theta) const {
        return amplitude * (1 + std::cos(periodicity * theta - theta0));
    }
    Real torque(Real theta) const {
        return periodicity * amplitude * std::sin(periodicity * theta - theta0);
    }

    std::ostream& generateSelfCode(std::ostream& os) const {
        os << ", " << periodicity;
        os << ", " << amplitude;
        os << ", " << theta0 * DuMM::Rad2Deg;

        return os;
    }

    int periodicity; // 1=360, 2=180, 3=120, etc.
    Real amplitude;  // energy units (kJ)
    Real theta0;     // radians
};


//-----------------------------------------------------------------------------
//                              BOND TORSION
//-----------------------------------------------------------------------------
class BondTorsion {
    public:
    // no default constructor
    BondTorsion(AtomClassIndexQuad key)
        : classes(key) {
        assert(classes.isValid());
    }

    // Clean up CustomBondTorsion terms.
    ~BondTorsion() {
        for (int i = 0; i < (int)customTerms.size(); ++i) {
            delete customTerms[i];
        }
    }

    bool hasBuiltinTerm() const {
        return !terms.empty();
    }
    bool hasCustomTerm() const {
        return !customTerms.empty();
    }

    void addBuiltinTerm(const TorsionTerm& tt) {
        assert(!hasTermWithPeriod(tt.periodicity));
        terms.push_back(tt);
    }

    int addCustomTerm(DuMM::CustomBondTorsion* custom) {
        const int index = (int)customTerms.size();
        customTerms.push_back(custom);
        return index;
    }

    bool isValid() const {
        return !(terms.empty() && customTerms.empty());
    }

    bool hasTermWithPeriod(int n) const {
        for (int i = 0; i < (int)terms.size(); ++i) {
            if (terms[i].periodicity == n) {
                return true;
            }
        }
        return false;
    }

    const TorsionTerm& getTermWithPeriod(int n) const {
        static const TorsionTerm dummy; // invalid
        for (int i = 0; i < (int)terms.size(); ++i) {
            if (terms[i].periodicity == n) {
                return terms[i];
            }
        }
        return dummy;
    }

    // equality operator to help handle case where user innocently
    // attempts to add the same torsion a second time
    // WARNING: this is very inefficient
    bool operator==(const BondTorsion& other) const {
        if (terms.size() != other.terms.size()) {
            return false;
        }
        Array_<TorsionTerm>::const_iterator iTerm;
        for (iTerm = terms.begin(); iTerm != terms.end(); ++iTerm) {
            const TorsionTerm& myTerm = *iTerm;
            if (!other.hasTermWithPeriod(myTerm.periodicity)) {
                return false;
            }

            Array_<TorsionTerm>::const_iterator iOtherTerm;
            for (iOtherTerm = other.terms.begin(); iOtherTerm != other.terms.end(); ++iOtherTerm) {
                const TorsionTerm& otherTerm = *iOtherTerm;
                if (otherTerm.periodicity == myTerm.periodicity) {
                    if (myTerm.amplitude != otherTerm.amplitude) {
                        return false;
                    }
                    if (myTerm.theta0 != otherTerm.theta0) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // Given atom locations r-x-y-s in the ground frame, calculate the
    // torsion angle, energy and a force on each atom so that the desired
    // pure torque is produced.
    // TODO: no need for all this if there is a mobility that
    // corresponds directly to this dihedral angle.
    void calculateAtomForces(const Vec3& rG,
                             const Vec3& xG,
                             const Vec3& yG,
                             const Vec3& sG,
                             const Real& builtinScale,
                             const Real& customScale,
                             Real& theta,
                             Real& pe,
                             Vec3& rf,
                             Vec3& xf,
                             Vec3& yf,
                             Vec3& sf) const;

    // Type 1 => normal torsion parameters
    // Type 2 => amber improper torsion parameters
    std::ostream& generateSelfCode(std::ostream& os, int torsionType = 1) const {
        if (hasBuiltinTerm()) {
            if (torsionType == 1) {
                os << "defineBondTorsion";
            } else {
                os << "defineAmberImproperTorsion";
            }
            os << "( (DuMM::AtomClassIndex)" << classes[0];
            os << ", (DuMM::AtomClassIndex)" << classes[1];
            os << ", (DuMM::AtomClassIndex)" << classes[2];
            os << ", (DuMM::AtomClassIndex)" << classes[3];
            Array_<TorsionTerm>::const_iterator term;
            for (term = terms.begin(); term != terms.end(); ++term) {
                term->generateSelfCode(os);
            }
            os << ");";
        }
        if (hasCustomTerm()) {
            if (hasBuiltinTerm()) {
                os << std::endl;
            }
            os << "// atom class quad " << classes << " has " << customTerms.size() << " custom term(s).";
        }

        return os;
    }

    AtomClassIndexQuad classes;
    Array_<TorsionTerm> terms;

    // We own the heap space here. Don't forget to clean up!
    Array_<DuMM::CustomBondTorsion*> customTerms;
};


//-----------------------------------------------------------------------------
//                             ATOM PLACEMENT
//-----------------------------------------------------------------------------
class AtomPlacement {
    public:
    AtomPlacement()
        : atomIndex(-1) {
    }
    AtomPlacement(DuMM::AtomIndex dAIx, const Vec3& station)
        : station(station)
        , atomIndex(dAIx) {
        assert(isValid());
    }
    [[nodiscard]] auto isValid() const -> bool {
        return atomIndex.isValid();
    }

    Vec3 station; // in nm
    DuMM::AtomIndex atomIndex;

    void setStation(const Vec3& newStation) {
        station = newStation;
    }
};

inline auto operator<(const AtomPlacement& ap1, const AtomPlacement& ap2) -> bool {
    return ap1.atomIndex < ap2.atomIndex;
}
inline auto operator==(const AtomPlacement& ap1, const AtomPlacement& ap2) -> bool {
    return ap1.atomIndex == ap2.atomIndex;
}


//-----------------------------------------------------------------------------
//                              CLUSTER PLACEMENT
//-----------------------------------------------------------------------------
class ClusterPlacement {
    public:
    ClusterPlacement()
        : clusterIndex(-1) {
    }
    ClusterPlacement(DuMM::ClusterIndex clusterIx, const Transform& transform)
        : placement(transform)
        , clusterIndex(clusterIx) {
        assert(isValid());
    }
    [[nodiscard]] auto isValid() const -> bool {
        return clusterIndex.isValid();
    }

    Transform placement; // translation in nm
    DuMM::ClusterIndex clusterIndex;
};
inline auto operator<(const ClusterPlacement& rhs, const ClusterPlacement& lhs) -> bool {
    return rhs.clusterIndex < lhs.clusterIndex;
}
inline auto operator==(const ClusterPlacement& rhs, const ClusterPlacement& lhs) -> bool {
    return rhs.clusterIndex == lhs.clusterIndex;
}

using AtomPlacementArray = std::vector<AtomPlacement>;
using ClusterPlacementArray = std::vector<ClusterPlacement>;
using ShortAtomArray = std::vector<DuMM::AtomIndex>;


//-----------------------------------------------------------------------------
//                               DuMM ATOM
//-----------------------------------------------------------------------------
// This class represents a unique atom in the set of modeled molecules, as it
// was defined during construction. We will keep an entry in the atoms array
// for every DuMMAtom but they are not necessarily all involved in force
// calculations; only "included atoms" are involved in force calculations.
// There is a separate IncludedAtom type containing useful precalculations
// to expedite runtime processing.
class DuMMAtom {
    public:
    DuMMAtom() {
    }
    DuMMAtom(DuMM::ChargedAtomTypeIndex t, DuMM::AtomIndex aIx)
        : atomIndex(aIx)
        , chargedAtomTypeIndex(t) {
        assert(isValid());
    }

    bool isValid() const {
        return atomIndex.isValid() && chargedAtomTypeIndex.isValid();
    }

    bool isAttachedToBody() const {
        return mobodIx.isValid();
    }
    MobilizedBodyIndex getMobodIndex() const {
        assert(isAttachedToBody());
        return mobodIx;
    }
    const Vec3& getStationOnBody() const // in body frame
    {
        assert(isAttachedToBody());
        return station_B;
    }

    bool isIncludedAtom() const {
        return inclAtomIndex.isValid();
    }
    DuMM::IncludedAtomIndex getIncludedAtomIndex() const {
        assert(isIncludedAtom());
        return inclAtomIndex;
    }

    bool isNonbondAtom() const {
        return nonbondAtomIndex.isValid();
    }
    DuMM::NonbondAtomIndex getNonbondAtomIndex() const {
        assert(isNonbondAtom());
        return nonbondAtomIndex;
    }

    bool isBondStarterAtom() const {
        return bondStarterIndex.isValid();
    }
    DuMMBondStarterIndex getBondStarterIndex() const {
        assert(isBondStarterAtom());
        return bondStarterIndex;
    }

    // EU COMMENT BEGIN
    void attachToBody(MobilizedBodyIndex mbx, const Vec3& stationInB) {
        // EU BEGIN
        // void attachToBody(MobilizedBodyIndex mbx, Vec3 stationInB) {
        // EU END
        assert(!isAttachedToBody());
        mobodIx = mbx;
        station_B = stationInB;
    }

    bool isBondedTo(DuMM::AtomIndex ax) const {
        for (int i = 0; i < (int)bond12.size(); ++i) {
            if (bond12[i] == ax) {
                return true;
            }
        }
        return false;
    }

    // desk_mass_related
    const SimTK::mdunits::Mass getMass() const {
        return mass;
    }

    // desk_mass_related
    void setMass(SimTK::mdunits::Mass inpMass) {
        this->mass = inpMass;
    }

    void dump() const {
        printf(" atomIndex=%d\n", (int)atomIndex);
        printf(" chargedAtomType=%d mobod=%d station=%g %g %g\n",
               (int)chargedAtomTypeIndex,
               (int)mobodIx,
               station_B[0],
               station_B[1],
               station_B[2]);
        printf(" includedAtomIx=%d nonbondAtomIx=%d bondStarterAtomIx=%d\n",
               (int)inclAtomIndex,
               (int)nonbondAtomIndex,
               (int)bondStarterIndex);
        printf(" includedBodyIx=%d\n", (int)inclBodyIndex);

        printf("          bond 1-2 (AtomIndex):");
        for (int i = 0; i < (int)bond12.size(); ++i) {
            printf(" %d", (int)bond12[i]);
        }
        printf("\n");
    }

    void invalidateTopologicalCache() {
        mobodIx.invalidate();
        station_B = NaN;
        inclAtomIndex.invalidate();
        nonbondAtomIndex.invalidate();
        bondStarterIndex.invalidate();
        mass = NaN;
    }

    public:
    // TOPOLOGICAL STATE VARIABLES
    //   Filled in during construction.

    DuMM::AtomIndex atomIndex;
    DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    ShortAtomArray bond12;

    // TOPOLOGICAL CACHE ENTRIES
    //   These are calculated in realizeTopology() from topological
    //   state variables (from here or others in the DuMM class).

    // After the atom or a containing cluster has been attached to a
    // body, we fill these in. The atoms's station fixed in body bodyIx's
    // frame is given in nm.
    MobilizedBodyIndex mobodIx;
    Vec3 station_B;

    // If this ends up being one of the included atoms, record the assigned
    // IncludedAtomIndex and IncludedBodyIndex here, otherwise they will be
    // invalid.
    DuMM::IncludedAtomIndex inclAtomIndex;
    DuMMIncludedBodyIndex inclBodyIndex;

    // If this is an included atom that is involved in nonbond interactions,
    // record its NonbondAtomIndex here, otherwise it will be invalid.
    DuMM::NonbondAtomIndex nonbondAtomIndex;

    // If this is an included atom that serves as atom 1 for some bond, then
    // record its BondStarterIndex here.
    DuMMBondStarterIndex bondStarterIndex;


    // GMolModel - all topology lists and indexes
    Vec3 station_B_All;
    DuMM::IncludedAtomIndex AllAtomIndex;
    DuMMIncludedBodyIndex AllBodyIndex;
    DuMM::NonbondAtomIndex AllnonbondAtomIndex;
    DuMMBondStarterIndex AllbondStarterIndex;


    mdunits::Mass mass;
};


//-----------------------------------------------------------------------------
//                             INCLUDED ATOM
//-----------------------------------------------------------------------------
// This class represents an atom for which it has been determined that it
// will participate in runtime force calculations. This includes precalculated
// information necessary for high-speed processing of bonded and nonbond
// calculations. These entries should be thought of as topological cache
// information; they are filled in during realizeTopology() using the
// topological state information to be found in atom, cluster, body, and
// force field information provided during construction.
class IncludedAtom {
    public:
    IncludedAtom() {
    }
    IncludedAtom(DuMM::IncludedAtomIndex inclAtomIx,
                 DuMMIncludedBodyIndex inclBodyIx,
                 DuMM::AtomIndex atomIx,
                 DuMM::ChargedAtomTypeIndex chargeTypeIx)
        : inclAtomIndex(inclAtomIx)
        , inclBodyIndex(inclBodyIx)
        , atomIndex(atomIx)
        , chargedAtomTypeIndex(chargeTypeIx) {
    }

    void setIndices(DuMM::IncludedAtomIndex inclAtomIx,
                    DuMMIncludedBodyIndex inclBodyIx,
                    DuMM::AtomIndex atomIx,
                    DuMM::ChargedAtomTypeIndex chargeTypeIx) {
        inclAtomIndex = inclAtomIx;
        inclBodyIndex = inclBodyIx;
        atomIndex = atomIx;
        chargedAtomTypeIndex = chargeTypeIx;
    }

    bool isValid() const {
        return inclAtomIndex.isValid() && inclBodyIndex.isValid() && atomIndex.isValid()
               && chargedAtomTypeIndex.isValid();
    }

    void clear() {
        inclAtomIndex.invalidate();
        inclBodyIndex.invalidate();
        atomIndex.invalidate();
        chargedAtomTypeIndex.invalidate();

        scale12.clear();
        scale13.clear();
        scale14.clear();
        scale15.clear();
        force12.clear();
        force13.clear();
        force14.clear();
        force15.clear();
        forceImproper14.clear();

        stretch.clear();
        bend.clear();
        torsion.clear();
        aImproperTorsion.clear();
    }

    void dump() const;

    // This is the included atom index of this atom -- redundant information
    // probably, but useful for debugging.
    DuMM::IncludedAtomIndex inclAtomIndex;

    // This is the included body to which this atom is fixed.
    DuMMIncludedBodyIndex inclBodyIndex;

    // This is the original index of the user-defined atom corresponding to
    // this included atom. This should not be referenced at run time except
    // possibly for debugging and error messages. This is the only index
    // likely to be recognizable to the DuMM user.
    DuMM::AtomIndex atomIndex;

    // This gets us to properties of this atom that are needed for nonbonded
    // calculations, such as charge and pre-mixed van der Waals properties.
    DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;

    // If this is a nonbond atom, then these are precalculated lists of nearby
    // nonbond atoms for which the force field may require us to scale Coulomb
    // and/or van der Waals interactions. These lists do not include
    // atoms that are on the same body with this one since we don't
    // need to calculate nonbond interactions with those atoms so it would
    // be a waste of time to scale them.
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale12;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale13;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale14;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale15;

    // And for GMolModel - all topology lists
    DuMM::IncludedAtomIndex AllAtomIndex;
    DuMMIncludedBodyIndex AllBodyIndex;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale12All;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale13All;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale14All;
    Array_<DuMM::NonbondAtomIndex, unsigned short> scale15All;


    // If this is a bondStarter atom (that is, it is atom 1 in some bonded
    // force term), then these are precalculated lists of the atoms that
    // participate in each of the bonded force terms that this atom is
    // required to calculate. Duplicates (e.g. bond bend a-b-c as seen both
    // from atom "a" and atom "c") will already have been eliminated;
    // every one of these force terms for every bond starter atom must be
    // evaluated.
    Array_<DuMM::IncludedAtomIndex, unsigned short> force12;
    Array_<IncludedAtomIndexPair, unsigned short> force13;
    Array_<IncludedAtomIndexTriple, unsigned short> force14;
    Array_<IncludedAtomIndexQuad, unsigned short> force15; // not used
    Array_<IncludedAtomIndexTriple, unsigned short> forceImproper14;

    // And for GMolModel - - all topology lists
    Array_<DuMM::IncludedAtomIndex, unsigned short> force12All;
    Array_<IncludedAtomIndexPair, unsigned short> force13All;
    Array_<IncludedAtomIndexTriple, unsigned short> force14All;
    Array_<IncludedAtomIndexQuad, unsigned short> force15All; // not used
    Array_<IncludedAtomIndexTriple, unsigned short> forceImproper14All;


    // These are pointers into the various bonded maps that provide instant
    // access to the coefficients for particular bonds. These correspond
    // elementwise to the atom sets in the force1X arrays above, so have
    // the same length as the indicated list.
    Array_<const BondStretch*, unsigned short> stretch; // matches force12
    Array_<const BondBend*, unsigned short> bend;       // matches force13
    Array_<const BondTorsion*, unsigned short> torsion; // matches force14
    // Currently there are no supported 1-5 bonded force terms.
    Array_<const BondTorsion*, unsigned short> aImproperTorsion;
    // matches forceImproper14

    // And for GMolModel - all topology lists
    Array_<const BondStretch*, unsigned short> stretchAll; // matches force12
    Array_<const BondBend*, unsigned short> bendAll;       // matches force13
    Array_<const BondTorsion*, unsigned short> torsionAll; // matches force14
    Array_<const BondTorsion*, unsigned short> aImproperTorsionAll;
};


//-----------------------------------------------------------------------------
//                                  BOND
//-----------------------------------------------------------------------------
class Bond {
    public:
    Bond() {
    }
    Bond(DuMM::AtomIndex atom1, DuMM::AtomIndex atom2)
        : atoms(atom1, atom2) {
        assert(isValid());
    }
    bool isValid() const {
        return atoms.isValid();
    }

    AtomIndexPair atoms;
};

// NOT USED YET
class ChargeProperties {
    public:
    Real netCharge;            // in proton charge units e
    Vec3 centerOfCharge;       // in nm
    Vec3 dipoleMoment;         // units?? TODO
    SymMat33 quadrupoleMoment; // units?? TODO
};

// NOT USED YET
class GeometricProperties {
    public:
    Transform obbFrame;
    Vec3 obbHalfLengths;       // nm
    Real boundingSphereRadius; // nm
    Vec3 boundingSphereCenter; // nm
};


//-----------------------------------------------------------------------------
//                                 CLUSTER
//-----------------------------------------------------------------------------
// This class is a rigid grouping of atoms. It can contain atoms directly, and
// subclusters which can contain atoms or sub-subclusters, etc. As we build
// up a cluster, we keep a running "flat" view of all the atoms and all the
// clusters contained anywhere deep within, already transformed to this
// cluster's reference frame.
//
class Cluster {
    public:
    Cluster() = default;
    Cluster(const char* name)
        : name(name) {
        // not valid yet -- still need index assigned
    }

    [[nodiscard]] auto isValid() const -> bool {
        return clusterIndex.isValid();
    }
    [[nodiscard]] auto isAttachedToBody() const -> bool {
        return mobodIx.isValid();
    }
    [[nodiscard]] auto isTopLevelCluster() const -> bool {
        return parentClusters.empty();
    }
    [[nodiscard]] auto getMobodIndex() const -> MobilizedBodyIndex {
        assert(isAttachedToBody());
        return mobodIx;
    }

    [[nodiscard]] auto getDirectlyContainedAtoms() const -> const AtomPlacementArray& {
        return directAtomPlacements;
    }
    [[nodiscard]] auto updDirectlyContainedAtoms() -> AtomPlacementArray& {
        return directAtomPlacements;
    }

    [[nodiscard]] auto getAllContainedAtoms() const -> const AtomPlacementArray& {
        return allAtomPlacements;
    }
    [[nodiscard]] auto updAllContainedAtoms() -> AtomPlacementArray& {
        return allAtomPlacements;
    }

    [[nodiscard]] auto getDirectlyContainedClusters() const -> const ClusterPlacementArray& {
        return directClusterPlacements;
    }
    [[nodiscard]] auto updDirectlyContainedClusters() -> ClusterPlacementArray& {
        return directClusterPlacements;
    }

    [[nodiscard]] auto getAllContainedClusters() const -> const ClusterPlacementArray& {
        return allClusterPlacements;
    }
    [[nodiscard]] auto updAllContainedClusters() -> ClusterPlacementArray& {
        return allClusterPlacements;
    }

    [[nodiscard]] auto containsAtom(DuMM::AtomIndex atomIndex) const -> bool {
        return std::find(allAtomPlacements.begin(),
                         allAtomPlacements.end(),
                         AtomPlacement(atomIndex, Vec3(0)))
               != allAtomPlacements.end();
    }

    [[nodiscard]] auto containsCluster(DuMM::ClusterIndex clusterIndex) const -> bool {
        return std::find(allClusterPlacements.begin(),
                         allClusterPlacements.end(),
                         ClusterPlacement(clusterIndex, Transform()))
               != allClusterPlacements.end();
    }

    // See if a cluster contains any atoms which are already in
    // any of the cluster trees to which this cluster is associated.
    // TODO: can only handle top-level cluster so we don't have to run up the
    //       ancestor branches.
    // If we find an atom common to both clusters we'll return it to permit
    // nice error messages, otherwise we return false and -1 for the atomIndex.
    auto overlapsWithCluster(const Cluster& test, DuMM::AtomIndex& anAtomIndexInBothClusters) const -> bool {
        assert(isTopLevelCluster());

        const AtomPlacementArray& testAtoms = test.getAllContainedAtoms();
        // const AtomPlacementArray& myAtoms   = getAllContainedAtoms();

        auto atomPlacement = testAtoms.begin();
        while (atomPlacement != testAtoms.end()) {
            if (containsAtom(atomPlacement->atomIndex)) {
                anAtomIndexInBothClusters = atomPlacement->atomIndex;
                return true;
            }
            ++atomPlacement;
        }
        anAtomIndexInBothClusters = DuMM::InvalidAtomIndex;
        return false;
    }

    // Return true if this cluster contains (directly or indirectly) any atom which has already
    // been attached to a body. If so return one of the attached atoms and its body, which can
    // be helpful in error messages.
    auto containsAnyAtomsAttachedToABody(DuMM::AtomIndex& atomIndex,
                                         MobilizedBodyIndex& bodyIx,
                                         const DuMMForceFieldSubsystemRep& mm) const -> bool;

    // Translation is in nm.
    void attachToBody(MobilizedBodyIndex bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm);

    // Place an atom in this cluster. To be valid, the atom must not
    // already be
    //   (a) in any of the trees of which this group is a part, or
    //   (b) attached to a body.
    // TODO: (c) at the moment we don't allow placing an atom in a group unless
    //           that group is a top-level group (i.e., it has no parents).
    // If this group is already attached to a body, then we will update
    // the atom entry to note that it is now attached to the body also.

    // EU COMMENT BEGIN
    void placeAtom(DuMM::AtomIndex atomIndex, const Vec3& stationInNm, DuMMForceFieldSubsystemRep& mm);
    // EU BEGIN
    // void placeAtom(DuMM::AtomIndex atomIndex, Vec3 stationInNm,
    //               DuMMForceFieldSubsystemRep& mm);
    // EU END

    // Place a child cluster in this parent cluster. To be valid, the child
    // must not
    //   (a) already be contained in the parent group or one of the parent's
    //       subgroups, or
    //   (b) contain any atoms which are already present in the parent or any
    //       of the parent's subgroups, or
    //   (c) already be attached to a body.
    // TODO: (d) at the moment we don't allow adding a child group unless
    //           the parent (this) group is a top-level group (i.e., it has no
    //           parents).
    // If the parent is already attached to a body, then we will update
    // the child to note that it is now attached to the body also (and it
    // will update its contained atoms).
    // (translation is in nm)
    void placeCluster(DuMM::ClusterIndex childClusterIndex,
                      const Transform& placement,
                      DuMMForceFieldSubsystemRep& mm);

    // Calculate the composite mass properties for this cluster, transformed
    // into the indicated frame. Translation part of the Transform is in nm,
    // returned mass proprties are in daltons and nm.
    [[nodiscard]] auto calcMassProperties(const Transform& tr, const DuMMForceFieldSubsystemRep& mm) const
        -> MassProperties;


    // Recursively calculate composite properties for this group and all the
    // groups it contains. All groups were marked "invalid" at the beginning
    // of this step.
    void invalidateTopologicalCache() { // TODO
    }
    void realizeTopologicalCache(DuMMForceFieldSubsystemRep& mm) {
    }

    static void dumpX(const Transform& X) {
        if (X == Transform()) {
            std::cout << "Identity\n";
        } else {
            std::cout << X;
        }
    }

    void dump() const {
        printf("    clusterIndex=%d(%s)\n", (int)clusterIndex, name.c_str());
        printf("      direct atom placements (nm): ");
        auto ap = directAtomPlacements.begin();
        while (ap != directAtomPlacements.end()) {
            std::cout << " " << ap->atomIndex << ":" << ap->station;
            ++ap;
        }
        printf("\n      all atom placements (nm): ");
        auto aap = allAtomPlacements.begin();
        while (aap != allAtomPlacements.end()) {
            std::cout << " " << aap->atomIndex << ":" << aap->station;
            ++aap;
        }
        printf("\n      direct cluster placements (nm):\n");
        auto cp = directClusterPlacements.begin();
        while (cp != directClusterPlacements.end()) {
            std::cout << "      " << cp->clusterIndex << ":";
            dumpX(cp->placement);
            ++cp;
        }
        printf("\n      all cluster placements (nm):\n");
        auto acp = allClusterPlacements.begin();
        while (acp != allClusterPlacements.end()) {
            std::cout << "      " << acp->clusterIndex << ":";
            dumpX(acp->placement);
            ++acp;
        }
        printf("\n      parent cluster placements (nm):\n");
        auto pp = parentClusters.begin();
        while (pp != parentClusters.end()) {
            std::cout << "      " << pp->clusterIndex << ":";
            dumpX(pp->placement);
            ++pp;
        }

        if (mobodIx.isValid()) {
            std::cout << "\n      attached to body " << mobodIx << " at (nm) ";
            dumpX(placement_B);
        } else {
            std::cout << "\n      NOT ATTACHED TO ANY BODY.";
        }
        std::cout << "\n";
    }

    void clearAllCalculatedData() {
        chargeProps = ChargeProperties();
        geometricProps = GeometricProperties();
    }

    private:
    // translation is in nm
    void noteNewChildCluster(DuMM::ClusterIndex childClusterIndex, const Transform& X_PC) {
        directClusterPlacements.emplace_back(childClusterIndex, X_PC);
        allClusterPlacements.emplace_back(childClusterIndex, X_PC);
    }

    // translation is in nm
    void noteNewParentCluster(DuMM::ClusterIndex parentClusterIndex, const Transform& X_PC) {
        parentClusters.emplace_back(parentClusterIndex, X_PC);
    }

    public:
    // TOPOLOGICAL STATE VARIABLES
    //   Filled in during construction.
    DuMM::ClusterIndex clusterIndex;
    std::string name;

    // These are the *directly* attached atoms and clusters.
    AtomPlacementArray directAtomPlacements;
    ClusterPlacementArray directClusterPlacements;

    // These sets are kept up to date as we add atoms and clusters.
    // 'allAtomPlacements' contains *all* the atoms in this cluster
    // or its descendents, transformed into this cluster's frame.
    // 'allClusterPlacements' contains *all* the clusters in this
    // cluster or its subclusters, transformed into this cluster's frame.
    AtomPlacementArray allAtomPlacements;
    ClusterPlacementArray allClusterPlacements;

    // This is a list of all the immediate parents of this cluster, if any.
    // This is updated whenever this cluster is placed in another one. The
    // body is *not* considered a parent cluster; it is handled separately
    // below. Note that whenever an atom or cluster is added to this cluster,
    // the atom or atoms involved [SHOULD BE: TODO] added to each ancestor.
    ClusterPlacementArray parentClusters;

    // After this cluster or a containing cluster has been attached to a
    // body, we can fill these in.
    MobilizedBodyIndex mobodIx;
    Transform placement_B; // cluster's placement fixed in body bodyIx's frame (nm)

    // TOPOLOGICAL CACHE ENTRIES
    //   These are calculated in realizeTopology() from topological
    //   state variables (from here or others in the DuMM class).

    // These reflect composite properties built from the allAtoms list.
    ChargeProperties chargeProps;
    GeometricProperties geometricProps;
};


//-----------------------------------------------------------------------------
//                                DuMM BODY
//-----------------------------------------------------------------------------
// A DuMMBody has a reference to a top-level Cluster.
class DuMMBody {
    public:
    DuMMBody()
        : clusterIndex(DuMM::InvalidClusterIndex)
        , mobilizedBodyIndex(InvalidMobilizedBodyIndex) {
    }

    explicit DuMMBody(DuMM::ClusterIndex cIx, MobilizedBodyIndex mIx)
        : clusterIndex(cIx)
        , mobilizedBodyIndex(mIx) {
        assert(isValid());
    }

    bool isValid() const {
        return (clusterIndex.isValid()) && (mobilizedBodyIndex.isValid());
    }

    void invalidateTopologicalCache() {
    }
    void realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm);

    DuMM::ClusterIndex getClusterIndex() const {
        assert(isValid());
        return clusterIndex;
    }
    MobilizedBodyIndex getMobilizedBodyIndex() const {
        return mobilizedBodyIndex;
    }

    void dump() const {
        printf("    clusterIndex=%d\n", (int)clusterIndex);
        printf("    mobodIndex=%d\n", (int)mobilizedBodyIndex);
    }

    static std::string createClusterNameForBody(int bnum) {
        char buf[100];
        std::sprintf(buf, "DuMMBody %d", bnum);
        return std::string(buf);
    }

    DuMM::ClusterIndex clusterIndex;
    SimTK::MobilizedBodyIndex mobilizedBodyIndex;
};


//-----------------------------------------------------------------------------
//                             INCLUDED BODY
//-----------------------------------------------------------------------------
// An array of these objects is kept for fast runtime access. These are just
// the bodies that have at least one included atom attached. This should be
// considered part of the topology-stage cache -- all the information is
// derived during realizeTopology() from user-provided data stored elsewhere.
class IncludedBody {
    public:
    IncludedBody() {
    }
    bool isValid() const {
        return mobodIx.isValid();
    }
    MobilizedBodyIndex mobodIx;

    // Defined as in std:: classes; end is one past the last atom.
    DuMM::IncludedAtomIndex beginIncludedAtoms, endIncludedAtoms;
    DuMM::NonbondAtomIndex beginNonbondAtoms, endNonbondAtoms;
    DuMMBondStarterIndex beginBondStarterAtoms, endBondStarterAtoms;

    // GMolModel - all topology indexes
    DuMM::IncludedAtomIndex beginAllAtoms, endAllAtoms;
    DuMM::NonbondAtomIndex beginAllNonbondAtoms, endAllNonbondAtoms;
    DuMMBondStarterIndex beginAllBondStarterAtoms, endAllBondStarterAtoms;

    void dump() const {
        printf("    mobodIndex=%d\n", (int)mobodIx);
        printf("    includedAtoms=[%d,%d)\n", (int)beginIncludedAtoms, (int)endIncludedAtoms);
        printf("    nonbondAtoms=[%d,%d)\n", (int)beginNonbondAtoms, (int)endNonbondAtoms);
        printf("    bondStarterAtoms=[%d,%d)\n", (int)beginBondStarterAtoms, (int)endBondStarterAtoms);
    }
};


//-----------------------------------------------------------------------------
//                          INCLUSION LIST SPEC
//-----------------------------------------------------------------------------
// An object of this class is used to collect up all the user's instructions
// about which atoms and bonds to include in force calculations. We will use
// this information during realizeTopology() to build lists of atoms for
// processing at run time. The set of all atoms that are involved in either
// a nonbond calculation or a bond calculation is called the includedAtoms set.
// The subset of includedAtoms that is involved in nonbond calculations
// (coulomb, van der Waals, and/or gbsa) is called the nonbondAtoms set.
// The subset of includedAtoms that serves as atom 1 of any bond is called
// the bondStartAtoms set. There are also includedAtoms that are involved in
// bonded force calculations but don't serve as atom 1 for any bond. Note
// that the subsets of includedAtoms are not necessarily disjoint.
class InclusionListSpec {
    public:
    InclusionListSpec()
        : useDefaultNonbondList(true)
        , useDefaultBondList(true) {
    }

    bool isNonbondAtom(DuMM::AtomIndex ax) const {
        assert(ax.isValid());
        return includedNonbondAtoms.find(ax) != includedNonbondAtoms.end();
    }
    bool isNonbondBody(MobodIndex mbx) const {
        assert(mbx.isValid());
        return includedNonbondBodies.find(mbx) != includedNonbondBodies.end();
    }

    bool isBondAtom(DuMM::AtomIndex ax) const {
        assert(ax.isValid());
        return atomsWhoseBondsAreIncluded.find(ax) != atomsWhoseBondsAreIncluded.end();
    }
    bool isBondBody(MobodIndex mbx) const {
        assert(mbx.isValid());
        return bodiesWhoseBondsAreIncluded.find(mbx) != bodiesWhoseBondsAreIncluded.end();
    }
    bool isBondAtomPair(DuMM::AtomIndex a1x, DuMM::AtomIndex a2x) const {
        assert(a1x.isValid() && a2x.isValid() && a1x != a2x);
        return atomPairsWhoseConnectingBondsAreIncluded.find(
                   AtomIndexPair(a1x, a2x, true)) // canonicalize order
               != atomPairsWhoseConnectingBondsAreIncluded.end();
    }
    bool isBondBodyPair(MobodIndex b1x, MobodIndex b2x) const {
        assert(b1x.isValid() && b2x.isValid());
        if (b1x == b2x) {
            return false;
        }
        return bodyPairsWhoseConnectingBondsAreIncluded.find(
                   MobodIndexPair(b1x, b2x, true)) // canonicalize order
               != bodyPairsWhoseConnectingBondsAreIncluded.end();
    }

    // Nonbond
    bool useDefaultNonbondList; // use all atoms; ignore other nonbond fields
    // These atoms are explicitly included.
    std::set<DuMM::AtomIndex> includedNonbondAtoms;
    // All atoms on this body are explicitly included.
    std::set<MobilizedBodyIndex> includedNonbondBodies;

    // Bonded
    bool useDefaultBondList; // all cross-body bonds; ignore other fields
    // Any cross-body bond that involves one of these atoms is included.
    std::set<DuMM::AtomIndex> atomsWhoseBondsAreIncluded;
    // Any cross-body bond that includes *both* of these atoms is included,
    // even if this is not a 1-2 connection. Pairs are always in
    // (low,high) order.
    std::set<AtomIndexPair> atomPairsWhoseConnectingBondsAreIncluded;
    // Any cross-body bond that includes an atom from any of these bodies
    // is included.
    std::set<MobilizedBodyIndex> bodiesWhoseBondsAreIncluded;
    // Any cross-body bond that contains an atom from both of these bodies
    // is included. Pairs are always in (low,high) order.
    std::set<MobodIndexPair> bodyPairsWhoseConnectingBondsAreIncluded;

    // Same for GMolModel - to be checked
    // These atoms are explicitly All.
    std::set<DuMM::AtomIndex> AllNonbondAtoms;
    // All atoms on this body are explicitly included.
    std::set<MobilizedBodyIndex> AllNonbondBodies;
};

// Unfortunately required by Value<T>.
static inline std::ostream& operator<<(std::ostream& o, const InclusionListSpec& ils) {
    o << "InclusionListSpec\n";
    return o;
}

//-----------------------------------------------------------------------------
//                       DuMM FORCE FIELD SUBSYSTEM REP
//-----------------------------------------------------------------------------
class SimTK::DuMMForceFieldSubsystemRep : public ForceSubsystem::Guts {
    friend class DuMMForceFieldSubsystem;
    static const char* ApiClassName; // "DuMMForceFieldSubsystem"
    public:
    DuMMForceFieldSubsystemRep()
        : ForceSubsystem::Guts("DuMMForceFieldSubsystem", "2.2.0")
        , forceEvaluationCount(0)
        , nextUnusedAtomClassIndex(0)
        , nextUnusedChargedAtomTypeIndex(0) {
        vdwGlobalScaleFactor = coulombGlobalScaleFactor = 1;
        bondStretchGlobalScaleFactor = bondBendGlobalScaleFactor = bondTorsionGlobalScaleFactor =
            amberImproperTorsionGlobalScaleFactor = 1;
        customBondStretchGlobalScaleFactor = customBondBendGlobalScaleFactor =
            customBondTorsionGlobalScaleFactor = 1;

        vdwScale12 = coulombScale12 = vdwScale13 = coulombScale13 = 0;
        vdwScale14 = coulombScale14 = vdwScale15 = coulombScale15 = 1;

        nonbondedCutoff = 2.0;
        nonbondedMethod = 0;

        tracing = true; // EU

        const DuMM::ClusterIndex gid = addCluster(Cluster("free atoms and groups"));
        SimTK_ASSERT_ALWAYS(gid == 0, "Free atoms and groups cluster should have been the first one.");

        internalListsRealized = false; // EU
    }

    ~DuMMForceFieldSubsystemRep() {
    }

    [[nodiscard]] auto getNonBondedMappings() const -> const NonBondedMappings& {
        return nonBondedMappings;
    }

    // common checks when defining improper and proper torsions
    void defineAnyTorsion(DuMM::AtomClassIndex class1,
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
                          const char* CallingMethodName) const;

    // common checks when defining improper and proper torsions
    void defineAnyTorsion(DuMM::AtomClassIndex class1,
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
                          const char* CallingMethodName) const;

    // Common checks when defining improper and proper torsions
    // using 5 parameters.
    void defineAnyTorsion(DuMM::AtomClassIndex class1,
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
                          const char* CallingMethodName) const;


    bool isValidElement(int atomicNumber) const {
        if (1 > atomicNumber) {
            return false;
        }
        try {
            Element::getByAtomicNumber(atomicNumber);
            return true;
        } catch (...) {
            return false;
        }
    }

    bool isValidAtom(DuMM::AtomIndex atomNum) const {
        return 0 <= atomNum && atomNum < atoms.size() && atoms[atomNum].isValid();
    }

    bool isValidBond(DuMM::BondIndex bondNum) const {
        return 0 <= bondNum && bondNum < bonds.size() && bonds[bondNum].isValid();
    }

    bool isValidCluster(DuMM::ClusterIndex clusterIndex) const {
        return 0 <= clusterIndex && clusterIndex < clusters.size() && clusters[clusterIndex].isValid();
    }

    bool isValidDuMMBody(DuMMBodyIndex bodyIx) const {
        return 0 <= bodyIx && bodyIx < duMMSubsetOfBodies.size() && duMMSubsetOfBodies[bodyIx].isValid();
    }

    bool isValidChargedAtomType(DuMM::ChargedAtomTypeIndex typeNum) const {
        return 0 <= typeNum && typeNum < chargedAtomTypes.size() && chargedAtomTypes[typeNum].isValid();
    }

    bool isValidAtomClass(DuMM::AtomClassIndex classNum) const {
        return 0 <= classNum && classNum < atomClasses.size() && atomClasses[classNum].isValid();
    }

    bool hasAtomClass(DuMM::AtomClassIndex classNum) const {
        if (classNum < 0) {
            return false;
        }
        if (classNum >= (DuMM::AtomClassIndex)atomClasses.size()) {
            return false;
        }
        if (!atomClasses[classNum].isValid()) {
            return false;
        }
        return true;
    }

    bool hasAtomClass(const String& atomClassName) const {
        if (atomClassIndicesByName.find(atomClassName) == atomClassIndicesByName.end()) {
            return false;
        }
        DuMM::AtomClassIndex classNum = atomClassIndicesByName.find(atomClassName)->second;
        return hasAtomClass(classNum);
    }

    DuMM::AtomClassIndex getAtomClassIndex(const String& atomClassName) const {
        if (!hasAtomClass(atomClassName)) {
            throw(std::range_error(String("no such atom class name: ") + atomClassName));
        }
        DuMM::AtomClassIndex classNum = atomClassIndicesByName.find(atomClassName)->second;
        return classNum;
    }

    DuMM::AtomClassIndex getNextUnusedAtomClassIndex() const {
        return nextUnusedAtomClassIndex;
    }

    bool hasChargedAtomType(DuMM::ChargedAtomTypeIndex chargedTypeIndex) const {
        if (chargedTypeIndex < 0) {
            return false;
        }
        if (chargedTypeIndex >= chargedAtomTypes.size()) {
            return false;
        }
        if (!chargedAtomTypes[chargedTypeIndex].isValid()) {
            return false;
        }
        return true;
    }
    bool hasChargedAtomType(const String& chargedTypeName) const {
        if (chargedAtomTypeIndicesByName.find(chargedTypeName) == chargedAtomTypeIndicesByName.end()) {
            return false;
        }
        DuMM::ChargedAtomTypeIndex typeIndex = chargedAtomTypeIndicesByName.find(chargedTypeName)->second;
        return hasChargedAtomType(typeIndex);
    }
    DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex(const String& chargedTypeName) const {
        if (!hasChargedAtomType(chargedTypeName)) {
            throw(std::range_error(String("no such charged atom type name: ") + chargedTypeName));
        }
        DuMM::ChargedAtomTypeIndex typeIndex = chargedAtomTypeIndicesByName.find(chargedTypeName)->second;
        return typeIndex;
    }
    DuMM::ChargedAtomTypeIndex getNextUnusedChargedAtomTypeIndex() const {
        return nextUnusedChargedAtomTypeIndex;
    }


    /*! <!--
    --> */
    void setAtomMass(SimTK::DuMM::AtomIndex dAIx, SimTK::mdunits::Mass atomicMass) {
        DuMMAtom dAtom = updAtom(dAIx);
        atoms[dAIx].setMass(atomicMass);
    }

    void insertNewChargedAtomType(const ChargedAtomType& chargedAtomType) {
        // Update nextUnusedChargedAtomType, if necessary
        if (chargedAtomType.chargedAtomTypeIndex >= nextUnusedChargedAtomTypeIndex) {
            nextUnusedChargedAtomTypeIndex = chargedAtomType.chargedAtomTypeIndex;
            ++nextUnusedChargedAtomTypeIndex;
        }

        chargedAtomTypeIndicesByName[chargedAtomType.name] = chargedAtomType.chargedAtomTypeIndex;

        // Define the new charged atom type.
        chargedAtomTypes[chargedAtomType.chargedAtomTypeIndex] = chargedAtomType;
    }

    void insertNewAtomClass(const AtomClass& atomClass) {
        if (atomClass.atomClassIx >= nextUnusedAtomClassIndex) {
            nextUnusedAtomClassIndex = atomClass.atomClassIx;
            ++nextUnusedAtomClassIndex;
        }

        atomClassIndicesByName[atomClass.name] = atomClass.atomClassIx;

        atomClasses[atomClass.atomClassIx] = atomClass;
    }

    auto addCluster(const Cluster& cluster) -> DuMM::ClusterIndex {
        invalidateSubsystemTopologyCache();

        const DuMM::ClusterIndex clusterIndex = (DuMM::ClusterIndex)clusters.size();
        clusters.push_back(cluster);
        clusters[clusterIndex].clusterIndex = clusterIndex;
        return clusterIndex;
    }
    auto updCluster(DuMM::ClusterIndex clusterIndex) -> Cluster& {
        assert(isValidCluster(clusterIndex));

        invalidateSubsystemTopologyCache();
        return clusters[clusterIndex];
    }
    auto getCluster(DuMM::ClusterIndex clusterIndex) const -> const Cluster& {
        assert(isValidCluster(clusterIndex));
        return clusters[clusterIndex];
    }
    auto getNumClusters() const -> std::size_t {
        return clusters.size();
    }

    auto updDuMMBody(DuMMBodyIndex bodyIx) -> DuMMBody& {
        assert(isValidDuMMBody(bodyIx));
        invalidateSubsystemTopologyCache();
        return duMMSubsetOfBodies[bodyIx];
    }
    auto getDuMMBody(DuMMBodyIndex duMMBodyIx) const -> const DuMMBody& {
        assert(isValidDuMMBody(duMMBodyIx));
        return duMMSubsetOfBodies[duMMBodyIx];
    }

    auto getNumBonds() const -> int {
        return (int)bonds.size();
    }

    // DuMM atoms (all atoms in the molecule)
    auto getNumAtoms() const -> int {
        return (int)atoms.size();
    }
    auto getAtom(DuMM::AtomIndex atomIndex) const -> const DuMMAtom& {
        assert(isValidAtom(atomIndex));
        return atoms[atomIndex];
    }
    auto updAtom(DuMM::AtomIndex atomIndex) -> DuMMAtom& {
        assert(isValidAtom(atomIndex));
        return atoms[atomIndex];
    }


    /*! <!-- desk_mass_related --> */
    const SimTK::mdunits::Mass getAtomMass(DuMM::AtomIndex atomIndex) const {
        assert(isValidAtom(atomIndex));
        return getAtom(atomIndex).getMass();
        // return atoms[atomIndex].getMass();
    }

    // Included bodies
    int getNumIncludedBodies() const {
        return (int)includedBodies.size();
    }
    const IncludedBody& getIncludedBody(DuMMIncludedBodyIndex incBodyIx) const {
        return includedBodies[incBodyIx];
    }
    IncludedBody& updIncludedBody(DuMMIncludedBodyIndex incBodyIx) {
        return includedBodies[incBodyIx];
    }

    // Included atoms
    int getNumIncludedAtoms() const {
        return (int)includedAtoms.size();
    }
    DuMM::AtomIndex getAtomIndexOfIncludedAtom(DuMM::IncludedAtomIndex incAtomIx) const {
        return includedAtoms[incAtomIx].atomIndex;
    }
    const IncludedAtom& getIncludedAtom(DuMM::IncludedAtomIndex incAtomIx) const {
        return includedAtoms[incAtomIx];
    }
    IncludedAtom& updIncludedAtom(DuMM::IncludedAtomIndex incAtomIx) {
        return includedAtoms[incAtomIx];
    }
    const Vec3& getIncludedAtomStation(DuMM::IncludedAtomIndex incAtomIx) const {
        return includedAtomStations[incAtomIx];
    }
    Vec3& updIncludedAtomStation(DuMM::IncludedAtomIndex incAtomIx) {
        return includedAtomStations[incAtomIx];
    }

    // Nonbond included atoms
    int getNumNonbondAtoms() const {
        return (int)nonbondAtoms.size();
    }
    DuMM::IncludedAtomIndex getIncludedAtomIndexOfNonbondAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return nonbondAtoms[nbAtomIx];
    }
    DuMM::AtomIndex getAtomIndexOfNonbondAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return getAtomIndexOfIncludedAtom(nonbondAtoms[nbAtomIx]);
    }
    const IncludedAtom& getNonbondIncludedAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return getIncludedAtom(nonbondAtoms[nbAtomIx]);
    }


    // for GMolModel - All

    // All bodies
    int getNumAllBodies() const {
        return (int)AllBodies.size();
    }
    const IncludedBody& getAllBody(DuMMIncludedBodyIndex incBodyIx) const {
        return AllBodies[incBodyIx];
    }
    IncludedBody& updAllBody(DuMMIncludedBodyIndex incBodyIx) {
        return AllBodies[incBodyIx];
    }

    // All atoms
    int getNumAllAtoms() const {
        return (int)AllAtoms.size();
    }
    DuMM::AtomIndex getAtomIndexOfAllAtom(DuMM::IncludedAtomIndex incAtomIx) const {
        return AllAtoms[incAtomIx].atomIndex;
    }
    const IncludedAtom& getAllAtom(DuMM::IncludedAtomIndex incAtomIx) const {
        return AllAtoms[incAtomIx];
    }
    IncludedAtom& updAllAtom(DuMM::IncludedAtomIndex incAtomIx) {
        return AllAtoms[incAtomIx];
    }
    const Vec3& getAllAtomStation(DuMM::IncludedAtomIndex incAtomIx) const {
        return AllAtomStations[incAtomIx];
    }
    Vec3& updAllAtomStation(DuMM::IncludedAtomIndex incAtomIx) {
        return AllAtomStations[incAtomIx];
    }

    // Nonbond All atoms
    int getNumAllNonbondAtoms() const {
        return (int)AllnonbondAtoms.size();
    }
    DuMM::IncludedAtomIndex getAllAtomIndexOfNonbondAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return AllnonbondAtoms[nbAtomIx];
    }
    DuMM::AtomIndex getAtomIndexOfAllNonbondAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return getAtomIndexOfAllAtom(AllnonbondAtoms[nbAtomIx]);
    }
    const IncludedAtom& getNonbondAllAtom(DuMM::NonbondAtomIndex nbAtomIx) const {
        return getAllAtom(AllnonbondAtoms[nbAtomIx]);
    }

    /*
    int getNumAllBondStarterAtoms()  const {return (int)AllbondStarterAtoms.size();}
    DuMM::IncludedAtomIndex getAllAtomIndexOfBondStarterAtom
       (DuMMBondStarterIndex bsAtomIx) const
    {   return bondStarterAtoms[bsAtomIx]; }
    DuMM::AtomIndex getAtomIndexOfBondStarterAtom
       (DuMMBondStarterIndex bsAtomIx) const
    {   return getAtomIndexOfAllAtom(bondStarterAtoms[bsAtomIx]); }
    const IncludedAtom& getBondStarterAllAtom
       (DuMMBondStarterIndex bsAtomIx) const
    {   return getAllAtom(bondStarterAtoms[bsAtomIx]); }
     */


    /////////////////////////


    int getNumBondStarterAtoms() const {
        return (int)bondStarterAtoms.size();
    }
    DuMM::IncludedAtomIndex getIncludedAtomIndexOfBondStarterAtom(DuMMBondStarterIndex bsAtomIx) const {
        return bondStarterAtoms[bsAtomIx];
    }
    DuMM::AtomIndex getAtomIndexOfBondStarterAtom(DuMMBondStarterIndex bsAtomIx) const {
        return getAtomIndexOfIncludedAtom(bondStarterAtoms[bsAtomIx]);
    }
    const IncludedAtom& getBondStarterIncludedAtom(DuMMBondStarterIndex bsAtomIx) const {
        return getIncludedAtom(bondStarterAtoms[bsAtomIx]);
    }

    DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex(DuMM::AtomIndex atomIx) const {
        return getAtom(atomIx).chargedAtomTypeIndex;
    }

    DuMM::AtomClassIndex getAtomClassIndex(DuMM::AtomIndex atomIndex) const {
        const ChargedAtomType& type = chargedAtomTypes[getChargedAtomTypeIndex(atomIndex)];
        return type.atomClassIx;
    }

    int getAtomElementNum(DuMM::AtomIndex atomIndex) const {
        const AtomClass& cl = atomClasses[getAtomClassIndex(atomIndex)];
        return cl.element;
    }

    Element getElement(int atomicNumber) const {
        assert(isValidElement(atomicNumber));
        return Element::getByAtomicNumber(atomicNumber);
    }


    // These return 0 if no entry can be found for the given series of atom classes.
    const BondStretch* getBondStretch(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2) const;
    const BondBend*
    getBondBend(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3) const;
    const BondTorsion* getBondTorsion(DuMM::AtomClassIndex class1,
                                      DuMM::AtomClassIndex class2,
                                      DuMM::AtomClassIndex class3,
                                      DuMM::AtomClassIndex class4) const;
    const BondTorsion* getAmberImproperTorsion(DuMM::AtomClassIndex class1,
                                               DuMM::AtomClassIndex class2,
                                               DuMM::AtomClassIndex class3,
                                               DuMM::AtomClassIndex class4) const;

    // We don't realize included atom velocities unless someone asks for them.
    // This won't do anything if the atom velocities are already realized.
    void realizeIncludedAtomVelocityCache(const State&) const;

    // We won't calculate forces and energy until Dynamics stage unless
    // someone asks for them earlier (they only depend on Positions).
    void realizeForcesAndEnergy(const State&) const;

    // Override virtual methods from Subsystem::Guts class.

    DuMMForceFieldSubsystemRep* cloneImpl() const {
        // return new DuMMForceFieldSubsystemRep(*this);
        // TODO fix this please
        assert(!"NOT IMPLEMENTED");
        return nullptr;
    }

    // Figure out the molecular topology and how the atoms are mapped
    // onto bodies.
    // EU BEGIN
    // void markInternalListsRealized();
    int realizeInternalLists(State& s) const;
    // EU END
    int realizeSubsystemTopologyImpl(State& s) const;

    int realizeSubsystemModelImpl(State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    // Calculate the ground-frame locations of every included atom.
    int realizeSubsystemPositionImpl(const State& s) const;

    // Nothing to do here since included atom velocities are lazy.
    int realizeSubsystemVelocityImpl(const State& state) const;

    // Calculate forces.
    int realizeSubsystemDynamicsImpl(const State& s) const;

    int realizeSubsystemAccelerationImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    // This will cause realization of at least the energy (probably
    // forces too) if it hasn't already been realized since the
    // last change to Position-stage state variables.
    Real calcPotentialEnergy(const State& state) const;

    void dump() const;

    // How many times has the forcefield been evaluated?
    long long getForceEvaluationCount() const {
        return forceEvaluationCount;
    }

    std::ostream& generateBiotypeChargedAtomTypeSelfCode(std::ostream& os) const;
    DuMM::ChargedAtomTypeIndex getBiotypeChargedAtomType(BiotypeIndex biotypeIx) const;

    SimTK_DOWNCAST(DuMMForceFieldSubsystemRep, ForceSubsystem::Guts);
    SimTK_DOWNCAST(DuMMForceFieldSubsystemRep, Subsystem::Guts);

    protected:
    std::ostream& generateBiotypeChargedAtomTypeSelfCode(std::ostream& os, BiotypeIndex biotypeIx) const;
    void setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, BiotypeIndex biotypeIx);

    private:
    DuMMBodyIndex ensureDuMMBodyEntryExists(MobilizedBodyIndex mobodIx) {
        DuMMBodyIndex duMMBodyIndex;

        if (dummBodyIndexByMobodIndex.find(mobodIx) == dummBodyIndexByMobodIndex.end()) {
            // Create a new DuMMBody for this MobilizedBody
            duMMBodyIndex = (DuMMBodyIndex)duMMSubsetOfBodies.size();

            const DuMM::ClusterIndex clusterIndex =
                addCluster(Cluster(DuMMBody::createClusterNameForBody(mobodIx).c_str()));
            clusters[clusterIndex].attachToBody(mobodIx, Transform(), *this);

            duMMSubsetOfBodies.push_back(DuMMBody(clusterIndex, mobodIx));
            dummBodyIndexByMobodIndex[mobodIx] = duMMBodyIndex;
        } else {
            // Sanity check of preexisting DuMMBody
            duMMBodyIndex = dummBodyIndexByMobodIndex.at(mobodIx);
        }

        assert(duMMSubsetOfBodies[duMMBodyIndex].isValid());

        return duMMBodyIndex;
    }

    void invalidateNecTopologicalCacheEntries() {
        // If any of these objects are invalid, the invalidateTopologicalCache() call does nothing (i.e., it
        // doesn't blow up!).
        inclAtomStationCacheIndex.invalidate();
        inclAtomPositionCacheIndex.invalidate();
        inclAtomVelocityCacheIndex.invalidate();
        inclAtomForceCacheIndex.invalidate();
        inclBodyForceCacheIndex.invalidate();
        energyCacheIndex.invalidate();
    }

    void invalidateAllTopologicalCacheEntries() {
        // If any of these objects are invalid, the invalidateTopologicalCache()
        // call does nothing (i.e., it doesn't blow up!).

        // molecule
        for (DuMM::AtomIndex i(0); i < atoms.size(); ++i) {
            atoms[i].invalidateTopologicalCache();
        }
        for (DuMM::ClusterIndex i(0); i < clusters.size(); ++i) {
            clusters[i].invalidateTopologicalCache();
        }
        for (DuMMBodyIndex i(0); i < duMMSubsetOfBodies.size(); ++i) {
            duMMSubsetOfBodies[i].invalidateTopologicalCache();
        }

        // force field
        for (DuMM::AtomClassIndex i(0); i < atomClasses.size(); ++i) {
            atomClasses[i].invalidateTopologicalCache();
        }

        inclAtomStationCacheIndex.invalidate();
        inclAtomPositionCacheIndex.invalidate();
        inclAtomVelocityCacheIndex.invalidate();
        inclAtomForceCacheIndex.invalidate();
        inclBodyForceCacheIndex.invalidate();
        energyCacheIndex.invalidate();
    }

    const Vector_<Vec3>& getIncludedAtomStationsInG(const State& s) const {
        return getIncludedAtomStationCache(s);
    }
    const Vector_<Vec3>& getIncludedAtomPositionsInG(const State& s) const {
        return getIncludedAtomPositionCache(s);
    }

    // Atom velocities won't be realized unless someone asks.
    const Vector_<Vec3>& getIncludedAtomVelocitiesInG(const State& s) const {
        realizeIncludedAtomVelocityCache(s);
        return getIncludedAtomVelocityCache(s);
    }

    // Included atom forces can be realized any time after position stage but
    // we normally delay until dynamics stage unless someone asks earlier.
    const Vector_<Vec3>& getIncludedAtomForcesInG(const State& s) const {
        realizeForcesAndEnergy(s);
        return getIncludedAtomForceCache(s);
    }

    // Included body forces are available at the same time as atom forces.
    const Vector_<SpatialVec>& getIncludedBodyForcesInG(const State& s) const {
        realizeForcesAndEnergy(s);
        return getIncludedBodyForceCache(s);
    }

    // Access to cache entries.

    Vector_<Vec3>& updIncludedAtomStationCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, inclAtomStationCacheIndex));
    }
    const Vector_<Vec3>& getIncludedAtomStationCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, inclAtomStationCacheIndex));
    }

    Vector_<Vec3>& updIncludedAtomPositionCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, inclAtomPositionCacheIndex));
    }
    const Vector_<Vec3>& getIncludedAtomPositionCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, inclAtomPositionCacheIndex));
    }


    // GMolModel
    Vector_<Vec3>& updAllAtomStationCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, AllAtomStationCacheIndex));
    }
    const Vector_<Vec3>& getAllAtomStationCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, AllAtomStationCacheIndex));
    }
    Vector_<Vec3>& updAllAtomPositionCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, AllAtomPositionCacheIndex));
    }
    const Vector_<Vec3>& getAllAtomPositionCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, AllAtomPositionCacheIndex));
    }


    // Atom velocities are lazy evaluated.
    Vector_<Vec3>& updIncludedAtomVelocityCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, inclAtomVelocityCacheIndex));
    }
    const Vector_<Vec3>& getIncludedAtomVelocityCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, inclAtomVelocityCacheIndex));
    }
    bool isIncludedAtomVelocityCacheRealized(const State& s) const {
        return isCacheValueRealized(s, inclAtomVelocityCacheIndex);
    }
    void markIncludedAtomVelocityCacheRealized(const State& s) const {
        markCacheValueRealized(s, inclAtomVelocityCacheIndex);
    }


    // Potential energy can be realized any time after Position stage
    // but won't be until someone asks for it.
    Real& updEnergyCache(const State& s) const {
        return Value<Real>::downcast(updCacheEntry(s, energyCacheIndex));
    }
    Real getEnergyCache(const State& s) const {
        return Value<Real>::downcast(getCacheEntry(s, energyCacheIndex));
    }
    bool isEnergyCacheRealized(const State& s) const {
        return isCacheValueRealized(s, energyCacheIndex);
    }
    void markEnergyCacheRealized(const State& s) const {
        markCacheValueRealized(s, energyCacheIndex);
    }

    // Forces can be realized any time after Position stage but won't be until
    // someone asks for them.

    // We cache the forces on each included atom.
    Vector_<Vec3>& updIncludedAtomForceCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(updCacheEntry(s, inclAtomForceCacheIndex));
    }
    const Vector_<Vec3>& getIncludedAtomForceCache(const State& s) const {
        return Value<Vector_<Vec3>>::downcast(getCacheEntry(s, inclAtomForceCacheIndex));
    }
    bool isIncludedAtomForceCacheRealized(const State& s) const {
        return isCacheValueRealized(s, inclAtomForceCacheIndex);
    }
    void markIncludedAtomForceCacheRealized(const State& s) const {
        markCacheValueRealized(s, inclAtomForceCacheIndex);
    }

    // We cache the accumulated rigid body force and torque on each of the
    // included bodies.
    Vector_<SpatialVec>& updIncludedBodyForceCache(const State& s) const {
        return Value<Vector_<SpatialVec>>::downcast(updCacheEntry(s, inclBodyForceCacheIndex));
    }
    const Vector_<SpatialVec>& getIncludedBodyForceCache(const State& s) const {
        return Value<Vector_<SpatialVec>>::downcast(getCacheEntry(s, inclBodyForceCacheIndex));
    }
    bool isIncludedBodyForceCacheRealized(const State& s) const {
        return isCacheValueRealized(s, inclBodyForceCacheIndex);
    }
    void markIncludedBodyForceCacheRealized(const State& s) const {
        markCacheValueRealized(s, inclBodyForceCacheIndex);
    }

    private:
    // TODO: I made this public so that the OpenMM Plugin could see these data
    // members. It should use an appropriate set of access methods instead.
    public:
    String forcefieldName;

    // keep track of how many forceEvaluations have been computed
    mutable long long forceEvaluationCount;

    // TOPOLOGICAL STATE VARIABLES
    //   Filled in during construction.

    // molecule

    Array_<DuMMAtom, DuMM::AtomIndex> atoms;
    Array_<Bond, DuMM::BondIndex> bonds;
    Array_<Cluster, DuMM::ClusterIndex> clusters;

    std::map<MobilizedBodyIndex, DuMMBodyIndex> dummBodyIndexByMobodIndex;

    // This defines the partitioning of atoms onto the matter subsystem's
    // bodies. The indices here correspond to the DuMM body numbers (not
    // MobilizedBodyIndices). Only entries for bodies on which our atoms have
    // been attached will be valid.
    Array_<DuMMBody, DuMMBodyIndex> duMMSubsetOfBodies;

    // This collects user instructions about what atoms should be involved in
    // force calculations (for Sam Flores' "physics where you want it" idea).
    InclusionListSpec inclList;

    // force field

    std::map<BiotypeIndex, DuMM::ChargedAtomTypeIndex> chargedAtomTypesByBiotype;

    DuMM::AtomClassIndex nextUnusedAtomClassIndex;
    std::map<String, DuMM::AtomClassIndex> atomClassIndicesByName;

    DuMM::ChargedAtomTypeIndex nextUnusedChargedAtomTypeIndex;
    std::map<String, DuMM::ChargedAtomTypeIndex> chargedAtomTypeIndicesByName;


    // Force field description. These are not necessarily fully populated;
    // check the "isValid()" method to see if anything is there.
    Array_<AtomClass, DuMM::AtomClassIndex> atomClasses;
    Array_<ChargedAtomType, DuMM::ChargedAtomTypeIndex> chargedAtomTypes;

    // These relate atom classes, not charged atom types.
    std::map<AtomClassIndexPair, BondStretch> bondStretch;
    std::map<AtomClassIndexTriple, BondBend> bondBend;
    std::map<AtomClassIndexQuad, BondTorsion> bondTorsion;
    std::map<AtomClassIndexQuad, BondTorsion> amberImproperTorsion;

    // GMolModel
    std::map<AtomClassIndexPair, BondStretch> bondStretchAll;
    std::map<AtomClassIndexTriple, BondBend> bondBendAll;
    std::map<AtomClassIndexQuad, BondTorsion> bondTorsionAll;
    std::map<AtomClassIndexQuad, BondTorsion> amberImproperTorsionAll;

    // Scale factors for nonbonded forces when applied to
    // atoms which are near in the graph formed by the bonds.
    Real vdwScale12, coulombScale12; // default 0,0
    Real vdwScale13, coulombScale13; // default 0,0
    Real vdwScale14, coulombScale14; // default 1,1
    Real vdwScale15, coulombScale15; // default 1,1

    // Global scale factors for non-physical disabling or fiddling with
    // individual force field terms.
    Real vdwGlobalScaleFactor, coulombGlobalScaleFactor, gbsaGlobalScaleFactor;
    Real bondStretchGlobalScaleFactor, bondBendGlobalScaleFactor, bondTorsionGlobalScaleFactor,
        amberImproperTorsionGlobalScaleFactor;
    Real customBondStretchGlobalScaleFactor, customBondBendGlobalScaleFactor,
        customBondTorsionGlobalScaleFactor;

    // These affect GBSA.
    bool gbsaIncludeAceApproximation;
    Real gbsaSolventDielectric; // typically 80 for water
    Real gbsaSoluteDielectric;  // typically 1 or 2 for protein

    Real nonbondedCutoff;
    int nonbondedMethod;

    bool tracing; // for debugging

    // TOPOLOGICAL CACHE ENTRIES
    //   These cache entries are allocated in realizeTopology().

    // Arrays set up for fast computation.

    // This is the list of all bodies on which we will be generating forces,
    // sorted in increasing order of MobilizedBodyIndex. The presence of an
    // entry here indicates the body has at least one includedAtom attached.
    Array_<IncludedBody, DuMMIncludedBodyIndex> includedBodies;
    Array_<IncludedBody, DuMMIncludedBodyIndex> AllBodies;
    // This is the list of all atoms that participate in *any* force
    // calculation, nonbonded or bonded. They are grouped in the same order
    // as the includedBodies list so that all the included atoms for the first
    // included body come first, then all included atoms for the second
    // included body, etc. Use these entries to index the atoms array.
    Array_<IncludedAtom, DuMM::IncludedAtomIndex> includedAtoms;
    Array_<IncludedAtom, DuMM::IncludedAtomIndex> AllAtoms;
    // These are the stations for each included atom on its included body. These
    // are kept separately since they are only needed during realizePosition()
    // when calculating the atom locations, but logically they are part of the
    // IncludedAtom information so this array is the same length as includedAtoms.
    Array_<Vec3, DuMM::IncludedAtomIndex> includedAtomStations;
    Array_<Vec3, DuMM::IncludedAtomIndex> AllAtomStations;

    Array_<DuMM::IncludedAtomIndex, DuMM::NonbondAtomIndex> nonbondAtoms;
    Array_<DuMM::IncludedAtomIndex, DuMMBondStarterIndex> bondStarterAtoms;

    // GMolModel
    Array_<DuMM::IncludedAtomIndex, DuMM::NonbondAtomIndex> AllnonbondAtoms;
    Array_<DuMM::IncludedAtomIndex, DuMMBondStarterIndex> AllbondStarterAtoms;

    CacheEntryIndex inclAtomStationCacheIndex;
    CacheEntryIndex inclAtomPositionCacheIndex;
    CacheEntryIndex inclAtomVelocityCacheIndex;
    CacheEntryIndex inclAtomForceCacheIndex;
    CacheEntryIndex inclBodyForceCacheIndex;
    CacheEntryIndex energyCacheIndex;

    // GMolModel
    CacheEntryIndex AllAtomStationCacheIndex;
    CacheEntryIndex AllAtomPositionCacheIndex;
    bool internalListsRealized; // EU

    mutable NonBondedMappings nonBondedMappings;
};


#endif // SimTK_MOLMODEL_DUMM_FORCE_FIELD_SUBSYSTEM_REP_H_
