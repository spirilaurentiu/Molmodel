#ifndef SimTK_COMPOUNDREP_H_
#define SimTK_COMPOUNDREP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

#include "molmodel/internal/bondGeometry.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/CompoundSystem.h"
#include "molmodel/internal/Superpose.h"

#include "CompoundAtom.h"
#include <fstream>

#include <vector>
#include <map>
#include <string>
#include <set>
using std::string;

namespace SimTK {

static const String InboardBondName = "inboard bond";


//void invalidateAtomFrameCache(std::vector<Transform>& atomFrameCache, int numAtoms);



class DihedralAngle {
public:
    DihedralAngle() 
        : nomenclatureOffset(0*Deg2Rad)
    {}

    DihedralAngle(Compound::BondCenterIndex bc1, Compound::BondCenterIndex bc2, Angle o = 0*Deg2Rad) 
        : nomenclatureOffset(o), bondCenter1(bc1), bondCenter2(bc2)
        // , internalOffset(0)
    {}

    Compound::BondCenterIndex getBondCenter1Id() const {return bondCenter1;}
    Compound::BondCenterIndex getBondCenter2Id() const {return bondCenter2;}
    Angle getNomenclatureOffset() const {return nomenclatureOffset;}

    // Angle getInternalOffset() const {return internalOffset;}
    // void setInternalOffset(Angle a) {internalOffset = a;}

private:
    // Angle internalOffset; // this dihedral may be offset from canonical dihedral for this bond
    Angle nomenclatureOffset; // internal + offset = nominal
    Compound::BondCenterIndex bondCenter1;
    Compound::BondCenterIndex bondCenter2;
};

    //////////////////
    // COMPOUND REP //
    //////////////////

class CompoundRep : public PIMPLImplementation<Compound,CompoundRep> 
{
public:
    friend class CompoundSystem;
    friend class ResidueInfo;

    explicit CompoundRep(const String& n="UnknownCompoundType", const Transform& transform = Transform()) :
        ownerSystem(0),
        topLevelTransform(transform),
        name(n),
        pdbResidueNumber(-111111),
        pdbResidueName("UNK"),
        pdbChainId(' ')
        //haveParentCompound(false)
    {
    	addCompoundSynonym(n);
    }

    CompoundRep&  setCompoundName(const Compound::Name& n) {
        name=n; 
        addCompoundSynonym(name);
        return *this;
    }
    CompoundRep& addCompoundSynonym(Compound::Name synonym) {
        synonyms.insert(synonym);
        return *this;
    }
    const Compound::Name& getCompoundName() const          {return name;}

    // This is concrete, but can be extended.
    virtual ~CompoundRep() { }
    virtual CompoundRep* clone() const {return new CompoundRep(*this);}

    // int getNumSubcompounds() const {return allSubcompounds.size();}

    void setMultibodySystem(MultibodySystem& system) 
    {
        ownerSystem = &system;

    }

    bool isOwnedBySystem() const {return ownerSystem != 0;}
    // const CompoundSystem& getOwnerCompoundSystem() const {assert(ownerSystem); return *ownerSystem;}
    const MultibodySystem& getOwnerMultibodySystem() const {assert(ownerSystem); return *ownerSystem;}
    // Compound::Index getIdWithinOwnerCompoundSystem() const {assert(ownerSystem); return ixWithinOwnerSystem;}

    // Add one simple atom unconnected to anything else
    CompoundRep& setBaseAtom(
        const Compound::AtomName& name, 
        const Element& element,
        const Transform& location);

    CompoundRep& setBaseAtom(
        const Compound::AtomName& name, 
        const Biotype& biotype,
        const Transform& location);

    // Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
    // This atom is not connected to anything else
    CompoundRep& setBaseAtom(
        const Compound::SingleAtom& compound,
        const Transform& location);

    // Add a subcompound without attaching it to anything
    Compound::BondCenterIndex setBaseCompound(
        const Compound::Name& n, 
        const Compound& c,
        const Transform& location);

    // Add a subcompound containing exactly one atom, so the Compound::AtomName can be reused for the Compound::Name
    // This atom is connected to existing material
    CompoundRep& bondAtom(
        const Compound::SingleAtom&   compound, 
        const Compound::BondCenterPathName& parentBondName, 
        mdunits::Length                      distance,
        Angle                         dihedral = 0,
        BondMobility::Mobility        mobility = BondMobility::Default
        );

    // Bond atom using default bond length and dihedral angle
    CompoundRep& bondAtom(
        const Compound::SingleAtom& compound, 
        const Compound::BondCenterPathName& parentBondName) 
    {
        // assert(! hasParentCompound());
        // assert(! compound.getImpl().hasParentCompound());

        // There are two choice for how to delegate this method.
        // 1) deduce the bond length and dihedral and call the other 
        //    bondAtom() that takes those parameters
        // 2) deduce the compound name and call the bondCompound()
        //    that does not take geometry parameters.
        // 
        // It is better to choose course 2), because it is less complex
        // (at the time of this writing)

        const CompoundRep& scRep = compound.getImpl();
        Compound::AtomName atomName = scRep.atomName_To_atomId.begin()->first;

        bondCompound(atomName, compound, parentBondName);
        inheritAtomNames(atomName);

        // assert(getSubcompound(atomName).getImpl().hasParentCompound());

        return *this;
    }

    // Add a subcompound attached by a bond to an existing atom
    // bondCompound("H1", MonovalentAtom(Element::Hydrogen()), "bond", "C/bond2", C_Hdistance );
    CompoundRep& bondCompound(
        const Compound::Name& name, 
        const Compound& subcompound, 
        const Compound::BondCenterPathName& parentBondName, 
        mdunits::Length distance,
        Angle dihedral = 0,
        BondMobility::Mobility mobility = BondMobility::Default
        );

    // Shorter version uses default bond length and dihedral angle
    CompoundRep& bondCompound(
        const Compound::Name& n, 
        const Compound& c, 
        const Compound::BondCenterPathName& parentBondName);
    // sam added polymorphism
  
    CompoundRep& bondCompound(
        const Compound::Name& n,
        const Compound& c,
        const Compound::BondCenterPathName& parentBondName,
        BondMobility::Mobility mobility
        );
    // deprecate removeSubcompound for now -- I'm not using it
    // CompoundRep& removeSubcompound(const Compound::Name& name);

    CompoundRep& setInboardBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName, 
        Angle zRotation,
        Angle oldXRotation);

    CompoundRep& setDefaultInboardBondLength(mdunits::Length d) {
        updInboardBondCenter().setDefaultBondLength(d);
        return *this;
    }

    CompoundRep& setDefaultInboardDihedralAngle(Angle a) {
        updInboardBondCenter().setDefaultDihedralAngle(a);
        return *this;
    }


    CompoundRep& addFirstBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName);

    CompoundRep& addSecondBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1
        );

    CompoundRep& addFirstTwoBondCenters(
            const Compound::BondCenterName& centerName1,
            const Compound::BondCenterName& centerName2,
            const Compound::AtomName& atomName,
            UnitVec3 dir1, UnitVec3 dir2
    );

    CompoundRep& addPlanarBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2);

    CompoundRep& addRightHandedBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2
        );

    CompoundRep& addLeftHandedBondCenter(
        const Compound::BondCenterName& centerName, 
        const Compound::AtomName& atomName,
        Angle bondAngle1,
        Angle bondAngle2
        );


    CompoundRep& addBondCenterInfo(
        const Compound::AtomIndex   atomId,
        const CompoundAtom::BondCenterIndex atomCenterIndex);

    CompoundRep& addRingClosingBond(
        const Compound::BondCenterName& centerName1, 
        const Compound::BondCenterName& centerName2 
        );
    CompoundRep& addRingClosingBond(
        const Compound::BondCenterName& centerName1, 
        const Compound::BondCenterName& centerName2,
        mdunits::Length bondLength,
        Angle dihedral,
        BondMobility::Mobility mobility
        );

    int getNumAtoms() const;

    const Compound::AtomName getAtomName(Compound::AtomIndex) const;

    const Element& getAtomElement(Compound::AtomIndex atomIndex) const {
        return getAtom(atomIndex).getElement();
    }

    const Element& getAtomElement(const Compound::AtomName& atomName) const {
        return getAtom(atomName).getElement();
    }


    BiotypeIndex getAtomBiotypeIndex(Compound::AtomIndex) const;
    void setDuMMAtomIndex(Compound::AtomIndex aid, DuMM::AtomIndex dummId) {
        updAtom(aid).setDuMMAtomIndex(dummId);
    }
    DuMM::AtomIndex getDuMMAtomIndex(Compound::AtomIndex aid) const {
        return getAtom(aid).getDuMMAtomIndex();
    }

    size_t getNumBondCenters() const;
	size_t getNumBondCenters(Compound::AtomIndex atomIndex) const;

    // const Compound::BondCenterName& getBondCenterName(Compound::BondCenterIndex bondCenterIndex) const;

    CompoundRep& nameAtom(const Compound::AtomName& newName, Compound::AtomIndex atomId);
    CompoundRep& nameAtom(const Compound::AtomName& newName, const Compound::AtomPathName& oldName);

    CompoundRep& nameAtom(
        const Compound::AtomName& newName, 
        const Compound::AtomPathName& oldName, 
        BiotypeIndex biotype);

    // setBiotype("C", Biotype::MethaneC);
    CompoundRep& setBiotypeIndex(const Compound::AtomName& atomName, BiotypeIndex biotype);

    // REX
    CompoundRep& setAtomMobilizedBodyIndex(const Compound::AtomIndex& atomIndex, const MobilizedBodyIndex mbx);
        
    CompoundRep& nameBondCenter(Compound::BondCenterName newName, Compound::BondCenterPathName oldName);

    // Use atoms names as found in subcompound
    CompoundRep& inheritAtomNames(const Compound::Name& scName);
    CompoundRep& inheritBondCenterNames(const Compound::Name& scName);

    bool hasDihedral(const Compound::DihedralName& angleName) const {
        return ( AtomName_To_dihedralAngles.find(angleName) != AtomName_To_dihedralAngles.end() );
    }

    bool atomsAreBonded(const AtomInfo& atom1, const AtomInfo& atom2) const 
    {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(atom1.getIndex(), atom2.getIndex());
        return (AIxPair_To_BondIx.find(key) != AIxPair_To_BondIx.end());
    }


    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const Compound::AtomName& atom1,
        const Compound::AtomName& atom2,
        const Compound::AtomName& atom3,
        const Compound::AtomName& atom4,
        Angle nomenclatureOffset
        ) 
    {
        assert( ! hasDihedral(angleName) );
        assert( atomsAreBonded(getAtomInfo(atom1), getAtomInfo(atom2)) );
        assert( atomsAreBonded(getAtomInfo(atom2), getAtomInfo(atom3)) );
        assert( atomsAreBonded(getAtomInfo(atom3), getAtomInfo(atom4)) );

        const BondCenterInfo& bond1 = getBondCenterInfo(atom2, atom1);
        const BondCenterInfo& bond2 = getBondCenterInfo(atom3, atom4);

        defineDihedralAngle( angleName, bond1, bond2, nomenclatureOffset );

        assert( hasDihedral(angleName) );

        return *this;
    }

    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const Compound::BondCenterName& bondName1,
        const Compound::BondCenterName& bondName2,
        Angle nomenclatureOffset
        ) 
    {
        // assert( ! hasDihedral(angleName) );

        const BondCenterInfo& bond1 = getBondCenterInfo(bondName1);
        const BondCenterInfo& bond2 = getBondCenterInfo(bondName2);

        defineDihedralAngle(angleName, bond1, bond2, nomenclatureOffset);

        assert( hasDihedral(angleName) ); 

        return *this;
    }

    CompoundRep& defineDihedralAngle(
        const Compound::DihedralName& angleName,
        const BondCenterInfo& bond1,
        const BondCenterInfo& bond2,
        Angle nomenclatureOffset
        )
    {
        assert(AtomName_To_dihedralAngles.find(angleName) == AtomName_To_dihedralAngles.end());
    
        AtomName_To_dihedralAngles[angleName] = DihedralAngle(bond1.getIndex(), bond2.getIndex(), nomenclatureOffset);

        assert(AtomName_To_dihedralAngles.find(angleName) != AtomName_To_dihedralAngles.end());

        //// Define internal offset
        //DihedralAngle& dihedral = dihedralAnglesByName.find(angleName)->second;
        //const BondCenterInfo& bc21 = getBondCenterInfo(dihedral.getBondCenter1Id());
        //const BondCenterInfo& bc34 = getBondCenterInfo(dihedral.getBondCenter2Id());

        //// Find bond axis to project onto
        //const AtomInfo& atom2 = getAtomInfo(bc21.getAtomIndex());
        //const AtomInfo& atom3 = getAtomInfo(bc34.getAtomIndex());
        //const BondCenterInfo& bondBondCenter = getBondCenterInfo(atom2, atom3);

        //assert(bondBondCenter.isBonded());
        //assert(bondBondCenter.getIndex() != bc21.getIndex());
        //assert(bondBondCenter.getIndex() != bc34.getIndex());
        //assert(bc21.getIndex() != bc34.getIndex());

        //UnitVec3 xAxis(1,0,0);

        //// vector v1: from atom 1 to atom 2
        //Transform C_X_A2 = calcDefaultAtomFrameInCompoundFrame(atom2.getIndex());
        //Transform A2_X_BC21 = calcDefaultBondCenterFrameInAtomFrame(bc21);
        //Transform C_X_BC21 = C_X_A2 * A2_X_BC21;
        //UnitVec3 v1(C_X_BC21 * -xAxis); // negative x-axis because want 1->2, not 2->1 vector

        //// vector v2: from atom 2 to atom 3
        //Transform A2_X_BCB = calcDefaultBondCenterFrameInAtomFrame(bondBondCenter);
        //Transform C_X_BCB = C_X_A2 * A2_X_BCB;
        //UnitVec3 v2(C_X_BCB * xAxis);

        //// vector v3: from atom 3 to atom 4
        //Transform C_X_A3 = calcDefaultAtomFrameInCompoundFrame(atom3.getIndex());
        //Transform A3_X_BC34 = calcDefaultBondCenterFrameInAtomFrame(bc34);
        //Transform C_X_BC34 = C_X_A3 * A3_X_BC34;
        //UnitVec3 v3(C_X_BC34 * xAxis);

        //Angle nominalDihedralAngle = calcDihedralAngle(v1, v2, v3);

        //const Bond& bond = getBond(getBondInfo(bondBondCenter.getBondIndex()));
        //Angle internalDihedralAngle = bond.getDefaultDihedralAngle();

        //// internal + offset = nominal
        //Angle offset = nominalDihedralAngle - internalDihedralAngle;
        //dihedral.setInternalOffset(offset);

        return *this;
    }

    Bond& updBondByDihedral(DihedralAngle& dihedral) 
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());

        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );

        BondInfo& bondInfo = updBondInfo(atom1, atom2);
        Bond& bond = updBond(bondInfo);

        return bond;
    }

    Bond& updBondByDihedralName(const String& bondName) 
    {
        assert( AtomName_To_dihedralAngles.find(bondName) != AtomName_To_dihedralAngles.end() );
        DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(bondName)->second;
        return  updBondByDihedral(dihedral);
    }

    const Bond& getBondByDihedral(const DihedralAngle& dihedral) const 
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());

        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );

        const BondInfo& bondInfo = getBondInfo(atom1, atom2);
        const Bond& bond = getBond(bondInfo);

        return bond;
    }

    const Bond& getBondByDihedralName(const String& dihedralName) const {
        assert( AtomName_To_dihedralAngles.find(dihedralName) != AtomName_To_dihedralAngles.end() );

        const DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(dihedralName)->second;

        return getBondByDihedral(dihedral);
    }


    /**
     * \brief Sets dihedral angle without modifying bond-length or bond-angles.
     *
     * Modifying bond-angles, on the other hand, can modify those dihedral angles that involve
     * BondCenters other than the first two BondCenters on each atom.
     */
    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::AtomIndex atomIndex1, 
            Compound::AtomIndex atomIndex2, 
            Compound::AtomIndex atomIndex3, 
            Compound::AtomIndex atomIndex4)
    {
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo( getAtomInfo(atomIndex2), getAtomInfo(atomIndex1) );
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo( getAtomInfo(atomIndex3), getAtomInfo(atomIndex4) );

        // for debugging
        //String atom1Name = getAtomName(atomIndex1);
        //String atom2Name = getAtomName(atomIndex2);
        //String atom3Name = getAtomName(atomIndex3);
        //String atom4Name = getAtomName(atomIndex4);
        //std::cout << atom1Name << "->" << atom2Name << "->" << atom3Name << "->" << atom4Name << std::endl;
        //std::cout << "RECONSTRUCT STEP 1.0.1 " << offsetAngle4 << std::endl << std::flush;

        return setDefaultDihedralAngle( angle, bondCenterInfo21.getIndex(), bondCenterInfo34.getIndex() );
    }


    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::AtomName atom1, 
            Compound::AtomName atom2, 
            Compound::AtomName atom3, 
            Compound::AtomName atom4)
    {
    	return setDefaultDihedralAngle(angle, 
    			getAtomInfo(atom1).getIndex(),
    			getAtomInfo(atom2).getIndex(),
    			getAtomInfo(atom3).getIndex(),
    			getAtomInfo(atom4).getIndex() );
    }
    
    // determine difference, in radians, between dihedral defined by these bond centers (nominal),
    // and dihedral defined by "canonical" bond centers (internal).
    // nominal = internal + offset => offset = nominal - internal
    Angle calcDefaultInternalDihedralOffsetAngle(
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34) const
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        //const BondInfo& bondInfo23 = getBondInfo(atomInfo2, atomInfo3);
        //const Bond& bond23 = getBond(bondInfo23);

        // Identify the bond centers associated with the atom2-atom3 bond
        const BondCenterInfo& bondCenterInfo23 = getBondCenterInfo(atomInfo2, atomInfo3);
        const BondCenterInfo& bondCenterInfo32 = getBondCenterInfo(atomInfo3, atomInfo2);
        // sanity check those central bond centers
        assert(bondCenterInfo23.getAtomIndex() == atomIndex2);
        assert(bondCenterInfo32.getAtomIndex() == atomIndex3);

        // 1) Identify canonical bond centers for internal dihedral angle
        // Usually bond-center number zero(0), unless zero participates in the atom2-atom3 bond
        CompoundAtom::BondCenterIndex canonicalCenterIndex2(0); // default to zero
        if (bondCenterInfo23.getAtomBondCenterIndex() == 0) // unless zero is used for 2->3 bond
            canonicalCenterIndex2 = CompoundAtom::BondCenterIndex(1);

        CompoundAtom::BondCenterIndex canonicalCenterIndex3(0); // default to zero
        if (bondCenterInfo32.getAtomBondCenterIndex() == 0) // unless zero is used for 2->3 bond
            canonicalCenterIndex3 = CompoundAtom::BondCenterIndex(1);

        // debug
        // Compound::AtomName n2 = getAtomName(atomIndex2);
        // Compound::AtomName n3 = getAtomName(atomIndex3);

        // 2) Compute offsets for actual bond centers
        // * offsetAngle1 is counter-clockwise angle from canonical bond center on atom2 to atom1, viewed
        // down the atom3-atom2 axis.
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo(bondCenterIndex21);
        Angle offsetAngle1;
        if (canonicalCenterIndex2 == bondCenterInfo21.getAtomBondCenterIndex())
            offsetAngle1 = 0.0;
        else
        {

            // trick the bond-vector version of calcDihedralAngle into giving the offset angle at the atom
            const CompoundAtom& atom2 = getAtom(atomIndex2);
            UnitVec3 dirAtom1    = -atom2.getBondCenterDirectionInAtomFrame(bondCenterInfo21.getAtomBondCenterIndex());
            UnitVec3 dirBond     = atom2.getBondCenterDirectionInAtomFrame(bondCenterInfo23.getAtomBondCenterIndex());
            UnitVec3 dirRefAtom1 = atom2.getBondCenterDirectionInAtomFrame(canonicalCenterIndex2);

            // Sometimes bond direction is colinear with atom direction, if chirality is hosed
            double problemCheck = std::abs(dot(dirBond, dirAtom1));
            if (problemCheck > 0.999)
                offsetAngle1 = 0.0;
            else
                offsetAngle1 = SimTK::calcDihedralAngle(dirRefAtom1, dirBond, dirAtom1);

            // assert(offsetAngle1 != 0);
        }

        // * offsetAngle4 is counter-clockwise angle from canonical bond center on atom3 to atom4, viewed
        // down the atom3-atom2 axis.
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo(bondCenterIndex34);
        Angle offsetAngle4 = std::numeric_limits<Angle>::max(); // TODO might use std::optionatl, should look into it
        if (canonicalCenterIndex3 == bondCenterInfo34.getAtomBondCenterIndex())
            offsetAngle4 = 0.0;
        else
        {

            // trick the bond-vector version of calcDihedralAngle into giving the offset angle at the atom
            const CompoundAtom& atom3 = getAtom(atomIndex3);
            UnitVec3 dirAtom4    = -atom3.getBondCenterDirectionInAtomFrame(bondCenterInfo34.getAtomBondCenterIndex());
            UnitVec3 dirBond     = -atom3.getBondCenterDirectionInAtomFrame(bondCenterInfo32.getAtomBondCenterIndex());
            UnitVec3 dirRefAtom4 = atom3.getBondCenterDirectionInAtomFrame(canonicalCenterIndex3);

            // Sometimes bond direction is colinear with atom direction, if chirality is hosed
            double problemCheck = std::abs(dot(dirBond, dirAtom4));
            if (problemCheck > 0.999)
                offsetAngle1 = 0.0;
            else
                offsetAngle4 = SimTK::calcDihedralAngle(dirRefAtom4, dirBond, dirAtom4);

            // assert(offsetAngle4 != 0);
        }

        if(offsetAngle4 == std::numeric_limits<Angle>::max()) {
            // Should never get here, but compiler keeps warning.
            // Se above for a more elegant solution.
            assert(false);
        }

        // nominal = internal + offset
        // offset = nominal - internal
        // internal = nominal - offset
        // Angle internalDihedralAngle = angle + offsetAngle1 - offsetAngle4;

        Angle offsetAngle = offsetAngle4 - offsetAngle1;

        if(offsetAngle4 != std::numeric_limits<Angle>::max()) {
            while ( -SimTK::Pi >= offsetAngle ) offsetAngle += 2 * SimTK::Pi;
            while ( SimTK::Pi < offsetAngle ) offsetAngle -= 2 * SimTK::Pi;
        }
        //std::cout << "RECONSTRUCT STEP 1.0.3 " << offsetAngle4 << std::endl << std::flush;

        // debugging
        //std::cout << "  total offset = " << offsetAngle * DuMM::Rad2Deg;
        //std::cout << "; offset1 = " << offsetAngle1 * DuMM::Rad2Deg;
        //std::cout << "; offset4 = " << offsetAngle4 * DuMM::Rad2Deg << std::endl;

        return offsetAngle;
    }


    /**
     * \brief Sets dihedral angle without modifying bond-length or bond-angles.
     *
     * Modifying bond-angles, on the other hand, can modify those dihedral angles that involve
     * BondCenters other than the first two BondCenters on each atom.
     */
    CompoundRep& setDefaultDihedralAngle( 
            Angle angle, 
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34)
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        BondInfo& bondInfo23 = updBondInfo(atomInfo2, atomInfo3);
        Bond& bond23 = updBond(bondInfo23);

        Angle internalDihedralOffsetAngle = calcDefaultInternalDihedralOffsetAngle(bondCenterIndex21, bondCenterIndex34);

        Angle internalDihedralAngle = angle - internalDihedralOffsetAngle;

        if(internalDihedralOffsetAngle != std::numeric_limits<Angle>::max()) {
            while ( -SimTK::Pi >= internalDihedralAngle ) internalDihedralAngle += 2 * SimTK::Pi;
            while ( SimTK::Pi < internalDihedralAngle ) internalDihedralAngle -= 2 * SimTK::Pi;
        }else{
            internalDihedralAngle = std::numeric_limits<Angle>::max();
        }
        //std::cout << "RECONSTRUCT STEP 1.0.2 " << offsetAngle4 << std::endl << std::flush;

        //std::cout << "old internal angle = " << bond23.getDefaultDihedralAngle() * DuMM::Rad2Deg << std::endl;
        //std::cout << "new internal angle = " << internalDihedralAngle * DuMM::Rad2Deg << std::endl;

		// debug - notice when angle changes
		//Real diff = internalDihedralAngle - bond23.getDefaultDihedral();
  //      while ( -SimTK::Pi >= diff ) diff += 2 * SimTK::Pi;
  //      while ( SimTK::Pi < diff ) diff -= 2 * SimTK::Pi;
		//diff = diff < 0 ? -diff : diff;
		//if (diff > 0.005) 
		//{
		//	int x = 1;
		//}

        bond23.setDefaultDihedralAngle(internalDihedralAngle);

        return *this;
    }

    // setDefaultDihedral changes no bond lengths or bond angles
    CompoundRep& setDefaultDihedralAngle(const String& dihedralName, Angle finalNominalAngle) 
    {
        //Bond& bond = updBondByDihedralName(dihedralName);
        DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(dihedralName)->second;

        // internal = nominal - offset
        Angle angle = finalNominalAngle - dihedral.getNomenclatureOffset();

        setDefaultDihedralAngle(angle, dihedral.getBondCenter1Id(), dihedral.getBondCenter2Id());
        // Angle internalAngle = finalNominalAngle - dihedral.getInternalOffset() - dihedral.getNomenclatureOffset();

        // bond.setDefaultDihedralAngle(internalAngle);

        return *this;
    }

// EU BEGIN
    Angle bgetDefaultDihedralAngle(Compound::BondIndex bondIx) const 
    {
        const BondInfo& bondInfo = getBondInfo(bondIx);
        const Bond& bond = getBond(bondInfo);
        Angle angle = -111111;
        if (bond.isRingClosingBond()){ // ring closing bonds cannot be part of tree structure
        return angle;
        }
        if (bond.getMobility() == BondMobility::Free){
        return angle;
        }
        else if (bond.getMobility() == BondMobility::Rigid){
        return angle;
        }
        else if (bond.getMobility() == BondMobility::Torsion){
        angle = bond.getDefaultDihedralAngle();
        }
        return angle;
    }


    Angle bgetDefaultInboardDihedralAngle(Compound::AtomIndex atomIx) const 
    {
        // Get atom
        const CompoundAtom& atom = getAtom(atomIx);

        // Get inboard bond index (in Compound not in Atom)
        CompoundAtom::BondCenterIndex inboardBondCenterIx = atom.getInboardBondCenterIndex();
        const BondCenterInfo& inboardBondCenterInfo = getBondCenterInfo(atomIx, inboardBondCenterIx);
        Compound::BondIndex inboardBondIndex = inboardBondCenterInfo.getBondIndex();
        //const BondInfo& inboardBondInfo = getBondInfo((getBondCenterInfo(atomIx, (atom.getInboardBondCenterIndex()))).getBondIndex());

        // Get the inboard bond
        //const BondInfo& inboardBondInfo = getBondInfo(inboardBondIndex);
        //const Bond& inboardBond = getBond(inboardBondInfo);
        return bgetDefaultDihedralAngle(inboardBondIndex);
    }

    const Transform& getFrameInMobilizedBodyFrame(Compound::AtomIndex atomIx) const
    {
        const CompoundAtom& atom = getAtom(atomIx);
        return atom.getFrameInMobilizedBodyFrame();
    }

    const Transform& bgetLocalTransform(Compound::AtomIndex atomIx) const
    {
        const CompoundAtom& atom = getAtom(atomIx);
        return atom.getDefaultFrameInCompoundFrame();
    }


    /*!
    * <!-- Print Vec3 -->
    */
    void PrintTransform(SimTK::Transform T, int decimal_places,
        std::string header = "", std::string rowPrefix = "")
    {
        std::cout << header << std::endl;
        const SimTK::Mat44 M = T.toMat44();

        for(int i = 0; i < 4; i++){
            std::cout << rowPrefix;
            for(int k = 0; k < 4; k++){
                std::cout
                    << std::setw(6 + decimal_places) << std::fixed
                    << std::setprecision(decimal_places)			
                    << M(i, k) << " ";
            }
            std::cout << std::endl;
        }
    }

    /*!
    * <!-- Print Transform -->
    */
    void PrintVec3(SimTK::Vec3 vec, int decimal_places,
        std::string header = "", std::string rowPrefix = "")
    {
        std::cout << header << std::endl;

        for(int i = 0; i < 3; i++){
            std::cout << rowPrefix
                << std::setw(6 + decimal_places) << std::fixed
                << std::setprecision(decimal_places)			
                << vec(i) << " ";
            std::cout << std::endl;
        }

    }


    /*!
    * <!-- Print Compound geometry (which is the most detailed) -->
    */
    CompoundRep& PrintCompoundGeometry(const Compound::AtomTargetLocations& atomTargets){

        // Iterate atoms
        std::vector< AtomIndexVector > atomRun = getBondedAtomRuns(1, atomTargets);
        std::cout << "CompoundRep::PrintCompoundGeometry atomTargets\n";
        for(const auto& atomRIx : atomRun) {
            const Compound::AtomIndex atomIx = atomRIx[0];
                
                SimTK::Vec3 loc = atomTargets.at(atomIx);

                std::cout << " cAIx " << atomIx
                    << " loc " << loc[0] <<" " << loc[1] <<" " << loc[2];

                std::cout << std::endl;

        }

        for(Compound::AtomIndex atomIx(0); atomIx < getNumAtoms(); atomIx++){    
            CompoundAtom& atom = updAtom(atomIx);
            const AtomInfo& atomInfo = getAtomInfo(atomIx);
           
            // Go through bond centers on atom.
            for (CompoundAtom::BondCenterIndex BCIx(0); BCIx < atom.getNumBonds(); ++BCIx) {
                BondCenter &BC = atom.updBondCenter(CompoundAtom::BondCenterIndex(BCIx));
                BondCenter& bondCenter = updBondCenter(Compound::BondCenterIndex(BCIx));
                SimTK::UnitVec3 dir = atom.getBondCenterDirectionInAtomFrame(BCIx);
                std::cout << "CompoundRep::PrintCompoundGeometry"
                    << " cAIx " << atomIx
                    << " BCIx " << BCIx
                    << " dir " << dir[0] << " " << dir[1] << " " << dir[2]
                    << " chirality " << BC.getChirality();

                std::cout << std::endl;

            }
        }

        return *this;

    }

    /*!
    <!-- Set atom frame in mobod frame -->
    */
    CompoundRep& bsetFrameInMobilizedBodyFrame(Compound::AtomIndex atomIx, Transform B_X_atom)
    {
        CompoundAtom& atom = updAtom(atomIx);
        atom.setFrameInMobilizedBodyFrame(B_X_atom);
        return *this;
    }

    /*!
    <!-- WIP Get the inboard atom index of a given atom implementation -->
    */
    Compound::AtomIndex getInboardAtomIndex(Compound::AtomIndex& atomIx) const
    {
        const CompoundAtom& atom = getAtom(atomIx);
        const CompoundAtom::BondCenterIndex inboardBondCenterIx =
            atom.getInboardBondCenterIndex();
        const BondCenterInfo& inboardBondCenterInfo =
            getBondCenterInfo(atomIx, inboardBondCenterIx);
        const Compound::BondIndex inboardBondIndex =
            inboardBondCenterInfo.getBondIndex();

        // Get the inboard bond
        // const BondInfo& inboardBondInfo = getBondInfo(inboardBondIndex);
        // const Bond& inboardBond = getBond(inboardBondInfo);
        // const Compound::BondIndex inboardBondIx = inboardBondInfo.getIndex();
        // const Compound::BondCenterIndex parentBCIx = 
        //     inboardBondInfo.getParentBondCenterIndex();
        
        // Aparently 0 is parent and 1 is child
        int pbc = 0;
        const Compound::AtomIndex inboardAIx = getBondAtomIndex(inboardBondIndex, pbc);
        return inboardAIx;

    }


// EU END

///* GMolModel Try other Mobilizers
  mdunits::Length bgetDefaultInboardBondLength(Compound::AtomIndex atomIx) const 
  {
      // Get atom
      const CompoundAtom& atom = getAtom(atomIx);

      // Get inboard bond index (in Compound not in Atom)
      CompoundAtom::BondCenterIndex inboardBondCenterIx = atom.getInboardBondCenterIndex();
      const BondCenterInfo& inboardBondCenterInfo = getBondCenterInfo(atomIx, inboardBondCenterIx);
      Compound::BondIndex inboardBondIndex = inboardBondCenterInfo.getBondIndex();

      // Get the inboard bond
      const BondInfo& inboardBondInfo = getBondInfo(inboardBondIndex);
      const Bond& inboardBond = getBond(inboardBondInfo);
      return inboardBond.getDefaultBondLength();

  }
// GMolmodel END */

    Angle calcDefaultDihedralAngle(const String& dihedralName) const 
    {
        assert( AtomName_To_dihedralAngles.find(dihedralName) != AtomName_To_dihedralAngles.end() );

        const DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(dihedralName)->second;

        return calcDefaultDihedralAngle(dihedral);
    }

    Angle calcDefaultDihedralAngle(            
            Compound::BondCenterIndex bondCenterIndex21, 
            Compound::BondCenterIndex bondCenterIndex34)
    {
        Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondCenterIndex21).getAtomIndex();
        Compound::AtomIndex atomIndex3 = getBondCenterInfo(bondCenterIndex34).getAtomIndex();

        const AtomInfo& atomInfo2 = getAtomInfo(atomIndex2);
        const AtomInfo& atomInfo3 = getAtomInfo(atomIndex3);

        // Sanity check topology
        assert( atomsAreBonded(atomInfo2, atomInfo3) ); // absolutely required

        // Find central bond
        const BondInfo& bondInfo23 = getBondInfo(atomInfo2, atomInfo3);
        const Bond& bond23 = getBond(bondInfo23);

        Angle internalDihedralAngle = bond23.getDefaultDihedralAngle();

        Angle nominalDihedralAngle = internalDihedralAngle + 
            calcDefaultInternalDihedralOffsetAngle(bondCenterIndex21, bondCenterIndex34);

        return nominalDihedralAngle;
    }

    Angle calcDefaultDihedralAngle( 
            Compound::AtomIndex atomIndex1, 
            Compound::AtomIndex atomIndex2, 
            Compound::AtomIndex atomIndex3, 
            Compound::AtomIndex atomIndex4)
    {
        // This belongs to atom2
        const BondCenterInfo& bondCenterInfo21 = getBondCenterInfo( getAtomInfo(atomIndex2), getAtomInfo(atomIndex1) );

        // This belong to atom3
        const BondCenterInfo& bondCenterInfo34 = getBondCenterInfo( getAtomInfo(atomIndex3), getAtomInfo(atomIndex4) );

        return calcDefaultDihedralAngle( bondCenterInfo21.getIndex(), bondCenterInfo34.getIndex() );
    }

    CompoundRep& setDihedralAngle(State& state, const String& dihedralName, Angle angleInRadians) 
    {
        assert(ownerSystem != NULL);
        Bond& bond = updBondByDihedralName(dihedralName);
        DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(dihedralName)->second;

        // case1 : Pin dihedral
        if (bond.getPinJointId().isValid()) {
            SimbodyMatterSubsystem &matter = ownerSystem->updMatterSubsystem();
            if(bond.getMobility() == BondMobility::Torsion) {
                MobilizedBody::Pin &body = (MobilizedBody::Pin &) matter.updMobilizedBody(bond.getPinJointId());
                // TODO - create calcDihedralOffset(State&...) method and use it here, instead of default
                Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(dihedral.getBondCenter1Id(),
                                                                              dihedral.getBondCenter2Id());
                // nominal = internal + offset
                Angle internalAngle = angleInRadians - internalOffset - dihedral.getNomenclatureOffset();
                body.setAngle(state, internalAngle);
            }else if(bond.getMobility() == BondMobility::BallF) { // Gmol
                MobilizedBody::Ball &ball = (MobilizedBody::Ball &) matter.updMobilizedBody(bond.getPinJointId());
                // TODO - create calcDihedralOffset(State&...) method and use it here, instead of default
                Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(dihedral.getBondCenter1Id(),
                                                                              dihedral.getBondCenter2Id());
                // nominal = internal + offset
                Angle internalAngle = angleInRadians - internalOffset - dihedral.getNomenclatureOffset();
                ball.setQ(state, SimTK::Rotation(internalAngle,
                        CoordinateAxis::ZCoordinateAxis()).convertRotationToQuaternion().asVec4());
            }
        }

        else  // TODO
        {
            assert(false);

            // Dihedral may be offset from "standard" dihedral for bond
            assert(ownerSystem);
            SimbodyMatterSubsystem& matter = ownerSystem->updMatterSubsystem();

            Angle previousInternalDihedral = bond.getDihedralAngle(state, matter);

            Angle previousNominalDihedral = calcDihedralAngle(state, dihedralName); // requires realizePosition
            // Nominal = internal + offset
            Angle offsetAngle = previousNominalDihedral - previousInternalDihedral;
            Angle newInternalDihedral = angleInRadians - offsetAngle;
            // Restrict to range +-Pi
            while (newInternalDihedral <= -Pi) newInternalDihedral += 2*Pi;
            while (newInternalDihedral > Pi) newInternalDihedral -= 2*Pi;
            bond.setDihedralAngle(state, matter, newInternalDihedral); // clears realizePosition

            Angle testAngle = calcDihedralAngle(state, dihedralName); // requires realizePosition
            Angle error = calcDihedralAngle(state, dihedralName) - testAngle;
            while (error <= -Pi) error += 2*Pi;
            while (error > Pi) error -= 2*Pi;
            error *= error;
            assert(error < 1e-6);
        }

        return *this;
    }

    CompoundRep& setDefaultBond1Angle(const String& bondName, Angle angle) {
        updBondCenter(bondName).setDefaultBond1Angle(angle);
        return *this;
    }
    CompoundRep& setDefaultBond2Angle(const String& bondName, Angle angle) {
        updBondCenter(bondName).setDefaultBond2Angle(angle);
        return *this;
    }

    CompoundRep& setDefaultBondAngle(
        Angle angle, 
        const Compound::AtomName& atom1Name, 
        const Compound::AtomName& atom2Name, 
        const Compound::AtomName& atom3Name) 
    {
        const Compound::AtomIndex atom1Id   = getAtomInfo(atom1Name).getIndex();
        const Compound::AtomIndex atom2Id   = getAtomInfo(atom2Name).getIndex();
        const Compound::AtomIndex atom3Id   = getAtomInfo(atom3Name).getIndex();
        
        return setDefaultBondAngle(angle, atom1Id, atom2Id, atom3Id);
    }
    
    /*! <!-- Set angles to BC0 and BC1
     * Version that takes IDs, to reduce string lookups
     * 
    --> */    
    CompoundRep& setDefaultBondAngle(
        Angle angle, 
        const Compound::AtomIndex cAIx_atom1, 
        const Compound::AtomIndex cAIx_atom2, 
        const Compound::AtomIndex cAIx_atom3) 
    {

        CompoundAtom& atom2 = updAtom(cAIx_atom2);
        const AtomInfo& atom2Info = getAtomInfo(cAIx_atom2);

        // go through bond centers on atom2
        // and establish which one is linked to atom1 and which to atom3
        CompoundAtom::BondCenterIndex atom2_to_atom1_BCIx;
        CompoundAtom::BondCenterIndex atom2_to_atom3_BCIx;
        
        for (CompoundAtom::BondCenterIndex atom2_BCIx(0); atom2_BCIx < atom2.getNumBonds(); ++atom2_BCIx) {
            // Get BC info
            const BondCenterInfo& atom2BondCenterInfo =
                getBondCenterInfo(atom2Info.getIndex(), atom2_BCIx);

            if (atom2BondCenterInfo.isBonded()) {
            
                const BondCenterInfo& partnerBondCenterInfo = 
                    getBondCenterInfo(atom2BondCenterInfo.getBondPartnerBondCenterIndex());
            
                if (partnerBondCenterInfo.getAtomIndex() == cAIx_atom1)
                    {atom2_to_atom1_BCIx = atom2_BCIx;}
                else if (partnerBondCenterInfo.getAtomIndex() == cAIx_atom3)
                    {atom2_to_atom3_BCIx = atom2_BCIx;}
            }
        }

        assert(atom2_to_atom1_BCIx.isValid());
        assert(atom2_to_atom3_BCIx.isValid());
        assert(atom2_to_atom1_BCIx != atom2_to_atom3_BCIx);

        // Get order relationship between atom1- and atom3- BC indexes
        CompoundAtom::BondCenterIndex largerId, smallerId;
        if (atom2_to_atom1_BCIx > atom2_to_atom3_BCIx) {
            largerId = atom2_to_atom1_BCIx;
            smallerId = atom2_to_atom3_BCIx;
        } else {
            largerId = atom2_to_atom3_BCIx;
            smallerId = atom2_to_atom1_BCIx;
        }

        // The smallest BCIx dictates the bond angles inside BCs 
        // one of the bond centers must be bond1 or bond2
        // assert(smallerId < 2);

        if (smallerId == 0) { // BC0 - atom2 - (BC1, BC2, BC3)
            atom2.updBondCenter(largerId).setDefaultBond1Angle(angle);
        }
        else if (smallerId == 1) { // BC1 - atom2 - (BC2, BC3)
            atom2.updBondCenter(largerId).setDefaultBond2Angle(angle);
        }
        else {
            // std::cout << "[WARNING]: setDefaultBondAngle " << atom1Id <<" "
            //     << atom2Id <<" " << atom3Id <<" "
            //     << "smallerId " << smallerId <<" "
            //     << "largerId " << largerId <<" "
            //     << std::endl;
            // assert(false); // TODO
        }

        return *this;
    }

    CompoundRep& setDefaultBondLength(mdunits::Length length, const Compound::AtomName& atom1Name, const Compound::AtomName& atom2Name) 
    {
        const CompoundAtom&             atom2       = getAtom(atom2Name);
        const AtomInfo&         atom2Info   = getAtomInfo(atom2Name);
        const Compound::AtomIndex  atom1Id     = getAtomInfo(atom1Name).getIndex();

        // go through bond centers on atom2
        CompoundAtom::BondCenterIndex center1Id;
        for (CompoundAtom::BondCenterIndex b(0); b < atom2.getNumBonds(); ++b) {
            const BondCenterInfo& bondCenterInfo = getBondCenterInfo(atom2Info.getIndex(), b);
            if (getBondCenter(bondCenterInfo).isBonded()) {
                const BondCenterInfo& partnerBondCenterInfo = getBondCenterInfo(bondCenterInfo.getBondPartnerBondCenterIndex());
                if (partnerBondCenterInfo.getAtomIndex() == atom1Id) 
                {
                    center1Id = b;

                    // set length from atom1 direction
                    updBondCenter(bondCenterInfo).setDefaultBondLength(length);

                    // for good measure, set length from atom2 direction
                    updBondCenter(partnerBondCenterInfo).setDefaultBondLength(length);

                    break;
                }
            }
        }

        assert(center1Id.isValid());

        return *this;
    }

    /*!
    * <!-- Helper for calcDefaultAtomFramesInCompoundFrame. It sets a NaN flag for
    Top to inboard bond center transforms passed. -->
    */
    void invalidateAtomFrameCache(
        std::vector<Transform>& atomFrameCache,
        int numAtoms) const;

    /*!
     * <!-- Calculate Top to inboard bond center transform for all atoms.
     * invalidateAtomFrameCache(atomFrameCache) must be called before this -->
     */
    // More efficient getting of all atoms at once
    // Compound::AtomTargetLocations calcDefaultAtomLocationsInCompoundFrame1() const;
    void calcDefaultAtomFramesInCompoundFrame(std::vector<Transform>& atomFrameCache) const;
    
    // Compute atom location in local compound frame
    Transform calcDefaultAtomFrameInCompoundFrame(const Compound::AtomName& name) const;
    Transform calcDefaultAtomFrameInGroundFrame(const Compound::AtomName& name) const;
    Vec3 calcDefaultAtomLocationInCompoundFrame(const Compound::AtomName& name) const;
    Vec3 calcDefaultAtomLocationInGroundFrame(const Compound::AtomName& name) const;

    MobilizedBodyIndex getAtomMobilizedBodyIndex(Compound::AtomIndex atomId) const 
    {
        const CompoundAtom& atom = getAtom(atomId);
        return atom.getMobilizedBodyIndex();
    }
    Vec3 getAtomLocationInMobilizedBodyFrame(Compound::AtomIndex atomId) const {
        const CompoundAtom& atom = getAtom(atomId);
        return atom.getLocationInMobilizedBodyFrame();
    }
    Vec3 calcAtomLocationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Position);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.getBodyTransform(state)*loc;
    }
    Vec3 calcAtomVelocityInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Velocity);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.findStationVelocityInGround(state, loc);
    }
    Vec3 calcAtomAccelerationInGroundFrame(const State& state, Compound::AtomIndex atomId) const {
        ownerSystem->realize(state, Stage::Acceleration);
        const CompoundAtom& atom = getAtom(atomId);
        Vec3 loc = atom.getLocationInMobilizedBodyFrame();
        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomId));
        return body.findStationAccelerationInGround(state, loc);
    }

    
    Transform calcDefaultBondCenterFrameInCompoundFrame(const BondCenterInfo& info) const;

    Transform calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId) const;

    /*! <!-- for O(n) version of all atom Frame computation
    * Version with caching for O(n) performance -->
    */ 
    const Transform& calcDefaultAtomFrameInCompoundFrame(Compound::AtomIndex atomId, std::vector<Transform>& atomFrameCache) const;

    Transform calcDefaultAtomFrameInGroundFrame(Compound::AtomIndex atomId) const;

    typedef std::vector<Compound::AtomIndex> AtomIndexVector;


    /*! <!-- get list of all runs of consecutive bonded atoms of run-length n from the atoms mentions in an AtomTargetLocations structure 
    * for example, to get a list of all bonded pairs, set run-length to 2.-->
    */ 
    std::vector< AtomIndexVector > getBondedAtomRuns(int atomRunCount, const Compound::AtomTargetLocations& atomTargets) const 
    {
        std::vector< AtomIndexVector >  answer;

        typedef std::map<Compound::AtomIndex, AtomIndexVector > BondMap;
        BondMap bondMap;

        // 1) Hash bonding data
        for (Compound::BondIndex bondCnt(0); bondCnt < getNumBonds(); ++bondCnt) 
        {
            const BondInfo& bondInfo = getBondInfo(bondCnt);

            // ignore bonds without known atom positions at both ends
            // i.e. keep the previous default bond lengths for those
            Compound::AtomIndex atomIndex1 = getBondCenterInfo(bondInfo.getParentBondCenterIndex()).getAtomIndex();
            Compound::AtomIndex atomIndex2 = getBondCenterInfo(bondInfo.getChildBondCenterIndex()).getAtomIndex();
            if (atomTargets.find(atomIndex1) == atomTargets.end()) continue;
            if (atomTargets.find(atomIndex2) == atomTargets.end()) continue;

            assert(atomIndex1 != atomIndex2);

            if (bondMap.find(atomIndex1) == bondMap.end()) bondMap[atomIndex1] = AtomIndexVector();
            bondMap[atomIndex1].push_back(atomIndex2);

            if (bondMap.find(atomIndex2) == bondMap.end()) bondMap[atomIndex2] = AtomIndexVector();
            bondMap[atomIndex2].push_back(atomIndex1);

            // Update n==2 version of answer
            answer.push_back(AtomIndexVector());
            answer.back().push_back(atomIndex1);
            answer.back().push_back(atomIndex2);

            // and in reverse order
            answer.push_back(AtomIndexVector());
            answer.back().push_back(atomIndex2);
            answer.back().push_back(atomIndex1);
        }

        // 2) Create bonded atom run lists
        for (int n = 3; n <= atomRunCount; ++n)
        {
            std::vector< AtomIndexVector >  newAnswer;
            std::vector< AtomIndexVector >::const_iterator oldRunIt;
            for (oldRunIt = answer.begin(); oldRunIt != answer.end(); ++oldRunIt)
            {
                // remember which atoms are already in this run
                std::set<Compound::AtomIndex> oldAtoms;
                AtomIndexVector::const_iterator oldAtomIt;
                for (oldAtomIt = oldRunIt->begin(); oldAtomIt != oldRunIt->end(); ++oldAtomIt)
                    oldAtoms.insert(*oldAtomIt);

                // look for new atoms bonded to end of old run
                Compound::AtomIndex oldTailIndex = oldRunIt->back(); // final atom of shorter run
                const AtomIndexVector& newTailCandidates = bondMap[oldTailIndex];
                AtomIndexVector::const_iterator newTailI;
                for (newTailI = newTailCandidates.begin(); newTailI != newTailCandidates.end(); ++newTailI)
                {
                    if (oldAtoms.find(*newTailI) != oldAtoms.end()) continue; // ignore all atoms already in this run

                    // first place a copy of the old run in the new set
                    newAnswer.push_back(*oldRunIt);
                    // then add the new atom to the end of the run
                    newAnswer.back().push_back(*newTailI);
                }
            }

            // overwrite the old set of shorter runs with the new one containing longer runs.
            answer = newAnswer;
        }

        return answer;
    }


    /*! <!-- From a vector of pairs (vector of two AtomIndexes) - atomPartners_cAIxs
     * build a map - atomPartners_Map - from an index to a set of partners.
     * similar to bSpecificAtom.neighborsIndex built by InternalCoordinates class --> */ 
    std::map<Compound::AtomIndex, std::set<Compound::AtomIndex>>
    buildAtomNeighbors_Map(const Compound::AtomTargetLocations& atomTargets){

        // Get list of consecutive bonded atoms from the atoms mentioned in an AtomTargetLocations
        std::vector<AtomIndexVector> atomPairs_cAIxs = getBondedAtomRuns(2, atomTargets);

        // Allocate map
        std::map<Compound::AtomIndex, std::set<Compound::AtomIndex>> atomNeighbors_Map;
        Compound::AtomTargetLocations::const_iterator atomTargetLocIt;
        for (atomTargetLocIt = atomTargets.begin(); atomTargetLocIt != atomTargets.end(); ++atomTargetLocIt){
            atomNeighbors_Map[atomTargetLocIt->first] = std::set<Compound::AtomIndex>();
        }

        // Insert bonds in map
        std::vector< AtomIndexVector >::const_iterator atomPairIt;
        for (atomPairIt = atomPairs_cAIxs.begin(); atomPairIt != atomPairs_cAIxs.end(); ++atomPairIt)
        {
            Compound::AtomIndex atomIndex1 = (*atomPairIt)[0];
            Compound::AtomIndex atomIndex2 = (*atomPairIt)[1];
            atomNeighbors_Map[atomIndex1].insert(atomIndex2);
            atomNeighbors_Map[atomIndex2].insert(atomIndex1);
        }

        return atomNeighbors_Map;
    }

    /*! <!-- Set BondCenters chirality by comparing with AtomTargetLocations
     * (A) Get bonds actual and target vectors in atom frame
     * (B) Permutations: get reference bond centers
     * (C) Chirality: s source, t target (old vs new)
     * If dot(cross(s1, s2),s3) and dot(cross(t1, t2),t3) have different
     * signs, then switch chirality
     *  -->
     */
    CompoundRep& matchDefaultAtomChirality(
        const Compound::AtomTargetLocations& atomTargets,
        Angle breakPlanarityThreshold,
        bool flipAll=true)
    {

        // Build a map from an AtomIndex to a set of bonded AtomIndexes
        std::map<Compound::AtomIndex, std::set<Compound::AtomIndex>>
        atomNeighbors_Map = buildAtomNeighbors_Map(atomTargets); 

        // Main loop over atoms: Check the chirality of each atom
        Compound::AtomTargetLocations::const_iterator atomTargetLocIt;
        for (atomTargetLocIt = atomTargets.begin(); atomTargetLocIt != atomTargets.end(); ++atomTargetLocIt){

            // Get atom (cAIx and CompoundAtom) and it's neighbors (atomNeighborsIxs)
            Compound::AtomIndex atomIndex = atomTargetLocIt->first;
            const SimTK::AtomInfo& atomInfo = getAtomInfo(atomIndex);
            CompoundAtom& atom = updAtom(atomIndex);

            const std::set<Compound::AtomIndex>& atomNeighborsIxs = atomNeighbors_Map[atomIndex];
            int numberOfBonds = atomNeighborsIxs.size();

            // No chirality for less than 3 partner atoms
            if (numberOfBonds < 3) {
                continue;
            }

            // Get atoms location
            Vec3 targetAtomLocation = atomTargetLocIt->second;

            // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // (A) Get source (actual) and target bond vectors
            // ----------------------------------------------------------------
            std::vector< std::pair<UnitVec3, UnitVec3> > bondVectors;
            std::vector< Compound::BondCenterIndex > compoundBCIxes;
            std::set<Compound::AtomIndex>::const_iterator neighborIt;

            for (neighborIt = atomNeighborsIxs.begin(); neighborIt != atomNeighborsIxs.end(); ++neighborIt) {

                // Get source bond direction
                Compound::AtomIndex neighborAtomIndex = atomTargets.find(*neighborIt)->first;                   // get neighbor cAIx
                const SimTK::AtomInfo& neighborAtomInfo = getAtomInfo(neighborAtomIndex);                       // get atomInfo
                const BondCenterInfo& neighborBCInfo = getBondCenterInfo(atomInfo, neighborAtomInfo);           // get BCInfo
                CompoundAtom::BondCenterIndex neighborBCIx_inAtom = neighborBCInfo.getAtomBondCenterIndex();    // get BCIx in atom !
                Transform neighborBCFrame = atom.calcDefaultBondCenterFrameInAtomFrame(neighborBCIx_inAtom);    // get BCFrame
                UnitVec3 sourceDirection(neighborBCFrame * UnitVec3(1, 0, 0));

                // Get target bond direction
                Vec3 neighborAtomLocation = atomTargets.find(*neighborIt)->second;
                UnitVec3 targetDirection(neighborAtomLocation - targetAtomLocation);

                // Store Compound bond center index
                Compound::BondCenterIndex neighborBCIx_inCompound = neighborBCInfo.getIndex();                  //  get BCIx in Compound
                compoundBCIxes.push_back( neighborBCIx_inCompound );

                // Store source and target bond directions (Horea)
                bondVectors.push_back(std::pair<UnitVec3, UnitVec3>(sourceDirection, targetDirection));

            } // every neighbor
            assert(numberOfBonds == bondVectors.size());

            // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // (B) Permutations: get reference bond centers: 0 and 1
            // ----------------------------------------------------------------
            // Identify a mapping between internal atom BondCenter indices and the recently constructed
            // neighbors atom indices
            // Because we need to distinguish left-handed from right-handed geometry in the atom frame,
            // we should use bondcenters number 0 and 1 to define the plane, so that the target structure
            // geometry matches that of the internal atom geometry
            int zeroBondCenterIndex = 0;
            int oneBondCenterIndex = 1;
            int twoBondCenterIndex = 2;

            for (int bondCnt = 0; bondCnt < (int) bondVectors.size(); ++bondCnt){

                // Get bond center index
                const BondCenterInfo& bondCenterInfo = getBondCenterInfo(compoundBCIxes[bondCnt]);
                const CompoundAtom::BondCenterIndex neighborBCIx_inAtom = bondCenterInfo.getAtomBondCenterIndex();

                if (neighborBCIx_inAtom == 0) {

                    // (one OR two) = zero
                    // zero = bondCnt
                    if ( oneBondCenterIndex == bondCnt ){
                        oneBondCenterIndex = zeroBondCenterIndex;
                    }
                    else if ( twoBondCenterIndex == bondCnt ){
                        twoBondCenterIndex = zeroBondCenterIndex;
                    }

                    zeroBondCenterIndex = bondCnt;
                
                }else if (neighborBCIx_inAtom == 1){

                    // (zero OR one) = one
                    // one = bondCnt
                    if ( zeroBondCenterIndex == bondCnt ){
                        zeroBondCenterIndex = oneBondCenterIndex;
                    }
                    else if ( twoBondCenterIndex == bondCnt ){
                        twoBondCenterIndex = oneBondCenterIndex;
                    }

                    oneBondCenterIndex = bondCnt;
                }

            } // every bondVector

            #ifdef DEBUG_MOLMODEL
                std::cout << "bondVectorIx";
                for (int bondCnt = 0; bondCnt < (int) bondVectors.size(); ++bondCnt){
                    std::cout <<" "<< getBondCenterInfo(bondCenterIndexes[bondCnt]).getAtomBondCenterIndex();
                }
                std::cout << std::endl;
                std::cout << "zeroOneTwo"
                    <<" "<< zeroBondCenterIndex <<" "<< oneBondCenterIndex <<" "<< twoBondCenterIndex
                    << std::endl;
            #endif

            // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // (C) Chirality: s source, t target (old vs new)
            // If dot(cross(s1, s2),s3) and dot(cross(t1, t2),t3) have different
            // signs, then switch chirality
            // ----------------------------------------------------------------
            // Use the first three atoms to detect chirality
            // This should work well for 3 and four atom case
            // With more than four bonded atoms, well... that's tricky.

            // Three ordered source vectors
            UnitVec3 s1 = bondVectors[zeroBondCenterIndex].first;
            UnitVec3 s2 = bondVectors[oneBondCenterIndex].first;
            UnitVec3 s3 = bondVectors[twoBondCenterIndex].first;

            // And three ordered target vectors
            UnitVec3 t1 = bondVectors[zeroBondCenterIndex].second;
            UnitVec3 t2 = bondVectors[oneBondCenterIndex].second;
            UnitVec3 t3 = bondVectors[twoBondCenterIndex].second;

            // Target reference plane
            UnitVec3 targetPlaneNormal(cross(t1, t2));
            UnitVec3 sourcePlaneNormal(cross(s1, s2));

            // If any Planar groups of the target structure is farther than <tolerance> from planar break them all
            bool doBreakPlanes = false;
            for (int bondCnt = 0; bondCnt < (int) bondVectors.size(); ++bondCnt) {

                const BondCenter& bondCenter = getBondCenter(compoundBCIxes[bondCnt]);

                if (bondCenter.getChirality() != BondCenter::Planar){ // not defined as Planar
                    continue;}

                Angle sinePlaneDeviation = dot(bondVectors[bondCnt].second, targetPlaneNormal);
                if ( std::abs(sinePlaneDeviation) < std::sin(breakPlanarityThreshold) ){ // close to Planar
                    continue;}

                doBreakPlanes = true;

                break;
            }

            // Set chiralities for former Planar bonds
            if (doBreakPlanes) {

                for (int bondCnt = 0; bondCnt < (int) bondVectors.size(); ++bondCnt) {

                    BondCenter& bondCenter = updBondCenter(compoundBCIxes[bondCnt]);

                    // Can't break planarity if it's not planar to begin with
                    if (bondCenter.getChirality() != BondCenter::Planar) 
                        {continue;}

                    // Don't try to break planarity of first and second bond centers
                    // because they define the plane
                    const BondCenterInfo& bondCenterInfo = getBondCenterInfo(compoundBCIxes[bondCnt]);
                    if (bondCenterInfo.getAtomBondCenterIndex() == 0) continue;
                    if (bondCenterInfo.getAtomBondCenterIndex() == 1) continue;

                    // OK, if we got this far, we must break planarity
                    #ifdef DEBUG_MOLMODEL
                        std::cerr <<__FILE__<<":"<<__LINE__<< " WARNING: matching out-of-plane atoms about atom ";
                        std::cerr << getAtomName(atomIndex);
                        std::cerr << ". Note that here we are using residue INDEX, not residue NUMBER. Residue indices start at 0.";
                        std::cerr << std::endl;
                    #endif

                    Angle sinePlaneDeviation = dot(bondVectors[bondCnt].second, targetPlaneNormal);

                    if (sinePlaneDeviation < 0){ 
                        bondCenter.setChirality(BondCenter::LeftHanded);
                    }
                    else{
                        bondCenter.setChirality(BondCenter::RightHanded);
                    }
                }
            }

            if (flipAll) { // flip entire atom or nothing
                Real sourceChirality = dot(cross(s1, s2),s3);
                Real targetChirality = dot(cross(t1, t2),t3);

                // Reverse chirality of bond centers that differ from those in target structure
                // same sign means same chirality
                if (sourceChirality * targetChirality < 0)
                { // mismatch
                    #ifdef DEBUG_MOLMODEL
                        std::cerr << "WARNING: Using unexpected chirality about atom ";
                        std::cerr << getAtomName(atomIndex);
                        std::cerr << std::endl;
                    #endif

                    // flip the chirality of every handed bond center in the target atom
                    for (CompoundAtom::BondCenterIndex bcIx(0); bcIx < atom.getNumBonds(); ++bcIx)
                    {
                        BondCenter& bondCenter = updBondCenter(getBondCenterInfo(atomIndex, bcIx));
                        BondCenter::Chirality newChirality = bondCenter.getChirality();
                        switch(bondCenter.getChirality()) {
                            case BondCenter::RightHanded:
                                newChirality = BondCenter::LeftHanded;
                                break;
                            case BondCenter::LeftHanded:
                                newChirality = BondCenter::RightHanded;
                                break;
                            default:
                                break;
                        }
                        bondCenter.setChirality(newChirality);
                    }
                } // end if chirality differs
            }
            else { // flip on a bondcenter by bondcenter basis
                // Flip chirality on a BondCenter by BondCenter basis
                for (int bondCnt = 0; bondCnt < (int) bondVectors.size(); ++bondCnt) 
                {
                    // First two bond centers cannot be chiral
                    if (bondCnt == zeroBondCenterIndex) continue;
                    if (bondCnt == oneBondCenterIndex) continue;

                    // Measure source and target chiralities
                    UnitVec3 sourceBondVec = bondVectors[bondCnt].first;
                    UnitVec3 targetBondVec = bondVectors[bondCnt].second;
                    Real sourceChirality = dot(cross(s1, s2),sourceBondVec);
                    Real targetChirality = dot(cross(t1, t2),targetBondVec);

                    if (sourceChirality * targetChirality < 0) // mismatched chirality
                    {
                        const BondCenterInfo& bondCenterInfo = getBondCenterInfo(compoundBCIxes[bondCnt]);

                        #ifdef DEBUG_MOLMODEL
                            Compound::AtomIndex partnerAtomIndex = 
                                getBondCenterInfo(bondCenterInfo.getBondPartnerBondCenterIndex())
                                .getAtomIndex();
                            std::cerr << "WARNING: Flipping chirality of bond from atom ";
                            std::cerr << getAtomName(atomIndex);
                            std::cerr << " to atom ";
                            std::cerr << getAtomName(partnerAtomIndex);
                            std::cerr << std::endl;
                        #endif
                        
                        BondCenter& bondCenter = updBondCenter(bondCenterInfo);
                        switch(bondCenter.getChirality()) {
                            case BondCenter::RightHanded:
                                bondCenter.setChirality(BondCenter::LeftHanded);
                                break;
                            case BondCenter::LeftHanded:
                                bondCenter.setChirality(BondCenter::RightHanded);
                                break;
                            default:
                                break; // TODO what should planar do here?
                        }
                    }
                }
            }


            
        } // every atom

        return *this;
    }

    /*!
    * <!-- Set bond length in vector<BondInfo> allBonds -->
    */
    CompoundRep& matchDefaultBondLengths(const Compound::AtomTargetLocations& atomTargets) 
    {
        // Get neighbour list: std::vector<std::vector<cAIx>>
        std::vector< AtomIndexVector > atomPairs = getBondedAtomRuns(2, atomTargets);

        // Loop over those pairs of atoms and set the bond length default to the target distances
        // This method is broken into two parts like this to serve as an example for the more
        // complex methods to follow.
        std::vector< AtomIndexVector >::const_iterator bonds12Ix;
        for (bonds12Ix = atomPairs.begin(); bonds12Ix != atomPairs.end();
        ++bonds12Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds12Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds12Ix)[1];

            // For efficiency, only set bonds lengths once per bond
            if (atomIndex2 > atomIndex1) continue;

            // compute distance
            Vec3 bondVector = atomTargets.find(atomIndex1)->second - atomTargets.find(atomIndex2)->second;
            Real distance = std::sqrt(dot(bondVector, bondVector));

            #ifdef DEBUG_MOLMODEL
                //std::cout<<__FILE__<<":"<<__LINE__<<" atomIndex1 "<<atomIndex1
                //  <<" atomIndex2 "<< atomIndex2 << " distance = "<< distance
                //  <<std::endl;
            #endif

            // Set bond length
            SimTK::Bond &  thisBond = updBond(updBondInfo(updAtomInfo(atomIndex1), updAtomInfo(atomIndex2)));
            thisBond.setDefaultBondLength(distance);

            // updBond(updBondInfo(
            //     updAtomInfo(atomIndex1),
            //     updAtomInfo(atomIndex2))).setDefaultBondLength(distance);


        } // every atom pair

        return *this;
    }

    /*!
    * <!-- atom2.updBondCenter(largerId).setDefaultBond1Angle(angle)
    * of the middle atom BC -->
    */
    CompoundRep& matchDefaultBondAngles(const Compound::AtomTargetLocations& atomTargets) 
    {
        //std::cout << "matchDefaultBondAngles" << std::endl;
        //std::cout<<"BEGIN  matchDefaultBondAngles"<<std::endl;
        std::vector< AtomIndexVector > atomTriples = getBondedAtomRuns(3, atomTargets);

        std::vector< AtomIndexVector >::const_iterator bonds13Ix;
        for (bonds13Ix = atomTriples.begin(); bonds13Ix != atomTriples.end(); ++bonds13Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds13Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds13Ix)[1];
            Compound::AtomIndex atomIndex3 = (*bonds13Ix)[2];

            // for efficiency, set each angle only once, not both 3->2->1 and 1->2->3
            if (atomIndex3 < atomIndex1) continue;

            // Calculate atomTargets angle
            UnitVec3 v1(atomTargets.find(atomIndex1)->second - atomTargets.find(atomIndex2)->second);
            UnitVec3 v2(atomTargets.find(atomIndex3)->second - atomTargets.find(atomIndex2)->second);

            Real dotProduct = dot(v1, v2);
            assert(dotProduct < 1.1);
            assert(dotProduct > -1.1);
            if (dotProduct > 1.0) dotProduct = 1.0;
            if (dotProduct < -1.0) dotProduct = -1.0;
            Real angle = std::acos(dotProduct);

            // Set the larger BC angle in relation to the smaller (0 or 1) BC
            // std::cerr << angle / SimTK::Deg2Rad << std::endl;
            //std::cout<<__FILE__<<":"<<__LINE__
            //  <<" angle, atomIndex1, atomIndex2, atomIndex3 "<<angle<<" , "
            //  << atomIndex1<<" , "<< atomIndex2<<" , "<< atomIndex3
            //  <<std::endl;
            setDefaultBondAngle(angle, atomIndex1, atomIndex2, atomIndex3);
        }

        return *this;
    }

    /*!
    * <!-- Continuation of the matchDefaultBondAngles: sets the direction of
    * the BCs >= 1 -->
    */
    CompoundRep& matchDefaultDirections(const Compound::AtomTargetLocations& atomTargets){
        
        std::vector< AtomIndexVector > atomRun = getBondedAtomRuns(1, atomTargets);
        
        for(const auto& atomRIx : atomRun) {

            const Compound::AtomIndex atomIx = atomRIx[0];
            CompoundAtom& atom = updAtom(atomIx);
            const AtomInfo& atomInfo = getAtomInfo(atomIx);

            const Compound::AtomIndex neighborAtomIx = atomRIx[1];
            const BondCenterInfo& bondCenterInfo = getBondCenterInfo(
                    getAtomInfo(atomIx),
                    getAtomInfo(neighborAtomIx) );
            CompoundAtom::BondCenterIndex atomBondCenterIndex =
                bondCenterInfo.getAtomBondCenterIndex();

            CompoundAtom::BondCenterIndex  BCIx = atomBondCenterIndex;

            // Go through bond centers on atom. Order counts.
            //for (CompoundAtom::BondCenterIndex BCIx(0); BCIx < atom.getNumBonds(); ++BCIx) { // RESTORE

                BondCenter &BC0 = atom.updBondCenter(CompoundAtom::BondCenterIndex(0));

                if(BCIx == 0) {
                    continue;

                }else if(BCIx == 1){ // Rotate with theta in the initial plane
                    BondCenter &BC1 = atom.updBondCenter(BCIx);
                    const UnitVec3& BC0_dir = BC0.updDirection();
                    const UnitVec3& BC1_dir = BC1.updDirection();

                    // const UnitVec3 rotAxis(BC1_dir % BC0_dir); // RESTORE
                    const UnitVec3 rotAxis(0, 0, -1);

                    const Angle rotAngle = BC1.getDefaultBond1Angle();
                    const Rotation rotMat = Rotation(rotAngle, rotAxis);
                    const UnitVec3 newDir = rotMat * BC0.getDirection();

                    // std::cout << "CompoundRep::matchDefaultDirections cAIx " << atomIx << " BCIx " << BCIx
                    //     << " BC0_dir " << BC0_dir
                    //     << " BC1_dir " << BC1_dir
                    //     << " rotAxis " << rotAxis
                    //     << " rotAngle " << rotAngle
                    //     <<" newDir " << newDir
                    //     << std::endl;

                    BC1.setDirection(newDir);

                }else if(BCIx > 1) { // Use Paul's method
                    BondCenter &BC_gt1 = atom.updBondCenter(BCIx);
                    const UnitVec3 a1 = atom.getBondCenterDirectionInAtomFrame(CompoundAtom::BondCenterIndex(0));
                    const UnitVec3 a2 = atom.getBondCenterDirectionInAtomFrame(CompoundAtom::BondCenterIndex(1));
                    const Angle theta1 = BC_gt1.getDefaultBond1Angle();
                    const Angle theta2 = BC_gt1.getDefaultBond2Angle();
                    const BondCenter::Chirality chirality = BC_gt1.getChirality();

                    // std::cout << "CompoundRep::matchDefaultDirections: bc a1 theta1 a2 theta2 chir dir"
                    //     <<" " << BCIx <<" "<< a1 <<" "<< theta1 <<" "<< a2 <<" "<< theta2 <<" "<< chirality
                    //     <<" "<< BondCenter::getBondDirection(a1, theta1, a2, theta2, chirality) << std::endl;

                    BC_gt1.setDirection(BondCenter::getBondDirection(a1, theta1, a2, theta2, chirality));
                }
            //} // every bond center RESTORE
        }

        return *this;
    }

    /*!
    * <!-- Helper method for matchDefaultDihedralAngles -->
    */
    bool isPlanarBond(
            Compound::AtomIndex atomIndex2,
            Compound::AtomIndex atomIndex3) 
    {
        const CompoundAtom& atom2 = getAtom(atomIndex2);
        const CompoundAtom& atom3 = getAtom(atomIndex3);

        // Three criteria for whether bond is planar

        // 1) both central atoms have three bonds
        if (atom2.getNumBondCenters() != 3) return false;
        if (atom3.getNumBondCenters() != 3) return false;

        // 2) third bond center on each of those atoms is planar
        if (atom2.getBondCenter(CompoundAtom::BondCenterIndex(2)).getChirality() != BondCenter::Planar)
            return false;
        if (atom3.getBondCenter(CompoundAtom::BondCenterIndex(2)).getChirality() != BondCenter::Planar)
            return false;

        // 3) initial dihedral angle is near 0 or 180 degrees
        const BondCenter& bondCenter23 = 
            getBondCenter(getBondCenterInfo(getAtomInfo(atomIndex2), getAtomInfo(atomIndex3)).getIndex());
        Angle initialAngle = bondCenter23.getDefaultDihedralAngle();
        // normalize to be near zero
        while (initialAngle < -90.0 * Deg2Rad) initialAngle += 180.0 * Deg2Rad;
        while (initialAngle > 90.0 * Deg2Rad) initialAngle -= 180.0 * Deg2Rad;
        if (std::abs(initialAngle) > 0.001) return false;

        // If we got this far, it must be planar
        return true;
    }

    /*!
    * <!--  -->
    */
    CompoundRep& matchDefaultDihedralAngles(
            const Compound::AtomTargetLocations& atomTargets, 
            Compound::PlanarBondMatchingPolicy policy) 
    {
        //std::cout << "matchDefaultDihedralAngles" << std::endl;
        //std::cout<<"BEGIN   matchDefaultDihedralAngles"<<std::endl;
        std::vector< AtomIndexVector > atomQuads = getBondedAtomRuns(4, atomTargets);

        std::vector< AtomIndexVector >::const_iterator bonds14Ix;
        for (bonds14Ix = atomQuads.begin(); bonds14Ix != atomQuads.end(); ++bonds14Ix) 
        {
            Compound::AtomIndex atomIndex1 = (*bonds14Ix)[0];
            Compound::AtomIndex atomIndex2 = (*bonds14Ix)[1];
            Compound::AtomIndex atomIndex3 = (*bonds14Ix)[2];
            Compound::AtomIndex atomIndex4 = (*bonds14Ix)[3];
            // for efficiency, set each dihedral only once
            if (atomIndex4 < atomIndex1) continue;

			// Don't set dihedrals involving ring-closing bonds, as these can damage "real" dihedrals
			if ( getBond(atomIndex2, atomIndex1).isRingClosingBond() ) continue;
			if ( getBond(atomIndex3, atomIndex4).isRingClosingBond() ) continue;

            // Compute and set dihedral angle
            UnitVec3 bond12(atomTargets.find(atomIndex2)->second - atomTargets.find(atomIndex1)->second);
            UnitVec3 bond23(atomTargets.find(atomIndex3)->second - atomTargets.find(atomIndex2)->second);
            UnitVec3 bond34(atomTargets.find(atomIndex4)->second - atomTargets.find(atomIndex3)->second);
            Angle angle = SimTK::calcDihedralAngle(bond12, bond23, bond34);

            // assert(false);  // need to implement general setDefaultDihedralAngle method

            // Don't set torsion for planar bonds, except maybe to Flip them
            if  (policy == Compound::KeepPlanarBonds)
            {
                if ( isPlanarBond(atomIndex2, atomIndex3) )
                    continue;
            }

            if (policy == Compound::FlipPlanarBonds)
                if ( isPlanarBond(atomIndex2, atomIndex3) ) {
                    // TODO - decide whether to flip the dihedral angle 180 degrees
                    Angle initialAngle = calcDefaultDihedralAngle(                
                        atomIndex1, 
                        atomIndex2, 
                        atomIndex3, 
                        atomIndex4);
                    Angle diffAngle = angle - initialAngle;
                    // normalize to range (-180 degrees, 180 degrees)
                    while ( -SimTK::Pi >= diffAngle ) diffAngle += 2 * SimTK::Pi;
                    while ( SimTK::Pi < diffAngle )   diffAngle -= 2 * SimTK::Pi;
                    // Either flip the dihedral 180 degrees...
                    if (std::abs(diffAngle) > 0.5 * SimTK::Pi)
                        angle = initialAngle + SimTK::Pi;          
                    else // ... or do nothing.
                        continue;
                }

/*            std::cout<<" angle, atomIndex1, atomIndex2, atomIndex3, atomIndex4 "<<angle<<" , "<< atomIndex1 ;

            //std::cout<<__FILE__<<":"<<__LINE__<<" angle, atomIndex1, atomIndex2, atomIndex3, atomIndex4 "<<angle<<" , "<< atomIndex1;
            std::cout << ":"<<getAtomName(atomIndex1);
            std::cout<<" , "<< atomIndex2;
            std::cout << ":"<<getAtomName(atomIndex2);
            std::cout<<" , "<< atomIndex3;
            std::cout << ":"<<getAtomName(atomIndex3)  <<" , "<< atomIndex4;
            std::cout << ":"<<getAtomName(atomIndex4)   <<std::endl;*/
            //std::cout << "RECONSTRUCT STEP 1.0.0 " << offsetAngle4 << std::endl << std::flush;
            setDefaultDihedralAngle( 
                angle, 
                atomIndex1, 
                atomIndex2, 
                atomIndex3, 
                atomIndex4);
        }

        return *this;
    }

    /*!
     * <!--  -->
     */
    CompoundRep& matchDefaultTopLevelTransform(const Compound::AtomTargetLocations& atomTargets) 
    {

        Transform adjustment = getTransformAndResidual(atomTargets).transform;

        //std::cout<<"adjustment:"<<std::endl<<adjustment<<std::endl;
        setTopLevelTransform( adjustment * getTopLevelTransform() );



        // // STUDY localTransform
        // Compound::AtomTargetLocations::const_iterator tI;
        // for (tI = atomTargets.begin(); tI != atomTargets.end(); ++tI) 
        // {
        //     Compound::AtomIndex cAIx = tI->first;
        //     const SimTK::AtomInfo & atomInfo = getAtomInfo(cAIx);
        //     const SimTK::CompoundAtom & atom = getAtom(atomInfo);

        //     std::cout << "STUDY CompoundRep::matchDefaultTopLevelTransform cAIx T " << cAIx <<" \n "<< atom.getLocalTransform() << std::endl;
        // }


        return *this;
    }


    /*!
     * <!--  -->
     */
    TransformAndResidual getTransformAndResidual(
        const Compound::AtomTargetLocations& atomTargets) const
    {
        // Declare a std::vector<Vec3Pair>
        Kabsch78::VectorSet vecPairs;

        // Try for more efficient calculation of atom starting locations
        std::vector<Transform> atomSourceFrames(getNumAtoms());
        invalidateAtomFrameCache(atomSourceFrames, getNumAtoms());
        calcDefaultAtomFramesInCompoundFrame(atomSourceFrames);
        Real weight = 1.0;
        
        Compound::AtomTargetLocations::const_iterator tI;
        for (tI = atomTargets.begin(); tI != atomTargets.end(); ++tI) 
        {
            Compound::AtomIndex atomIndex = tI->first;
            const Vec3& target = tI->second;
           
            // slow
            // Vec3 source = calcDefaultAtomLocationInGroundFrame(getAtomName(atomIndex));

            // faster
            const Vec3 source = getTopLevelTransform() * atomSourceFrames[atomIndex].T();

            vecPairs.push_back(Vec3Pair(source, target, weight));

        }

        return Kabsch78::superpose(vecPairs);
    }


    const std::set<Compound::AtomName>& getAtomSynonyms(Compound::AtomIndex a) const
    {
        const AtomInfo& atomInfo = getAtomInfo(a);
        return atomInfo.getNames();
    }
    // scf added a new parameter .. when guessCoordinates is true, default atom
    // positions from the biopolymer are pushed into the returned AtomTargetLocations.
    virtual Compound::AtomTargetLocations createAtomTargets
       (const PdbStructure& targetStructure, bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            int residueNumber = getPdbResidueNumber();
            char insertionCode = ' '; // TODO - make this an attribute of the residue?
            String chainId = getPdbChainId();

            String atomName = getAtomName(a);
            ////std::cout<<__FILE__<<":"<<__LINE__<<" "<<a<<" "<<atomName<<" "<<residueNumber<<std::endl;
            // search synonyms if we cannot find this atom in the structure
            if (! targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) )
                        break;
                }
            }

            if ( targetStructure.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode), chainId) ) {
                const PdbAtom& pdbAtom = targetStructure.getAtom( atomName, PdbResidueId(residueNumber, insertionCode), chainId );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else {
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

    virtual Compound::AtomTargetLocations createAtomTargets
       (const PdbChain& targetChain, bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            int residueNumber = getPdbResidueNumber();
            char insertionCode = ' '; // TODO - make this an attribute of the residue?
            String chainId = getPdbChainId();

            String atomName = getAtomName(a);
            // search synonyms if we cannot find this atom in the structure
            if (! targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) )
                        break;
                }
            }

            if ( targetChain.hasAtom(atomName, PdbResidueId(residueNumber, insertionCode)) ) {
                const PdbAtom& pdbAtom = targetChain.getAtom( atomName, PdbResidueId(residueNumber, insertionCode) );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else {
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

	// TODO - too much copy paste in these createAtomTargets methods
    virtual Compound::AtomTargetLocations createAtomTargets(const PdbResidue& targetResidue, const bool guessCoordinates = false) const 
    {
        Compound::AtomTargetLocations answer;
		if (getPdbResidueNumber() != targetResidue.getPdbResidueNumber()) return answer;

        for (Compound::AtomIndex a(0); a < getNumAtoms(); ++a) 
        {
            //int residueNumber = getPdbResidueNumber();
            //char insertionCode = ' '; // TODO - make this an attribute of the residue?
            //String chainId = getPdbChainId();

            String atomName = getAtomName(a);
            // search synonyms if we cannot find this atom in the structure
            if (! targetResidue.hasAtom(atomName) )
            {
                const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(a);
                std::set<Compound::AtomName>::const_iterator nameIx;
                for (nameIx = atomNames.begin(); nameIx != atomNames.end(); ++nameIx)
                {
                    atomName = *nameIx;
                    if ( targetResidue.hasAtom(atomName) )
                        break;
                }
            }

            if ( targetResidue.hasAtom(atomName) ) {
                const PdbAtom& pdbAtom = targetResidue.getAtom( atomName );
                if (pdbAtom.hasLocation())
                    answer[a] = pdbAtom.getLocation();
            }
            else { // skip atoms not found in pdbResidue
                // std::cerr << atomName << std::endl;
            }
        }

        return answer;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populateDefaultPdbChain(
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
		Transform myTransform = transform;
            // if (!hasParentCompound()) {
                myTransform = myTransform * getTopLevelTransform();
            // }

        int residueNumber = getPdbResidueNumber();

        // try to guess when to use internal PdbResidueNumber vs. defaultNextResidueNumber
        if (-9999 > getPdbResidueNumber()) 
            residueNumber = defaultNextResidueNumber;

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        pdbChain.appendResidue( PdbResidue(getOwnerHandle(), residueNumber, myTransform) );

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populatePdbChain(
        const State& state, 
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {

		Transform myTransform = transform;
		// Don't apply top-level transform for state-taking methods!!!
        //    if (!hasParentCompound()) {
        //        myTransform = myTransform * getTopLevelTransform();
        //    }

        int residueNumber = getPdbResidueNumber();

        // try to guess when to use internal PdbResidueNumber vs. defaultNextResidueNumber
        if (-9999 > getPdbResidueNumber()) 
            residueNumber = defaultNextResidueNumber;

        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber)))
        {
            ++defaultNextResidueNumber;
            residueNumber = defaultNextResidueNumber;        
        }

        pdbChain.appendResidue( PdbResidue(state, getOwnerHandle(), residueNumber, myTransform) );

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    // One argument version of writeDefaultPdb begins numbering atoms at 1
    std::ostream& writeDefaultPdb(std::ostream& os, const Transform& transform) const;
    std::ostream& writeDefaultPdb(std::ostream& os, int& nextSerialNumber, const Transform& transform) const;

    //std::ostream& writeDefaultAtomPdb(
    //    const Compound::AtomName& name, 
    //    std::ostream& os, 
    //    int& nextSerialNumber,
    //    const Transform& transform
    //    ) const;

    std::ostream& writePdb(
        const State& state, 
        std::ostream& os, 
        const Transform& transform) const;

    std::ostream& writePdb(
        const State& state, 
        std::ostream& os, 
        int& nextSerialNumber,
        const Transform& transform) const;

    //std::ostream& writeAtomPdb(
    //    const State& state, 
    //    const Compound::AtomName&   name, 
    //    std::ostream& os, 
    //    int& nextSerialNumber, 
    //    const Transform& transform) const;

    //std::ostream& writeAtomPdb(
    //    const Compound::AtomName&   name, 
    //    std::ostream&               os, 
    //    int&                        nextSerialNumber,
    //    const Vec3&                 location
    //    ) const;
// protected:

    bool hasInboardBondCenter() const;

    CompoundRep& convertInboardBondCenterToOutboard();

    const BondCenter& getInboardBondCenter() const;
    BondCenter& updInboardBondCenter();
    const BondCenterInfo& getInboardBondCenterInfo() const;
    BondCenterInfo& updInboardBondCenterInfo();

    CompoundRep& setInboardBondCenter(const Compound::BondCenterName& n);
    CompoundRep& setInboardBondCenter(Compound::BondCenterIndex id);

    Compound::BondCenterIndex addLocalCompound(
        const Compound::Name& scName, 
        const Compound& subcompound,
        const Transform& location = Transform());

    // Copy atoms etc.
    // Returns new bond center index of absorbed inboard bond center
    Compound::BondCenterIndex absorbSubcompound(const Compound::Name& scName, const Compound& subcompound, bool isBase);

    const BondInfo& getBondInfo(Compound::BondIndex bi) const {
        return allBonds[bi];
    }
    BondInfo& updBondInfo(Compound::BondIndex bi) {
        return allBonds[bi];
    }

    const BondInfo& getBondInfo(const AtomInfo& a1, const AtomInfo& a2) const {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(a1.getIndex(), a2.getIndex());
        Compound::BondIndex bi = AIxPair_To_BondIx.find(key)->second;

        return allBonds[bi];
    }
    BondInfo& updBondInfo(const AtomInfo& a1, const AtomInfo& a2) {
        std::pair<Compound::AtomIndex, Compound::AtomIndex> key(a1.getIndex(), a2.getIndex());
        Compound::BondIndex bi = AIxPair_To_BondIx.find(key)->second;

        return allBonds[bi];
    }

	const Bond& getBond(Compound::AtomIndex atom1, Compound::AtomIndex atom2) {
		return getBond( getBondInfo(getAtomInfo(atom1), getAtomInfo(atom2)) );
	}

    const Bond& getBond(const BondInfo& bondInfo) const {
        return bondInfo.getBond();
        //if ( bondInfo.isLocalBond() || bondInfo.isRingClosingBond() );
        //else {
        //    assert(false);
        //}
        //    assert(bondInfo.isSubcompoundBond());
        //    Compound::Index subcompoundId = bondInfo.getSubcompoundId();
        //    // const CompoundRep& scRep = bondInfo.getSubcompound().getImpl();
        //    const CompoundRep& scRep = getSubcompound(subcompoundId).getImpl();
        //    return scRep.getBond(scRep.getBondInfo(bondInfo.getSubcompoundBondIndex()));
        //}
    }
    Bond& updBond(BondInfo& bondInfo) {
        return bondInfo.updBond();
        //if ( bondInfo.isLocalBond() || bondInfo.isRingClosingBond() );
        //else {
        //    assert(false);
        //}
        //    assert(bondInfo.isSubcompoundBond());
        //    CompoundRep& scRep = updSubcompound(bondInfo.getSubcompoundId()).updImpl();
        //    // CompoundRep& scRep = bondInfo.updSubcompound().updImpl();
        //    return scRep.updBond(scRep.updBondInfo(bondInfo.getSubcompoundBondIndex()));
        //}
    }

    Transform calcDefaultBondCenterFrameInAtomFrame(const BondCenterInfo& info) const;
    const Transform calcDefaultBondCenterFrameInCompoundFrame(const Compound::BondCenterName name) const;

    /*!
    * <!-- Cache method used in O(n) all atom Frame computation --> 
    */
    const Transform calcDefaultBondCenterFrameInCompoundFrame(
        const BondCenterInfo& info,
        std::vector<Transform>& atomFrameCache) const;

    Compound::BondCenterIndex getBondCenterIndex(const Compound::BondCenterName& name) const;

    // return the bond center on atom1 that is attached to atom2
    const BondCenterInfo& getBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2) const;
    BondCenterInfo& updBondCenterInfo(const Compound::AtomName& atom1, const Compound::AtomName& atom2);
    const BondCenterInfo& getBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) const;
    BondCenterInfo& updBondCenterInfo(const AtomInfo& atom1, const AtomInfo& atom2) ;
    BondCenterInfo& updBondCenterInfo(const Compound::BondCenterName&);
    const BondCenterInfo& getBondCenterInfo(const Compound::BondCenterName&) const;
    BondCenterInfo& updBondCenterInfo(Compound::BondCenterIndex);
    const BondCenterInfo& getBondCenterInfo(Compound::BondCenterIndex) const;
    BondCenterInfo& updBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex);
    const BondCenterInfo& getBondCenterInfo(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) const;
    BondCenterInfo& updBondCenterInfo(BondCenterInfo::AtomKey key);
    const BondCenterInfo& getBondCenterInfo(BondCenterInfo::AtomKey key) const;

    bool hasBondCenter(const Compound::BondCenterName&) const;
    bool hasBondCenter(Compound::AtomIndex atomId, CompoundAtom::BondCenterIndex atomBondCenterIndex) const {
        return hasBondCenter(BondCenterInfo::AtomKey(atomId, atomBondCenterIndex));
    }
    bool hasBondCenter(Compound::BondCenterIndex id) const;
    bool hasBondCenter(const BondCenterInfo::AtomKey& key) const {
        return bondCenterIndicesByAtomKey.find(key) != bondCenterIndicesByAtomKey.end();
    }

    BondCenter& updBondCenter(const Compound::BondCenterName& name);
    const BondCenter& getBondCenter(const Compound::BondCenterName& name) const;
    BondCenter& updBondCenter(Compound::BondCenterIndex id);

    const BondCenter& getBondCenter(Compound::BondCenterIndex id) const;

    const BondCenter& getBondCenter(const BondCenterInfo& info) const;

    BondCenter& updBondCenter(const BondCenterInfo& info);

    Compound::AtomIndex getAtomIndex(const Compound::AtomName& atomName) const;

    AtomInfo& updAtomInfo(const Compound::AtomName&);

    const AtomInfo& getAtomInfo(const Compound::AtomName& name) const {
        // assert(CompoundPathName::isValidAtomName(name));
        assert(hasAtom(name));
        return getAtomInfo(atomName_To_atomId.find(name)->second);
    }

    AtomInfo& updAtomInfo(Compound::AtomIndex);
    const AtomInfo& getAtomInfo(Compound::AtomIndex id) const {
        assert((Compound::AtomIndex)allAtoms.size() > id);
        assert(0 <= id);
        return allAtoms[id];
    }
    //AtomInfo& updAtomInfo(Compound::Index subcompoundId, Compound::AtomIndex subAtomIndex) {
    //    const CompoundRep&      scRep           = getSubcompound(subcompoundId).getImpl();
    //    const AtomInfo&         scAtomInfo      = scRep.getAtomInfo(subAtomIndex);
    //    const Compound::AtomIndex  parentAtomIndex    = scAtomInfo.getParentCompoundAtomIndex();
    //    return updAtomInfo(parentAtomIndex);
    //}
    //const AtomInfo& getAtomInfo(Compound::Index subcompoundId, Compound::AtomIndex subAtomIndex) const {
    //    const CompoundRep&      scRep           = getSubcompound(subcompoundId).getImpl();
    //    const AtomInfo&         scAtomInfo      = scRep.getAtomInfo(subAtomIndex);
    //    const Compound::AtomIndex  parentAtomIndex    = scAtomInfo.getParentCompoundAtomIndex();
    //    return getAtomInfo(parentAtomIndex);
    //}

    // desk_mass_related

    const SimTK::mdunits::Mass getAtomMass(Compound::AtomIndex id) const
    {
        return getAtom(id).getMass();
    }

    void setAtomMass(Compound::AtomIndex id, const SimTK::mdunits::Mass& mass) {
        updAtom(id).setMass(mass);
    }

    void updAtomMass(Compound::AtomIndex id, const SimTK::mdunits::Mass& mass) {
        updAtom(id).updateMass(mass);
    }

    // _end_ desk_mass_related 

    bool hasAtom(const Compound::AtomName& name) const;

    bool hasAtom(Compound::AtomIndex atomId) const {
        if (atomId < 0) return false;
        if (atomId >= (Compound::AtomIndex) allAtoms.size()) return false;

        return true;
    }

    const CompoundAtom& getAtom(const Compound::AtomName& name) const {
        return getAtom(getAtomInfo(name));
    }
    CompoundAtom& updAtom(const Compound::AtomName& name) {
        return updAtom(updAtomInfo(name));
    }
    const CompoundAtom& getAtom(Compound::AtomIndex id) const {
        return getAtom(getAtomInfo(id));
    }
    CompoundAtom& updAtom(Compound::AtomIndex id) {
        return updAtom(updAtomInfo(id));
    }
    CompoundAtom& updAtom(AtomInfo& info) {
        return info.updAtom();
    }
    const CompoundAtom& getAtom(const AtomInfo& info) const {
        return info.getAtom();
    }

    Compound::BondIndex getNumBonds() const {return Compound::BondIndex(allBonds.size());}
    Compound::AtomIndex getBondAtomIndex(Compound::BondIndex bid, int which) const;

    //const CompoundInfo& getSubcompoundInfo(const Compound::Name& name) const 
    //{
    //    assert( hasSubcompound(name) );

    //    // TODO - parse "X/Y" indirect subcompound identifiers
    //    Compound::Index subcompoundId;

    //    // First check for simple subcompound name without any "/" separators
    //    if (CompoundPathName::isValidSubcompoundName(name)) 
    //    {
    //        subcompoundId = subcompoundIdsByName.find(name)->second;
    //    }
    //    else // parse "X/Y" path type subcompound names
    //    {
    //        std::vector<String> tokens;
    //        if (CompoundPathName::isValidSubcompoundPathName(name, &tokens)) 
    //        {
    //            String subcompoundName = tokens[0];
    //            if (! hasSubcompound(subcompoundName)) 
    //            {
    //                assert(false); // TODO - raise exception - hasSubcompound() check should have caught this
    //            }

    //            Compound::Index topSubcompoundId = getSubcompoundInfo(subcompoundName).getIndex();

    //            const CompoundRep& scRep = getSubcompound(subcompoundName).getImpl();
    //            Compound::Index childSubcompoundId = scRep.getSubcompoundInfo(CompoundPathName::shiftLeftPathName(name)).getIndex();

    //            // TODO find CompoundInfo that matches topSubcompoundId and childSubcompoundId
    //            // TODO this is not efficient, checking every subcompound
    //            std::vector<CompoundInfo>::const_iterator scI;
    //            for (scI = allSubcompounds.begin(); scI != allSubcompounds.end(); ++scI) 
    //            {
    //                if (scI->isLocal()) continue;
    //                if (scI->isBonded()) continue;
    //                if (scI->getIntermediateSubcompoundId() != topSubcompoundId) continue;
    //                if (scI->getIntermediateSubcompoundSubcompoundId() != childSubcompoundId) continue;

    //                // If we get this far, we have found the correct subcompound
    //                subcompoundId = scI->getIndex();
    //                break;
    //            }

    //            Compound::Index invalidCompoundId;
    //            assert(subcompoundId != invalidCompoundId);
    //        }
    //        else { // string not well formed
    //            assert(false);
    //            // TODO raise exception - hasSubcompound() check should have caught this
    //        }
    //    }

    //    return getSubcompoundInfo(subcompoundId);
    //}
    //CompoundInfo& updSubcompoundInfo(const Compound::Name& name) {
    //    assert( hasSubcompound(name) );

    //    const Compound::Index id = subcompoundIdsByName.find(name)->second;
    //    return updSubcompoundInfo(id);
    //}
    //const CompoundInfo& getSubcompoundInfo(Compound::Index id) const {
    //    assert (0 <= id);
    //    assert ((Compound::Index)allSubcompounds.size() > id);

    //    return allSubcompounds[id];
    //}
    //CompoundInfo& updSubcompoundInfo(Compound::Index id) {
    //    assert (0 <= id);
    //    assert ((Compound::Index)allSubcompounds.size() > id);

    //    return allSubcompounds[id];
    //}

    //Compound& updSubcompound(const Compound::Name& name) {
    //    CompoundInfo& info = updSubcompoundInfo(name);
    //    return updSubcompound(info);
    //}
    //const Compound& getSubcompound(const Compound::Name& name) const {
    //    const CompoundInfo& info = getSubcompoundInfo(name);
    //    return getSubcompound(info);
    //}
    //Compound& updSubcompound(Compound::Index id) {
    //    CompoundInfo& info = updSubcompoundInfo(id);
    //    return updSubcompound(info);
    //}
    //const Compound& getSubcompound(Compound::Index id) const {
    //    const CompoundInfo& info = getSubcompoundInfo(id);
    //    return getSubcompound(info);
    //}

    //const Compound& getSubcompound(const CompoundInfo& info) const {
    //    if (info.isLocal()) {
    //        return info.getCompound();
    //    }
    //    else if (info.isBonded()) {
    //        return info.getCompound();
    //    }
    //    else { // subcompound of subcompound
    //        const CompoundRep& sc1 = getSubcompound(info.getIntermediateSubcompoundId()).getImpl();
    //        return sc1.getSubcompound(info.getIntermediateSubcompoundSubcompoundId());
    //    }
    //}

    //Compound& updSubcompound(CompoundInfo& info) {
    //    if (info.isLocal()) {
    //        return info.updCompound();
    //    }
    //    else if (info.isBonded()) {
    //        return info.updCompound();
    //    }
    //    else { // subcompound of subcompound
    //        CompoundRep& sc1 = updSubcompound(info.getIntermediateSubcompoundId()).updImpl();
    //        return sc1.updSubcompound(info.getIntermediateSubcompoundSubcompoundId());
    //    }
    //}


    //// const Compound& getSubcompound(int) const;
    //// Compound& getSubcompound(int);

    //CompoundRep& nameSubcompound(Compound::Name newName, Compound::Name olderName) 
    //{
    //    assert(hasSubcompound(olderName));
    //    assert(!hasSubcompound(newName));
    //    Compound::Index subcompoundId = getSubcompoundInfo(olderName).getIndex();
    //    nameSubcompound(newName, subcompoundId);
    //    assert(hasSubcompound(newName));
    //    return *this;
    //}

    //CompoundRep& nameSubcompound(Compound::Name newName, Compound::Index subcompoundId) 
    //{
    //    assert(hasSubcompound(subcompoundId));
    //    assert(! hasSubcompound(newName));
    //    subcompoundIdsByName[newName] = subcompoundId;
    //    assert(hasSubcompound(newName));
    //    return *this;
    //}


    //bool hasSubcompound(const Compound::Name& name) const 
    //{
    //    // First check for simple subcompound name without any "/" separators
    //    if (CompoundPathName::isValidSubcompoundName(name)) {
    //        return subcompoundIdsByName.find(name) != subcompoundIdsByName.end();
    //    }
    //    else // parse "X/Y" path type subcompound names
    //    {
    //        std::vector<String> tokens;
    //        if (CompoundPathName::isValidSubcompoundPathName(name, &tokens)) 
    //        {
    //            String subcompoundName = tokens[0];
    //            if (! hasSubcompound(subcompoundName)) 
    //                return false;

    //            const CompoundRep& scRep = getSubcompound(subcompoundName).getImpl();
    //            return scRep.hasSubcompound(CompoundPathName::shiftLeftPathName(name));
    //        }
    //        else { // string not well formed
    //            return false;
    //        }
    //    }
    //}


    //bool hasSubcompound(const Compound::Index cId) const {
    //    if (cId < 0) return false;
    //    if (cId >= (Compound::Index)allSubcompounds.size()) return false;
    //    return true;
    //}

    static mdunits::Length getConsensusBondLength(const BondCenter& c1, const BondCenter& c2) 
    {
        mdunits::Length d1 = c1.getDefaultBondLength();
        mdunits::Length d2 = c2.getDefaultBondLength();

        mdunits::Length answer = d1;

        // Need to address the following cases:
        // 1) neither center has distance defined -> raise error
        // 2) one or other has distance defined -> OK, use that distance
        // 3) both have distance defined and agree -> OK, use that distance
        // 4) both have distance defined and disagree -> raise error

        // No information available to determine bond length
        if (isNaN(d1) && isNaN(d2)) assert(false); // case 1

        else if (isNaN(d1)) answer = d2; // case 2a
        else if (isNaN(d2)) answer = d1; // case 2b
        else if (d1 == d2) answer = d1; // case 3
        else assert(false); // case 4

        return answer;
    }

    static Angle getConsensusDihedralAngle(const BondCenter& c1, const BondCenter& c2) 
    {
        Angle a1 = c1.getDefaultDihedralAngle();
        Angle a2 = c2.getDefaultDihedralAngle();

        Angle answer = a1;

        // The same rules for dihedral angle as for bond length,
        // except that if both are NaN, default to 180 degrees
        if (isNaN(a1) && isNaN(a2)) answer = 180*Deg2Rad; // case 1
        else if (isNaN(a1)) answer = a2; // case 2a
        else if (isNaN(a2)) answer = a1; // case 2b
        else if (a1 == a2) answer = a1; // case 3
        else assert(false); // case 4

        return answer;
    }

    //Compound::BondIndex bondBondCenters(
    //    Compound::BondCenterIndex outboardId, 
    //    Compound::BondCenterIndex inboardId
    //    ) 
    //{
    //    mdunits::Length bondLength = getConsensusBondLength   (getBondCenter(outboardId), getBondCenter(inboardId));
    //    Angle    dihedral   = getConsensusDihedralAngle(getBondCenter(outboardId), getBondCenter(inboardId));

    //    return bondBondCenters(outboardId, inboardId, bondLength, dihedral);
    //}

    // TODO When a new bond is created, always call indexNewBond() to establish cross references
    //Compound::BondIndex bondBondCenters(
    //    Compound::BondCenterIndex outboardId, 
    //    Compound::BondCenterIndex inboardId,
    //    mdunits::Length               distance,
    //    Angle                  dihedral
    //    ) 

    void indexNewBond(const BondInfo& newBondInfo)
    {
        const Compound::BondCenterIndex outboardId = newBondInfo.getChildBondCenterIndex();
        const Compound::BondCenterIndex inboardId = newBondInfo.getParentBondCenterIndex();
        const Compound::BondIndex bondIndex = newBondInfo.getIndex();

        BondCenter&     outboardBondCenter      = updBondCenter(outboardId);
        BondCenter&     inboardBondCenter       = updBondCenter(inboardId);
        BondCenterInfo& outboardBondCenterInfo  = updBondCenterInfo(outboardId);
        BondCenterInfo& inboardBondCenterInfo   = updBondCenterInfo(inboardId);

        assert(! outboardBondCenter.isBonded() );
        assert(! inboardBondCenter.isBonded() );
        assert(! outboardBondCenterInfo.isBonded() );
        assert(! inboardBondCenterInfo.isBonded() );

        // Add a new BondInfo to this compound.
        // const Compound::BondIndex bondIndex = Compound::BondIndex(allBonds.size());

        // TODO - BondInfo constructor should depend on type of bond
        // allBonds.push_back( BondInfo(bondIndex, outboardId, inboardId, distance, dihedral) );

        // Mark the two newly-connected BondCenterInfos so they know they're connected.
        outboardBondCenterInfo.setBondPartnerBondCenterIndex(inboardBondCenterInfo.getIndex());
        inboardBondCenterInfo.setBondPartnerBondCenterIndex(outboardBondCenterInfo.getIndex());
        outboardBondCenterInfo.setBondIndex( bondIndex );
        inboardBondCenterInfo.setBondIndex( bondIndex );

        // Reach down to the Atoms and mark the physical bonds as in use, although they
        // don't know to whom they are connected.
        outboardBondCenter.setBonded(true);
        inboardBondCenter.setBonded(true);

        // Build a map so that we can find the connecting BondCenterInfo given the
        // AtomInfo indexes in either order.
        std::pair<Compound::AtomIndex, Compound::AtomIndex> 
            key1(outboardBondCenterInfo.getAtomIndex(), inboardBondCenterInfo.getAtomIndex());
        std::pair<Compound::AtomIndex, Compound::AtomIndex> 
            key2(inboardBondCenterInfo.getAtomIndex(), outboardBondCenterInfo.getAtomIndex());
        AIxPair_To_BondIx[key1] = bondIndex;
        AIxPair_To_BondIx[key2] = bondIndex;

        assert( outboardBondCenter.isBonded() );
        assert( inboardBondCenter.isBonded() );
        assert( outboardBondCenterInfo.isBonded() );
        assert( inboardBondCenterInfo.isBonded() );

        // return bondIndex;
    }

    //bool hasLocalSubcompound(const Compound::Name& name) const {
    //    return localSubcompoundIdsByName.find(name) != localSubcompoundIdsByName.end();
    //}

    //bool hasBondedSubcompound(const Compound::Name& name) const {
    //    return bondCenterIndexesByCompoundName.find(name) != bondCenterIndexesByCompoundName.end();
    //}

    CompoundRep& setBondMobility(BondMobility::Mobility mobility, const Compound::AtomName& atom1, const Compound::AtomName& atom2) 
    {
        AtomInfo& atomInfo1 = updAtomInfo(atom1);
        AtomInfo& atomInfo2 = updAtomInfo(atom2);
        BondInfo& bondInfo = updBondInfo(atomInfo1, atomInfo2);
        Bond& bond = updBond(bondInfo);

        bond.setMobility(mobility);
        // bond.setRotatable(isRotatable);
        
        return *this;
    }

    CompoundRep& setBondMobility(BondMobility::Mobility mobility, const Compound::BondIndex bondIndex) 
    {
        BondInfo& bondInfo = updBondInfo(bondIndex);
        Bond& bond = updBond(bondInfo);

        bond.setMobility(mobility);
        // bond.setRotatable(isRotatable);
        
        return *this;
    }


    CompoundRep& setPdbResidueNumber(int);
    int getPdbResidueNumber() const;

    CompoundRep& setPdbResidueName(const String&);
    const String& getPdbResidueName() const;

    CompoundRep& setPdbChainId(String);
    String getPdbChainId() const;


    // const Bond& getBond(const String& name) const;
    // Bond& getBond(const String& name);
    // const Bond& getBond(int id) const;
    // Bond& getBond(int id);
    // bool hasBond(const String& name) const;

     
    std::ostream& dumpCompoundRepToStream(std::ostream& o, int level=0) const;

    // bool hasParentCompound() const {return haveParentCompound;}


    CompoundRep& setTopLevelTransform(const Transform& transform) {
        // assert(!hasParentCompound());
        topLevelTransform = transform;

        return *this;
    }

    const Transform& getTopLevelTransform() const {
        // assert(!hasParentCompound());
        return topLevelTransform;
    }

protected:

    // offset + internal = nominal; offset + Bond = Dihedral
    Angle calcDefaultDihedralAngle(const DihedralAngle& dihedral) const
    {
        const BondCenterInfo& bc1 = getBondCenterInfo(dihedral.getBondCenter1Id());
        const BondCenterInfo& bc2 = getBondCenterInfo(dihedral.getBondCenter2Id());
        const AtomInfo& atom1 = getAtomInfo(bc1.getAtomIndex());
        const AtomInfo& atom2 = getAtomInfo(bc2.getAtomIndex());
        assert( atomsAreBonded(atom1, atom2) );
        const BondInfo& bondInfo = getBondInfo(atom1, atom2);
        const Bond& bond = getBond(bondInfo);

        Angle internalAngle = bond.getDefaultDihedralAngle();

        Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(bc1.getIndex(), bc2.getIndex());

        Angle nominalAngle = internalAngle + internalOffset + dihedral.getNomenclatureOffset();

        return nominalAngle;
    }


    Angle calcDihedralAngle(const State& state, const String& dihedralName) const
    {
        assert( AtomName_To_dihedralAngles.find(dihedralName) != AtomName_To_dihedralAngles.end() );

        const DihedralAngle& dihedral = AtomName_To_dihedralAngles.find(dihedralName)->second;

        return calcDihedralAngle(state, dihedral);
    }

    Transform calcAtomFrameInGroundFrame(const State& state, Compound::AtomIndex atomId) const
    {
        const CompoundAtom& atom = getAtom(atomId);

        // Frame of parent body
        //DuMM::AtomIndex dummAtomIndex = atom.getDuMMAtomIndex();
        //const DuMMForceFieldSubsystem& dumm = ownerSystem->getMolecularMechanicsForceSubsystem();

        MobilizedBodyIndex bodyId = atom.getMobilizedBodyIndex();

        const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();
        const MobilizedBody& body = matter.getMobilizedBody(bodyId);
        const Transform& G_X_B = body.getBodyTransform(state);

        // Frame of atom in body;
        Transform B_X_A = atom.getFrameInMobilizedBodyFrame();

        return G_X_B * B_X_A;
    }
    
    Angle calcDihedralAngle(const State& state, const DihedralAngle& dihedral) const
    {
        const Bond& bond = getBondByDihedral(dihedral);
        MobilizedBodyIndex bodyId = bond.getPinJointId();
        if (bodyId.isValid())
        {
            assert(ownerSystem != NULL);
            const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();

            Angle internalAngle = 0;
            if(bond.getMobility() == BondMobility::Torsion) {
                const MobilizedBody::Pin &body = (const MobilizedBody::Pin &) matter.getMobilizedBody(bodyId);
                internalAngle = body.getAngle(state);
            }else if(bond.getMobility() == BondMobility::BallF){// Gmol
                const MobilizedBody::Ball &ball = (const MobilizedBody::Ball &) matter.getMobilizedBody(bodyId);
                // Return psi Euler angle and ignore phi and theta
                Vec4 q = SimTK::Quaternion(ball.getQ(state));
                double psi;
                // Deal with singularity
                if( (std::abs((q[1] * q[2]) + (q[3] * q[0])) - 0.5) < 0.01 ){
                    psi = 0.0;
                }else {
                    double q0q3 = q[0] * q[3];
                    double q1q2 = q[1] * q[2];
                    double q2sq = q[2] * q[2];
                    double q3sq = q[3] * q[3];
                    psi = atan2(2 * (q0q3 + q1q2),
                                1 - 2 * (q2sq + q3sq));
                    internalAngle = psi;
                }
            } else {
                // Should never get here, but compiler keeps complaining.
                assert(false);
            }

            // TODO - use simtime offset, not default

            Angle internalOffset = calcDefaultInternalDihedralOffsetAngle(dihedral.getBondCenter1Id(), dihedral.getBondCenter2Id());

            Angle nominalAngle = internalAngle + internalOffset + dihedral.getNomenclatureOffset();

            return nominalAngle;
        }
        else
        {
            assert(false); // TODO

            assert(ownerSystem != NULL);
            // Requires that system be realized to Position stage
            ownerSystem->realize(state, Stage::Position);

            const BondCenterInfo& bc21 = getBondCenterInfo(dihedral.getBondCenter1Id());
            const BondCenterInfo& bc34 = getBondCenterInfo(dihedral.getBondCenter2Id());

            // Find bond axis to project onto
            const AtomInfo& atom2 = getAtomInfo(bc21.getAtomIndex());
            const AtomInfo& atom3 = getAtomInfo(bc34.getAtomIndex());
            const BondCenterInfo& bondBondCenter = getBondCenterInfo(atom2, atom3);

            // Debug TODO - remove this stanza
            const Bond& bond = getBond(getBondInfo(bondBondCenter.getBondIndex()));
            MobilizedBodyIndex bodyId = bond.getPinJointId();
            std::cout << "Pin joint id = " << bodyId;
            std::cout << "\t";
            const SimbodyMatterSubsystem& matter = ownerSystem->getMatterSubsystem();

            if(bond.getMobility() == BondMobility::Torsion) {
                const MobilizedBody::Pin &body = (const MobilizedBody::Pin &) matter.getMobilizedBody(bodyId);
                std::cout << "angle = " << body.getAngle(state) * DuMM::Rad2Deg << " degrees";
                std::cout << std::endl;
            }else if(bond.getMobility() == BondMobility::BallF) { // Gmol
                const MobilizedBody::Ball &body = (const MobilizedBody::Ball &) matter.getMobilizedBody(bodyId);
                std::cout << "angle = ";

                // Gregory G. Slabaugh description
                SimTK::Rotation R;
                R.setRotationFromQuaternion(Quaternion(body.getQ(state)));
                Angle theta1 = -1.0 * std::asin(R[2][0]);
                //Angle theta2 = SimTK::Pi - theta1;
                double cosTheta1 = std::cos(theta1);
                //double cosTheta2 = std::cos(theta2);
                Angle psi1 =  std::atan2(R[2][1] / cosTheta1, R[2][2] / cosTheta1);
                //Angle psi2 =  std::atan2(R[2][1] / cosTheta2, R[2][2] / cosTheta2);

                std::cout << psi1 * DuMM::Rad2Deg;
                std::cout << " degrees";
                std::cout << std::endl;
            }

            assert(bondBondCenter.isBonded());
            assert(bondBondCenter.getIndex() != bc21.getIndex());
            assert(bondBondCenter.getIndex() != bc34.getIndex());
            assert(bc21.getIndex() != bc34.getIndex());

            UnitVec3 xAxis(1,0,0);

            // vector v1: from atom 1 to atom 2
            Transform G_X_A2 = calcDefaultAtomFrameInCompoundFrame(atom2.getIndex());
            Transform A2_X_BC21 = calcDefaultBondCenterFrameInAtomFrame(bc21);
            Transform G_X_BC21 = G_X_A2 * A2_X_BC21;
            UnitVec3 v1(G_X_BC21 * -xAxis); // negative x-axis because want 1->2, not 2->1 vector

            // vector v2: from atom 2 to atom 3
            Transform A2_X_BCB = calcDefaultBondCenterFrameInAtomFrame(bondBondCenter);
            Transform G_X_BCB = G_X_A2 * A2_X_BCB;
            UnitVec3 v2(G_X_BCB * xAxis);

            // vector v3: from atom 3 to atom 4
            Transform G_X_A3 = calcDefaultAtomFrameInCompoundFrame(atom3.getIndex());
            Transform A3_X_BC34 = calcDefaultBondCenterFrameInAtomFrame(bc34);
            Transform G_X_BC34 = G_X_A3 * A3_X_BC34;
            UnitVec3 v3(G_X_BC34 * xAxis);

            Angle nominalDihedralAngle = SimTK::calcDihedralAngle(v1, v2, v3) + dihedral.getNomenclatureOffset();

            return nominalDihedralAngle;
        }
    }

    // synonyms are alternate names for this compound type
    // synonyms should include the primary name of the compound
    // one use of synonyms is to help resolve biotypes atoms found in tinker parameter files
    std::set<Compound::Name> synonyms;

    // AtomInfo references
    std::vector<AtomInfo>                          allAtoms;    // [Compound::AtomIndex]
    // std::vector<CompoundAtom>                              localAtoms;  // [Compound::LocalAtomIndex]
    std::map<Compound::AtomName, Compound::AtomIndex> atomName_To_atomId;

    // bool haveParentCompound;

private:
    friend class Compound;
    friend class Bond;

    // ownerSystem is being used in two ways:
    // 1) ownerSystem plus ixWithinOwnerSystem represent handle for compounds directly owned by a CompoundSystem
    // 2) ownerSystem with invalid ixWithinOwnerSystem represents a handle to the system for subcompounds of
    //    a compound that is in turn directly owned by a CompoundSystem
    MultibodySystem* ownerSystem;
    // CompoundSystem* ownerSystem;
    // Compound::Index          ixWithinOwnerSystem;

    // Global transform that applies only to top level compound
    Transform topLevelTransform;

    Compound::Name name; // set on construction; means whatever you like

    // The following comment may be wrong cmb Feb 2009
    // local subcompounds placed directly - NOT those placed by bonds
    // std::vector<CompoundInfo> allSubcompounds; // [Compound::Index]
    // std::vector<Compound> localSubcompounds;
    // std::map<Compound::Name, Compound::Index> subcompoundIdsByName;
    // std::map<int, Transform> localSubcompoundTransformsById;

    // BondCenters
    std::vector<BondCenterInfo>              allBondCenters; // [Compound::BondCenterIndex]
    std::map<String, Compound::BondCenterIndex> BCName_To_BCIx;
    BondCenterInfo::AtomKeyMap               bondCenterIndicesByAtomKey;

    // Bonds
    // bonds do not contain subcompounds
    std::vector<BondInfo> allBonds; // [Compound::BondIndex]
    std::map< std::pair<Compound::AtomIndex, Compound::AtomIndex>, Compound::BondIndex > 
                          AIxPair_To_BondIx;

    // Dihedral Angles
    std::map<String, DihedralAngle> AtomName_To_dihedralAngles;

    int    pdbResidueNumber;
    String pdbResidueName;
    String   pdbChainId;
    
    class MemberForDebuggingCopyCtor {
    public:
        MemberForDebuggingCopyCtor() {}
        MemberForDebuggingCopyCtor(const MemberForDebuggingCopyCtor&) {
            // Put a breakpoint here to notice when CompoundRep copyCtor is called
            //int x = 5;
        }
    };
    // MemberForDebuggingCopyCtor testMember;
};


BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const Compound::Name& resName, 
                        const Compound::AtomName& atomName, 
                        Ordinality::Residue ordinality = Ordinality::Any);

// Return a biotype index matching any of the residue/atom names supplied
// Returns invalid index if not found
BiotypeIndex SimTK_MOLMODEL_EXPORT getBiotypeIndex(
                        const std::set<Compound::Name>& resNames, 
                        const std::set<Compound::AtomName>& atomNames, 
                        Ordinality::Residue ordinality = Ordinality::Any);

class BiopolymerResidueRep : public CompoundRep {
public:
    /*virtual*/ ~BiopolymerResidueRep() { }
    /*virtual*/ BiopolymerResidueRep* clone() const {return new BiopolymerResidueRep(*this);}

    BiopolymerResidueRep(String name, String tlc = "Unk", char olc = '?')
        : residueName(name), threeLetterCode(tlc), oneLetterCode(olc)
        {}

    BiopolymerResidueRep& setOneLetterCode(char olc) {
        oneLetterCode = olc;
        return *this;
    }
    BiopolymerResidueRep& setThreeLetterCode(const String& tlc) {
        threeLetterCode = tlc;
        return *this;
    }
    BiopolymerResidueRep& setResidueTypeName(const String& name) {
        residueName = name;
        return *this;
    }

    const String& getResidueTypeName() const {return residueName;}
    const String& getThreeLetterCode() const {return threeLetterCode;}
    char getOneLetterCode() const {return oneLetterCode;}

    // Attempt to deduce correct biotypes from global biotypes database
    /// @return true if all atoms have a valid biotype assigned, false otherwise
    bool assignBiotypes(Ordinality::Residue ordinality = Ordinality::Any) 
    {
        bool answer = true; // start optimistic

        // for each named atom, look up resname, atomname, ordinality
        // if biotype is still undefined, raise exception
        std::vector<AtomInfo>::iterator atomI;
        for (atomI = allAtoms.begin(); atomI != allAtoms.end(); ++atomI) 
        {
            // debugging
            // int atomIndex = atomI->getIndex();
            // std::cout << "biotype for atom " << atomI->getIndex() << std::endl;
            
            CompoundAtom& atom = updAtom(*atomI);
            
			// if the atom already has a valid biotype, keep it.
			if (atom.getBiotypeIndex().isValid()) continue;

            // Examine all possible residue names
            const std::set<Compound::Name>& residueNames = synonyms;

            // Create a container to hold variations of atom name
            const std::set<Compound::AtomName>& atomNames = getAtomSynonyms(atomI->getIndex());

            // Loop over residue names and atom names until a match is found
            bool foundBiotype = false;
            BiotypeIndex index = getBiotypeIndex(residueNames, atomNames, ordinality);
            if (index.isValid()) {
                atom.setBiotypeIndex(index);
                foundBiotype = true;
            }
            else {
                foundBiotype = false;
            }

            if (!foundBiotype) {
                answer = false;
                std::set<Compound::AtomName>::const_iterator atomNamesIterator; 
                for (atomNamesIterator = atomNames.begin(); atomNamesIterator != atomNames.end(); atomNamesIterator++){
                    //std::cout<<__FILE__<<":"<<__LINE__<<" atomNames = "<<string(*atomNamesIterator)<<std::endl;
                }
            }
            // Perhaps the atom already had a usable biotype...
            /// assert(atom.getBiotypeIndex().isValid());
            assert(foundBiotype);

        }

        return answer;
    }

private:
    String residueName;
    String threeLetterCode;
    char   oneLetterCode;
};



class BiopolymerRep : public CompoundRep {
public:
    friend class Biopolymer;

    BiopolymerRep* clone() const {return new BiopolymerRep(*this);}
    //const std::vector<String>& getResidueNames() const {
    //    return residueNames;
    //}

    //std::vector<String>& updResidueNames() {
    //    return residueNames;
    //}

/*BiopolymerRep& fitDefaultConfiguration(
        const Compound::AtomTargetLocations& atomTargets,
        SimTK::Real targetRms,
        bool useObservedPointFitter,
        Real minimizerTolerance//,
        //Compound compoundCopy //= *this;
        )
{
    // this is a pointer, *this is its value
    Compound compoundCopy((*this));
    // TODO - DuMM should not be required
    CompoundSystem matchingSystem;
    SimbodyMatterSubsystem matchingMatter(matchingSystem);
    DuMMForceFieldSubsystem dumm(matchingSystem);
    dumm.loadAmber99Parameters();
    dumm.setAllGlobalScaleFactors(0);
    GeneralForceSubsystem forces(matchingSystem);
    matchingSystem.adoptCompound(compoundCopy);
    matchingSystem.modelCompounds();
    matchingSystem.realizeTopology();
    State& state = matchingSystem.updDefaultState();
    matchingSystem.realize(state, Stage::Position);
    // cout << "Number of atom matches(2) = " << optimizationAtomTargets.size() << endl;
    std::map<MobilizedBodyIndex, std::vector<Vec3> > stations;
    std::map<MobilizedBodyIndex, std::vector<Vec3> > targetLocations;
    for (Compound::AtomTargetLocations::const_iterator targetIx = atomTargets.begin();
         targetIx != atomTargets.end();
         ++targetIx)
    {
        Compound::AtomIndex atomId = targetIx->first;
        MobilizedBodyIndex bodyId = compoundCopy.getAtomMobilizedBodyIndex(atomId);
        stations[bodyId].push_back(compoundCopy.getAtomLocationInMobilizedBodyFrame(atomId));
        targetLocations[bodyId].push_back(targetIx->second);
    }
   
    // Use ObservedPointFitter to optimize geometry
    std::vector<MobilizedBodyIndex> bodyList;
    std::vector<std::vector<Vec3> > stationList;
    std::vector<std::vector<Vec3> > targetList;
    for (std::map<MobilizedBodyIndex, std::vector<Vec3> >::const_iterator iter = stations.begin(); iter != stations.end(); iter++) {
        bodyList.push_back(iter->first);
        stationList.push_back(iter->second);
        targetList.push_back(targetLocations.find(iter->first)->second);
    }


    // ObservedPointFitter takes a while, and occasionally aborts with line search trouble,
    // So lets try a minimization using custom forces
    //bool useObservedPointFitter = true;
    if (useObservedPointFitter) {
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
        ObservedPointFitter::findBestFit(matchingSystem, state, bodyList, stationList, targetList, tolerance);
    }
    else {
        const MobilizedBody& groundBody = matchingMatter.getGround();
        for (int b = 0; b < (int)bodyList.size(); ++b)
        {
            const MobilizedBody& atomBody = matchingMatter.getMobilizedBody(bodyList[b]);
            for (int s = 0; s < (int)stationList[b].size(); ++s)
            {
                const Vec3& atomLocation = stationList[b][s];
                const Vec3& targetLocation = targetList[b][s];
                Force::TwoPointLinearSpring(forces, atomBody, atomLocation, groundBody, targetLocation, 1000000.0, 0.0);
            }
        }

        state = matchingSystem.realizeTopology();
        matchingSystem.realize(state, Stage::Position);
        matchingSystem.realize(state, Stage::Position);
            LocalEnergyMinimizer::minimizeEnergy(matchingSystem, state,  minimizerTolerance);

            // Stuff optimized coordinates into a string
            std::ostringstream optimizedPdbStringOut;
            matchingSystem.realize(state, Stage::Position);
            compoundCopy.writePdb(state, optimizedPdbStringOut);
            // Create another PdbStructure, and match the dihedral angles to that
            std::istringstream optimizedPdbStringIn(optimizedPdbStringOut.str());
            PdbStructure optimizedStructure(optimizedPdbStringIn);
            //Compound::AtomTargetLocations optimizedAtomTargets =
            //        createAtomTargets(optimizedStructure,false); // scf set guessCoordinates to false here to make sure it's done exactly as before
        std::ofstream myofstream("match1e.pdb");
        optimizedStructure.write(myofstream, SimTK::Transform(Vec3(0)));
    }
    // Stuff optimized coordinates into a string
    std::ostringstream optimizedPdbStringOut;
    matchingSystem.realize(state, Stage::Position);
    compoundCopy.writePdb(state, optimizedPdbStringOut);
    // Create another PdbStructure, and match the dihedral angles to that
    std::istringstream optimizedPdbStringIn(optimizedPdbStringOut.str());
    PdbStructure optimizedStructure(optimizedPdbStringIn);
    Compound::AtomTargetLocations optimizedAtomTargets =
            createAtomTargets(optimizedStructure,false); // scf set guessCoordinates to false here to make sure it's done exactly as before
    bool matchHydrogenAtomLocations = false;
    if (! matchHydrogenAtomLocations) {
        std::map<Compound::AtomIndex, Vec3>::iterator it;
        std::map<Compound::AtomIndex, Vec3>::iterator next;
        next = optimizedAtomTargets.begin();
        while (next != optimizedAtomTargets.end()) {
            it = next;
            Compound::AtomIndex m = (*it).first;
            Element myAtomElement = getAtomElement(m);
            next++;
            if  ((myAtomElement.getName()).compare("hydrogen") == 0) {
                //cout<<__FILE__<<":"<<__LINE__<<" erasing "<<m<<endl;
                optimizedAtomTargets.erase(it);
            }
            //cout<<__FILE__<<":"<<__LINE__<<" "<<m<<","<<(getAtomName(m))<<endl;
        }
    }
    
    //not sure that this was ever needed.  In any event, cutting out from the BiopolymerRep version:
    matchDefaultBondLengths(optimizedAtomTargets);
    matchDefaultBondAngles(optimizedAtomTargets);
    matchDefaultDihedralAngles(optimizedAtomTargets);
    // Use original atom locations for top level transform
    matchDefaultTopLevelTransform(optimizedAtomTargets);
    
    return *this;
} */



    int getNumResidues() const;
    const String& getResidueName(int residueIndex) const;

    const ResidueInfo& getResidue(ResidueInfo::Index residueIndex) const;
    ResidueInfo& updResidue(ResidueInfo::Index residueIndex);

    const ResidueInfo& getResidue(const Compound::Name& residueName) const;
    ResidueInfo& updResidue(const Compound::Name& residueName);

    Transform calcDefaultResidueFrameInBiopolymerFrame(ResidueInfo::Index r, const std::vector<Transform>& atomFrameCache) const {
        // frame of inboard atom
        Compound::AtomIndex inboardAtomIndex = getResidue(r).getAtomIndex(ResidueInfo::AtomIndex(0));
        return atomFrameCache[inboardAtomIndex];
        // return calcDefaultAtomFrameInCompoundFrame(inboardAtomIndex);
    }

    virtual const CompoundRep& populateResidueDefaultPdbChain(
        ResidueInfo::Index r,
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform,
        const std::vector<Transform>& atomFrameCache) const 
    {
        //Transform myTransform = transform * calcDefaultResidueFrameInBiopolymerFrame(r, atomFrameCache);
        const ResidueInfo& residue = getResidue(r);

        int residueNumber = residue.getPdbResidueNumber();
        char insertionCode = residue.getPdbInsertionCode(); // For some reason, insertion code was not used before. SCF: correcting this now.
        if (-9999 > residueNumber) {
            residueNumber = defaultNextResidueNumber;
        }
        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber,insertionCode))) {
            std::cout<<__FILE__<<":"<<__LINE__<<" The residue ID passed on by the compound, "<< residueNumber<<insertionCode<<" conflicts with an existing residue ID in pdbChain.  In the past we would invent a new residueNumber, but this probably indicates a deeper underlying problem. "    <<std::endl; exit(1);
            ++residueNumber;
            insertionCode = ' '; // make sure insertion code goes back to the default of ' '.
        }
        PdbResidue pdbResidue(residue.getPdbResidueName(), PdbResidueId(residueNumber,insertionCode));
        for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) 
        {
            Compound::AtomIndex atomIx = residue.getAtomIndex(a);

            PdbAtom pdbAtom(residue.getAtomName(a), getAtomElement(atomIx));
            pdbAtom.setLocation(PdbAtomLocation(transform * atomFrameCache[atomIx].p()));
            pdbResidue.addAtom(pdbAtom);
        }

        pdbChain.appendResidue(pdbResidue);

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }

    /// New way to do PDB writing: create intermediate PdbChain object
    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populateDefaultPdbChain(
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        std::vector<Transform> atomFrameCache(getNumAtoms());
        invalidateAtomFrameCache(atomFrameCache, getNumAtoms());
        calcDefaultAtomFramesInCompoundFrame(atomFrameCache);

        for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) {
            populateResidueDefaultPdbChain(r, pdbChain, defaultNextResidueNumber, transform, atomFrameCache);
        }

        return *this;
    }

    /// Write current default(initial) Compound configuration into a PdbChain object
    virtual const CompoundRep& populatePdbChain(
        const State& state, 
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        for (ResidueInfo::Index r(0); r < getNumResidues(); ++r) 
            populateResiduePdbChain(state, r, pdbChain, defaultNextResidueNumber, transform);

        return *this;
    }

    const CompoundRep& populateResiduePdbChain(
        const State& state,
        ResidueInfo::Index r,
        class PdbChain& pdbChain, 
        int& defaultNextResidueNumber,
        const Transform& transform) const 
    {
        // Transform myTransform = transform * calcResidueFrameInBiopolymerFrame(r);
        const ResidueInfo& residue = getResidue(r);

        int residueNumber = residue.getPdbResidueNumber();
        int insertionCode = residue.getPdbInsertionCode();
        if (-9999 > residueNumber) residueNumber = defaultNextResidueNumber;
        // In case of residue number conflicts, find a new number
        if (pdbChain.hasResidue(PdbResidueId(residueNumber,insertionCode)))
            residueNumber = defaultNextResidueNumber;
        while (pdbChain.hasResidue(PdbResidueId(residueNumber,insertionCode)))
            ++residueNumber;

        PdbResidue pdbResidue(residue.getPdbResidueName(), PdbResidueId(residueNumber,insertionCode ));
        for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a) 
        {
            Compound::AtomIndex atomIx = residue.getAtomIndex(a);
            //const CompoundAtom& atom = getAtom(atomIx);

            PdbAtom pdbAtom(residue.getAtomName(a), getAtomElement(atomIx));
            pdbAtom.setLocation(
                PdbAtomLocation(transform * calcAtomLocationInGroundFrame(state, atomIx))
            );
            pdbResidue.addAtom(pdbAtom);
        }
        pdbChain.appendResidue(pdbResidue);

        defaultNextResidueNumber = residueNumber + 1;

        return *this;
    }


    //virtual std::ostream& writeDefaultPdb(
    //    std::ostream& os, 
    //    int& nextSerialNumber, 
    //    const Transform& transform) const;

    //virtual std::ostream& writePdb(
    //    const State& state, 
    //    std::ostream& os, 
    //    int& nextSerialNumber,
    //    const Transform& transform) const;

private:
    // std::vector<String> residueNames;
    std::vector<ResidueInfo> residues;
    std::map<Compound::Name, ResidueInfo::Index> residueIdsByName;
};

} // namespace SimTK

#endif // SimTK_COMPOUNDREP_H_
