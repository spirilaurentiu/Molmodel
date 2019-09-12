#include "molmodel/internal/CompoundSystem.h"

#include "CompoundRep.h"
#include "SimTKmath.h"

#include "molmodel/internal/RiboseMobilizer.h"

#include <set>

using namespace std;

namespace SimTK {

void SimTK::CompoundSystem::modelCompounds(String mobilizedBodyType) 
{
    // Turn off default decorations, since we'll make our own decorations.
    updMatterSubsystem().setShowDefaultGeometry(false); 
    for (CompoundSystem::CompoundIndex c(0); c < getNumCompounds(); ++c) 
        modelOneCompound(c, mobilizedBodyType);
}

// RigidUnit data structure is for use in modelCompounds() method
class RigidUnit {
public:
    RigidUnit() {}
    RigidUnit(DuMM::ClusterIndex id) 
        : clusterIx(id), hasChild(false) {}

    DuMM::ClusterIndex clusterIx; // primary key
    MobilizedBodyIndex bodyId; // populated toward the end
    DuMM::ClusterIndex parentId; // InvalidId implies parented to Ground
    bool hasChild; // Whether a child body is tethered to this one

    Transform frameInTopCompoundFrame; // useful intermediate computation
    Transform frameInParentFrame; // what we ultimately want

    Compound::BondCenterIndex inboardBondCenterIndex;
    Angle inboardBondDihedralAngle;
    Compound::BondIndex inboardBondIndex;

    std::set<Compound::AtomIndex> clusterAtoms;
};

// AtomBonding data structure represents one Atom in the modelCompounds() method
class AtomBonding {
public:
    AtomBonding() {}
    AtomBonding(Compound::AtomIndex id)
        : atomId(id) {}

    Compound::AtomIndex atomId; // primary key
    Compound::AtomIndex parentAtomIndex; // for finding rootiest atom in cluster
    DuMM::AtomIndex dummAtomIndex;
    DuMM::ClusterIndex clusterIx;

    // Store only those bonds that are part of the multibody tree structure
    std::set<Compound::AtomIndex> treeBonds;
    std::set<Compound::AtomIndex> freeTreeBonds;

    Vec3 locationInBodyFrame;
};

// Recursive method to create rigid bodies from a seed atom
void buildUpRigidBody(Compound::AtomIndex atomId, 
                      DuMM::ClusterIndex clusterIx,
                      std::set<Compound::AtomIndex>& clusterAtoms, 
                      std::map<Compound::AtomIndex, AtomBonding>& atomBondings,
                      CompoundRep& compoundRep
                      ) 
{
    clusterAtoms.insert(atomId);

    // Store DuMM::Cluster id in Compound::Atom object
    CompoundAtom& atom = compoundRep.updAtom(atomId);
    assert(!atom.getDuMMPrimaryClusterIndex().isValid());
    atom.setDuMMPrimaryClusterIndex(clusterIx);

    assert(atomBondings.find(atomId) != atomBondings.end());
    AtomBonding& atomBonding = atomBondings.find(atomId)->second;

    atomBonding.clusterIx = clusterIx;
    assert(atomBonding.clusterIx.isValid());


    // Unbonded atoms always represent a lone rigid body
    int numberOfNonBendyBondsOnThisAtom = atomBonding.treeBonds.size() - atomBonding.freeTreeBonds.size();
    if (numberOfNonBendyBondsOnThisAtom == 0) return;

    // recursively check all atoms bonded to this one for membership in group
    std::set<Compound::AtomIndex>::const_iterator b;
    for (b = atomBonding.treeBonds.begin(); b != atomBonding.treeBonds.end(); ++b)
    {
        // Ignore atoms already assigned to a rigid cluster
        if (clusterAtoms.find(*b) != clusterAtoms.end()) 
             continue; // already assigned

        // Bond archaeology
        const AtomInfo& atomInfo1 = compoundRep.getAtomInfo(atomId);
        const AtomInfo& atomInfo2 = compoundRep.getAtomInfo(*b);
        const BondInfo& bondInfo = compoundRep.getBondInfo(atomInfo1, atomInfo2);
        const Bond& bond = compoundRep.getBond(bondInfo);

        if (bond.getMobility() == BondMobility::Free) 
        {
            // continue;  // cannot restrict rigid body if free
        }
        else if (bond.getMobility() == BondMobility::Torsion)
        {
            // TORSION-only bonds to atoms with just one bond are included in rigid body
            // for inertian reasons

            const AtomBonding& bondedAtomBonding = atomBondings.find(*b)->second;
            int numberOfNonBendyBondsOnBondedAtom = bondedAtomBonding.treeBonds.size() - bondedAtomBonding.freeTreeBonds.size();
            assert(numberOfNonBendyBondsOnBondedAtom > 0);

            // TODO - Also join atoms with just one bond among Torsion or Rigid (or possibly Angle) bond
            
            // atoms with just one bond are not rotatable for inertia reasons
            // and thus get included in rigid body
            if ( (numberOfNonBendyBondsOnBondedAtom == 1) || (numberOfNonBendyBondsOnThisAtom == 1) )
            {
                buildUpRigidBody(*b, clusterIx, clusterAtoms, atomBondings, compoundRep);
            }
        }

        else if (bond.getMobility() == BondMobility::Rigid) // definitely add to rigid body
        {
            buildUpRigidBody(*b, clusterIx, clusterAtoms, atomBondings, compoundRep);
        }
        else assert(false); // unexpected mobility
    }

}


void CompoundSystem::modelOneCompound(CompoundIndex compoundId, String mobilizedBodyType ) 
{
    bool showDebugMessages = false;
    if (showDebugMessages) cout << "modelOneCompound" << endl;

    // Turn off default decorations, since we'll make our own decorations.
    updMatterSubsystem().setShowDefaultGeometry(false); 

    Compound& compound = updCompound(compoundId);
    CompoundRep& compoundRep = compound.updImpl();

    // const Transform& compoundTransform = compound.getTopLevelTransform();
    
    // Cache default atom frames for performance
    std::vector<Transform> defaultAtomFrames(compoundRep.getNumAtoms());
    invalidateAtomFrameCache(defaultAtomFrames, compoundRep.getNumAtoms());
    compoundRep.calcDefaultAtomFramesInCompoundFrame(defaultAtomFrames);

    DuMMForceFieldSubsystem& dumm = (DuMMForceFieldSubsystem&) updMolecularMechanicsForceSubsystem();
    SimbodyMatterSubsystem&  matter  = updMatterSubsystem();

    std::map<DuMM::ClusterIndex, RigidUnit> rigidUnits;
    std::map<Compound::AtomIndex, AtomBonding> atomBonds;

    if (showDebugMessages) cout << "Step 1 create atomBonds" << endl;
    // 1) Create initial atomBonds data structure for each atom (no bonds are cached in this loop)
    for (Compound::AtomIndex a(0); a < compound.getNumAtoms(); ++a) 
    {
        // Create AtomBonding object
        atomBonds[a] = AtomBonding(a);
        AtomBonding& atomBonding = atomBonds[a];

        // Assign DuMM::AtomIndex for linking to simbody
        BiotypeIndex biotypeIx = compound.getAtomBiotypeIndex(a);
        assert(biotypeIx.isValid());
        DuMM::ChargedAtomTypeIndex chargedTypeId = dumm.getBiotypeChargedAtomType(biotypeIx);
        assert(chargedTypeId.isValid());
        atomBonding.dummAtomIndex = dumm.addAtom(chargedTypeId);
        assert(atomBonding.dummAtomIndex.isValid());

        // Store DuMMAtomIndex in Compound::Atom object
        CompoundAtom& atom = compoundRep.updAtom(a);

        assert(!atom.getDuMMPrimaryClusterIndex().isValid());
        assert(!atom.getDuMMAtomIndex().isValid());

        atom.setDuMMAtomIndex(atomBonding.dummAtomIndex);
    }

    if (showDebugMessages) cout << "Step 2 analyze bonding structure" << endl;
    // 2) Analyze atom bonding structure, and populate dumm bond structure
    for (Compound::BondIndex i(0); i < compound.getNumBonds(); ++i) 
    {
        const BondInfo& bondInfo = compoundRep.getBondInfo(i);
        const BondCenterInfo& parentBondCenterInfo = compoundRep.getBondCenterInfo(bondInfo.getParentBondCenterIndex());
        const BondCenterInfo& childBondCenterInfo = compoundRep.getBondCenterInfo(bondInfo.getChildBondCenterIndex());
        const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentBondCenterInfo.getAtomIndex());
        const AtomInfo& childAtomInfo = compoundRep.getAtomInfo(childBondCenterInfo.getAtomIndex());
        const Bond& bond = compoundRep.getBond(bondInfo);

        Compound::AtomIndex a0 = parentAtomInfo.getIndex();
        Compound::AtomIndex a1 = childAtomInfo.getIndex();
        assert(a0 != a1);

        assert(atomBonds.find(a0) != atomBonds.end());
        assert(atomBonds.find(a1) != atomBonds.end());

        // Tell subsystem about the bond (including ring closing bonds)
        DuMM::BondIndex dummId = dumm.addBond(atomBonds[a0].dummAtomIndex, atomBonds[a1].dummAtomIndex);

        // store bond info on each atom
        // only store those bonds that are part of the tree structure
        if (! bond.isRingClosingBond())
        {
            atomBonds[a0].treeBonds.insert(a1);
            atomBonds[a1].treeBonds.insert(a0);
            
            if (bond.getMobility() == BondMobility::Free) {
                atomBonds[a0].freeTreeBonds.insert(a1);
                atomBonds[a1].freeTreeBonds.insert(a0);            	
            }
            
            // store parent-child relationship
            atomBonds[a1].parentAtomIndex = a0;
        }
    }

    if (showDebugMessages) cout << "Step 3 distribute atoms" << endl;
    // 3) Distribute atoms to bodies, using DuMM::ClusterIndex as proxy for body for the present
    std::map<Compound::AtomIndex, AtomBonding>::iterator atomBondI;
    for (atomBondI = atomBonds.begin(); atomBondI != atomBonds.end(); ++atomBondI) 
    {
        AtomBonding& atomBonding = atomBondI->second;

        // Ignore atoms already in a cluster
        if (atomBonding.clusterIx.isValid()) continue;

        Compound::AtomIndex atomId = atomBonding.atomId;

        // Start a new body
        // Create new clusterIx as primary key for new rigid body
        DuMM::ClusterIndex clusterIx = dumm.createCluster(String(atomId));
        assert( rigidUnits.find(clusterIx) == rigidUnits.end() );
        rigidUnits[clusterIx] = RigidUnit(clusterIx);
        RigidUnit& rigidUnit = rigidUnits[clusterIx];
        assert( rigidUnits.find(clusterIx) != rigidUnits.end() );

        // use recursive method to find all of the atoms in the cluster
        buildUpRigidBody(atomId, clusterIx, rigidUnit.clusterAtoms, atomBonds, compoundRep);
    }


    if (showDebugMessages) cout << "Step 4 assign rigid body parents" << endl;
    // 4) Assign rigid body parents and set body frames relative to top compound
    for (Compound::BondIndex i(0); i < compound.getNumBonds(); ++i) 
    {
        const BondInfo& bondInfo = compoundRep.getBondInfo(i);

        // Don't use ring closing bonds to assign parent child relationships
        // because we only want to use the tree structure of bonds to make a multibody tree structure
        const Bond& bond = compoundRep.getBond(bondInfo);

        if (bond.isRingClosingBond()) // ring closing bonds cannot be part of tree structure
            continue;

        switch (bond.getMobility()) {
            case BondMobility::Rigid:
                break; // same body => ignore
            case BondMobility::Free:
                break; // no multibody parent/child relationship => ignore
            case BondMobility::Torsion:
                {
                    // This might represent a parent/child relationship
                    const BondCenterInfo& parentBondCenterInfo = compoundRep.getBondCenterInfo(bondInfo.getParentBondCenterIndex());
                    const BondCenterInfo& childBondCenterInfo = compoundRep.getBondCenterInfo(bondInfo.getChildBondCenterIndex());
                    const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentBondCenterInfo.getAtomIndex());
                    const AtomInfo& childAtomInfo = compoundRep.getAtomInfo(childBondCenterInfo.getAtomIndex());

                    Compound::AtomIndex a0 = parentAtomInfo.getIndex();
                    Compound::AtomIndex a1 = childAtomInfo.getIndex();
                    assert(a0 != a1);

                    assert(atomBonds.find(a0) != atomBonds.end());
                    assert(atomBonds.find(a1) != atomBonds.end());

                    DuMM::ClusterIndex parentClusterIndex = atomBonds[a0].clusterIx;
                    DuMM::ClusterIndex childClusterIndex = atomBonds[a1].clusterIx;

                    if (parentClusterIndex == childClusterIndex) break;  // same body, no action

                    RigidUnit& childUnit = rigidUnits.find(childClusterIndex)->second;
                    // parent should be same or undefined
                    assert(!childUnit.parentId.isValid());

                    childUnit.parentId = parentClusterIndex;
                    if (parentClusterIndex.isValid()) {
                        RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;
                        parentUnit.hasChild = true;
                    }

                    childUnit.inboardBondCenterIndex = childBondCenterInfo.getIndex();

                    childUnit.inboardBondDihedralAngle = bond.getDefaultDihedralAngle();

                    childUnit.inboardBondIndex = bondInfo.getIndex();

                    // based child frame on most inboard atom/bondCenter of the cluster

                    // child body frame in top compound frame including default dihedral rotation
                    Transform T_X_Mr = 
                        compoundRep.calcDefaultBondCenterFrameInCompoundFrame(compoundRep.getBondCenterInfo(childUnit.inboardBondCenterIndex), defaultAtomFrames);

                    childUnit.frameInTopCompoundFrame = T_X_Mr;
                }

                break;
            default:
                assert(false); // Uh oh, unrecognized bond mobility
                break;
        }

    }

    // 4.5) - temporary debugging status
    std::map<DuMM::ClusterIndex, RigidUnit>::iterator rigidUnitI;
    static bool doPrintRigidUnitDebugInfo = false;
    if (doPrintRigidUnitDebugInfo) 
    {
        for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
        {
            const RigidUnit& unit = rigidUnitI->second;

            cout << "Cluster number " << unit.clusterIx;
            cout << "\t";
            cout << "parent = " << unit.parentId;
            cout << "\t";
            cout << "# atoms = " << unit.clusterAtoms.size();
            std::set<Compound::AtomIndex>::const_iterator atomI;
            for (atomI = unit.clusterAtoms.begin(); atomI != unit.clusterAtoms.end(); ++atomI)
            {
                if (atomI == unit.clusterAtoms.begin()) 
                    cout << "(";
                else 
                    cout << ", ";

                cout << compound.getAtomName(*atomI);
            }
            cout << ")";

            // cout << "\t";
            // cout << "inboard bond dihedral angle = " << unit.inboardBondDihedralAngle;

            cout << endl;
        }
    }

    if (showDebugMessages) cout << "Step 5 set ground frames" << endl;
    // 5) Set body frames lacking parent bodies relative to Ground
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies having parent bodies
        if (rigidUnit.parentId.isValid()) continue;

        // Use root atom as basis for 
        // Follow atom tree to find root atom
        Compound::AtomIndex seedAtomIndex = *(rigidUnit.clusterAtoms.begin());
        const AtomBonding* rootAtomPtr = &atomBonds.find( seedAtomIndex )->second;
        while (rootAtomPtr->parentAtomIndex.isValid())
        {
            AtomBonding& candidateAtom = atomBonds[rootAtomPtr->parentAtomIndex];
            if (candidateAtom.clusterIx == rigidUnit.clusterIx)
                rootAtomPtr = &candidateAtom;
            else {
                break;
            }
        }
        assert(rootAtomPtr->clusterIx == rigidUnit.clusterIx);
       
        // possibly slow
        // Transform T_X_atom = compoundRep.calcDefaultAtomFrameInCompoundFrame(rootAtomPtr->atomId);
        
        // possibly faster
        const Transform& T_X_atom = defaultAtomFrames[rootAtomPtr->atomId];
        
        const Transform& G_X_T = compoundRep.getTopLevelTransform();
        Transform G_X_atom = G_X_T * T_X_atom;

        rigidUnit.frameInTopCompoundFrame = T_X_atom;
        rigidUnit.frameInParentFrame = G_X_atom;
    }

    if (showDebugMessages) cout << "Step 6 set child frames" << endl;
    // 6) Set body frames relative to parent frame for bodies that have a parent
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies without a parent body
        if (!rigidUnit.parentId.isValid()) continue;

        // reset child frame to zero dihedral angle
        Transform Mr_X_M0 = Rotation( rigidUnit.inboardBondDihedralAngle, XAxis);

        const Transform& T_X_Mr = rigidUnit.frameInTopCompoundFrame;
        Transform T_X_M0 = T_X_Mr * Mr_X_M0;

        // Express body frame relative to parent frame
        RigidUnit& parentUnit = rigidUnits.find(rigidUnit.parentId)->second;
        const Transform& T_X_Fr = parentUnit.frameInTopCompoundFrame;
        Transform Fr_X_T = ~T_X_Fr;

        // unrotated mobile frame in parent body frame
        Transform Fr_X_M0 = Fr_X_T * T_X_M0;
        rigidUnit.frameInParentFrame = Fr_X_M0;
    }

    if (showDebugMessages) cout << "Step 7 populate dumm clusters" << endl;
    // 7) Populate DuMMClusters
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        const RigidUnit& unit = rigidUnitI->second;

        std::set<Compound::AtomIndex>::const_iterator atomI;
        for (atomI = unit.clusterAtoms.begin(); atomI != unit.clusterAtoms.end(); ++atomI) 
        {
            Compound::AtomIndex atomId(*atomI);
            AtomBonding& atomBonding = atomBonds.find(atomId)->second;

            // possibly slow
            // Transform T_X_atom = compoundRep.calcDefaultAtomFrameInCompoundFrame(atomId);
            
            // possibly faster
            const Transform& T_X_atom = defaultAtomFrames[atomId];
            
            Transform T_X_B = unit.frameInTopCompoundFrame;
            Transform B_X_T = ~T_X_B;
            Transform B_X_atom = B_X_T * T_X_atom;

            Vec3 atomLocationInBody = B_X_atom.p();

            dumm.placeAtomInCluster(atomBonding.dummAtomIndex, unit.clusterIx, atomLocationInBody);
            atomBonding.locationInBodyFrame = atomLocationInBody;

            // Store atom location in Compound::Atom object
            CompoundAtom& atom = compoundRep.updAtom(atomId);
            atom.setFrameInMobilizedBodyFrame(B_X_atom);
        }
    }


    if (showDebugMessages) cout << "Step 8 create mobilized bodies" << endl;
    // 8) create MobilizedBodies
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        RigidUnit& unit = rigidUnitI->second;

        // skip this body if MobilizedBodyIndex is already defined
        if (unit.bodyId.isValid()) continue;

        // Case A: body to be attached to Ground
        if (!unit.parentId.isValid())
        {
            assert (unit.clusterAtoms.size() > 0);

            // body frame in ground frame
            const Transform& G_X_T = compoundRep.getTopLevelTransform();
            const Transform& T_X_B = unit.frameInTopCompoundFrame;

            // Anything with a child body attached gets a 6-dof Free Mobilizer,
            // because it can handle six degrees of freedom.  (assuming all children
            // are attached by Pin Mobilizers.)  Also anything with three or more
            // non-colinear atoms can be a Free Mobilizer.
            if ( (unit.clusterAtoms.size() > 2) || (unit.hasChild) ) 
            {   
                if (mobilizedBodyType.compare("Free") == 0) {
                    MobilizedBody::Free freeBody
                       (matter.Ground(), 
				        G_X_T * T_X_B, 
				        dumm.calcClusterMassProperties(unit.clusterIx), 
				        Transform());
		            unit.bodyId = freeBody.getMobilizedBodyIndex();
                } else if (mobilizedBodyType.compare("Weld") == 0) {
		            MobilizedBody::Weld weldBody
                       (matter.Ground(), 
					    G_X_T * T_X_B, 
					    dumm.calcClusterMassProperties(unit.clusterIx), 
					    Transform());
		            unit.bodyId = weldBody.getMobilizedBodyIndex();
                }
                dumm.attachClusterToBody(unit.clusterIx, unit.bodyId);
            }
            else if (unit.clusterAtoms.size() == 1) // One atom, no children => Cartesian mobility
            {
                MobilizedBody::Cartesian particleBody(matter.Ground(), 
                                                      G_X_T * T_X_B, 
                                                      dumm.calcClusterMassProperties(unit.clusterIx), 
                                                      Transform());

                unit.bodyId = particleBody.getMobilizedBodyIndex();
                dumm.attachClusterToBody(unit.clusterIx, unit.bodyId);
            }
            else // Two atoms, no children => FreeLine mobility
            {
                assert(unit.clusterAtoms.size() == 2);

                // Compute a rotation that places the atoms along the ZAxis
                // It does not matter the order of the atoms
                std::set<Compound::AtomIndex>::const_iterator atomI = unit.clusterAtoms.begin();
                Compound::AtomIndex atomId1 = *atomI;
                ++atomI;
                Compound::AtomIndex atomId2 = *atomI;

                Vec3 atom1Location = atomBonds[atomId1].locationInBodyFrame;
                Vec3 atom2Location = atomBonds[atomId2].locationInBodyFrame;

                UnitVec3 zAxis(0,0,1); // Required direction of freeline
                UnitVec3 bondDirection(atom2Location - atom1Location);

                // Compute transform that positions bond direction on z-axis
                // Br is body frame rotated so that atom axis is parallel to ZAxis
                Transform B_X_Br;
                Real rotAngle = std::acos( ~zAxis * bondDirection );
                if (rotAngle != 0) { // TODO - numerical precision
                    UnitVec3 rotAxis( bondDirection % zAxis );
                    B_X_Br = Rotation(rotAngle, rotAxis);
                }

                MobilizedBody::FreeLine freeBody(matter.Ground(), 
                                                 G_X_T * B_X_Br, // TODO - is this right?
                                                 dumm.calcClusterMassProperties(unit.clusterIx), 
                                                 T_X_B * B_X_Br);

                unit.bodyId = freeBody.getMobilizedBodyIndex();
                dumm.attachClusterToBody(unit.clusterIx, unit.bodyId);
            }
        }

        // Case B: rigid unit has a parent rigid unit
        else // unit has a parent body
        {
            DuMM::ClusterIndex parentClusterIndex = unit.parentId;
            RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;

            // If this assertion fails I may need to create a recursive method --CMB
            assert(parentUnit.bodyId.isValid());

            Transform P_X_M = unit.frameInParentFrame;

            // Move rotation axis (x) to z-axis
            Transform M_X_pin = Rotation(-90*Deg2Rad, YAxis);

            Bond& bond = compoundRep.updBond(compoundRep.updBondInfo(unit.inboardBondIndex));

            // this is just a hack for testing the unfinished ribose mobilizer
            bool testRiboseMobilizer = false;
            
            // CMB -- temporarily comment out Pin mobilizer while we test 
            // function based mobilizer for ribose pseudorotation
            if (testRiboseMobilizer) {
	            
	            RiboseNu3Mobilizer torsionBody(
	                                           matter.updMobilizedBody(parentUnit.bodyId),
	                                           P_X_M * M_X_pin,
	                                           dumm.calcClusterMassProperties(unit.clusterIx),
	                                           M_X_pin
	                                           );
	            
	            // Save a pointer to the pin joint in the bond object
	            // (ensure that the default angle of the MobilizedBody::Pin matches that of 
	            // the bond, in Atom.h)
	            // NOTE - setPinBody automatically sets the torsionBody default torsion angle
	            bond.setRiboseBody(torsionBody);
	            unit.bodyId = torsionBody.getMobilizedBodyIndex();
            }
            else {
            	
	            MobilizedBody::Pin torsionBody(
	                                           matter.updMobilizedBody(parentUnit.bodyId),
	                                           P_X_M * M_X_pin,
	                                           dumm.calcClusterMassProperties(unit.clusterIx),
	                                           M_X_pin);
	            // Save a pointer to the pin joint in the bond object
	            // (ensure that the default angle of the MobilizedBody::Pin matches that of 
	            // the bond, in Atom.h)
	            // NOTE - setPinBody automatically sets the torsionBody default torsion angle
	            bond.setPinBody(torsionBody);
	            unit.bodyId = torsionBody.getMobilizedBodyIndex();
            }
            
            dumm.attachClusterToBody(unit.clusterIx, unit.bodyId);

        }

        // Set mobilized body index of compound atoms
        std::set<Compound::AtomIndex>::const_iterator atomI;
        for (atomI = unit.clusterAtoms.begin(); atomI != unit.clusterAtoms.end(); ++atomI) {
            compoundRep.updAtom(*atomI).setMobilizedBodyIndex(unit.bodyId);
        }

    }

    if (showDebugMessages) cout << "Step 9 create decorations" << endl;
    // 9) Create nice visualization geometry
    if (hasDecorationSubsystem()) 
    {
        DecorationSubsystem&     artwork = updDecorationSubsystem();
        DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);

        for (DuMM::BondIndex i(0); i < dumm.getNumBonds(); ++i) {
            const DuMM::AtomIndex    a1 = dumm.getBondAtom(i,0), a2 = dumm.getBondAtom(i,1);
            const MobilizedBodyIndex b1 = dumm.getAtomBody(a1),  b2 = dumm.getAtomBody(a2);
            if (b1==b2)
                artwork.addBodyFixedDecoration(b1, Transform(),
                                               DecorativeLine(dumm.getAtomStationOnBody(a1), dumm.getAtomStationOnBody(a2))
                                                 .setColor(Gray).setLineThickness(3));
            else
                artwork.addRubberBandLine(b1, dumm.getAtomStationOnBody(a1),
                                          b2, dumm.getAtomStationOnBody(a2), crossBodyBond);
        }

        for (DuMM::AtomIndex anum(0); anum < dumm.getNumAtoms(); ++anum) {
            Real shrink = 0.25 /* 1 */, opacity = dumm.getAtomElement(anum)==1?0.5:1;
            Real r = dumm.getAtomRadius(anum);
            if (r<.001) r=0.1; //nm
            //opacity=0.5;//XXX
            artwork.addBodyFixedDecoration(dumm.getAtomBody(anum), dumm.getAtomStationOnBody(anum),
                DecorativeSphere(shrink*r)
                    .setColor(dumm.getAtomDefaultColor(anum)).setOpacity(opacity).setResolution(3));
        }
    }
    if (showDebugMessages) cout << "Finished modelOneCompound" << endl;
}

} // namespace SimTK

