#include "molmodel/internal/CompoundSystem.h"

#include "CompoundRep.h"
#include "SimTKmath.h"

#include "molmodel/internal/RiboseMobilizer.h"

#include <set>

using namespace std;

namespace SimTK {

void SimTK::CompoundSystem::modelCompounds(String mobilizedBodyType) 
{
    std::cout << "CompoundSystem::modelCompounds" << std::endl;
    // Turn off default decorations, since we'll make our own decorations.
    updMatterSubsystem().setShowDefaultGeometry(false); 
    for (CompoundSystem::CompoundIndex c(0); c < getNumCompounds(); ++c){
        
        // Get Compound
        Compound& compound = updCompound(c);
        
        std::vector<Transform> atomFrameCache(compound.getNAtoms());

        modelOneCompound(c, atomFrameCache, mobilizedBodyType);
    }
}

// RigidUnit data structure is for use in modelCompounds() method
class RigidUnit {
public:
    RigidUnit() {}
    RigidUnit(DuMM::ClusterIndex id) 
        : clusterIx(id), hasChild(false) {}

    DuMM::ClusterIndex clusterIx; // primary key
    MobilizedBodyIndex mbx; // populated toward the end
    DuMM::ClusterIndex parentId; // InvalidId implies parented to Ground
    bool hasChild; // Whether a child body is tethered to this one

    Transform frameInTopCompoundFrame; // useful intermediate computation
    Transform frameInParentFrame; // what we ultimately want

    Compound::BondCenterIndex inboardBondCenterIndex;
    Angle inboardBondDihedralAngle;
    Compound::BondIndex inboardBondIndex;

    std::set<Compound::AtomIndex> clusterAtoms;

    const void Print(void) const {
        
        std::cout << "unit" <<" "
            << "clustIx " << clusterIx <<" "
            << "parClustIx " << parentId <<" "
            << "mbx " << mbx <<" "
            << "inBCIx " << inboardBondCenterIndex <<" "
            << "inBoIx " << inboardBondIndex <<" "
            << "inBoDih " << inboardBondDihedralAngle <<" "
            << std::endl;

        std::cout << "unit aIxs" <<" ";
        for (const auto& clustAtom : clusterAtoms) {
            std::cout << clustAtom << " ";
        } std::cout << std::endl;
    }

    const void PrintTransforms(void) const {

        std::cout << "unit frameInTopCompoundFrame ";
        std::cout << frameInTopCompoundFrame;
        std::cout << "unit frameInParentFrame ";
        std::cout << frameInParentFrame;
    
    }

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
    // assert(!atom.getDuMMPrimaryClusterIndex().isValid()); // TODO why assert?
    atom.setDuMMPrimaryClusterIndex(clusterIx);

    assert(atomBondings.find(atomId) != atomBondings.end());
    AtomBonding& atomBonding = atomBondings.find(atomId)->second;

    atomBonding.clusterIx = clusterIx;
    assert(atomBonding.clusterIx.isValid());


    // Unbonded atoms always represent a lone rigid body
    int numberOfNonBendyBondsOnThisAtom = atomBonding.treeBonds.size() - atomBonding.freeTreeBonds.size();
    if (numberOfNonBendyBondsOnThisAtom == 0) return;

    // recursively check all atoms bonded to this one for membership in group
    std::set<Compound::AtomIndex>::const_iterator atomBondingIt;
    for (atomBondingIt = atomBonding.treeBonds.begin();
        atomBondingIt != atomBonding.treeBonds.end();
        ++atomBondingIt)
    {
        // Ignore atoms already assigned to a rigid cluster
        if (clusterAtoms.find(*atomBondingIt) != clusterAtoms.end()) 
             continue; // already assigned

        // Bond archaeology
        const AtomInfo& atomInfo1 = compoundRep.getAtomInfo(atomId);
        const AtomInfo& atomInfo2 = compoundRep.getAtomInfo(*atomBondingIt);
        const BondInfo& bondInfo = compoundRep.getBondInfo(atomInfo1, atomInfo2);
        const Bond& bond = compoundRep.getBond(bondInfo);

        //if (bond.getMobility() == BondMobility::Free) // OLDMOB
        if(false)
        {
            // continue;  // cannot restrict rigid body if free
        }
        else if ( (bond.getMobility() == BondMobility::Translation) // NEWMOB
                || (bond.getMobility() == BondMobility::FreeLine) // NEWMOB
                || (bond.getMobility() == BondMobility::Free) // NEWMOB
                || (bond.getMobility() == BondMobility::LineOrientationF) // NEWMOB
                || (bond.getMobility() == BondMobility::LineOrientationM)
                || (bond.getMobility() == BondMobility::UniversalM) // NEWMOB
                || (bond.getMobility() == BondMobility::Torsion)
                || (bond.getMobility() == BondMobility::AnglePin) // NEWMOB
                || (bond.getMobility() == BondMobility::BendStretch) // NEWMOB
                || (bond.getMobility() == BondMobility::Slider) // NEWMOB
                || (bond.getMobility() == BondMobility::Cylinder)
                || (bond.getMobility() == BondMobility::Spherical) // NEWMOB
                || (bond.getMobility() == BondMobility::BallF) // Gmol
                || (bond.getMobility() == BondMobility::BallM) ) // NEWMOB
        {
            // TORSION-only bonds to atoms with just one bond are included in rigid body
            // for inertian reasons

            const AtomBonding& bondedAtomBonding = atomBondings.find(*atomBondingIt)->second;
            int numberOfNonBendyBondsOnBondedAtom = bondedAtomBonding.treeBonds.size() - bondedAtomBonding.freeTreeBonds.size();
            assert(numberOfNonBendyBondsOnBondedAtom > 0);

            // TODO - Also join atoms with just one bond among Torsion or Rigid (or possibly Angle) bond
            
            // atoms with just one bond are not rotatable for inertia reasons
            // and thus get included in rigid body
            if ( (numberOfNonBendyBondsOnBondedAtom == 1) || (numberOfNonBendyBondsOnThisAtom == 1) )
            {
                //buildUpRigidBody(*b, clusterIx, clusterAtoms, atomBondings, compoundRep); // OLDMOB
                ; // NEWMOB
            }
        }

        else if (bond.getMobility() == BondMobility::Rigid) // definitely add to rigid body
        {
            buildUpRigidBody(*atomBondingIt, clusterIx, clusterAtoms, atomBondings, compoundRep);
        }
        else assert(false); // unexpected mobility
    }

}

/*!
 * <!-- Build Mobilized bodies based on informations in Compound, which in turn
  * contains CompoundAtoms, Bonds and BondCenters. -->
*/
void CompoundSystem::modelOneCompound(
    CompoundIndex compoundId,
    std::vector<Transform>& atomFrameCache,
    String mobilizedBodyType
)
{
    bool showDebugMessages = true;

    // ------------------------------------------------------------------------
    // (0) calc default Compound atom frames in Top
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages){cout << "modelOneCompound" << endl;}

    // Turn off default decorations, since we'll make our own decorations.
    updMatterSubsystem().setShowDefaultGeometry(false); 

    // Get Compound
    Compound& compound = updCompound(compoundId);
    CompoundRep& compoundRep = compound.updImpl();
    
    // Cache default atom frames for performance
    //std::vector<Transform> atomFrameCache(compoundRep.getNumAtoms());
    // if(atomFrameCache.size() == 0){
    //     atomFrameCache.resize(compoundRep.getNumAtoms());
    // }

    compoundRep.invalidateAtomFrameCache(atomFrameCache, compoundRep.getNumAtoms());
    compoundRep.calcDefaultAtomFramesInCompoundFrame(atomFrameCache);

    // Get the force field and Simbody's MatterSubsystem
    DuMMForceFieldSubsystem& dumm =
        (DuMMForceFieldSubsystem&) updMolecularMechanicsForceSubsystem();
    SimbodyMatterSubsystem&  matter  = updMatterSubsystem();

    // Convenient maps needed to manage rigid rigid units
    std::map<DuMM::ClusterIndex, RigidUnit> rigidUnits;
    std::map<Compound::AtomIndex, AtomBonding> atomBonds;


    // ------------------------------------------------------------------------
    // (1) Create initial atomBonds data structure for each atom.
    // No bonds are cached in this loop.
    // Here we also add DuMM atoms
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 1 create atomBonds" << endl;}

    // Iterate Compound atoms and assign DuMM::AtomIndex
    for (Compound::AtomIndex cACnt(0); cACnt < compound.getNumAtoms(); ++cACnt) 
    {
        // Create AtomBonding object
        atomBonds[cACnt] = AtomBonding(cACnt);
        AtomBonding& atomBonding = atomBonds[cACnt];

        // Get biotype index
        BiotypeIndex biotypeIx = compound.getAtomBiotypeIndex(cACnt);
        assert(biotypeIx.isValid());

        // Get ChargedAtomTypeIndex
        DuMM::ChargedAtomTypeIndex chargedTypeId = dumm.getBiotypeChargedAtomType(biotypeIx);
        assert(chargedTypeId.isValid());

        // Add Dumm atom
        atomBonding.dummAtomIndex = dumm.addAtom(chargedTypeId);
        assert(atomBonding.dummAtomIndex.isValid());

        // Store DuMMAtomIndex in Compound::Atom object
        CompoundAtom& atom = compoundRep.updAtom(cACnt);
        atom.setDuMMAtomIndex(atomBonding.dummAtomIndex);

        // Check
        if (showDebugMessages){
            std::cout << "SP_NEW_LAB cAIx dAIx bioIx chATIx "
                << cACnt << " " << atomBonding.dummAtomIndex <<" " 
                << biotypeIx <<" " << chargedTypeId <<" "
                << std::endl;
        }

    } // every atom


    // ------------------------------------------------------------------------
    // (2) Analyze atom bonding structure, and populate dumm bond structure
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 2 analyze bonding structure" << endl;}

    // Iterate Compound bonds and add to DuMM
    for (Compound::BondIndex bCnt(0); bCnt < compound.getNumBonds(); ++bCnt) 
    {
        // Get bond and bond info
        const BondInfo& bondInfo = compoundRep.getBondInfo(bCnt);
        const Bond& bond = compoundRep.getBond(bondInfo);

        // Get parent and child bond center infos
        const BondCenterInfo& parentBondCenterInfo =
            compoundRep.getBondCenterInfo(bondInfo.getParentBondCenterIndex());
        const BondCenterInfo&  childBondCenterInfo =
            compoundRep.getBondCenterInfo( bondInfo.getChildBondCenterIndex());
        
        // Get parent and child atom infos
        const AtomInfo& parentAtomInfo =
            compoundRep.getAtomInfo(parentBondCenterInfo.getAtomIndex());
        const AtomInfo& childAtomInfo =
            compoundRep.getAtomInfo( childBondCenterInfo.getAtomIndex());
        
        // Get parent and child Compound atom index
        Compound::AtomIndex a0 = parentAtomInfo.getIndex();
        Compound::AtomIndex a1 = childAtomInfo.getIndex();
        assert(a0 != a1);

        assert(atomBonds.find(a0) != atomBonds.end());
        assert(atomBonds.find(a1) != atomBonds.end());

        // Add bond to DuMM
        dumm.addBond(
            atomBonds[a0].dummAtomIndex,
            atomBonds[a1].dummAtomIndex);

        // Store bond info on each atom
        // only store those bonds that are part of the tree structure
        if (! bond.isRingClosingBond())
        {
            atomBonds[a0].treeBonds.insert(a1);
            atomBonds[a1].treeBonds.insert(a0);
            
            if (bond.getMobility() == BondMobility::Free) {
                //std::cout
                // << "CompoundSystem: Step2: Found Free bond between atoms "
                // << a0 << " " << a1 << std::endl;
                atomBonds[a0].freeTreeBonds.insert(a1);
                atomBonds[a1].freeTreeBonds.insert(a0);            	
            }
            
            // Store parent-child relationship
            atomBonds[a1].parentAtomIndex = a0;
        }
    } // every bond


    // ------------------------------------------------------------------------
    // (3) Distribute atoms to bodies, using DuMM::ClusterIndex as proxy for 
    // body for the present.
    // We use recursive function buildUpRigidBody to construct rigid units and
    // clusters. Clusters are primarly used in DuMM
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 3 distribute atoms" << endl;}

    // AtomBonding stores cAIx, dAIx, clusterIx, adjacency and its station in
    // Body
    std::map<Compound::AtomIndex, AtomBonding>::iterator atomBondIt;

    // Iterate AtomBondings
    for (atomBondIt = atomBonds.begin(); atomBondIt != atomBonds.end(); ++atomBondIt) 
    {
        AtomBonding& atomBonding = atomBondIt->second;

        // Ignore atoms already in a cluster
        if (atomBonding.clusterIx.isValid()) continue;

        Compound::AtomIndex cAIx = atomBonding.atomId;

        // Create new cluster and use clusterIx as primary key for new rigid body
        DuMM::ClusterIndex clusterIx = dumm.createCluster(String(cAIx));
        assert( rigidUnits.find(clusterIx) == rigidUnits.end() );

        // Insert new rigid unit into the cluster rigid units map
        rigidUnits[clusterIx] = RigidUnit(clusterIx);
        RigidUnit& rigidUnit = rigidUnits[clusterIx];
        assert( rigidUnits.find(clusterIx) != rigidUnits.end() );

        // Use recursive method to find all of the atoms in the cluster
        buildUpRigidBody(
            cAIx, clusterIx, rigidUnit.clusterAtoms,
            atomBonds, compoundRep);

    }


    // ------------------------------------------------------------------------
    // (4) Assign rigid body parents and set body frames relative to top  
    // compound
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 4 assign rigid body parents" << endl;
    
    // Iterate Compound bonds
    for (Compound::BondIndex bCnt(0); bCnt < compound.getNumBonds(); ++bCnt) 
    {
        //std::cout << "CompoundSystem: BondIndex " << i << std::endl;
        const BondInfo& bondInfo = compoundRep.getBondInfo(bCnt);

        // Get bond
        const Bond& bond = compoundRep.getBond(bondInfo);

        // Don't use ring closing bonds to assign parent child relationships
        // because we only want to use the tree structure of bonds to make a
        // multibody tree structure
        if (bond.isRingClosingBond())
            {continue;}

        switch (bond.getMobility()) {
            case BondMobility::Rigid:
                break; // same body => ignore
            case BondMobility::Free:
            case BondMobility::Translation:
            case BondMobility::FreeLine:
            case BondMobility::LineOrientationF:
            case BondMobility::LineOrientationM:
            { // NEWMOB introduce parent-child here too
                            
                // Get parent and child bond center infos                            
                const BondCenterInfo& parentBondCenterInfo =
                    compoundRep.getBondCenterInfo(bondInfo.getParentBondCenterIndex());
                const BondCenterInfo& childBondCenterInfo =
                    compoundRep.getBondCenterInfo(bondInfo.getChildBondCenterIndex());
                
                // Get parent and child atom infos
                const AtomInfo& parentAtomInfo =
                    compoundRep.getAtomInfo(parentBondCenterInfo.getAtomIndex());
                const AtomInfo& childAtomInfo =
                    compoundRep.getAtomInfo(childBondCenterInfo.getAtomIndex());

                // Get atoms
                Compound::AtomIndex a0 = parentAtomInfo.getIndex();
                Compound::AtomIndex a1 = childAtomInfo.getIndex();
                assert(a0 != a1);
                assert(atomBonds.find(a0) != atomBonds.end());
                assert(atomBonds.find(a1) != atomBonds.end());

                // Get parent and child cluster indexes
                DuMM::ClusterIndex parentClusterIndex = atomBonds[a0].clusterIx;
                DuMM::ClusterIndex childClusterIndex = atomBonds[a1].clusterIx;

                // This doesn't seem necessary
                if (parentClusterIndex == childClusterIndex){
                    break;  // same body, no action
                }

                // Get child rigid unit
                RigidUnit& childUnit = rigidUnits.find(childClusterIndex)->second;
                assert(!childUnit.parentId.isValid());

                // Get parent rigid unit
                childUnit.parentId = parentClusterIndex;
                if (parentClusterIndex.isValid()) {
                    RigidUnit& parentUnit =
                        rigidUnits.find(parentClusterIndex)->second;
                    parentUnit.hasChild = true;
                }

                // Set child rigid unit:
                //      inboardBCIx, inboardDihedral and inboardBondIx
                childUnit.inboardBondCenterIndex = childBondCenterInfo.getIndex();
                childUnit.inboardBondDihedralAngle = bond.getDefaultDihedralAngle();
                childUnit.inboardBondIndex = bondInfo.getIndex();

                // Set child body frame in top compound frame including
                // default dihedral rotation
                Transform T_X_Mr =
                    compoundRep.calcDefaultBondCenterFrameInCompoundFrame(
                        compoundRep.getBondCenterInfo(
                            childUnit.inboardBondCenterIndex),
                            atomFrameCache
                    );

                childUnit.frameInTopCompoundFrame = T_X_Mr;
            } break; // Gmol
            case BondMobility::Torsion:
            case BondMobility::Cylinder:
            case BondMobility::BallF:
            case BondMobility::BallM:
            case BondMobility::UniversalM:
            case BondMobility::Spherical:
            case BondMobility::AnglePin:
            case BondMobility::BendStretch:
            case BondMobility::Slider: {

                // Get parent and child bond center infos                            
                const BondCenterInfo& parentBondCenterInfo =
                    compoundRep.getBondCenterInfo(bondInfo.getParentBondCenterIndex());
                const BondCenterInfo& childBondCenterInfo =
                    compoundRep.getBondCenterInfo(bondInfo.getChildBondCenterIndex());

                // Get parent and child atom infos
                const AtomInfo& parentAtomInfo =
                    compoundRep.getAtomInfo(parentBondCenterInfo.getAtomIndex());
                const AtomInfo& childAtomInfo =
                    compoundRep.getAtomInfo(childBondCenterInfo.getAtomIndex());

                // Get atoms
                Compound::AtomIndex a0 = parentAtomInfo.getIndex();
                Compound::AtomIndex a1 = childAtomInfo.getIndex();
                assert(a0 != a1);
                assert(atomBonds.find(a0) != atomBonds.end());
                assert(atomBonds.find(a1) != atomBonds.end());

                // Get parent and child cluster indexes
                DuMM::ClusterIndex parentClusterIndex = atomBonds[a0].clusterIx;
                DuMM::ClusterIndex childClusterIndex = atomBonds[a1].clusterIx;

                // This doesn't seem necessary
                if (parentClusterIndex == childClusterIndex){
                    break; // same body, no action
                }

                // Get child rigid unit
                RigidUnit& childUnit = rigidUnits.find(childClusterIndex)->second;
                assert(!childUnit.parentId.isValid());

                // Get parent rigid unit
                childUnit.parentId = parentClusterIndex;
                if (parentClusterIndex.isValid()) {
                    RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;
                    parentUnit.hasChild = true;
                }

                // Set child rigid unit:
                //      inboardBCIx, inboardDihedral and inboardBondIx
                childUnit.inboardBondCenterIndex = childBondCenterInfo.getIndex();
                childUnit.inboardBondDihedralAngle = bond.getDefaultDihedralAngle();
                childUnit.inboardBondIndex = bondInfo.getIndex();

                // Set child body frame in top compound frame including
                // default dihedral rotation
                Transform T_X_Mr = 
                    compoundRep.calcDefaultBondCenterFrameInCompoundFrame(
                        compoundRep.getBondCenterInfo(childUnit.inboardBondCenterIndex), atomFrameCache
                    );

                childUnit.frameInTopCompoundFrame = T_X_Mr;

            } break;

            default:
                assert(false); // Uh oh, unrecognized bond mobility
                break;
        }

    } // every bond

    // 4.5) - temporary debugging status
    std::map<DuMM::ClusterIndex, RigidUnit>::iterator rigidUnitI;
    static bool doPrintRigidUnitDebugInfo = true;
    if (doPrintRigidUnitDebugInfo) 
    {
        for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
        {
            const RigidUnit& unit = rigidUnitI->second;

            //cout << "Cluster number " << unit.clusterIx <<"\t"
            // << "parent = " << unit.parentId <<"\t"
            //<< "# atoms = " << unit.clusterAtoms.size();

            std::set<Compound::AtomIndex>::const_iterator atomI;
            for (atomI = unit.clusterAtoms.begin();
            atomI != unit.clusterAtoms.end();
            ++atomI){}
            //cout << endl;
        }
    }

    // ------------------------------------------------------------------------
    // (5) Set rigid units frames lacking parent bodies relative to Ground 
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages){std::cout << "Step 5 set ground frames" << std::endl;}

    // Iterate rigid units
    for (rigidUnitI = rigidUnits.begin();
    rigidUnitI != rigidUnits.end();
    ++rigidUnitI)
    {
        // Get rigid unit
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies having parent bodies
        if (rigidUnit.parentId.isValid()) continue;

        // Start from the first atom in the rigid unit
        Compound::AtomIndex seedAtomIndex = *(rigidUnit.clusterAtoms.begin());
        const AtomBonding* rootAtomPtr = &atomBonds.find( seedAtomIndex )->second;

        // And follow atom tree to find root atom
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
        
        // Get it's Top frame from the already calculated frames at step 0
        const Transform& T_X_rootAtom = atomFrameCache[rootAtomPtr->atomId];
        
        const Transform& G_X_T = compoundRep.getTopLevelTransform();
        Transform G_X_rootAtom = G_X_T * T_X_rootAtom;

        rigidUnit.frameInTopCompoundFrame = T_X_rootAtom;
        rigidUnit.frameInParentFrame = G_X_rootAtom;

    } // every rigid unit


    // ------------------------------------------------------------------------
    // (6) Set body frames relative to parent frame for bodies that have a parent 
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 6 set child frames" << endl;
    
    // Iterate rigid units
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        // Get rigid units
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies without a parent body
        if (!rigidUnit.parentId.isValid()) continue;

        // Reset child frame to zero dihedral angle
        Transform Mr_X_M0 = Rotation(rigidUnit.inboardBondDihedralAngle, XAxis);

        // Get frame in Top Compound frame
        const Transform& T_X_Mr = rigidUnit.frameInTopCompoundFrame;
        Transform T_X_M0 = T_X_Mr * Mr_X_M0;

        // Express body frame relative to parent frame
        RigidUnit& parentUnit = rigidUnits.find(rigidUnit.parentId)->second;
        const Transform& T_X_Fr = parentUnit.frameInTopCompoundFrame;
        Transform Fr_X_T = ~T_X_Fr;

        // Unrotated mobile frame in parent body frame
        Transform Fr_X_M0 = Fr_X_T * T_X_M0;
        rigidUnit.frameInParentFrame = Fr_X_M0;

    } // every rigid unit


    // ------------------------------------------------------------------------
    // (7) Populate DuMMClusters (atom station in cluster and atom transform 
    // in body)
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 7 populate dumm clusters" << endl;
    
    // Iterate rigid units
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end();
    ++rigidUnitI)
    {
        // Get rigid unit
        const RigidUnit& unit = rigidUnitI->second;
        
        // Iterate atoms inside the unit
        std::set<Compound::AtomIndex>::const_iterator atomI;
        for (atomI = unit.clusterAtoms.begin(); atomI != unit.clusterAtoms.end();
        ++atomI) 
        {
            // Get atom
            Compound::AtomIndex atomId(*atomI);
            AtomBonding& atomBonding = atomBonds.find(atomId)->second;
            
            // Get it's Top frame from the already calculated frames at step 0
            const Transform& T_X_atom = atomFrameCache[atomId];
            
            // Calculate frame in mobod frame
            Transform T_X_B = unit.frameInTopCompoundFrame;
            Transform B_X_T = ~T_X_B;
            Transform B_X_atom = B_X_T * T_X_atom;

            // Also calculate station
            Vec3 atomLocationInBody = B_X_atom.p();
            dumm.placeAtomInCluster(
                atomBonding.dummAtomIndex,
                unit.clusterIx,
                atomLocationInBody);
            atomBonding.locationInBodyFrame = atomLocationInBody;

            // Store atom location in Compound::Atom object
            CompoundAtom& atom = compoundRep.updAtom(atomId);
            atom.setFrameInMobilizedBodyFrame(B_X_atom);
        }

    } // every unit


    // ------------------------------------------------------------------------
    // (8) create MobilizedBodies 
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 8 create mobilized bodies" << endl;
    
    // Iterate rigid units
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end();
    ++rigidUnitI)
    {
        RigidUnit& unit = rigidUnitI->second;

        // skip this body if MobilizedBodyIndex is already defined
        if (unit.mbx.isValid()) continue;

        ///////////////////////////////////////////////////////////////////////
        ////////////////// Case A: body to be attached to Ground //////////////
        ///////////////////////////////////////////////////////////////////////
        if (!unit.parentId.isValid())
        {
            assert (unit.clusterAtoms.size() > 0);

            // body frame in ground frame
            const Transform& G_X_T = compoundRep.getTopLevelTransform();
            const Transform& T_X_B = unit.frameInTopCompoundFrame;

            //if ( (unit.clusterAtoms.size() > 2) || (unit.hasChild) )
            if ( true ) // TEST
            {   
                if (mobilizedBodyType.compare("Free") == 0) {
                    MobilizedBody::Free freeBody
                            (matter.Ground(),
                             G_X_T * T_X_B,
                             dumm.calcClusterMassProperties(unit.clusterIx),
                             Transform());
                    unit.mbx = freeBody.getMobilizedBodyIndex();
                    std::cout << "First body Free mobodIx " << unit.mbx << std::endl;
                }else if (mobilizedBodyType.compare("Cartesian") == 0) {
                        MobilizedBody::Translation cartesianBody
                                (matter.Ground(),
                                 G_X_T * T_X_B,
                                 dumm.calcClusterMassProperties(unit.clusterIx),
                                 Transform());
                        unit.mbx = cartesianBody.getMobilizedBodyIndex();
                        std::cout << "First body Cartesian mobodIx " << unit.mbx << std::endl;
                } else if (mobilizedBodyType.compare("Weld") == 0) {
		            MobilizedBody::Weld weldBody
                       (matter.Ground(), 
					    G_X_T * T_X_B, 
					    dumm.calcClusterMassProperties(unit.clusterIx), 
					    Transform());
		            unit.mbx = weldBody.getMobilizedBodyIndex();
                    std::cout << "First body Weld mobodIx " << unit.mbx << std::endl;
                }else if (mobilizedBodyType.compare("FreeLine") == 0) {
                    MobilizedBody::FreeLine freeLineBody
                            (matter.Ground(),
                             G_X_T * T_X_B,
                             dumm.calcClusterMassProperties(unit.clusterIx),
                             Transform());
                    unit.mbx = freeLineBody.getMobilizedBodyIndex();
                    std::cout << "First body FreeLine mobodIx " << unit.mbx << std::endl;
                }else if (mobilizedBodyType.compare("Ball") == 0) {
                    MobilizedBody::Ball ballBody
                            (matter.Ground(),
                             G_X_T * T_X_B,
                             dumm.calcClusterMassProperties(unit.clusterIx),
                             Transform());
                    unit.mbx = ballBody.getMobilizedBodyIndex();
                    std::cout << "First body Ball mobodIx " << unit.mbx << std::endl;
                }else if (mobilizedBodyType.compare("Pin") == 0) {
                    MobilizedBody::Pin pinBody
                            (matter.Ground(),
                             G_X_T * T_X_B,
                             dumm.calcClusterMassProperties(unit.clusterIx),
                             Transform());
                    unit.mbx = pinBody.getMobilizedBodyIndex();
                    std::cout << "First body Pin mobodIx " << unit.mbx << std::endl;
                }else{
		            MobilizedBody::Weld weldBody
                       (matter.Ground(), 
					    G_X_T * T_X_B, 
					    dumm.calcClusterMassProperties(unit.clusterIx), 
					    Transform());
		            unit.mbx = weldBody.getMobilizedBodyIndex();
                    std::cout << "First body implicitly Weld mobodIx " << unit.mbx << std::endl;                    
                }
                
                dumm.attachClusterToBody(unit.clusterIx, unit.mbx);
            }
            /* else if (unit.clusterAtoms.size() == 1) // One atom, no children => Cartesian mobility
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
                // std::cout << " got FreeLine mobodIx " << unit.bodyId << std::endl;
            } */
        }

        ///////////////////////////////////////////////////////////////////////
        //////////////// Case B: rigid unit has a parent rigid unit ///////////
        ///////////////////////////////////////////////////////////////////////
        else // Unit has a parent body
        {

            // This is just a hack for testing the unfinished ribose mobilizer
            bool testRiboseMobilizer = false;

            // Get the parent rigid unit
            // If assertion fails I may need to create a recursive method --CMB
            DuMM::ClusterIndex parentClusterIndex = unit.parentId;
            RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;
            assert(parentUnit.mbx.isValid());


            // Get rigid unit inboard bond
            Bond& bond = compoundRep.updBond(compoundRep.updBondInfo(
                unit.inboardBondIndex));

            // Get cAIx, AtomBonding and AtomInfo at the origin of the Unit
            Compound::AtomIndex originAtomId(*unit.clusterAtoms.begin());
            AtomBonding& originAtomBonding = atomBonds.find(originAtomId)->second;
            const AtomInfo& originAtomInfo = compoundRep.getAtomInfo(originAtomId);

            // Get the parent atom cAIx and AtomInfo
            Compound::AtomIndex parentAtomId = originAtomBonding.parentAtomIndex;
            const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentAtomId);

            // Get BondInfo
            const BondInfo& bondInfo =
                compoundRep.getBondInfo(originAtomInfo, parentAtomInfo);

            // Get frame in parent rigid unit frame = Fr_X_M0
            Transform Fr_X_M0 = unit.frameInParentFrame;

            // Axis switching Rotations
            Transform XAxis_To_ZAxis = Rotation(-90*Deg2Rad, YAxis);
            Transform YAxis_To_ZAxis = Rotation(-90*Deg2Rad, XAxis);
            Transform XAxis_To_YAxis = Rotation(-90*Deg2Rad, ZAxis);

            // Get parent BC to child BC transforms:
            //     1) rotate about x-axis by dihedral angle
            //     2) translate along x-axis by bond length
            //     3) rotate 180 degrees about y-axis to face the parent bond center
            Transform X_parentBC_childBC =
                bondInfo.getBond().getDefaultBondCenterFrameInOtherBondCenterFrame();
            Transform X_childBC_parentBC = ~X_parentBC_childBC;

            // -------------- Old mobod transforms:
            //    - X_PF is parent rigid unit inboard_BC to child rigid unit
            //       inboard_BC without the default dihedral and with the X
            //        axis switched to Z axis
            //     - X_MB switches the Z axis back to X axis
            Transform oldX_PF = Fr_X_M0 * XAxis_To_ZAxis;
            Transform oldX_BM = XAxis_To_ZAxis;
            Transform oldX_MB = ~oldX_BM; // ZAxis_To_XAxis
            Transform oldX_FM = Rotation(bond.getDefaultDihedral(), ZAxis);
            Transform oldX_PB = (oldX_PF * oldX_FM * oldX_MB);

            // -------------- New mobod transforms:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC) with the X axis switched to
            //       Z axis
            //    - X_MB is the inboard bond parent BC to child BC transform
            //      with the Z axis switched to X
            Transform newX_BM = X_childBC_parentBC * XAxis_To_ZAxis;
            Transform newX_PF = oldX_PB * newX_BM;

            // -------------- Universal transforms:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC) with the Y axis switched to
            //       Z axis
            //    - X_MB is the inboard bond parent BC to child BC transform
            //      with the Z axis switched to Y
            Transform universalX_BM = X_childBC_parentBC * YAxis_To_ZAxis;
            Transform universalX_PF = oldX_PB * universalX_BM;

            // -------------- BendStretch, AnglePin and Slider transforms:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC)
            //    - X_MB is the inboard bond parent BC to child BC transform
            Transform bendStretchX_BM = X_childBC_parentBC;
            Transform bendStretchX_PF = oldX_PB * bendStretchX_BM;

            // -------------- AnglePin alias:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC)
            //    - X_MB is the inboard bond parent BC to child BC transform
            Transform anglePinX_BM = bendStretchX_BM;
            Transform anglePinX_PF = bendStretchX_PF;

            // -------------- Slider alias:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC)
            //    - X_MB is the inboard bond parent BC to child BC transform
            Transform sliderX_BM = bendStretchX_BM;
            Transform sliderX_PF = bendStretchX_PF;

            // -------------- Spherical transforms:
            //    - X_PF is parent rigid unit inboard_BC to the outboard
            //       atom's BC (parent BC) with the Y axis switched to
            //       Z axis
            //    - X_MB is the inboard bond parent BC to child BC transform
            //      with the Z axis switched to Y
            Transform sphericalX_BM = X_childBC_parentBC * XAxis_To_ZAxis * XAxis_To_YAxis;
            Transform sphericalX_PF = oldX_PB * sphericalX_BM;

	        // -------------- Special transform for free line
            Transform newX_BM_FreeLine;
            std::set<Compound::AtomIndex>::const_iterator atomI =
                unit.clusterAtoms.begin();
            Compound::AtomIndex atomId1 = *atomI;
            ++atomI;
            Compound::AtomIndex atomId2 = *atomI;

            Vec3 atom1Location = atomBonds[atomId1].locationInBodyFrame;
            Vec3 atom2Location = atomBonds[atomId2].locationInBodyFrame;
	    
            UnitVec3 zAxis(0,0,1); // Required direction of freeline
            UnitVec3 bondDirection(atom2Location - atom1Location);

            // Compute transform that positions bond direction on z-axis
            // Br is body frame rotated so that atom axis is parallel to ZAxis
            Real rotAngle = std::acos( ~zAxis * bondDirection );
            if (rotAngle != 0) { // TODO - numerical precision
                 UnitVec3 rotAxis( bondDirection % zAxis );
                 newX_BM_FreeLine = newX_BM * Transform(Rotation(rotAngle, rotAxis));
            }

            // GMOL END

            // CMB -- temporarily comment out Pin mobilizer while we test 
            // function based mobilizer for ribose pseudorotation
            if (testRiboseMobilizer) {
	            
	            RiboseNu3Mobilizer torsionBody(
                    matter.updMobilizedBody(parentUnit.mbx),
                    Fr_X_M0 * XAxis_To_ZAxis,
                    dumm.calcClusterMassProperties(unit.clusterIx),
                    XAxis_To_ZAxis);
	            
	            // Save a pointer to the pin joint in the bond object
	            // (ensure that the default angle of the MobilizedBody::Pin matches that of 
	            // the bond, in Atom.h)
	            // NOTE - setPinBody automatically sets the torsionBody default torsion angle
	            bond.setRiboseBody(torsionBody);
	            unit.mbx = torsionBody.getMobilizedBodyIndex();
            }
            else {
                
               if(bond.getMobility() == BondMobility::Torsion) {
                    // GMOL BIG RB ======
                    MobilizedBody::Pin torsionBody(
                            matter.updMobilizedBody(parentUnit.mbx),
                            newX_PF,
                            dumm.calcClusterMassProperties(unit.clusterIx),
                            newX_BM);

                    bond.setPinBody(torsionBody, 0);
                    torsionBody.setDefaultAngle(0); // no chem change to 0 for chem
                    unit.mbx = torsionBody.getMobilizedBodyIndex();
                    //std::cout << " got Pin mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Pin";

                }else if(bond.getMobility() == BondMobility::AnglePin) {

                   MobilizedBody::Torsion anglePinBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           anglePinX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           anglePinX_BM
                   );

                   bond.setAnglePinBody(anglePinBody, 0);
                   unit.mbx = anglePinBody.getMobilizedBodyIndex();
                   //std::cout << " got AnglePin mobodIx " << unit.bodyId << std::endl;
                    std::cout << " aPin";

                }else if(bond.getMobility() == BondMobility::BendStretch) {

                   MobilizedBody::BendStretch bendStretchBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           bendStretchX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           bendStretchX_BM
                   );

                   bond.setBendStretchBody(bendStretchBody, 0, 0);
                   unit.mbx = bendStretchBody.getMobilizedBodyIndex();
                   //std::cout << " got BendStretch mobodIx " << unit.bodyId << std::endl;
                    std::cout << " BeSt";

               }else if(bond.getMobility() == BondMobility::Slider) {

                   MobilizedBody::Slider sliderBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           sliderX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           sliderX_BM
                   );

                   bond.setSliderBody(sliderBody, 0);
                   unit.mbx = sliderBody.getMobilizedBodyIndex();
                   //std::cout << " got Slider mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Sli";

               } else if(bond.getMobility() == BondMobility::BallF) {

                    MobilizedBody::Ball ballBody(
                            matter.updMobilizedBody(parentUnit.mbx),
                            oldX_PF,
                            dumm.calcClusterMassProperties(unit.clusterIx),
                            oldX_BM
                    );

                    bond.setBallFBody(ballBody, 0);
                    unit.mbx = ballBody.getMobilizedBodyIndex();
                    //std::cout << " got BallF mobodIx " << unit.bodyId << std::endl;
                    std::cout << " BalF";

                }else if(bond.getMobility() == BondMobility::BallM) {

                   MobilizedBody::Ball ballBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           newX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           newX_BM
                   );

                   bond.setBallMBody(ballBody, 0);
                   unit.mbx = ballBody.getMobilizedBodyIndex();
                   //std::cout << " got BallM mobodIx " << unit.bodyId << std::endl;
                    std::cout << " BalM";

               }else if(bond.getMobility() == BondMobility::Spherical) {

                   MobilizedBody::SphericalCoords sphereBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           sphericalX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           sphericalX_BM
                   );

                   //sphereBody.setRadialAxis(XAxis);
                   bond.setSphericalBody(sphereBody, 0, 0, 0); // BAT coordinates
                   unit.mbx = sphereBody.getMobilizedBodyIndex();
                   //std::cout << " got Spherical mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Sphe";

               }else if(bond.getMobility() == BondMobility::UniversalM) {

                   MobilizedBody::Universal universalMBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           universalX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           universalX_BM
                   );

                   bond.setUniversalMBody(universalMBody, 0);
                   unit.mbx = universalMBody.getMobilizedBodyIndex();
                   //std::cout << " got UniversalM mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Univ";

               }else if(bond.getMobility() == BondMobility::Translation) {

                   MobilizedBody::Translation transBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           newX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           newX_BM
                   );

                   bond.setTransBody(transBody, 0);
                   unit.mbx = transBody.getMobilizedBodyIndex();
                   //std::cout << " got Trans mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Trans";

               }else if(bond.getMobility() == BondMobility::FreeLine) {
                   MobilizedBody::FreeLine freeLineBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           newX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           newX_BM
                   );

                   bond.setFreeLineBody(freeLineBody, 0, 0);
                   unit.mbx = freeLineBody.getMobilizedBodyIndex();
                   //std::cout << " got FreeLine mobodIx " << unit.bodyId << std::endl;
                    std::cout << " FreeL";

               }else if(bond.getMobility() == BondMobility::LineOrientationF) {
                   MobilizedBody::LineOrientation lineOrientationFBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           oldX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           oldX_BM
                   );

                   bond.setLineOrientationFBody(lineOrientationFBody, 0);
                   unit.mbx = lineOrientationFBody.getMobilizedBodyIndex();
                   //std::cout << " got LineOrientationF mobodIx " << unit.bodyId << std::endl;
                    std::cout << " LiOF";

               }else if(bond.getMobility() == BondMobility::LineOrientationM) {
                   MobilizedBody::LineOrientation lineOrientationMBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           newX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           newX_BM
                   );

                   bond.setLineOrientationMBody(lineOrientationMBody, 0);
                   unit.mbx = lineOrientationMBody.getMobilizedBodyIndex();
                   //std::cout << " got LineOrientationM mobodIx " << unit.bodyId << std::endl;
                    std::cout << " LiOM";

               }else if(bond.getMobility() == BondMobility::Free) {

                   MobilizedBody::Free freeBody(
                           matter.updMobilizedBody(parentUnit.mbx),
                           newX_PF,
                           dumm.calcClusterMassProperties(unit.clusterIx),
                           newX_BM
                   );

                   bond.setFreeBody(freeBody, 0, 0);
                   unit.mbx = freeBody.getMobilizedBodyIndex();
                   //std::cout << " got Free mobodIx " << unit.bodyId << std::endl;
                    std::cout << " Free";

               }else if(bond.getMobility() == BondMobility::Cylinder) {
                    MobilizedBody::Cylinder cylinderBody(
                               matter.updMobilizedBody(parentUnit.mbx),
                               newX_PF,
                               dumm.calcClusterMassProperties(unit.clusterIx),
                               newX_BM);
                    // Save a pointer to the pin joint in the bond object
                    // (ensure that the default angle of the MobilizedBody::Pin matches that of
                    // the bond, in Atom.h)
                    // NOTE - setPinBody automatically sets the torsionBody default torsion angle
                    bond.setCylinderBody(cylinderBody, 0, 0);
                    unit.mbx = cylinderBody.getMobilizedBodyIndex();
                    //std::cout << " got Cylinder mobodIx " << unit.bodyId << std::endl;
                     std::cout << " Cyl";

               }
            }
            
            dumm.attachClusterToBody(unit.clusterIx, unit.mbx);

        }

        // Set mobilized body index of compound atoms
        std::set<Compound::AtomIndex>::const_iterator atomI;
        for (atomI = unit.clusterAtoms.begin(); atomI != unit.clusterAtoms.end(); ++atomI) {
            compoundRep.updAtom(*atomI).setMobilizedBodyIndex(unit.mbx);
        }

    } // every rigid unit
    std::cout << std::endl;



    // Print Rigid units
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end();
    ++rigidUnitI)
    {
        RigidUnit& unit = rigidUnitI->second;

        unit.Print();
        unit.PrintTransforms();
    }    


    if (showDebugMessages) cout << "Step 9 create decorations" << endl;
    // 9) Create nice visualization geometry
 //   /*
    if (hasDecorationSubsystem()) 
    {
        //DecorationSubsystem&     artwork = updDecorationSubsystem();
        DecorativeLine crossBodyBond; crossBodyBond.setColor(Orange).setLineThickness(5);
        DecorativeLine sameBodyBond; sameBodyBond.setColor(Gray).setLineThickness(3);

        for (DuMM::BondIndex i(0); i < dumm.getNumBonds(); ++i) {
            const DuMM::AtomIndex    a1 = dumm.getBondAtom(i,0), a2 = dumm.getBondAtom(i,1);
            const MobilizedBodyIndex b1 = dumm.getAtomBody(a1),  b2 = dumm.getAtomBody(a2);
            if (b1==b2) {
                //artwork.addBodyFixedDecoration(b1, Transform(),
                //                               DecorativeLine(dumm.getAtomStationOnBody(a1),
                //                                              dumm.getAtomStationOnBody(a2))
                //                                       .setColor(Gray).setLineThickness(3));
/*                artwork.addRubberBandLine(b1, dumm.getAtomStationOnBody(a1),
                                          b2, dumm.getAtomStationOnBody(a2), sameBodyBond);*/
            }else {
/*                artwork.addRubberBandLine(b1, dumm.getAtomStationOnBody(a1),
                                          b2, dumm.getAtomStationOnBody(a2), crossBodyBond);*/
            }
        }

/*        for (DuMM::AtomIndex anum(0); anum < dumm.getNumAtoms(); ++anum) {
            Real shrink = 0.25 , opacity = dumm.getAtomElement(anum)==1?0.5:1;
            Real r = dumm.getAtomRadius(anum);
            if (r<.001) r=0.1; //nm
            //opacity=0.5;//XXX
            artwork.addBodyFixedDecoration(dumm.getAtomBody(anum), dumm.getAtomStationOnBody(anum),
                DecorativeSphere(shrink*r)
                    .setColor(dumm.getAtomDefaultColor(anum)).setOpacity(opacity).setResolution(3));
        }*/
	
    }

    if (showDebugMessages) cout << "Finished modelOneCompound" << endl;
}

} // namespace SimTK

