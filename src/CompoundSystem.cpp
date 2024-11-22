#include "molmodel/internal/CompoundSystem.h"

#include "CompoundRep.h"
#include "SimTKmath.h"

#include "molmodel/internal/RiboseMobilizer.h"

#include <set>

using namespace std;

namespace SimTK {

/*! <!-- Model all Compounds. 
 * --> */
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

/*! <!-- RigidUnit data structure is for use in modelCompounds() method. 
 * --> */
class RigidUnit {
public:
    RigidUnit() {}
    RigidUnit(DuMM::ClusterIndex id) 
        : clusterIx(id), hasChild(false) {}

    DuMM::ClusterIndex clusterIx; // primary key
    MobilizedBodyIndex mbx; // populated toward the end
    DuMM::ClusterIndex parentClusterIx; // InvalidId implies parented to Ground
    bool hasChild; // Whether a child body is tethered to this one

    Transform frameInTopCompoundFrame; // useful intermediate computation
    Transform frameInParentFrameRotated; // added for tracking
    Transform frameInParentFrame; // what we ultimately want

    Compound::BondCenterIndex inboardBCIx;
    Angle inboardBondDihedralAngle;
    Compound::BondIndex inboardBondIx;

    std::set<Compound::AtomIndex> clusterCAIxs; // atoms included in this unit

    const void Print(void) const {
        
        std::cout << "unit" <<" "
            << "clustIx " << clusterIx <<" "
            << "parClustIx " << parentClusterIx <<" "
            << "mbx " << mbx <<" "
            << "inBCIx " << inboardBCIx <<" "
            << "inBoIx " << inboardBondIx <<" "
            << "inBoDih " << inboardBondDihedralAngle <<" "
            << std::endl;

        std::cout << "unit aIxs" <<" ";
        for (const auto& clustAtom : clusterCAIxs) {
            std::cout << clustAtom << " ";
        } std::cout << std::endl;
    }

    const void PrintTransforms(void) const {

        std::cout << "unit frameInTopCompoundFrame ";
        //std::cout << frameInTopCompoundFrame;
        SimTK::Test::PrintTransform(frameInTopCompoundFrame, 6, "frameInTopCompoundFrame", "X_TopUnit");

        std::cout << "unit frameInParentFrame ";
        //std::cout << frameInParentFrame;
        SimTK::Test::PrintTransform(frameInParentFrame, 6, "frameInParentFrame", "X_ParUnit");
        SimTK::Test::PrintTransform(frameInParentFrameRotated, 6, "frameInParentFrameRotated", "frameInParentFrameRotated");
    
    }

};

/*! <!-- AtomBonding data structure represents one Atom in the modelCompounds() method.
 * --> */
class AtomBonding {
public:
    AtomBonding() {}
    AtomBonding(Compound::AtomIndex id)
        : compundAIx(id) {}

    Compound::AtomIndex compundAIx; // primary key
    Compound::AtomIndex parentCAIx; // for finding rootiest atom in cluster
    DuMM::AtomIndex dummAIx;
    DuMM::ClusterIndex clusterIx;

    // Store only those bonds that are part of the multibody tree structure
    std::set<Compound::AtomIndex> treeBondedCAIxs;
    std::set<Compound::AtomIndex> freeTreeBondedCAIxs;

    Vec3 locationInBodyFrame;
};

/*! <!-- Recursive method to create rigid bodies from a seed atom 
 * --> */
void buildUpRigidBody(Compound::AtomIndex cAIx,
                      DuMM::ClusterIndex clusterIx,
                      std::set<Compound::AtomIndex>& clusterAIxs, 
                      std::map<Compound::AtomIndex, AtomBonding>& atomBondings,
                      CompoundRep& compoundRep
                      ) 
{
    clusterAIxs.insert(cAIx);

    // Store DuMM::Cluster id in Compound::Atom object
    CompoundAtom& atom = compoundRep.updAtom(cAIx);
    // assert(!atom.getDuMMPrimaryClusterIndex().isValid()); // TODO why assert?
    atom.setDuMMPrimaryClusterIndex(clusterIx);

    assert(atomBondings.find(cAIx) != atomBondings.end());
    AtomBonding& atomBonding = atomBondings.find(cAIx)->second;

    atomBonding.clusterIx = clusterIx;
    assert(atomBonding.clusterIx.isValid());


    // Unbonded atoms always represent a lone rigid body
    int numberOfNonBendyBondsOnThisAtom = atomBonding.treeBondedCAIxs.size() - atomBonding.freeTreeBondedCAIxs.size();
    if (numberOfNonBendyBondsOnThisAtom == 0) return;

    // recursively check all atoms bonded to this one for membership in group
    std::set<Compound::AtomIndex>::const_iterator atomBondingIt;
    for (atomBondingIt = atomBonding.treeBondedCAIxs.begin();
        atomBondingIt != atomBonding.treeBondedCAIxs.end();
        ++atomBondingIt)
    {

        // Ignore atoms already assigned to a rigid cluster
        if (clusterAIxs.find(*atomBondingIt) != clusterAIxs.end()) 
        {
            continue; // already assigned
        }

        // Bond archaeology
        const AtomInfo& atomInfo1 = compoundRep.getAtomInfo(cAIx);
        const AtomInfo& atomInfo2 = compoundRep.getAtomInfo(*atomBondingIt);
        const BondInfo& bondInfo = compoundRep.getBondInfo(atomInfo1, atomInfo2);
        const Bond& bond = compoundRep.getBond(bondInfo);

        if((bond.getMobility() == BondMobility::Translation)
        || (bond.getMobility() == BondMobility::FreeLine)
        || (bond.getMobility() == BondMobility::Free)
        || (bond.getMobility() == BondMobility::LineOrientationF)
        || (bond.getMobility() == BondMobility::LineOrientationM)
        || (bond.getMobility() == BondMobility::UniversalM)
        || (bond.getMobility() == BondMobility::Torsion)
        || (bond.getMobility() == BondMobility::AnglePin)
        || (bond.getMobility() == BondMobility::BendStretch)
        || (bond.getMobility() == BondMobility::Slider)
        || (bond.getMobility() == BondMobility::Cylinder)
        || (bond.getMobility() == BondMobility::Spherical)
        || (bond.getMobility() == BondMobility::BallF)
        || (bond.getMobility() == BondMobility::BallM))
        {
            // TORSION-only bonds to atoms with just one bond are included in rigid body
            // for inertian reasons

            const AtomBonding& bondedAtomBonding = atomBondings.find(*atomBondingIt)->second;
            int numberOfNonBendyBondsOnBondedAtom = bondedAtomBonding.treeBondedCAIxs.size() - bondedAtomBonding.freeTreeBondedCAIxs.size();
            assert(numberOfNonBendyBondsOnBondedAtom > 0);

            // TODO - Also join atoms with just one bond among Torsion or Rigid (or possibly Angle) bond
            
            // atoms with just one bond are not rotatable for inertia reasons
            // and thus get included in rigid body
            if ( (numberOfNonBendyBondsOnBondedAtom == 1) || (numberOfNonBendyBondsOnThisAtom == 1) )
            {
                //buildUpRigidBody(*b, clusterIx, clusterAtoms, atomBondings, compoundRep); // OLDMOB
                ; // NEWMOB
            }
        
        // definitely add to rigid body
        }else if (bond.getMobility() == BondMobility::Rigid)
        {
            buildUpRigidBody(*atomBondingIt,
                              clusterIx,
                              clusterAIxs,
                              atomBondings,
                              compoundRep);
        }

        // 
        else {
            assert(false); // unexpected mobility
        }
    }

}


/*! <!--
 * Once matchDefaults are done, we should be able to calculate mobod
 * transforms.
 * --> */
CompoundSystem&
CompoundSystem::calc_XPF_XBM(
    SimTK::Compound& compound,
    SimTK::Compound::AtomIndex originAtomIx,
    SimTK::Compound::AtomIndex parentAtomIx,
    BondMobility::Mobility bondMobility,
    SimTK::Transform& Fr_X_M0,
    SimTK::Angle rigidUnitInboardDihedral,
    std::vector<SimTK::Transform>& PFBM)
{

    CompoundRep compoundRep = compound.updImpl();

    const AtomInfo& originAtomInfo = compoundRep.getAtomInfo(originAtomIx);

    const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentAtomIx);

    const BondInfo& chemBondInfo = compoundRep.getBondInfo(originAtomInfo, parentAtomInfo);

    const Bond& chemBond = compoundRep.getBond(chemBondInfo);

    // Axis switching Rotations
    Transform XAxis_To_ZAxis = Rotation(-90*Deg2Rad, YAxis);
    Transform YAxis_To_ZAxis = Rotation(-90*Deg2Rad, XAxis);
    Transform XAxis_To_YAxis = Rotation(-90*Deg2Rad, ZAxis);

    // Get parent BC to child BC transforms:
    //     1) rotate about x-axis by dihedral angle
    //     2) translate along x-axis by bond length
    //     3) rotate 180 degrees about y-axis to face the parent bond center
    Transform X_parentBC_childBC =
        chemBondInfo.getBond().getDefaultBondCenterFrameInOtherBondCenterFrame();
    Transform X_childBC_parentBC = ~X_parentBC_childBC;

    // -------------- Old mobod transforms:
    //    - X_PF is parent rigid unit inboard_BC to child rigid unit
    //       inboard_BC without the default dihedral and with the X
    //        axis switched to Z axis
    //     - X_MB switches the Z axis back to X axis
    Transform oldX_PF = Fr_X_M0 * XAxis_To_ZAxis;
    Transform oldX_BM = XAxis_To_ZAxis;
    Transform oldX_MB = ~oldX_BM; // ZAxis_To_XAxis
    Transform oldX_FM = Rotation(rigidUnitInboardDihedral, ZAxis);
    Transform oldX_PB = (oldX_PF * oldX_FM * oldX_MB);

    // std::cout << "STUDY_calc_XPF_XBM rigidUnitInboardDihedral " << rigidUnitInboardDihedral << std::endl;
    // SimTK::Test::PrintTransform(XAxis_To_ZAxis * oldX_FM * oldX_MB, 3, "STUDY xz_ZPhi_zx", "STUDY xz_ZPhi_zx");
    // SimTK::Test::PrintTransform(oldX_FM, 3, "STUDY oldX_FM", "STUDY oldX_FM");

    // -------------- New mobod transforms:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC) with the X axis switched to
    //       Z axis
    //    - X_MB is the inboard bond parent BC to child BC transform
    //      with the Z axis switched to X
    if(bondMobility == BondMobility::Torsion){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::Translation){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::Free){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::Cylinder){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::BallM){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::LineOrientationM){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    else if(bondMobility == BondMobility::FreeLine){
        PFBM[1] = X_childBC_parentBC * XAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- Universal transforms:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC) with the Y axis switched to
    //       Z axis
    //    - X_MB is the inboard bond parent BC to child BC transform
    //      with the Z axis switched to Y
    else if(bondMobility == BondMobility::UniversalM){
        PFBM[1] = X_childBC_parentBC * YAxis_To_ZAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- BendStretch, AnglePin and Slider transforms:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC)
    //    - X_MB is the inboard bond parent BC to child BC transform
    else if(bondMobility == BondMobility::BendStretch){
        PFBM[1] = X_childBC_parentBC;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- AnglePin alias:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC)
    //    - X_MB is the inboard bond parent BC to child BC transform
    else if(bondMobility == BondMobility::AnglePin){
        PFBM[1] = X_childBC_parentBC;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- Slider alias:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC)
    //    - X_MB is the inboard bond parent BC to child BC transform
    else if(bondMobility == BondMobility::Slider){
        PFBM[1] = X_childBC_parentBC;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- Spherical transforms:
    //    - X_PF is parent rigid unit inboard_BC to the outboard
    //       atom's BC (parent BC) with the Y axis switched to
    //       Z axis
    //    - X_MB is the inboard bond parent BC to child BC transform
    //      with the Z axis switched to Y
    else if(bondMobility == BondMobility::Spherical){
        PFBM[1] = X_childBC_parentBC 
            * XAxis_To_ZAxis * XAxis_To_YAxis;
        PFBM[0] = oldX_PB * PFBM[1];
    }

    // -------------- Old transforms:
    else if(bondMobility == BondMobility::BallF){
        PFBM[1] = XAxis_To_ZAxis;
        PFBM[0] = Fr_X_M0 * XAxis_To_ZAxis;
    }

    else if(bondMobility == BondMobility::LineOrientationF){
        PFBM[1] = XAxis_To_ZAxis;
        PFBM[0] = Fr_X_M0 * XAxis_To_ZAxis;
    }


    return *this;

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
    bool showDebugMessages = false;

    // ------------------------------------------------------------------------
    // (0) Calc default Compound atom frames in Top
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
    std::map<Compound::AtomIndex, AtomBonding> atomBondings;


    // ------------------------------------------------------------------------
    // (1) - Create initial AtomBonding data for each atom.
    //     - Add DuMM atoms
    // No bonds are cached in this loop.
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 1 create atomBonds" << endl;}

    // Iterate Compound atoms and assign DuMM::AtomIndex
    for (Compound::AtomIndex cACnt(0); cACnt < compound.getNumAtoms(); ++cACnt) 
    {
        // Create AtomBonding object
        atomBondings[cACnt] = AtomBonding(cACnt);
        AtomBonding& atomBonding = atomBondings[cACnt];

        // Get biotype index
        BiotypeIndex biotypeIx = compound.getAtomBiotypeIndex(cACnt);
        assert(biotypeIx.isValid());

        // Get ChargedAtomTypeIndex
        DuMM::ChargedAtomTypeIndex chargedTypeId = dumm.getBiotypeChargedAtomType(biotypeIx);
        assert(chargedTypeId.isValid());

        // Add Dumm atom
        atomBonding.dummAIx = dumm.addAtom(chargedTypeId);
        assert(atomBonding.dummAIx.isValid());

        // Store DuMMAtomIndex in Compound::Atom object
        CompoundAtom& atom = compoundRep.updAtom(cACnt);
        atom.setDuMMAtomIndex(atomBonding.dummAIx);

        // // Check
        // if (showDebugMessages){
        //     std::cout << "SP_NEW_LAB cAIx dAIx bioIx chATIx "
        //         << cACnt << " " << atomBonding.dummAtomIndex <<" " 
        //         << biotypeIx <<" " << chargedTypeId <<" "
        //         << std::endl;
        // }

    } // every atom


    // ------------------------------------------------------------------------
    // (2) - Add DuMM bonds
    //     - Update AtomBondings
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 2 analyze bonding structure" << endl;}

    // Iterate Compound bonds and add to DuMM
    for (Compound::BondIndex bIxCnt(0); bIxCnt < compound.getNumBonds(); ++bIxCnt) 
    {
        // Get bond and bond info
        const BondInfo& bondInfo = compoundRep.getBondInfo(bIxCnt);
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
        Compound::AtomIndex cAIx_0 = parentAtomInfo.getIndex();
        Compound::AtomIndex cAIx_1 = childAtomInfo.getIndex();
        assert(cAIx_0 != cAIx_1);

        assert(atomBondings.find(cAIx_0) != atomBondings.end());
        assert(atomBondings.find(cAIx_1) != atomBondings.end());

        // Add bond to DuMM
        dumm.addBond(
            atomBondings[cAIx_0].dummAIx,
            atomBondings[cAIx_1].dummAIx);

        // Store bond info on each atom
        // only store those bonds that are part of the tree structure
        if (! bond.isRingClosingBond())
        {
            atomBondings[cAIx_0].treeBondedCAIxs.insert(cAIx_1);
            atomBondings[cAIx_1].treeBondedCAIxs.insert(cAIx_0);
            
            if (bond.getMobility() == BondMobility::Free) {

                atomBondings[cAIx_0].freeTreeBondedCAIxs.insert(cAIx_1);
                atomBondings[cAIx_1].freeTreeBondedCAIxs.insert(cAIx_0);            	
            }
            
            // Store parent-child relationship
            atomBondings[cAIx_1].parentCAIx = cAIx_0;
        }
    } // every bond


    // ------------------------------------------------------------------------
    // (3) - Create DuMM clusters
    //     - Create RigidUnits
    //     - Call buildUpRigidBody to populate RigidUnits
    // Distribute atoms to bodies, using DuMM::ClusterIndex as proxy for 
    // body for the present.
    // We use recursive function buildUpRigidBody to construct rigid units and
    // clusters. Clusters are primarly used in DuMM
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) {cout << "Step 3 distribute atoms" << endl;}

    // AtomBonding stores cAIx, dAIx, clusterIx, adjacency and its station in Body
    std::map<Compound::AtomIndex, AtomBonding>::iterator atomBondIt;

    // Iterate AtomBondings
    for (atomBondIt = atomBondings.begin(); atomBondIt != atomBondings.end(); ++atomBondIt) 
    {
        AtomBonding& atomBonding = atomBondIt->second;

        // Ignore atoms already in a cluster
        if (atomBonding.clusterIx.isValid()) continue;

        Compound::AtomIndex cAIx = atomBonding.compundAIx;

        // Create new cluster and use clusterIx as primary key for new rigid body
        DuMM::ClusterIndex clusterIx = dumm.createCluster(String(cAIx));
        assert( rigidUnits.find(clusterIx) == rigidUnits.end() );

        // Insert new rigid unit into the cluster rigid units map
        rigidUnits[clusterIx] = RigidUnit(clusterIx);
        RigidUnit& rigidUnit = rigidUnits[clusterIx];
        assert( rigidUnits.find(clusterIx) != rigidUnits.end() );

        // Use recursive method to find all of the atoms in the cluster
        buildUpRigidBody(
            cAIx, clusterIx, rigidUnit.clusterCAIxs,
            atomBondings, compoundRep);

    }


    // ------------------------------------------------------------------------
    // (4)  - Assign RigidUnits parent - child relationships (Input = BondInfo)
    //      - Set child rigid unit: inboardBCIx, inboardDihedral and inboardBondIx
    //      - Set RigidUnits frameInTopCompoundFrame = inboard BCFrameInCompoundFrame
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 4 assign rigid body parents" << endl;
    
    // Iterate Compound bonds
    for (Compound::BondIndex bIxCnt(0); bIxCnt < compound.getNumBonds(); ++bIxCnt) 
    {

        // Get bond
        const BondInfo& bondInfo = compoundRep.getBondInfo(bIxCnt);
        const Bond& bond = compoundRep.getBond(bondInfo);

        // Don't use ring closing bonds to assign parent child relationships
        // because we only want to use the tree structure of bonds to make a
        // multibody tree structure
        if (bond.isRingClosingBond()){
            continue;
        }
        
        if (bond.getMobility() == BondMobility::Rigid){

            // Do nothing (same body)

        } else { // joint

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
            Compound::AtomIndex cAIx_0 = parentAtomInfo.getIndex();
            Compound::AtomIndex cAIx_1 = childAtomInfo.getIndex();
            assert(cAIx_0 != cAIx_1);
            assert(atomBondings.find(cAIx_0) != atomBondings.end());
            assert(atomBondings.find(cAIx_1) != atomBondings.end());

            // Get parent and child cluster indexes
            DuMM::ClusterIndex parentClusterIndex = atomBondings[cAIx_0].clusterIx;
            DuMM::ClusterIndex childClusterIndex = atomBondings[cAIx_1].clusterIx;

            // This doesn't seem necessary
            if (parentClusterIndex == childClusterIndex){
                break; // same body, no action
            }

            // Get child rigid unit
            RigidUnit& childUnit = rigidUnits.find(childClusterIndex)->second;
            assert(!childUnit.parentClusterIx.isValid());

            // Get parent rigid unit
            childUnit.parentClusterIx = parentClusterIndex;
            if (parentClusterIndex.isValid()) {
                RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;
                parentUnit.hasChild = true;
            }

            // Set child rigid unit:
            //      inboardBCIx, inboardDihedral and inboardBondIx
            childUnit.inboardBCIx = childBondCenterInfo.getIndex();
            childUnit.inboardBondDihedralAngle = bond.getDefaultDihedralAngle();
            childUnit.inboardBondIx = bondInfo.getIndex();

            // Set child body frame in top compound frame including
            // default dihedral rotation
            Transform T_X_Mr = 
                compoundRep.calcDefaultBondCenterFrameInCompoundFrame(
                    compoundRep.getBondCenterInfo(
                        childUnit.inboardBCIx),
                    atomFrameCache
                );

            childUnit.frameInTopCompoundFrame = T_X_Mr;

        }

    } // every bond

    // 4.5) - temporary debugging status
    std::map<DuMM::ClusterIndex, RigidUnit>::iterator rigidUnitI;
    static bool doPrintRigidUnitDebugInfo = false;
    if (doPrintRigidUnitDebugInfo) 
    {
        for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
        {
            const RigidUnit& unit = rigidUnitI->second;

            //cout << "Cluster number " << unit.clusterIx <<"\t"
            // << "parent = " << unit.parentId <<"\t"
            //<< "# atoms = " << unit.clusterAtoms.size();

            std::set<Compound::AtomIndex>::const_iterator atomI;
            for (atomI = unit.clusterCAIxs.begin(); atomI != unit.clusterCAIxs.end(); ++atomI){

            }
            //cout << endl;
        }
    }

    // ------------------------------------------------------------------------
    // (5) Set rigid units lacking parent frameInParentFrame relative to Ground 
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    // if (showDebugMessages){std::cout << "Step 5 set ground frames" << std::endl;}

    // Iterate rigid units
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end(); ++rigidUnitI)
    {
        // Get rigid unit
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies having parent bodies
        if (rigidUnit.parentClusterIx.isValid()) continue;

        // Start from the first atom in the rigid unit
        Compound::AtomIndex seedAtomIndex = *(rigidUnit.clusterCAIxs.begin());
        const AtomBonding* rootAtomPtr = &atomBondings.find( seedAtomIndex )->second;

        // And follow atom tree to find root atom
        while (rootAtomPtr->parentCAIx.isValid())
        {
            AtomBonding& candidateAtom = atomBondings[rootAtomPtr->parentCAIx];
            if (candidateAtom.clusterIx == rigidUnit.clusterIx)
                {rootAtomPtr = &candidateAtom;}
            else {
                break;
            }
        }
        assert(rootAtomPtr->clusterIx == rigidUnit.clusterIx);
        
        // Get it's Top frame from the already calculated frames at step 0
        const Transform& T_X_rootAtom = atomFrameCache[rootAtomPtr->compundAIx];
        
        const Transform& G_X_T = compoundRep.getTopLevelTransform();
        Transform G_X_rootAtom = G_X_T * T_X_rootAtom;

        rigidUnit.frameInTopCompoundFrame = T_X_rootAtom;
        rigidUnit.frameInParentFrame = G_X_rootAtom;

        rigidUnit.frameInParentFrameRotated = G_X_rootAtom;

    } // every rigid unit


    // ------------------------------------------------------------------------
    // (6) Set body frames relative to parent frame for bodies that have a parent 
    // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&$$$$

    if (showDebugMessages) cout << "Step 6 set child frames" << endl;
    
    // Iterate rigid units
    int ruIx = -1;
    for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end();
    ++rigidUnitI)
    {
        ruIx++;

        // Get rigid units
        RigidUnit& rigidUnit = rigidUnitI->second;

        // Skip bodies without a parent body
        if (!rigidUnit.parentClusterIx.isValid()) continue;

        // Get parent frame in Top Compound frame
        RigidUnit& parentUnit = rigidUnits.find(rigidUnit.parentClusterIx)->second;
        const Transform& T_X_Fr = parentUnit.frameInTopCompoundFrame;
        Transform Fr_X_T = ~T_X_Fr;

        // Reset child frame to zero dihedral angle
        Transform Mr_X_M0 = Rotation(rigidUnit.inboardBondDihedralAngle, XAxis);

        // Get frame in Top Compound frame
        const Transform& T_X_Mr = rigidUnit.frameInTopCompoundFrame;
        Transform T_X_M0 = T_X_Mr * Mr_X_M0;

        // Unrotated mobile frame in parent body frame
        Transform Fr_X_M0 = Fr_X_T * T_X_M0;
        rigidUnit.frameInParentFrame = Fr_X_M0;

        rigidUnit.frameInParentFrameRotated = Fr_X_T * T_X_Mr;

        // std::cout << "STUDY RU_parentChild ruIx Fr_X_Mr Fr_X_M0 " << ruIx << std::endl;
        // SimTK::Test::PrintTransform(Fr_X_T * T_X_Mr, 3, "Fr_X_Mr", "Fr_X_Mr");
        // SimTK::Test::PrintTransform(Fr_X_M0, 3, "Fr_X_M0", "Fr_X_M0");      

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
        for (atomI = unit.clusterCAIxs.begin(); atomI != unit.clusterCAIxs.end();
        ++atomI)
        {
            // Get atom
            Compound::AtomIndex atomId(*atomI);
            AtomBonding& atomBonding = atomBondings.find(atomId)->second;
            
            // Get it's Top frame from the already calculated frames at step 0
            const Transform& T_X_atom = atomFrameCache[atomId];
            
            // Calculate frame in mobod frame
            Transform T_X_B = unit.frameInTopCompoundFrame;
            Transform B_X_T = ~T_X_B;
            Transform B_X_atom = B_X_T * T_X_atom;

            // Also calculate station
            dumm.placeAtomInCluster(atomBonding.dummAIx, unit.clusterIx, B_X_atom.p());
            
            atomBonding.locationInBodyFrame = B_X_atom.p();

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
        if (!unit.parentClusterIx.isValid())
        {
            assert (unit.clusterCAIxs.size() > 0);

            // Cluster mass properties
            SimTK::MassProperties massProps =
                dumm.calcClusterMassProperties(unit.clusterIx);

            // std::cout << "[INERTIA CompoundSystem::modelOneCompound] " << massProps.getInertia() << std::endl;

            // body frame in ground frame
            const Transform& G_X_T = compoundRep.getTopLevelTransform();
            const Transform& T_X_B = unit.frameInTopCompoundFrame;
            Transform G_X_B = G_X_T * T_X_B;
            SimTK::MobilizedBody::Ground & parentMobod = matter.Ground();

            //if ( (unit.clusterAtoms.size() > 2) || (unit.hasChild) )
            if ( true ) // TEST
            {   
                if (mobilizedBodyType.compare("Free") == 0) {
                    MobilizedBody::Free freeBody
                        (parentMobod, G_X_B, massProps, Transform());
                    unit.mbx = freeBody.getMobilizedBodyIndex();
                    // std::cout << "First body Free mobodIx " << unit.mbx << std::endl;

                }else if (mobilizedBodyType.compare("Cartesian") == 0) {
                    MobilizedBody::Translation cartesianBody
                        (parentMobod, G_X_B, massProps, Transform());
                    unit.mbx = cartesianBody.getMobilizedBodyIndex();
                    // std::cout << "First body Cartesian mobodIx " << unit.mbx << std::endl;

                }else if (mobilizedBodyType.compare("Weld") == 0) {
		            MobilizedBody::Weld weldBody
                        (parentMobod, G_X_B, massProps, Transform());
		            unit.mbx = weldBody.getMobilizedBodyIndex();
                    // std::cout << "First body Weld mobodIx " << unit.mbx << std::endl;

                }else if (mobilizedBodyType.compare("FreeLine") == 0) {
                    MobilizedBody::FreeLine freeLineBody
                        (parentMobod, G_X_B, massProps, Transform());
                    unit.mbx = freeLineBody.getMobilizedBodyIndex();
                    // std::cout << "First body FreeLine mobodIx " << unit.mbx << std::endl;

                }else if (mobilizedBodyType.compare("Ball") == 0) {
                    MobilizedBody::Ball ballBody
                        (parentMobod, G_X_B, massProps, Transform());
                    unit.mbx = ballBody.getMobilizedBodyIndex();
                    // std::cout << "First body Ball mobodIx " << unit.mbx << std::endl;

                }else if (mobilizedBodyType.compare("Pin") == 0) {
                    MobilizedBody::Pin pinBody
                        (parentMobod, G_X_B, massProps, Transform());
                    unit.mbx = pinBody.getMobilizedBodyIndex();
                    // std::cout << "First body Pin mobodIx " << unit.mbx << std::endl;

                }else{
		            MobilizedBody::Weld weldBody
                        (parentMobod, G_X_B, massProps, Transform());
		            unit.mbx = weldBody.getMobilizedBodyIndex();
                    // std::cout << "First body implicitly Weld mobodIx " << unit.mbx << std::endl;   

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
            DuMM::ClusterIndex parentClusterIndex = unit.parentClusterIx;
            RigidUnit& parentUnit = rigidUnits.find(parentClusterIndex)->second;
            assert(parentUnit.mbx.isValid());

            // Get rigid unit inboard bond (not chemical parent)
            Bond& unitInboardBond = compoundRep.updBond(compoundRep.updBondInfo(
                unit.inboardBondIx));

            // Get cAIx, AtomBonding and AtomInfo at the origin of the Unit
            Compound::AtomIndex originAtomId(*unit.clusterCAIxs.begin());
            AtomBonding& originAtomBonding = atomBondings.find(originAtomId)->second;
            const AtomInfo& originAtomInfo = compoundRep.getAtomInfo(originAtomId);

            // Get the parent atom cAIx and AtomInfo
            Compound::AtomIndex parentAtomId = originAtomBonding.parentCAIx;
            const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentAtomId);

            // Get BondInfo
            const BondInfo& bondInfo =
                compoundRep.getBondInfo(originAtomInfo, parentAtomInfo);

            // Get frame in parent rigid unit frame = Fr_X_M0
            Transform Fr_X_M0 = unit.frameInParentFrame;

            std::vector<SimTK::Transform> PFBM(2);
            calc_XPF_XBM(compound,
                originAtomId, parentAtomId,
                unitInboardBond.getMobility(),
                Fr_X_M0, unitInboardBond.getDefaultDihedral(),
                PFBM);

            if(false){ // STUDY

                CompoundRep compoundRep = compound.updImpl();

                const AtomInfo& originAtomInfo = compoundRep.getAtomInfo(originAtomId);
                const AtomInfo& parentAtomInfo = compoundRep.getAtomInfo(parentAtomId);

                const BondInfo& chemBondInfo = compoundRep.getBondInfo(originAtomInfo, parentAtomInfo);
                const Bond& chemBond = compoundRep.getBond(chemBondInfo);

                // Axis switching Rotations
                Transform XAxis_To_ZAxis = Rotation(-90*Deg2Rad, YAxis);
                Transform YAxis_To_ZAxis = Rotation(-90*Deg2Rad, XAxis);
                Transform XAxis_To_YAxis = Rotation(-90*Deg2Rad, ZAxis);

                // Get parent BC to child BC transforms:
                //     1) rotate about x-axis by dihedral angle
                //     2) translate along x-axis by bond length
                //     3) rotate 180 degrees about y-axis to face the parent bond center
                Transform X_parentBC_childBC = chemBondInfo.getBond().getDefaultBondCenterFrameInOtherBondCenterFrame();
                Transform X_childBC_parentBC = ~X_parentBC_childBC;

                // -------------- Old mobod transforms:
                //    - X_PF is parent rigid unit inboard_BC to child rigid unit
                //       inboard_BC without the default dihedral and with the X axis switched to Z axis
                //     - X_MB switches the Z axis back to X axis
                Transform oldX_PF = Fr_X_M0 * XAxis_To_ZAxis;
                Transform oldX_BM = XAxis_To_ZAxis;
                Transform oldX_MB = ~oldX_BM; // ZAxis_To_XAxis
                Transform oldX_FM = Rotation(unitInboardBond.getDefaultDihedral(), ZAxis);
                Transform oldX_PB = (oldX_PF * oldX_FM * oldX_MB);

                // std::cout << "STUDY_calc_XPF_XBM rigidUnitInboardDihedral " << rigidUnitInboardDihedral << std::endl;
                // SimTK::Test::PrintTransform(XAxis_To_ZAxis * oldX_FM * oldX_MB, 3, "STUDY xz_ZPhi_zx", "STUDY xz_ZPhi_zx");
                // SimTK::Test::PrintTransform(oldX_FM, 3, "STUDY oldX_FM", "STUDY oldX_FM");

                // -------------- New mobod transforms:
                //    - X_PF is parent rigid unit inboard_BC to the outboard
                //       atom's BC (parent BC) with the X axis switched to Z axis
                //    - X_MB is the inboard bond parent BC to child BC transform
                //      with the Z axis switched to X
                Transform PFBM_1 = X_childBC_parentBC * XAxis_To_ZAxis;
                Transform PFBM_0 = oldX_PB * PFBM_1;

            }

            // CMB -- temporarily comment out Pin mobilizer while we test
            // function based mobilizer for ribose pseudorotation
            if (testRiboseMobilizer) {
	            
	            // RiboseNu3Mobilizer torsionBody(
                //     matter.updMobilizedBody(parentUnit.mbx),
                //     Fr_X_M0 * Transform(Rotation(-90*Deg2Rad, YAxis)),
                //     dumm.calcClusterMassProperties(unit.clusterIx),
                //     Transform(Rotation(-90*Deg2Rad, YAxis)));
	            // // Save a pointer to the pin joint in the bond object
	            // // (ensure that the default angle of the MobilizedBody::Pin matches that of 
	            // // the bond, in Atom.h)
	            // // NOTE - setPinBody automatically sets the torsionBody default torsion angle
	            // unitInboardBond.setRiboseBody(torsionBody);
	            // unit.mbx = torsionBody.getMobilizedBodyIndex();
            }
            else {

                // Cluster mass properties
                SimTK::MassProperties massProps =
                    dumm.calcClusterMassProperties(unit.clusterIx);

                // Parent mobodd
                SimTK::MobilizedBody & parentMobod =
                    matter.updMobilizedBody(parentUnit.mbx);

                if(unitInboardBond.getMobility() == BondMobility::Torsion) {

                    MobilizedBody::Pin torsionBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setPinBody(torsionBody, 0);
                    torsionBody.setDefaultAngle(0); // no chem change to 0 for chem
                    unit.mbx = torsionBody.getMobilizedBodyIndex();
                    // std::cout << " Pin";

                }else if(unitInboardBond.getMobility() == BondMobility::AnglePin) {

                    MobilizedBody::Torsion anglePinBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setAnglePinBody(anglePinBody, 0);
                    unit.mbx = anglePinBody.getMobilizedBodyIndex();
                    // std::cout << " aPin";

                }else if(unitInboardBond.getMobility() == BondMobility::BendStretch) {

                    MobilizedBody::BendStretch bendStretchBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setBendStretchBody(bendStretchBody, 0, 0);
                    unit.mbx = bendStretchBody.getMobilizedBodyIndex();
                    // std::cout << " BeSt";

               }else if(unitInboardBond.getMobility() == BondMobility::Slider) {

                   MobilizedBody::Slider sliderBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setSliderBody(sliderBody, 0);
                    unit.mbx = sliderBody.getMobilizedBodyIndex();
                    // std::cout << " Sli";

               } else if(unitInboardBond.getMobility() == BondMobility::BallF) {

                    MobilizedBody::Ball ballBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setBallFBody(ballBody, 0);
                    unit.mbx = ballBody.getMobilizedBodyIndex();
                    // std::cout << " BalF";

                }else if(unitInboardBond.getMobility() == BondMobility::BallM) {

                    MobilizedBody::Ball ballBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setBallMBody(ballBody, 0);
                    unit.mbx = ballBody.getMobilizedBodyIndex();
                    // std::cout << " BalM";

               }else if(unitInboardBond.getMobility() == BondMobility::Spherical) {

                    MobilizedBody::SphericalCoords sphereBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setSphericalBody(sphereBody, 0, 0, 0); // BAT coordinates
                    unit.mbx = sphereBody.getMobilizedBodyIndex();
                    // std::cout << " Sphe";

               }else if(unitInboardBond.getMobility() == BondMobility::UniversalM) {

                    MobilizedBody::Universal universalMBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setUniversalMBody(universalMBody, 0);
                    unit.mbx = universalMBody.getMobilizedBodyIndex();
                    // std::cout << " Univ";

               }else if(unitInboardBond.getMobility() == BondMobility::Translation) {

                    MobilizedBody::Translation transBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setTransBody(transBody, 0);
                    unit.mbx = transBody.getMobilizedBodyIndex();
                    // std::cout << " Trans";

               }else if(unitInboardBond.getMobility() == BondMobility::FreeLine) {

                    MobilizedBody::FreeLine freeLineBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setFreeLineBody(freeLineBody, 0, 0);
                    unit.mbx = freeLineBody.getMobilizedBodyIndex();
                    // std::cout << " FreeL";

               }else if(unitInboardBond.getMobility() == BondMobility::LineOrientationF) {

                    MobilizedBody::LineOrientation lineOrientationFBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setLineOrientationFBody(lineOrientationFBody, 0);
                    unit.mbx = lineOrientationFBody.getMobilizedBodyIndex();
                    // std::cout << " LiOF";

               }else if(unitInboardBond.getMobility() == BondMobility::LineOrientationM) {

                    MobilizedBody::LineOrientation lineOrientationMBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setLineOrientationMBody(lineOrientationMBody, 0);
                    unit.mbx = lineOrientationMBody.getMobilizedBodyIndex();
                    // std::cout << " LiOM";

               }else if(unitInboardBond.getMobility() == BondMobility::Free) {

                    MobilizedBody::Free freeBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    unitInboardBond.setFreeBody(freeBody, 0, 0);
                    unit.mbx = freeBody.getMobilizedBodyIndex();
                    // std::cout << " Free";

               }else if(unitInboardBond.getMobility() == BondMobility::Cylinder) {

                    MobilizedBody::Cylinder cylinderBody(
                        parentMobod, PFBM[0], massProps, PFBM[1]);

                    // Save a pointer to the pin joint in the bond object
                    // (ensure that the default angle of the MobilizedBody::Pin matches that of
                    // the bond, in Atom.h)
                    // NOTE - setPinBody automatically sets the torsionBody default torsion angle
                    unitInboardBond.setCylinderBody(cylinderBody, 0, 0);
                    unit.mbx = cylinderBody.getMobilizedBodyIndex();
                    //std::cout << " got Cylinder mobodIx " << unit.bodyId << std::endl;
                    //  std::cout << " Cyl";

               }
            }
            
            dumm.attachClusterToBody(unit.clusterIx, unit.mbx);

        }

        // Set mobilized body index of compound atoms
        std::set<Compound::AtomIndex>::const_iterator atomI;
        for (atomI = unit.clusterCAIxs.begin(); atomI != unit.clusterCAIxs.end(); ++atomI) {
            compoundRep.updAtom(*atomI).setMobilizedBodyIndex(unit.mbx);
        }

    } // every rigid unit
    // std::cout << std::endl;

    // // Print Rigid units
    // for (rigidUnitI = rigidUnits.begin(); rigidUnitI != rigidUnits.end();
    // ++rigidUnitI)
    // {
    //     RigidUnit& unit = rigidUnitI->second;
    //     unit.Print();
    //     unit.PrintTransforms();
    // }    


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

