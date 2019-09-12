#include "molmodel/internal/AtomSubsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

using namespace SimTK::units::md;

namespace SimTK {

class AtomSubsystem::Impl : public Subsystem::Guts 
{
protected:
    MultibodySystem& system;

    // Topological set of atoms for this system
    std::vector<AtomSubsystem::Atom> atoms;
    std::vector<AtomSubsystem::Bond> bonds;

    // Relationship of atoms to all mobilized bodies in the matter subsystem
    std::map<MobilizedBodyIndex, AtomSubsystem::LocalBodyIndex> systemBodyIxToLocalBodyIx;
    std::map<AtomSubsystem::LocalBodyIndex, MobilizedBodyIndex> localBodyIxToSystemBodyIx;
    AtomSubsystem::AtomicBodies atomsByBody;

    mutable CacheEntryIndex atomPositionCacheIndex;
    mutable CacheEntryIndex atomBodyStationCacheIndex;
    mutable CacheEntryIndex atomVelocityCacheIndex;

public:
    // need to traverse bodies in the order they were inserted
    typedef std::vector< std::vector<AtomSubsystem::AtomIndex> > AtomAtomicBodies;

    friend class AtomSubsystem;
    friend class AtomSubsystem::PairIterator;

    // Only implement the realizeWhatever methods in this Impl() class
    // Other methods should be defined in the AtomSubsystem class
    // Keep the data in the Impl, the methods in the non-Impl
    
    Impl(MultibodySystem& system) 
        : system(system)
    {}

    Subsystem::Guts* cloneImpl() const {
        return new Impl(*this);
    }

    int realizeSubsystemTopologyImpl(State& state) const 
    {
        atomPositionCacheIndex = state.allocateCacheEntry(
                getMySubsystemIndex(), 
                Stage::Position, 
                new Value< Vector_<Vec3> >()
        );

        atomBodyStationCacheIndex = state.allocateCacheEntry(
                getMySubsystemIndex(), 
                Stage::Position, 
                new Value< Vector_<Vec3> >()
        );

        atomVelocityCacheIndex = state.allocateCacheEntry(
                getMySubsystemIndex(), 
                Stage::Velocity, 
                new Value< Vector_<Vec3> >()
        );

        return 0;
    }

    int realizeSubsystemPositionImpl(const State& state) const 
    {
        // Compute atom positions from body positions
        Vector_<Vec3>& atomPositionCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
            getMySubsystemIndex(), 
            atomPositionCacheIndex
            )).upd();

        atomPositionCache.resize( atoms.size() );

        Vector_<Vec3>& atomBodyStationCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
            getMySubsystemIndex(), 
            atomBodyStationCacheIndex
            )).upd();

        atomBodyStationCache.resize( atoms.size() );

        for (AtomSubsystem::LocalBodyIndex b(0); b < atomsByBody.size(); ++b) {
            const MobilizedBodyIndex& bodyIx = localBodyIxToSystemBodyIx.find(b)->second;
            const MobilizedBody& body = system.getMatterSubsystem().getMobilizedBody(bodyIx);
            const Transform& X_GB = body.getBodyTransform(state);

            const std::vector<AtomSubsystem::AtomIndex>& bodyAtoms = atomsByBody[b].atoms;
            for (size_t a(0); a < bodyAtoms.size(); ++a) {
                AtomSubsystem::AtomIndex atomIx = bodyAtoms[a];
                const AtomSubsystem::Atom& atom = atoms[atomIx];
                const Vec3& v_B_atom = atom.getStationInBodyFrame();

                atomBodyStationCache[atomIx] = X_GB.R() * v_B_atom;
                atomPositionCache[atomIx] = X_GB.T() + atomBodyStationCache[atomIx];
            }
        }

        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& state) const 
    {
        // Compute atom positions from body positions
        Vector_<Vec3>& atomVelocityCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
            getMySubsystemIndex(), 
            atomVelocityCacheIndex
            )).upd();

        atomVelocityCache.resize( atoms.size() );

        const Vector_<Vec3>& atomBodyStationCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
            getMySubsystemIndex(), 
            atomBodyStationCacheIndex
            )).get();

        for (AtomSubsystem::LocalBodyIndex b(0); b < atomsByBody.size(); ++b) {
            const MobilizedBodyIndex& bodyIx = localBodyIxToSystemBodyIx.find(b)->second;
            const MobilizedBody& body = system.getMatterSubsystem().getMobilizedBody(bodyIx);

            const Vec3& w = body.getBodyAngularVelocity(state); // in G
            const Vec3& v = body.getBodyOriginVelocity(state);  // in G

            const std::vector<AtomSubsystem::AtomIndex>& bodyAtoms = atomsByBody[b].atoms;
            for (size_t a(0); a < bodyAtoms.size(); ++a) {
                AtomSubsystem::AtomIndex atomIx = bodyAtoms[a];
                const AtomSubsystem::Atom& atom = atoms[atomIx];

                atomVelocityCache[atomIx] = v + w % atomBodyStationCache[a];
            }
        }

        return 0;
    }

};

// Subystem for managing atom locations.
// Depends on matter subsystem

AtomSubsystem::AtomSubsystem(MultibodySystem& system) 
{
    adoptSubsystemGuts( new AtomSubsystem::Impl(system) );
    system.adoptSubsystem(*this);
}

MultibodySystem& AtomSubsystem::updMultibodySystem() {
    return updRep().system;
}
const MultibodySystem& AtomSubsystem::getMultibodySystem() const {
    return getRep().system;
}

// Fetch atom list by local body ID in this subsystem
const AtomSubsystem::AtomIndexList& AtomSubsystem::getBodyAtoms( LocalBodyIndex bodyIx ) const {
    return getRep().atomsByBody[bodyIx].atoms;
}
// Fetch atom list by body ID in MultibodySystem
const AtomSubsystem::AtomIndexList& AtomSubsystem::getBodyAtoms( MobilizedBodyIndex bodyIx ) const {
    return getBodyAtoms(getRep().systemBodyIxToLocalBodyIx.find(bodyIx)->second);
}

SimbodyMatterSubsystem& AtomSubsystem::updMatterSubsystem() {
    return updRep().system.updMatterSubsystem();
}

AtomSubsystem::AtomIndex AtomSubsystem::addAtom( mass_t mass ) 
{
    AtomSubsystem::AtomIndex atomIndex( getRep().atoms.size() );
    updRep().atoms.push_back( Atom(atomIndex, mass) );

    return atomIndex;
}

AtomSubsystem& 
AtomSubsystem::setAtomMobilizedBodyIndex(AtomSubsystem::AtomIndex atomIx, MobilizedBodyIndex bodyIx) 
{
    MobilizedBodyIndex oldBodyIx = getAtom(atomIx).getMobilizedBodyIndex();
    assert( ! oldBodyIx.isValid() );
    
    if (getRep().systemBodyIxToLocalBodyIx.find(bodyIx) == getRep().systemBodyIxToLocalBodyIx.end()) {
        LocalBodyIndex localIndex( getRep().atomsByBody.size() );
        updRep().systemBodyIxToLocalBodyIx[bodyIx] = localIndex;
        updRep().localBodyIxToSystemBodyIx[localIndex] = bodyIx;
        updRep().atomsByBody.push_back( AtomsByBody(bodyIx) );
    }
    updRep().atomsByBody[ updRep().systemBodyIxToLocalBodyIx[bodyIx] ].atoms.push_back(atomIx);
    
    updAtom(atomIx).setMobilizedBodyIndex(bodyIx);
    
    return *this;
}

size_t AtomSubsystem::getNumAtoms() const { return getRep().atoms.size(); }

const AtomSubsystem::Atom& AtomSubsystem::getAtom( AtomSubsystem::AtomIndex ix ) const {
    return getRep().atoms[ix];
}

AtomSubsystem::Atom& AtomSubsystem::updAtom( AtomSubsystem::AtomIndex ix ) {
    return updRep().atoms[ix];
}
   
AtomSubsystem::BondIndex 
    AtomSubsystem::addBond(AtomSubsystem::AtomIndex a1, AtomSubsystem::AtomIndex a2)
{
    BondIndex bondIx( getRep().bonds.size() );
    updRep().bonds.push_back( Bond(a1,a2) );
    
    updAtom(a1).addBond(a2);
    updAtom(a2).addBond(a1);
    
    return bondIx;
}

size_t AtomSubsystem::getNumBonds() const {
    return getRep().bonds.size();
}

const AtomSubsystem::Bond& 
AtomSubsystem::getBond(AtomSubsystem::BondIndex ix) const
{
    return getRep().bonds[ix];
}

AtomSubsystem::Bond& 
AtomSubsystem::updBond(AtomSubsystem::BondIndex ix)
{
    return updRep().bonds[ix];
}

size_t AtomSubsystem::getNumBodies() const {return getRep().atomsByBody.size();}

const Vec3& AtomSubsystem::getAtomLocationInGround(const SimTK::State& state, AtomSubsystem::AtomIndex atomIx) const 
{
    // Must be realized to position stage
    const Vector_<Vec3>& atomPositionCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
        getMySubsystemIndex(), 
        getRep().atomPositionCacheIndex
        )).get();

    return atomPositionCache[atomIx];
}

const Vec3& AtomSubsystem::getAtomBodyStationInGround(const SimTK::State& state, AtomSubsystem::AtomIndex atomIx) const 
{
    // Must be realized to position stage
    const Vector_<Vec3>& atomBodyStationCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
        getMySubsystemIndex(), 
        getRep().atomBodyStationCacheIndex
        )).get();

    return atomBodyStationCache[atomIx];
}

const Vec3& AtomSubsystem::getAtomVelocityInGround(const SimTK::State& state, AtomSubsystem::AtomIndex atomIx) const 
{
    // Must be realized to velocity stage
    const Vector_<Vec3>& atomVelocityCache = Value< Vector_<Vec3> >::downcast(state.updCacheEntry(
        getMySubsystemIndex(), 
        getRep().atomVelocityCacheIndex
        )).get();

    return atomVelocityCache[atomIx];
}

const AtomSubsystem::Impl& AtomSubsystem::getRep() const {
    return static_cast<const AtomSubsystem::Impl&>( getSubsystemGuts() );
}

AtomSubsystem::Impl& AtomSubsystem::updRep() {
    return static_cast<AtomSubsystem::Impl&>( updSubsystemGuts() );
}


AtomSubsystem::PairIterator AtomSubsystem::pairBegin( 
        const State& state, 
        length_t cutoff, 
        AtomSubsystem::NeighborAlgorithm algorithm ) const 
{
    PairIterator answer( *this, state, cutoff, algorithm );
    return answer;
    // return PairIterator( *this, state, cutoff, algorithm );
}

AtomSubsystem::PairIterator AtomSubsystem::pairEnd() const {
    return PairIterator();
}

void AtomSubsystem::addBodyExclusion(MobilizedBodyIndex bodyIx1, MobilizedBodyIndex bodyIx2) 
{
    assert( bodyIx1 != bodyIx2 ); // self ignoring is automatic
    assert( getRep().systemBodyIxToLocalBodyIx.find(bodyIx1) != getRep().systemBodyIxToLocalBodyIx.end() );
    assert( getRep().systemBodyIxToLocalBodyIx.find(bodyIx2) != getRep().systemBodyIxToLocalBodyIx.end() );

    LocalBodyIndex localIx1 = getRep().systemBodyIxToLocalBodyIx.find(bodyIx1)->second;
    LocalBodyIndex localIx2 = getRep().systemBodyIxToLocalBodyIx.find(bodyIx2)->second;
    
    AtomsByBody& body1 = updRep().atomsByBody[localIx1];
    AtomsByBody& body2 = updRep().atomsByBody[localIx2];

    body1.bodyExclusions.insert(bodyIx2);
    body2.bodyExclusions.insert(bodyIx1);
}
    
///////////
// Atom ///
///////////



//////////////////
// PairIterator //
//////////////////


// Automatically use voxel hash for neighbor list if numAtoms/cutoff^3 > VOXEL_HASH_CUTOFF
// TODO tune this parameter
#define VOXEL_HASH_CUTOFF (inverse_volume_t(1000.0 / cubic_nanometer))

// 0.0 means no cutoff applied
AtomSubsystem::PairIterator::PairIterator(
        const AtomSubsystem& a, 
        const State& state, 
        length_t cutoff,
        NeighborAlgorithm algorithm
    ) :
    atomSubsystem(&a),
    state(&state),
    cutoffSquared(cutoff * cutoff),
    cutoff(cutoff),
    atomsByBody( &(a.getRep().atomsByBody) ), 
    voxelHash( cutoff, a.getNumAtoms() ),
    algorithm(algorithm)
{
    // Decide which algorithm to use; populate bUseVoxelHash
    if (algorithm == N_SQUARED) bUseVoxelHash = false;
    else if (algorithm == VOXEL_HASH) bUseVoxelHash = true;
    else if (algorithm == AUTOMATIC) {
        // Automatically choose algorithm based on number of atoms and cutoff

         // no cutoff? => always use n-squared
        if (cutoff <= 0.0) 
            bUseVoxelHash = false;

        else {
            // Cleverly estimate which algorithm to use
            // TODO - tune the VOXEL_HASH_CUTOFF parameter
            inverse_volume_t measuredBlah = dimensionless_t(a.getNumAtoms()) / (cutoff*cutoff*cutoff);
            inverse_volume_t cutoffBlah = VOXEL_HASH_CUTOFF;

            if ( measuredBlah > cutoffBlah )
                bUseVoxelHash = true;

            else 
                bUseVoxelHash = false;
        }
    }
    else throw std::string("Error parsing neighbor list algorithm");

    if (bUseVoxelHash) { // O(n) method
        // pairs are disallowed within the first body, so
        // set state to just before second body, then increment
        body1 = atomsByBody->begin();
        if (atEnd()) return; // there are no bodies

        pendingBodies.insert(body1);
        atom1 = body1->atoms.end();
        currentNeighborAtoms.clear();
        atom2Index = 0;
        ++(*this);
    } 
    else { // O(n^2) method
        // Initialize to begin state
        body1 = atomsByBody->begin();
        if (body1 == atomsByBody->end())
            return; // there are no bodies

        body2 = body1; ++body2;
        if (body2 == atomsByBody->end()) {
            ++body1;
            assert (body1 == atomsByBody->end());
            return; // there is only one body
        }

        atom1 = body1->atoms.begin();
        assert(atom1 != body1->atoms.end()); // all bodies must have atoms
        currentPair.first = *atom1;

        atom2Index = 0;
        assert( atom2Index < body2->atoms.size() ); // all bodies must have atoms
        currentPair.second = body2->atoms[atom2Index];
    }
}

AtomSubsystem::PairIterator::PairIterator() 
: atomsByBody(NULL), voxelHash(1.0 * nanometers, 0)
{}

bool AtomSubsystem::PairIterator::operator!=(const PairIterator& rhs) const 
{
    if (atEnd() && rhs.atEnd()) return false;
    else if (rhs.atEnd()) return true;
    else if (atEnd()) return true;

    else if (body1 != rhs.body1) return true;
    else if (body2 != rhs.body2) return true;
    else if (atom1 != rhs.atom1) return true;
    else if (atom2Index != rhs.atom2Index) return true;

    else return false;
}

bool AtomSubsystem::PairIterator::operator==(const PairIterator& rhs) const {
    return !(*this != rhs);
}

bool AtomSubsystem::PairIterator::atEnd() const {
    if (atomsByBody == NULL) return true;
    if (body1 == atomsByBody->end()) return true;

    else return false;
}

// postfix increment, defined in terms of prefix increment
AtomSubsystem::PairIterator AtomSubsystem::PairIterator::operator++(int) {
    PairIterator tmp(*this);
    ++(*this);
    return(tmp);
}


// prefix increment
// Step to the next pair of atoms in the series.
// The meat of the neighbor list algorithm goes in here
AtomSubsystem::PairIterator& AtomSubsystem::PairIterator::operator++() 
{
    area_t zeroArea(0.0 * square_nanometers);

    if ( atEnd() ) return *this; // no more pairs

    // O(n) voxel hash method...
    if (bUseVoxelHash) 
    {
        do {
            trialIncrementAtoms();

            if ( atEnd() ) return *this;

        } while ( !atEnd() && (atom2Index >= currentNeighborAtoms.size()) );
    }
    else 
    {
        // We can either return all possible atom pairs, or just those
        // pairs that lie within a cutoff distance of one another.
        // cutoff==0.0 means no cutoff is applied at all.
        if (cutoffSquared > zeroArea) {
            // Check cutoff, if set
            area_t dSquared;
            do {
                trialIncrementAtoms();

                if ( atEnd() ) return *this;

                const Vec3& pos1 = atomSubsystem->getAtomLocationInGround(*state, *atom1);
                const Vec3& pos2 = atomSubsystem->getAtomLocationInGround(*state, body2->atoms[atom2Index]);
                dSquared = area_t( (pos1 - pos2).normSqr() * square_nanometers );
            } while ( (!atEnd()) && (dSquared >= cutoffSquared) );
        } 
        else {
            trialIncrementAtoms();
        }
    }
        
    if ( atEnd() ) return *this; // no more pairs    

    currentPair.first = *atom1;
    if (bUseVoxelHash)
        currentPair.second = currentNeighborAtoms[atom2Index];
    else
        currentPair.second = body2->atoms[atom2Index];

    return *this;
}

void AtomSubsystem::PairIterator::trialIncrementAtoms()
{
    if (bUseVoxelHash) // O(n) algorithm
    {
        // consult local list of neighbor atoms
        if ( atom2Index < currentNeighborAtoms.size() ) 
            ++atom2Index;

        if ( atom2Index >= currentNeighborAtoms.size() ) // at end
        {
            if ( atom1 != body1->atoms.end() )
                ++atom1;

            if ( atom1 == body1->atoms.end() ) {
                trialIncrementBodies();

                if (atEnd()) 
                    return;
                else 
                    atom1 = body1->atoms.begin();
            }

            currentNeighborAtoms.clear();
            voxelHash.findNeighbors(
                currentNeighborAtoms, 
                atomSubsystem->getAtomLocationInGround(*state, *atom1),
                cutoff);

            atom2Index = 0;
        }
    }
    else // O(n^2) algorithm
    {
        // Turn non-iterator for-loop inside out
        if ( atom2Index < body2->atoms.size() )
            ++atom2Index;

        if (atom2Index >= body2->atoms.size())
        {
            if ( atom1 != body1->atoms.end() )
                ++atom1;

            if (atom1 == body1->atoms.end())
            {
                trialIncrementBodies();
                if (atEnd()) return; // end of series of pairs
                
                // Check for explicit body exclusions
                while ( !atEnd() &&
                        body1->bodyExclusions.find(body2->bodyIx) != body1->bodyExclusions.end()) 
                {
                    // These two bodies are excluded from interacting, try again.
                    trialIncrementBodies();
                }
                
                if (atEnd()) return; // end of series of atom pairs

                atom1 = body1->atoms.begin();
            }

            atom2Index = 0;
        }
    }
}

// Step to the next pair of MobilizedBodies, without checking
// whether they are excluded from comparison.
void AtomSubsystem::PairIterator::trialIncrementBodies()
{
    if (bUseVoxelHash) // O(n) method
    {
        if (body1 == atomsByBody->end()) return;
        ++body1;
        if (body1 == atomsByBody->end()) return;

        // Insert previous bodies into voxel hash, unless they are in exclusion list
        std::set<AtomicBodies::const_iterator>::const_iterator b2;
        // Remember which bodies get inserted into voxel hash, 
        // and remove them from the pending list, 
        // but only AFTER the pending list iterator has gone through its motions.
        std::vector<AtomicBodies::const_iterator> unPendingBodies;
        for (b2 = pendingBodies.begin(); b2 != pendingBodies.end(); ++b2) 
        {
            MobilizedBodyIndex body2Ix = (*b2)->bodyIx;
            // If old body is NOT excluded, put it in the voxelHash
            if (body1->bodyExclusions.find(body2Ix) == body1->bodyExclusions.end()) {
                populateVoxelHash(*b2);
                unPendingBodies.push_back(*b2);
                // pendingBodies.erase(b2); // oops! I erased from container while using an iterator
            }
        }
        // Now that b2 iterator is done, remove elements from pending
        std::vector<AtomicBodies::const_iterator>::const_iterator b3;
        for (b3 = unPendingBodies.begin(); b3 != unPendingBodies.end(); ++b3)
            pendingBodies.erase(*b3);

        pendingBodies.insert(body1);
    }
    else { // O(n^2) method
        if (body2 != atomsByBody->end())
            ++body2;

        if (body2 == atomsByBody->end())
        {
            if ( body1 != atomsByBody->end() )
                ++body1;

            if ( body1 == atomsByBody->end() ) 
                return; // end of pairs

            body2 = body1;
            ++body2;

            if (body2 == atomsByBody->end()) {
                // what if body2 is at end?
                body1 = atomsByBody->end();
                return; // end of pairs
            }
        }
    }
}

void AtomSubsystem::PairIterator::populateVoxelHash(AtomicBodies::const_iterator bodyIter) {
    AtomIndexList::const_iterator a;
    for (a = bodyIter->atoms.begin(); a != bodyIter->atoms.end(); ++a) {
        voxelHash.insert(*a, atomSubsystem->getAtomLocationInGround(*state, *a) );
    }
}

const AtomSubsystem::PairIterator::Pair& AtomSubsystem::PairIterator::operator*() const {
    assert(!atEnd());
    return currentPair;
}

const AtomSubsystem::PairIterator::Pair* AtomSubsystem::PairIterator::operator->() const {
    assert(!atEnd());
    return &currentPair;
}

} // namespace SimTK
