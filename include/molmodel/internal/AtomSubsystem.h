#ifndef ATOM_SUBSYSTEM_HPP_
#define ATOM_SUBSYSTEM_HPP_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "VoxelHash.h"
#include <iterator>
#include <set>

// #include "md_units.hpp"

// Subystem for managing atom locations.
// Depends on matter subsystem

namespace SimTK {

using SimTK::units::md::nanometers;
using SimTK::units::md::daltons;

class SimTK_MOLMODEL_EXPORT AtomSubsystem : public SimTK::Subsystem 
{
protected:
    class Impl;
    const Impl& getRep() const;
    Impl& updRep();

public:
    enum NeighborAlgorithm {N_SQUARED, VOXEL_HASH, AUTOMATIC};

    typedef SimTK::units::md::mass_t mass_t;
    typedef SimTK::units::md::length_t length_t;
    typedef SimTK::units::md::location_t location_t;
    typedef SimTK::units::md::area_t area_t;

    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(AtomSubsystem, AtomIndex); // index within this subsystem
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(AtomSubsystem, BodyAtomIndex); // index within one body
    // Mobilized body index within this subsystem, not the same as MultibodySystem MobilizedBodyIndex
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(AtomSubsystem, LocalBodyIndex);
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(AtomSubsystem, BondIndex);
 
    typedef std::vector<AtomIndex> AtomIndexList;
    
    class Atom 
    {
        friend class AtomSubsystem;
        
    private:
        // required static "Topology" stage properties
        AtomIndex atomIx;
        mass_t mass; // optional, unless used to set mass of bodies
        SimTK::MobilizedBodyIndex bodyIx;
        location_t stationInBody;
        std::set<AtomIndex> bonds;

        // optional identifying information of possible use to force subsystems
        char atomType[5];
        char residueType[5];
        int residueNumber;
        bool isFirstResidue;
        bool isFinalResidue;

        // private addBond method, so users will use AtomSubsystem::addBond method instead
        Atom& addBond(AtomIndex partnerIx) {
            bonds.insert(partnerIx);

            return *this;
        }
        
        // use AtomSubsystem::setAtomMobilizedBodyIndex() instead
        Atom& setMobilizedBodyIndex(MobilizedBodyIndex ix) {
            bodyIx = ix;
            return *this;
        }
        
    public:
        Atom( 
                AtomIndex atomIx, 
                mass_t mass )
            : atomIx(atomIx), mass(mass)
        {}

        AtomIndex getAtomIndex() const {return atomIx;}

        MobilizedBodyIndex getMobilizedBodyIndex() const {return bodyIx;} 

        const Vec3& getStationInBodyFrame() const {return stationInBody;}
        Atom& setStationInBodyFrame(const Vec3& s) {
            stationInBody = s;
            return *this;
        }

        mass_t getMass() const {return mass;}
        Atom& setMass(mass_t m) {
            mass = m;
            return *this;
        }
    };

    class AtomsByBody {
    public:
        explicit AtomsByBody(MobilizedBodyIndex bodyIx) 
                : bodyIx(bodyIx)
        {}

        MobilizedBodyIndex bodyIx;
        AtomIndexList atoms;
        std::set<MobilizedBodyIndex> bodyExclusions;
    };

    typedef std::vector<AtomsByBody> AtomicBodies;

    class Bond {
    public:
        Bond(AtomIndex a1, AtomIndex a2) 
        {
            assert(a1 != a2);
            if (a1 < a2) {
                atom1Ix = a1;
                atom2Ix = a2;
            }
            else {
                atom1Ix = a2;
                atom2Ix = a1;
            }
        }
        
    private:
        AtomIndex atom1Ix;
        AtomIndex atom2Ix;
    };
    
    
    class PairIterator;
    PairIterator pairBegin(
            const State& state, 
            length_t cutoff = 0.0 * nanometers, 
            NeighborAlgorithm algorithm = AUTOMATIC) const;
    PairIterator pairEnd() const;

    explicit AtomSubsystem(SimTK::MultibodySystem& system);

    MultibodySystem& updMultibodySystem();
    const MultibodySystem& getMultibodySystem() const;

    AtomSubsystem& setAtomMobilizedBodyIndex(AtomIndex atomIx, MobilizedBodyIndex bodyIx);
    
    // Fetch atom list by local body ID in this subsystem
    const AtomIndexList& getBodyAtoms( LocalBodyIndex bodyIx ) const;
    // Fetch atom list by body ID in MultibodySystem
    const AtomIndexList& getBodyAtoms( SimTK::MobilizedBodyIndex bodyIx ) const;

    AtomIndex addAtom( mass_t mass = 0.0 * daltons );

    size_t getNumAtoms() const;

    const Vec3& getAtomLocationInGround(const SimTK::State& state, AtomIndex atom) const;
    const Vec3& getAtomVelocityInGround(const SimTK::State& state, AtomIndex atom) const;

    // useful for calculating body forces
    const Vec3& getAtomBodyStationInGround(const SimTK::State& state, AtomIndex atom) const;

    const Atom& getAtom( AtomIndex ix ) const;
    Atom& updAtom( AtomIndex ix ) ;

    BondIndex addBond(AtomIndex a1, AtomIndex a2);
    size_t getNumBonds() const;
    const Bond& getBond(BondIndex ix) const;
    Bond& updBond(BondIndex ix);
    
    size_t getNumBodies() const;

    void addBodyExclusion(MobilizedBodyIndex bodyIx1, MobilizedBodyIndex bodyIx2);

    SimbodyMatterSubsystem& updMatterSubsystem();
};

// Iterates over pairs of atoms not in the same body
// Intended to use neighbor list if enough atoms and cutoff is supplied
class AtomSubsystem::PairIterator 
        : public std::iterator<std::forward_iterator_tag, std::pair<AtomSubsystem::AtomIndex, AtomSubsystem::AtomIndex> >
{
public:
    typedef std::pair<AtomIndex, AtomIndex> Pair;

private:
    Pair currentPair;
    const AtomSubsystem* atomSubsystem;
    const State* state;
    length_t cutoff;
    area_t cutoffSquared;
    const AtomicBodies* atomsByBody;
    AtomicBodies::const_iterator body1;
    AtomicBodies::const_iterator body2;
    AtomIndexList::const_iterator atom1;
    // atom2 iterator might point to a local container, when voxel hash is used.
    // which could become invalid when this iterator is copied.
    // So use an index for atom2, instead of an iterator
    // (body1, body2, and atom1 iterators always refer to atomsByBody pointee, and thus remain valid)
    size_t atom2Index;

    NeighborAlgorithm algorithm;
    bool bUseVoxelHash;
    VoxelHash<AtomSubsystem::AtomIndex> voxelHash; // for use in O(n) neighbor list
    std::set<AtomicBodies::const_iterator> pendingBodies;
    AtomIndexList currentNeighborAtoms;
    
    void trialIncrementAtoms();
    void trialIncrementBodies();
    void populateVoxelHash(AtomicBodies::const_iterator bodyIter);
    
public:
    // 0.0 means no cutoff applied
    SimTK_MOLMODEL_EXPORT PairIterator(
            const AtomSubsystem& a, 
            const State& state, 
            length_t cutoff = 0.0 * nanometers,
            NeighborAlgorithm algorithm = AUTOMATIC
            );
    SimTK_MOLMODEL_EXPORT PairIterator(); // need default constructor for what?

    SimTK_MOLMODEL_EXPORT bool operator!=(const PairIterator& rhs) const;
    SimTK_MOLMODEL_EXPORT bool operator==(const PairIterator& rhs) const;
    SimTK_MOLMODEL_EXPORT bool atEnd() const;
    
    // postfix increment, defined in terms of prefix increment
    SimTK_MOLMODEL_EXPORT PairIterator operator++(int);
    
    // prefix increment
    // TODO - the meat of the neighbor list algorithm must go in here
    SimTK_MOLMODEL_EXPORT PairIterator& operator++();
    
    SimTK_MOLMODEL_EXPORT const Pair& operator*() const;
    SimTK_MOLMODEL_EXPORT const Pair* operator->() const;
};

} // namespace SimTK

#endif /* ATOM_SUBSYSTEM_HPP_ */
