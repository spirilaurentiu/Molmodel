/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
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
 * Implementation of DuMMForceFieldSubsystemRep and its helper classes.
 */

#include "molmodel/internal/common.h"
#include "molmodel/internal/Biotype.h"
#include "DuMMForceFieldSubsystemRep.h"

using namespace SimTK;

// #ifndef DEBUG
// #define DEBUG 1
// #endif

//#ifdef DEBUG
#define TRACE(STR) printf("%s", STR);
//#else
//#define TRACE(STR)
//#endif

// Optimize for Robosample
#include <sstream>
#include <fstream>
#include <sys/resource.h> // memory

/*! <!-- Execute a command from within (Linux free) -->
*/
std::string exec_molmodel(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
	printf ("No pipe with error: %s\n",strerror(errno));
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

/*! <!-- Get Linux memory usage -->
*/
std::size_t getLinuxMemoryUsageFromProc_m() {
    std::ifstream file("/proc/self/status");
    std::string line;
    std::size_t memoryUsage = 0;

    while (std::getline(file, line)) {
        if (line.find("VmRSS:") != std::string::npos) {
            std::istringstream iss(line);
            std::string ignore;
            iss >> ignore >> memoryUsage;
            break;
        }
    }

    return memoryUsage; // Value in kB
}

/*! <!-- Get memory usage with getrusage -->
*/
long getResourceUsage_m() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

#ifndef MEMDEBUG
#define MEMDEBUG 0
#endif

// This is Coulomb's constant 1/(4*pi*e0) in units which convert
// e^2/nm to kJ/mol.
static const Real CoulombFac = (Real)SimTK_COULOMB_CONSTANT_IN_MD;

    ////////////////////////////////////
    // DUMM FORCE FIELD SUBSYSTEM REP //
    ////////////////////////////////////

/*static*/ const char* DuMMForceFieldSubsystemRep::ApiClassName
    = "DuMMForceFieldSubsystem";

/*! <!-- Helper -->
*/
const BondStretch* DuMMForceFieldSubsystemRep::getBondStretch
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2) const
{   const AtomClassIndexPair key(class1,class2,true);
    std::map<AtomClassIndexPair,BondStretch>::const_iterator
        bs = bondStretch.find(key);
    return (bs != bondStretch.end()) ? &bs->second : 0; }

/*! <!-- Helper -->
*/
const BondBend* DuMMForceFieldSubsystemRep::getBondBend
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
    DuMM::AtomClassIndex class3) const
{   const AtomClassIndexTriple key(class1, class2, class3, true);
    std::map<AtomClassIndexTriple,BondBend>::const_iterator
        bb = bondBend.find(key);
    return (bb != bondBend.end()) ? &bb->second : 0; }

/*! <!-- Helper -->
*/
const BondTorsion* DuMMForceFieldSubsystemRep::getBondTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4) const
{   const AtomClassIndexQuad key(class1, class2, class3, class4, true);
    std::map<AtomClassIndexQuad,BondTorsion>::const_iterator
        bt = bondTorsion.find(key);
    return (bt != bondTorsion.end()) ? &bt->second : 0; }

/*! <!-- Helper -->
*/
const BondTorsion* DuMMForceFieldSubsystemRep::getAmberImproperTorsion
   (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
    DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4) const
{
//xxx -> Randy's warning flag
    bool printCrapToTheScreen = false;
    if (printCrapToTheScreen)
    {
        printf("aImp--classes: %d-%d-%d-%d\n", (int)class1,
                                      (int)class2,
                                      (int)class3,
                                      (int)class4);
        std::map<AtomClassIndexQuad,BondTorsion>::const_iterator i;
        for (i=amberImproperTorsion.begin(); i!=amberImproperTorsion.end(); i++) {
            printf( "aImp-matches: %d-%d-%d-%d\n", (int)i->first[0],
                                          (int)(i->first[1]),
                                          (int)(i->first[2]),
                                          (int)(i->first[3]) );
        }
    }

    const AtomClassIndexQuad key(class1, class2, class3, class4, false);
    std::map<AtomClassIndexQuad,BondTorsion>::const_iterator bt = amberImproperTorsion.find(key);
    return (bt != amberImproperTorsion.end()) ? &bt->second : 0;
}


//==============================================================================
//                           CLASS CrossBodyBondInfo
//==============================================================================
// This class is used locally in realizeSubsystemTopologyImpl() below to
// temporarily accumulate for each atom lists of relevant bonded connections
// that include atoms from at least two bodies, so that forces applied by
// bonded force terms can produce motion.
// Note that all indices here are AtomIndex; we will use these to construct
// the final lists that will use IncludedAtomIndex and then throw these
// away.
struct CrossBodyBondInfo {
    Array_<DuMM::AtomIndex, unsigned short>     xbond12;
    Array_<AtomIndexPair,   unsigned short>     xbond13;
    Array_<AtomIndexTriple, unsigned short>     xbond14;
    Array_<AtomIndexQuad,   unsigned short>     xbond15;

    Array_<DuMM::AtomIndex, unsigned short>     xshortPath12;
    Array_<AtomIndexPair,   unsigned short>     xshortPath13;
    Array_<AtomIndexTriple, unsigned short>     xshortPath14;
    Array_<AtomIndexQuad,   unsigned short>     xshortPath15;

    // This is even less likely to be valid than bonds3Atoms above. It will
    // be valid iff (a) bonds3Atoms is valid, and (b) at least one of the
    // three atoms is on a different body from this one.
    AtomIndexTriple                             xbonds3Atoms; 


    // GMolModel - extra 
    Array_<DuMM::AtomIndex, unsigned short>     xbond12All;
    Array_<AtomIndexPair,   unsigned short>     xbond13All;
    Array_<AtomIndexTriple, unsigned short>     xbond14All;
    Array_<AtomIndexQuad,   unsigned short>     xbond15All;

    Array_<DuMM::AtomIndex, unsigned short>     xshortPath12All;
    Array_<AtomIndexPair,   unsigned short>     xshortPath13All;
    Array_<AtomIndexTriple, unsigned short>     xshortPath14All;
    Array_<AtomIndexQuad,   unsigned short>     xshortPath15All;

    AtomIndexTriple                             xbonds3AtomsAll;

};


//------------------------------------------------------------------------------
//                             REALIZE TOPOLOGY
//------------------------------------------------------------------------------

int DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl(State& s) const
{

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory 0.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory 0.\n" << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory 0.\n" << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    //std::cout << "LAB DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl BEGIN" << std::endl;
    if(includedAtomStations.size()){

        // At realization time, we need to verify that every atom has a valid atom
        // class id. TODO: should apply only to included atoms.
        for (DuMM::AtomIndex anum(0); anum < atoms.size(); ++anum) {
            if (!isValidChargedAtomType(atoms[anum].chargedAtomTypeIndex)) {
                throw Exception::Base("Atom must have valid charged atom type before realizing topology");
            }
        }
    
        // We need to write once onto the 'cache' portion of the object once
        // the topology is known.
        DuMMForceFieldSubsystemRep* mutableThis =
            const_cast<DuMMForceFieldSubsystemRep*>(this);
    
        mutableThis->invalidateNecTopologicalCacheEntries();
    
/*  MASS PROPERTIES 
        // Set mass properties for clusters bodies. We need write access
        // to mobilized bodies.
        const MultibodySystem&        mbs    = getMultibodySystem();
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
        for (DuMM::ClusterIndex cnum(0); cnum < clusters.size(); ++cnum) {
            Cluster& c = mutableThis->clusters[cnum];
            assert(c.isValid()); // Shouldn't be any unused cluster numbers.
            if(c.isTopLevelCluster()){
                MobilizedBodyIndex mbx = c.mobodIx;
                MobilizedBody&     mobod   = matter.updMobilizedBody(mbx);
                const Transform tr;
                const MassProperties& massProperties = c.calcMassProperties(tr, *mutableThis);
                mobod.setDefaultMassProperties(massProperties);
            }
        }
*/

        // Create cache entries for storing position info and forces for included
        // atoms and included bodies.
    
        // Included atom position information is realized unconditionally when
        // realizeSubsystemPosition() is called.
        mutableThis->inclAtomStationCacheIndex  = allocateCacheEntry
           (s, Stage::Position, new Value<Vector_<Vec3> >());
        mutableThis->inclAtomPositionCacheIndex = allocateCacheEntry
           (s, Stage::Position, new Value<Vector_<Vec3> >());
    
        // Included atom velocity information is "lazy evaluated" because it
        // usually isn't need for anything. We'll realize it if someone
        // asks for it.
        mutableThis->inclAtomVelocityCacheIndex = allocateLazyCacheEntry
           (s, Stage::Velocity, // no earlier than this stage
            new Value<Vector_<Vec3> >()); // don't allocate yet
    
        // Forces and potential energy here can be calculated any time
        // after Position stage has been realized, but we won't calculate
        // them until Dynamics stage unless someone asks for them earlier.
        mutableThis->inclAtomForceCacheIndex = allocateCacheEntry
           (s, Stage::Position, Stage::Dynamics,
            new Value<Vector_<Vec3> >());
    
        // Allocate cache entry to hold rigid body moments and forces
        // on included bodies, resulting from the included atom forces above.
        mutableThis->inclBodyForceCacheIndex  = allocateCacheEntry
           (s, Stage::Position, Stage::Dynamics,
            new Value<Vector_<SpatialVec> >());
    
        mutableThis->energyCacheIndex = allocateCacheEntry
           (s, Stage::Position, Stage::Dynamics, new Value<Real>());
    
        mutableThis->internalListsRealized = true;
    
        return 0;
    }else{
        int realizeInternalListsResult = realizeInternalLists(s);
        //std::cout << "DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl END" << std::endl;
        return realizeInternalListsResult;
    }
    
    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory .\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory .\n" << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeSubsystemTopologyImpl memory .\n" << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }
}

// All the force field and molecule parameters have been set, as well as
// instructions regarding which atoms should be allowed to participate in force
// calculations. Here we precalculate everything we can that derives from these
// parameters and write it to the Topology-stage cache; i.e. write-once
// data members of this object.
// EU RESTORE
//int DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl(State& s) const
// EU BEGIN
int DuMMForceFieldSubsystemRep::realizeInternalLists(State& s) const
// EU END
{

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 0.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.\n" << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.\n" << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    // At realization time, we need to verify that every atom has a valid atom
    // class id. TODO: should apply only to included atoms.
    for (DuMM::AtomIndex anum(0); anum < atoms.size(); ++anum) {
        if (!isValidChargedAtomType(atoms[anum].chargedAtomTypeIndex)) {
            throw Exception::Base("Atom must have valid charged atom type before realizing topology");
        }
    }

    // We need to write once onto the 'cache' portion of the object once
    // the topology is known.
    DuMMForceFieldSubsystemRep* mutableThis =
        const_cast<DuMMForceFieldSubsystemRep*>(this);

    mutableThis->invalidateAllTopologicalCacheEntries();

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 0.1.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.1. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.1. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

        // force field

    // Calculate effective van der Waals parameters for all pairs of atom
    // classes. We only fill in the diagonal and upper triangle; that is, each
    // class contains parameters for like classes and all classes whose
    // (arbitrary) class number is higher.
    for (DuMM::AtomClassIndex i(0); i < atomClasses.size(); ++i) {
        if (!atomClasses[i].isValid()) continue;
        if (!atomClasses[i].isComplete()) continue;

        AtomClass& iclass = mutableThis->atomClasses[i];
        iclass.vdwDij.resize((int)atomClasses.size()-i, NaN);
        iclass.vdwEij.resize((int)atomClasses.size()-i, NaN);
        for (DuMM::AtomClassIndex j=i; j < atomClasses.size(); ++j) {
            const AtomClass& jclass = atomClasses[j];
            if (jclass.isValid() && jclass.isComplete())
                applyMixingRule(iclass.vdwRadius,    jclass.vdwRadius,
                                iclass.vdwWellDepth, jclass.vdwWellDepth,
                                iclass.vdwDij[j-i],  iclass.vdwEij[j-i]);
        }
    }

    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 0.1.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.2. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 0.2. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

        // molecule

    // Process clusters & bodies (bodies are treated as top-level clusters)

    // We process clusters recursively, so we need to allow the clusters
    // writable access to the main DuMM object (i.e., "this").
    for (DuMM::ClusterIndex cnum(0); cnum < clusters.size(); ++cnum) {
        Cluster& c = mutableThis->clusters[cnum];
        assert(c.isValid()); // Shouldn't be any unused cluster numbers.
        c.realizeTopologicalCache(*mutableThis);
    }

    // Bodies, on the other hand, are always top level clusters and the
    // calculation here assumes that all the clusters have been processed.
    // Thus bodies need only read access to the main DuMM object,
    // although we're passing the mutable one in so we can use the
    // same routine (TODO).
    for (DuMMBodyIndex bnum(0); bnum < duMMSubsetOfBodies.size(); ++bnum) {
        DuMMBody& b = mutableThis->duMMSubsetOfBodies[bnum];
        if (!b.isValid())
            continue; // OK for these to be unused.
        b.realizeTopologicalCache(*mutableThis);
    }

    // Assign body & station to every atom that has been assigned to a body.
    // At the same time we can determine whether each atom will be involved
    // in nonbond force calculations. We hold off making assignments of nonbond
    // and included atom indices because we're going to hand them out in
    // order of included bodies in the hope of getting nicer memory cache
    // behavior at run time, when we'll be processing atoms in body order.
    // We can't know all the included bodies until we're done with both nonbond
    // and bonded processing so index assignments come last.

    // Construct a local set that will eventually contain all the mobilized
    // bodies that have any included atom attached. We'll also mark atoms as
    // included or nonbonded as we process them.

    std::set<MobilizedBodyIndex> allIncludedMobods;

// for GMolModel
    std::set<MobilizedBodyIndex> allAllMobods;

    for (DuMMBodyIndex bnum(0); bnum < duMMSubsetOfBodies.size(); ++bnum) {
        const DuMMBody& b = duMMSubsetOfBodies[bnum];
        if (!b.isValid())
            continue;   // Unused body numbers are OK.
        const MobodIndex mbx = b.getMobilizedBodyIndex();
        const Cluster& cluster = getCluster(b.getClusterIndex());
        for (AtomPlacementSet::const_iterator
             app = cluster.getAllContainedAtoms().begin();
             app != cluster.getAllContainedAtoms().end();
             ++app)
        {
            const AtomPlacement& ap = *app; assert(ap.isValid());
            DuMMAtom& a = mutableThis->atoms[ap.atomIndex]; assert(a.isValid());
            assert(!a.isAttachedToBody()); // Can only be on one body!!
            a.attachToBody(mbx, ap.station);

            if (inclList.useDefaultNonbondList
             || inclList.isNonbondBody(mbx)
             || inclList.isNonbondAtom(ap.atomIndex))
            {
                allIncludedMobods.insert(mbx);
	
                // These will be replaced by the actual assignments later.
                a.inclAtomIndex    = DuMM::IncludedAtomIndex(1); // mark "included"
                a.nonbondAtomIndex = DuMM::NonbondAtomIndex(1);  // mark "nonbond"

            }

	    // for GMolModel
	    allAllMobods.insert(mbx);
	    a.AllAtomIndex    = DuMM::IncludedAtomIndex(1); // mark "All"
	    a.AllnonbondAtomIndex = DuMM::NonbondAtomIndex(1);  // mark "All nonbond"

        }
    }
    for (DuMM::AtomIndex ax(0); ax < atoms.size(); ++ax) {
        assert(getAtom(ax).isAttachedToBody()); // TODO catch unassigned atoms
    }

    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 1.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 1.\n" << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 1.\n" << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    //------- Process bonds -------
    // Now we're going to look at each atom again and find all interesting
    // bonded connections for which an atom serves as atom 1 in a 1-2, 1-2-3,
    // 1-2-3-4, 1-2-3-4-5 sequence or 2-3-1-4 improper torsion term. These
    // connections may be used for two purposes: generation of bonded forces,
    // and scaling of nearby atoms when performing nonbond calculations.
    // Normally we are only interested in bonds that cross bodies; neither
    // bonded nor nonbonded forces are interesting if they involve only atoms
    // from the same body.
    //
    // There are several subtleties here. Bonded force terms are generated
    // for *all* possible connections between atoms, while scaling is
    // based only on the *shortest* possible connection. So if atoms A and B are
    // connected both directly by bond A-B and torsion A-c-d-B we want to
    // include both for force terms, but when we're processing atom A's nonbond
    // interactions B should be scaled based only on the 1-2 connection.
    // So we construct separate lists for all bonds and shortest bonds only,
    // using the former for forces and the latter for scaling.
    //
    // We don't know yet which bonds we're going to keep. For scaling, we
    // only need those shortest, cross-body connections for which our
    // atom 1 is a nonbond atom and so is the atom at the other end of
    // the bonded connection, that is both 1-2, 1-3, 1-4, or 1-5 atoms are
    // nonbond atoms (and we just determined above which atoms are nonbond).
    // So we don't need to compute shortest connections lists at all for an
    // atom that isn't involved in nonbond calculations.
    //
    // For bond force connections, we need keep only those cross-body bonds
    // which have been designated as included by one of the following
    // conditions:
    //    - any one if the bond's atoms has been marked as one all of whose
    //      bonds are to be included
    //    - any pair of atoms in the bond has been marked as a pair all of
    //      whose bonded interactions are to be included (these might not be
    //      adjacent in the bonded sequence)
    //    - similar conditions for the bodies to which the atoms are attached
    // Note that this can drag in more atoms than were explicitly included.
    // For example say we have 1-2-3-4 torsion bond. If only atom 3 has been
    // designated "include my bonds", then atoms 1,2, and 4 also become
    // included, though only for the purpose of bonds including atom 3.
    //
    // Finally, any of the sequences 1-2 through 1-...-5 will come up twice,
    // once when we process the first atom and once for the last. We don't
    // want to double-count so we only keep the one in which the atom index
    // of the first atom is smaller than the atom index of the last atom.
    // (This does not apply to improper torsions.)
    //
    // After having decided which bonds to keep for atom i, if there are
    // any remaining for which atom i serves as atom 1 then we add i to
    // the list of bondStartAtoms; those are the only atoms for which we
    // need to do bond processing.


    // This is a group of lists which identify atoms nearby in the molecule's
    // bond structure. Each atom's bond12 list already contains the directly
    // bonded (1-2) atoms; the 13 list below has the 1-(2)-3 bonded atoms (that
    // is, it includes the path to the "3" atom), etc. The current Atom is
    // always atom "1" so it isn't stored.
    //
    // Note that the shortPath and xshortPath arrays give the shortest path
    // between two atoms, while the bond and xbond arrays give *all*
    // connection paths, with bonds3Atoms giving at most one. Forces must
    // be calculated for all paths between two atoms, but scaling is only
    // done according to the shortest path.

    Array_<AtomIndexPair,  unsigned short>      bond13;
    Array_<AtomIndexTriple,unsigned short>      bond14;
    Array_<AtomIndexQuad,  unsigned short>      bond15;
    Array_<AtomIndexPair,  unsigned short>      shortPath13;
    Array_<AtomIndexTriple,unsigned short>      shortPath14;
    Array_<AtomIndexQuad,  unsigned short>      shortPath15;

    // This will be invalid unless we find that the current atom is directly
    // bonded to exactly three other atoms, in which case their atom indices
    // will be stored here and isValid() will return true.
    AtomIndexTriple                             bonds3Atoms;

    // This set is used to avoid duplicate paths in the shortestPath
    // calculation.
    std::set<DuMM::AtomIndex>                   allBondedSoFar;

    // We have to save these for each atom during a first pass, then
    // we'll use them to fill in the per-atom bond force and scaling
    // arrays in a second pass where the included atom indices are known.
    Array_<CrossBodyBondInfo, DuMM::AtomIndex>  crossBodyBondInfo(atoms.size());

    // need to chase bonds to fill in the bonded data
    // Be sure to distinguish the *shortest* path between two atoms from
    // the set of all paths between atoms.
    for (DuMM::AtomIndex anum(0); anum < atoms.size(); ++anum) {
        DuMMAtom&          a = mutableThis->atoms[anum];
        CrossBodyBondInfo& x = crossBodyBondInfo[anum];

        // Make sure all the reusable temporary lists defined above are
        // cleared before their first use for the current atom a.

        // Only the bond12 list has been filled in so far. We'll sort
        // all the lists when they're done for good hygiene.
        std::sort(a.bond12.begin(), a.bond12.end());

        allBondedSoFar.clear();
        // Add this atom and its direct (1-2) bonds to the list of all bonded
        // atoms.
        allBondedSoFar.insert(anum);
        allBondedSoFar.insert(a.bond12.begin(), a.bond12.end());

        // Find longer bond paths by building each list in turn from
        // the direct bonds of the atoms in the previous list.

        // build the bond13 and shortPath13 lists
        // - bond1x list gives *all* paths between bonded atoms where all the
        // atoms are distinct (i.e., no fair retracing one of the bonds or
        // running around a short loop to get back to the first atom again).
        // - shortPath1x list gives *shortest* path between bonded atoms
        bond13.clear();
        shortPath13.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j) {
            const DuMMAtom&       a12 = atoms[a.bond12[j]];
            const ShortAtomArray& a12_12 = a12.bond12;
            for (int k=0; k < (int)a12_12.size(); ++k) {
                const DuMM::AtomIndex newAtom = a12_12[k];
                assert(newAtom != a.bond12[j]);
                if (newAtom == anum)
                    continue; // no loop backs!
                bond13.push_back(AtomIndexPair(a.bond12[j], newAtom));

                // if no shorter path, note this short route
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    shortPath13.push_back(AtomIndexPair(a.bond12[j], newAtom));
                }
            }
        }
        std::sort(bond13.begin(), bond13.end());
        std::sort(shortPath13.begin(), shortPath13.end());

        // Randy was too big of a sissy to combine the bond14 and shortPath14
        // computations! Or, discretion is sometimes the better part of valor.

        // build the bond14 list (all non-overlapping, non-looped paths)
        bond14.clear();
        for (int j=0; j < (int)bond13.size(); ++j) {
            const DuMMAtom&       a13 = atoms[bond13[j][1]];
            const ShortAtomArray& a13_12 = a13.bond12;
            for (int k=0; k < (int)a13_12.size(); ++k) {
                const DuMM::AtomIndex newAtom = a13_12[k];
                assert(newAtom != bond13[j][1]);
                // avoid repeated atoms (loop back)
                if (newAtom!=anum && newAtom!=bond13[j][0]) {
                    bond14.push_back(AtomIndexTriple(bond13[j][0],
                                                 bond13[j][1], newAtom));
                }
            }
        }
        std::sort(bond14.begin(), bond14.end());

        // build the shortPath14 list
        shortPath14.clear();
        for (int j=0; j < (int)shortPath13.size(); ++j) {
            const DuMMAtom&       a13 = atoms[shortPath13[j][1]];
            const ShortAtomArray& a13_12 = a13.bond12;
            for (int k=0; k < (int)a13_12.size(); ++k) {
                const DuMM::AtomIndex newAtom = a13_12[k];

                 // check if there was already a shorter path
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    shortPath14.push_back(AtomIndexTriple(shortPath13[j][0],
                                                shortPath13[j][1], newAtom));
                }
            }
        }
        std::sort(shortPath14.begin(), shortPath14.end());


        // build the bond15 list
        bond15.clear();
        for (int j=0; j < (int)bond14.size(); ++j) {
            const DuMMAtom&       a14    = atoms[bond14[j][2]];
            const ShortAtomArray& a14_12 = a14.bond12;
            for (int k=0; k < (int)a14_12.size(); ++k) {
                const DuMM::AtomIndex newAtom = a14_12[k];
                assert(newAtom != bond14[j][2]);

                // avoid repeats and loop back
                if (newAtom!=anum && newAtom!=bond14[j][0] && newAtom!=bond14[j][1]) {
                    bond15.push_back(AtomIndexQuad(bond14[j][0],
                                               bond14[j][1],
                                               bond14[j][2], newAtom));
                }
            }
        }
        std::sort(bond15.begin(), bond15.end());

        // build the shortPath15 list
        shortPath15.clear();
        for (int j=0; j < (int)shortPath14.size(); ++j) {
            const DuMMAtom&       a14    = atoms[shortPath14[j][2]];
            const ShortAtomArray& a14_12 = a14.bond12;
            for (int k=0; k < (int)a14_12.size(); ++k) {
                const DuMM::AtomIndex newAtom = a14_12[k];

                // check if there was already a shorter path
                if (allBondedSoFar.find(newAtom) == allBondedSoFar.end()) {
                    allBondedSoFar.insert(newAtom);
                    shortPath15.push_back(AtomIndexQuad(shortPath14[j][0],
                                                    shortPath14[j][1],
                                                    shortPath14[j][2], newAtom));
                }
            }
        }
        std::sort(shortPath15.begin(), shortPath15.end());

        // Find all atom that are connected to three (and only three) other
        // atoms, then add all orderings of this to the improper torsion list.
        bonds3Atoms.invalidate();
        if (a.bond12.size() == 3) {
            bonds3Atoms = AtomIndexTriple(a.bond12[0], a.bond12[1], a.bond12[2]);
        }

        // Fill in the cross-body bond lists. We only keep bonds that involve
        // atoms which are not all attached to the same body. Also, we throw
        // away any bond sequences for which the atom index of atom 1 is
        // greater than the atom index of the last atom. That prevents
        // double counting since the bond sequence will show up in both orders
        // eventually. (This doesn't apply to improper torsions.)

        // TODO: need a way to weed out cross-body bonded terms when the
        // mobilities available prevent any use of that term. For example, if
        // there is only a torsion dof between two bodies, and it is aligned
        // with atoms A and B, then a bond stretch term between A and B
        // can't do anything and shouldn't be kept on the list below.

        x.xbond12.clear(); x.xbond12All.clear();
        for (int j=0; j < (int)a.bond12.size(); ++j) {
            if (anum > a.bond12[j]) continue;
            if (atoms[a.bond12[j]].mobodIx != a.mobodIx)
                x.xbond12.push_back(a.bond12[j]); 
            x.xbond12All.push_back( a.bond12[j] );
        }

        x.xbond13.clear(); x.xbond13All.clear();
        for (int j=0; j < (int)bond13.size(); ++j) {
            if (anum > bond13[j][1]) continue;
            if (   atoms[bond13[j][0]].mobodIx != a.mobodIx
                || atoms[bond13[j][1]].mobodIx != a.mobodIx)
                x.xbond13.push_back(bond13[j]);
            x.xbond13All.push_back( bond13[j] );
        }

        x.xbond14.clear(); x.xbond14All.clear();
        for (int j=0; j < (int)bond14.size(); ++j) {
            if (anum > bond14[j][2]) continue;
            if (   atoms[bond14[j][0]].mobodIx != a.mobodIx
                || atoms[bond14[j][1]].mobodIx != a.mobodIx
                || atoms[bond14[j][2]].mobodIx != a.mobodIx)
                x.xbond14.push_back(bond14[j]);
            x.xbond14All.push_back( bond14[j] );
        }

        x.xbond15.clear(); x.xbond15All.clear();
        for (int j=0; j < (int)bond15.size(); ++j) {
            if (anum > bond15[j][3]) continue;
            if (   atoms[bond15[j][0]].mobodIx != a.mobodIx
                || atoms[bond15[j][1]].mobodIx != a.mobodIx
                || atoms[bond15[j][2]].mobodIx != a.mobodIx
                || atoms[bond15[j][3]].mobodIx != a.mobodIx)
                x.xbond15.push_back(bond15[j]);
	    x.xbond15All.push_back( bond15[j] );

        }

        x.xbonds3Atoms.invalidate(); x.xbonds3AtomsAll.invalidate();
        // If there were exactly 3 bonds, and at least one of them is
        // on a different body, then we win!
        if (bonds3Atoms.isValid() &&
            (   atoms[bonds3Atoms[0]].mobodIx != a.mobodIx
             || atoms[bonds3Atoms[1]].mobodIx != a.mobodIx
             || atoms[bonds3Atoms[2]].mobodIx != a.mobodIx)) {
                 x.xbonds3Atoms = bonds3Atoms;
             }
	    if (bonds3Atoms.isValid()) {
            x.xbonds3AtomsAll = bonds3Atoms;
        }

        // By default, or if this atom or its body are on the "must include"
        // list then we have to keep all the cross-body bonds we just
        // discovered. Otherwise, we only keep the ones for which some other
        // atoms or bodies provides the necessity.
        if (!(   inclList.useDefaultBondList
              || inclList.isBondAtom(anum)
              || inclList.isBondBody(a.mobodIx)))
        {
            for (int j=0; j < (int)x.xbond12.size();) {
                bool keep = false;
                const DuMM::AtomIndex bnum = x.xbond12[j];
                const MobodIndex      mbx  = atoms[bnum].mobodIx;
                     if (inclList.isBondAtom(bnum)) keep=true;
                else if (inclList.isBondBody(mbx))  keep=true;
                else if (inclList.isBondAtomPair(anum, bnum))    keep=true;
                else if (inclList.isBondBodyPair(a.mobodIx, mbx)) keep=true;
                if (keep) ++j;
                else {
			x.xbond12.erase(&x.xbond12[j]); // don't increment j		
		    }
            }

            for (int j=0; j < (int)x.xbond13.size();) {
                bool keep = false;
                for (int k=0; k < 2 && !keep; ++k) {
                    const DuMM::AtomIndex bnum = x.xbond13[j][k];
                    const MobodIndex      mbx  = atoms[bnum].mobodIx;
                         if (inclList.isBondAtom(bnum)) keep=true;
                    else if (inclList.isBondBody(mbx))  keep=true;
                    else if (inclList.isBondAtomPair(anum, bnum))    keep=true;
                    else if (inclList.isBondBodyPair(a.mobodIx, mbx)) keep=true;
                    for (int c=k+1; c < 2 && !keep; ++c) {
                        const DuMM::AtomIndex cnum = x.xbond13[j][c];
                        const MobodIndex      mcx  = atoms[cnum].mobodIx;
                             if (inclList.isBondAtomPair(bnum, cnum)) keep=true;
                        else if (inclList.isBondBodyPair(mbx, mcx))   keep=true;
                    }
                }
                if (keep) ++j;
                else x.xbond13.erase(&x.xbond13[j]); // don't increment j
            }

            for (int j=0; j < (int)x.xbond14.size();) {
                bool keep = false;
                for (int k=0; k < 3 && !keep; ++k) {
                    const DuMM::AtomIndex bnum = x.xbond14[j][k];
                    const MobodIndex      mbx  = atoms[bnum].mobodIx;
                         if (inclList.isBondAtom(bnum)) keep=true;
                    else if (inclList.isBondBody(mbx))  keep=true;
                    else if (inclList.isBondAtomPair(anum, bnum))    keep=true;
                    else if (inclList.isBondBodyPair(a.mobodIx, mbx)) keep=true;
                    for (int c=k+1; c < 3 && !keep; ++c) {
                        const DuMM::AtomIndex cnum = x.xbond14[j][c];
                        const MobodIndex      mcx  = atoms[cnum].mobodIx;
                             if (inclList.isBondAtomPair(bnum, cnum)) keep=true;
                        else if (inclList.isBondBodyPair(mbx, mcx))   keep=true;
                    }
                }
                if (keep) ++j;
                else x.xbond14.erase(&x.xbond14[j]); // don't increment j
            }

            for (int j=0; j < (int)x.xbond15.size();) {
                bool keep = false;
                for (int k=0; k < 4 && !keep; ++k) {
                    const DuMM::AtomIndex bnum = x.xbond15[j][k];
                    const MobodIndex      mbx  = atoms[bnum].mobodIx;
                         if (inclList.isBondAtom(bnum)) keep=true;
                    else if (inclList.isBondBody(mbx))  keep=true;
                    else if (inclList.isBondAtomPair(anum, bnum))    keep=true;
                    else if (inclList.isBondBodyPair(a.mobodIx, mbx)) keep=true;
                    for (int c=k+1; c < 4 && !keep; ++c) {
                        const DuMM::AtomIndex cnum = x.xbond15[j][c];
                        const MobodIndex      mcx  = atoms[cnum].mobodIx;
                             if (inclList.isBondAtomPair(bnum, cnum)) keep=true;
                        else if (inclList.isBondBodyPair(mbx, mcx))   keep=true;
                    }
                }
                if (keep) ++j;
                else x.xbond15.erase(&x.xbond15[j]); // don't increment j
            }

            if (x.xbonds3Atoms.isValid()) { // there is only one of these
                bool keep = false;
                for (int k=0; k < 3 && !keep; ++k) {
                    const DuMM::AtomIndex bnum = x.xbonds3Atoms[k];
                    const MobodIndex      mbx  = atoms[bnum].mobodIx;
                         if (inclList.isBondAtom(bnum)) keep=true;
                    else if (inclList.isBondBody(mbx))  keep=true;
                    else if (inclList.isBondAtomPair(anum, bnum))    keep=true;
                    else if (inclList.isBondBodyPair(a.mobodIx, mbx)) keep=true;
                    for (int c=k+1; c < 3 && !keep; ++c) {
                        const DuMM::AtomIndex cnum = x.xbonds3Atoms[c];
                        const MobodIndex      mcx  = atoms[cnum].mobodIx;
                             if (inclList.isBondAtomPair(bnum, cnum)) keep=true;
                        else if (inclList.isBondBodyPair(mbx, mcx))   keep=true;
                    }
                }
                if (!keep) x.xbonds3Atoms.invalidate();
            }
        }

        // At this point the cross-body bond arrays contain only bonds that we
        // are keeping. Every atom mentioned in any of the kept bonds should be
        // made an included atom, and each of those atoms' mobilized bodies
        // should be made an included body. The present atom "a" and its body
        // are included if we kept any bonds at all, and we also note that
        // atom "a" is a "bond starter" atom, that is, it serves as atom 1
        // for some bond.
        if (x.xbond12.size() || x.xbond13.size() || x.xbond14.size()
         || x.xbond15.size() || x.xbonds3Atoms.isValid())
        {
            // We kept at least one bond.
            a.inclAtomIndex = DuMM::IncludedAtomIndex(1); // just mark it
            allIncludedMobods.insert(a.mobodIx);
            a.bondStarterIndex = DuMMBondStarterIndex(1); // mark

            for (int j=0; j < (int)x.xbond12.size(); ++j) {
                DuMMAtom& b = mutableThis->atoms[x.xbond12[j]];
                b.inclAtomIndex = DuMM::IncludedAtomIndex(1);
                allIncludedMobods.insert(b.mobodIx);
            }

            for (int j=0; j < (int)x.xbond13.size(); ++j)
                for (int k=0; k < 2; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond13[j][k]];
                    b.inclAtomIndex = DuMM::IncludedAtomIndex(1);
                    allIncludedMobods.insert(b.mobodIx);
                }

            for (int j=0; j < (int)x.xbond14.size(); ++j)
                for (int k=0; k < 3; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond14[j][k]];
                    b.inclAtomIndex = DuMM::IncludedAtomIndex(1);
                    allIncludedMobods.insert(b.mobodIx);
                }

            for (int j=0; j < (int)x.xbond15.size(); ++j)
                for (int k=0; k < 4; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond15[j][k]];
                    b.inclAtomIndex = DuMM::IncludedAtomIndex(1);
                    allIncludedMobods.insert(b.mobodIx);
                }

            if (x.xbonds3Atoms.isValid())
                for (int k=0; k < 3; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbonds3Atoms[k]];
                    b.inclAtomIndex = DuMM::IncludedAtomIndex(1);
                    allIncludedMobods.insert(b.mobodIx);
                }
        }


	    // Adding the same for GMolModel for All bonds.
        a.AllAtomIndex = DuMM::IncludedAtomIndex(1); // just mark it
        allAllMobods.insert(a.mobodIx);
        a.AllbondStarterIndex = DuMMBondStarterIndex(1); // mark

	    for (int j=0; j < (int)x.xbond12All.size(); ++j) {
                DuMMAtom& b = mutableThis->atoms[x.xbond12All[j]];
                b.AllAtomIndex = DuMM::IncludedAtomIndex(1);
                allAllMobods.insert(b.mobodIx);
        	}

        for (int j=0; j < (int)x.xbond13All.size(); ++j)
                for (int k=0; k < 2; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond13All[j][k]];
                    b.AllAtomIndex = DuMM::IncludedAtomIndex(1);
                    allAllMobods.insert(b.mobodIx);
        	}

        for (int j=0; j < (int)x.xbond14All.size(); ++j)
                for (int k=0; k < 3; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond14All[j][k]];
                    b.AllAtomIndex = DuMM::IncludedAtomIndex(1);
                    allAllMobods.insert(b.mobodIx);
                }

        for (int j=0; j < (int)x.xbond15All.size(); ++j)
                for (int k=0; k < 4; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbond15All[j][k]];
                    b.AllAtomIndex = DuMM::IncludedAtomIndex(1);
                    allAllMobods.insert(b.mobodIx);
                }

        if (x.xbonds3AtomsAll.isValid())
                for (int k=0; k < 3; ++k) {
                    DuMMAtom& b = mutableThis->atoms[x.xbonds3AtomsAll[k]];
                    b.AllAtomIndex = DuMM::IncludedAtomIndex(1);
                    allAllMobods.insert(b.mobodIx);
                }





        // Next we're going to work on the shortest-path interconnections
        // that will be used for nonbond scaling. We'll include
        // any shortest path connections that cross a body, provided that
        // both the first (i.e. current atom "a") and last atom are nonbond
        // atoms. Note that we want these bond paths to show up
        // in both orders so that we can do nonbond scaling from either end.

        x.xshortPath12.clear(); x.xshortPath13.clear();
        x.xshortPath14.clear(); x.xshortPath15.clear();

        x.xshortPath12All.clear(); x.xshortPath13All.clear();
        x.xshortPath14All.clear(); x.xshortPath15All.clear();

        if (a.isNonbondAtom()) {
            for (int j=0; j < (int)a.bond12.size(); ++j) {
                if (!atoms[a.bond12[j]].isNonbondAtom()) continue;
                if (atoms[a.bond12[j]].mobodIx != a.mobodIx)
                    x.xshortPath12.push_back(a.bond12[j]);
		        x.xshortPath12All.push_back( a.bond12[j] );
            }

            for (int j=0; j < (int)shortPath13.size(); ++j) {
                if (!atoms[shortPath13[j][1]].isNonbondAtom()) continue;
                if (   atoms[shortPath13[j][0]].mobodIx != a.mobodIx
                    || atoms[shortPath13[j][1]].mobodIx != a.mobodIx) {
                        x.xshortPath13.push_back(shortPath13[j]);
                    }
	            x.xshortPath13All.push_back( shortPath13[j] );
            }


            for (int j=0; j < (int)shortPath14.size(); ++j) {
                if (!atoms[shortPath14[j][2]].isNonbondAtom()) continue;
                if (   atoms[shortPath14[j][0]].mobodIx != a.mobodIx
                    || atoms[shortPath14[j][1]].mobodIx != a.mobodIx
                    || atoms[shortPath14[j][2]].mobodIx != a.mobodIx) {
                        x.xshortPath14.push_back(shortPath14[j]);
                    }
	            x.xshortPath14All.push_back(shortPath14[j]);
            }

            for (int j=0; j < (int)shortPath15.size(); ++j) {
                if (!atoms[shortPath15[j][3]].isNonbondAtom()) continue;
                if (   atoms[shortPath15[j][0]].mobodIx != a.mobodIx
                    || atoms[shortPath15[j][1]].mobodIx != a.mobodIx
                    || atoms[shortPath15[j][2]].mobodIx != a.mobodIx
                    || atoms[shortPath15[j][3]].mobodIx != a.mobodIx)
                    x.xshortPath15.push_back(shortPath15[j]);
                x.xshortPath15All.push_back(shortPath15[j]);
            }
        }
    }

    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 5.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 5. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 5. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    // We have processed all the atoms and marked them included if they
    // will appear in any nonbonded or bonded force calculation. The
    // nonbond atoms have been separately marked. We have also created a
    // set allIncludedMobods that contains the mobilized body index of
    // every body that has at least one included atom attached to it.
    // Now we want to build data structures suitable for high-speed
    // processing at run time. We're going to run through the included bodies
    // in order, look at the attached atoms, and assign included atom and
    // nonbond atom indices in order so that they are grouped by included body.
    // Then we can just keep (start,length) pairs with each included body
    // to locate its subset of included, nonbond, and bondstart atoms.

  // Build a temporary map of included bodies to their included atoms.
    std::map<MobodIndex, Array_<DuMM::AtomIndex> > inclBods2Atoms;
    unsigned numIncludedAtoms = 0;
    for (DuMM::AtomIndex ax(0); ax < atoms.size(); ++ax) {
        const DuMMAtom& atom = atoms[ax];
        if (!atom.isIncludedAtom()) continue;
        assert(atom.isAttachedToBody());
        inclBods2Atoms[atom.getMobodIndex()].push_back(ax);
        ++numIncludedAtoms;
    }


  // same for GMolModel 
    std::map<MobodIndex, Array_<DuMM::AtomIndex> > AllBods2Atoms;
    unsigned numAllAtoms = 0;
    for (DuMM::AtomIndex ax(0); ax < atoms.size(); ++ax) {
        const DuMMAtom& atom = atoms[ax];
        AllBods2Atoms[atom.getMobodIndex()].push_back(ax);
        ++numAllAtoms;
    }


    SimTK_ASSERT_ALWAYS(inclBods2Atoms.size() == allIncludedMobods.size(),
        "DuMMForceFieldSubsystem::realizeTopology(): inconsistent lists "
        "of included bodies generated.");
  
    // Allocate with default construction the included bodies and included
    // atoms arrays. We'll fill them in below.
    mutableThis->includedBodies.resize((unsigned)allIncludedMobods.size());
    mutableThis->includedAtoms.resize(numIncludedAtoms);
    mutableThis->includedAtomStations.resize(numIncludedAtoms);

    DuMMIncludedBodyIndex   nextInclBodyIx(0);
    DuMM::IncludedAtomIndex nextInclAtomIx(0);

    for (std::set<MobodIndex>::const_iterator
    allIncludedMobodsIt = allIncludedMobods.begin();
    allIncludedMobodsIt != allIncludedMobods.end();
    ++allIncludedMobodsIt, ++nextInclBodyIx){

        // Set included bodies from allIncludedMobods
        IncludedBody& inclBody = mutableThis->includedBodies[nextInclBodyIx];
        inclBody.mobodIx = *allIncludedMobodsIt;

        inclBody.beginIncludedAtoms    = inclBody.endIncludedAtoms =
            DuMM::IncludedAtomIndex(nextInclAtomIx);
        inclBody.beginNonbondAtoms     = inclBody.endNonbondAtoms =
            DuMM::NonbondAtomIndex(nonbondAtoms.size());
        inclBody.beginBondStarterAtoms = inclBody.endBondStarterAtoms =
            DuMMBondStarterIndex(bondStarterAtoms.size());

        const Array_<DuMM::AtomIndex>& inclBodAtoms = inclBods2Atoms[*allIncludedMobodsIt];

        for (unsigned i=0; i < inclBodAtoms.size(); ++i, ++nextInclAtomIx) {

            DuMMAtom& dummAtom = mutableThis->atoms[inclBodAtoms[i]];
            assert(dummAtom.atomIndex == inclBodAtoms[i]);
            assert(dummAtom.inclAtomIndex.isValid()); // should have been marked "1"

            dummAtom.inclAtomIndex = nextInclAtomIx;
            dummAtom.inclBodyIndex = nextInclBodyIx;
            IncludedAtom& includedAtom = mutableThis->includedAtoms[nextInclAtomIx];
            Vec3&         station_B = mutableThis->includedAtomStations[nextInclAtomIx];
            includedAtom.setIndices(dummAtom.inclAtomIndex, dummAtom.inclBodyIndex,
                          dummAtom.atomIndex, dummAtom.chargedAtomTypeIndex);
            station_B = dummAtom.station_B;
            ++inclBody.endIncludedAtoms;

            if (dummAtom.nonbondAtomIndex.isValid()) {
                dummAtom.nonbondAtomIndex = DuMM::NonbondAtomIndex(nonbondAtoms.size());
                mutableThis->nonbondAtoms.push_back(dummAtom.inclAtomIndex);
                ++inclBody.endNonbondAtoms;
            }

            if (dummAtom.bondStarterIndex.isValid()) {
                dummAtom.bondStarterIndex = DuMMBondStarterIndex(bondStarterAtoms.size());
                mutableThis->bondStarterAtoms.push_back(dummAtom.inclAtomIndex);
                ++inclBody.endBondStarterAtoms;
            }

            // std::cout << "SP_NEW_LAB i dAIx inclDAIx nbDAIx " << i <<" "
            //     << dummAtom.atomIndex <<" " << dummAtom.inclAtomIndex <<" "
            //     << dummAtom.nonbondAtomIndex << std::endl;

        }


    }

  // same for GMolModel
    mutableThis->AllBodies.resize((unsigned)allAllMobods.size());
    mutableThis->AllAtoms.resize(numAllAtoms);
    mutableThis->AllAtomStations.resize(numAllAtoms);

    DuMMIncludedBodyIndex   nextAllBodyIx(0);
    DuMM::IncludedAtomIndex nextAllAtomIx(0);
    for (std::set<MobodIndex>::const_iterator p = allAllMobods.begin();
         p != allAllMobods.end(); ++p, ++nextAllBodyIx)
    {
        IncludedBody& AllBody = mutableThis->AllBodies[nextAllBodyIx];
        AllBody.mobodIx = *p;
        AllBody.beginAllAtoms    = AllBody.endAllAtoms =
            DuMM::IncludedAtomIndex(nextAllAtomIx);
        AllBody.beginAllNonbondAtoms     = AllBody.endAllNonbondAtoms =
            DuMM::NonbondAtomIndex(AllnonbondAtoms.size());
        AllBody.beginAllBondStarterAtoms = AllBody.endAllBondStarterAtoms =
            DuMMBondStarterIndex(AllbondStarterAtoms.size());

        const Array_<DuMM::AtomIndex>& AllBodAtoms = AllBods2Atoms[*p];
        for (unsigned i=0; i < AllBodAtoms.size(); ++i, ++nextAllAtomIx) {
            DuMMAtom& a = mutableThis->atoms[AllBodAtoms[i]];
            assert(a.atomIndex == AllBodAtoms[i]);
            assert(a.AllAtomIndex.isValid()); // should have been marked "1"
            a.AllAtomIndex = nextAllAtomIx;
            a.AllBodyIndex = nextAllBodyIx;

            IncludedAtom& ia = mutableThis->AllAtoms[nextAllAtomIx];

            Vec3&         station_B_All = mutableThis->AllAtomStations[nextAllAtomIx];
            ia.setIndices(a.AllAtomIndex, a.AllBodyIndex,
                          a.atomIndex, a.chargedAtomTypeIndex);
            station_B_All = a.station_B;

            ++AllBody.endAllAtoms;

            if (a.AllnonbondAtomIndex.isValid()) {
                a.AllnonbondAtomIndex = DuMM::NonbondAtomIndex(AllnonbondAtoms.size());
                mutableThis->AllnonbondAtoms.push_back(a.AllAtomIndex);
                ++AllBody.endAllNonbondAtoms;

            }
            if (a.AllbondStarterIndex.isValid()) {
                a.AllbondStarterIndex = DuMMBondStarterIndex(AllbondStarterAtoms.size());
                mutableThis->AllbondStarterAtoms.push_back(a.AllAtomIndex);
                ++AllBody.endAllBondStarterAtoms;
            }
        }

    }


    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 6.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 6. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 6. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }


    // Now that included atom index assignments have been made, we can
    // allocate the includedAtoms array and fill each IncludedAtom with
    // bonded force arrays that use included atom indices rather than full
    // atom indices. Similarly, we can create nonbond scaling arrays that
    // use nonbond atom indices.

    for (DuMM::IncludedAtomIndex iax(0); iax < includedAtoms.size(); ++iax) {
        const DuMM::AtomIndex ax = getAtomIndexOfIncludedAtom(iax);
        const CrossBodyBondInfo& x = crossBodyBondInfo[ax]; // computed above
        IncludedAtom& ia = mutableThis->updIncludedAtom(iax);

        ia.scale12.clear(); ia.scale13.clear();
        ia.scale14.clear(); ia.scale15.clear();

        for (int j=0; j < (int)x.xshortPath12.size(); ++j) {
            const DuMM::AtomIndex a2x = x.xshortPath12[j];
            assert(atoms[a2x].nonbondAtomIndex.isValid());
            ia.scale12.push_back(atoms[a2x].nonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath13.size(); ++j) {
            const DuMM::AtomIndex a3x = x.xshortPath13[j][1];
            assert(atoms[a3x].nonbondAtomIndex.isValid());
            ia.scale13.push_back(atoms[a3x].nonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath14.size(); ++j) {
            const DuMM::AtomIndex a4x = x.xshortPath14[j][2];
            assert(atoms[a4x].nonbondAtomIndex.isValid());
            ia.scale14.push_back(atoms[a4x].nonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath15.size(); ++j) {
            const DuMM::AtomIndex a4x = x.xshortPath15[j][3];
            assert(atoms[a4x].nonbondAtomIndex.isValid());
            ia.scale15.push_back(atoms[a4x].nonbondAtomIndex);
        }

        const DuMM::AtomClassIndex c1 = getAtomClassIndex(ax);

        // Save a BondStretch entry for each cross-body 1-2 bond
        ia.force12.resize(x.xbond12.size());
        ia.stretch.resize(x.xbond12.size());


        for (int b12=0; b12 < (int)x.xbond12.size(); ++b12) {
            const DuMM::AtomIndex bx = x.xbond12[b12];
            ia.force12[b12] = atoms[bx].inclAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx);
            ia.stretch[b12] = getBondStretch(c1, c2);

            SimTK_REALIZECHECK2_ALWAYS(ia.stretch[b12],
                Stage::Topology, getMySubsystemIndex(), getName(),
                "Couldn't find bond stretch parameters for included "
                "cross-body atom class pair (%d,%d).", (int)c1, (int)c2);
        }

        // Save a BondBend entry for each cross-body 1-3 bond
        ia.force13.resize(x.xbond13.size());
        ia.bend.resize(x.xbond13.size());


        for (int b13=0; b13 < (int)x.xbond13.size(); ++b13) {
            const AtomIndexPair& bx = x.xbond13[b13];
            ia.force13[b13][0] = atoms[bx[0]].inclAtomIndex;
            ia.force13[b13][1] = atoms[bx[1]].inclAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx[0]);
            const DuMM::AtomClassIndex c3 = getAtomClassIndex(bx[1]);
            ia.bend[b13] = getBondBend(c1, c2, c3);

            SimTK_REALIZECHECK3_ALWAYS(ia.bend[b13],
                Stage::Topology, getMySubsystemIndex(), getName(),
                "Couldn't find bond bend parameters for included "
                "cross-body atom class triple (%d,%d,%d).",
                (int)c1, (int)c2, (int)c3);
        }

        // Save a BondTorsion entry for each cross-body 1-4 bond
        ia.force14.resize(x.xbond14.size());
        ia.torsion.resize(x.xbond14.size());

        for (int b14=0; b14 < (int)x.xbond14.size(); ++b14) {
            const AtomIndexTriple& bx = x.xbond14[b14];
            ia.force14[b14][0] = atoms[bx[0]].inclAtomIndex;
            ia.force14[b14][1] = atoms[bx[1]].inclAtomIndex;
            ia.force14[b14][2] = atoms[bx[2]].inclAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx[0]);
            const DuMM::AtomClassIndex c3 = getAtomClassIndex(bx[1]);
            const DuMM::AtomClassIndex c4 = getAtomClassIndex(bx[2]);
            ia.torsion[b14] = getBondTorsion(c1, c2, c3, c4);

            SimTK_REALIZECHECK4_ALWAYS(ia.torsion[b14],
                Stage::Topology, getMySubsystemIndex(), getName(),
                "Couldn't find bond torsion parameters for included "
                "cross-body atom class quad (%d,%d,%d,%d).",
                (int)c1, (int)c2, (int)c3, (int)c4);
        }


        // Save *all* Amber improper torsion entries if this atom is bonded to
        // three, and only three other atoms, *and* a matching Amber improper
        // torsion term is found in the amberImproperTorsion array. Note that
        // by convention, the center atom is in the third position. Also note
        // that unlike Amber, which keeps only *one* match, we keep *all*.
        // To correct for this we also scale my the total number of matches.
        // This is how Tinker implements Amber's improper torsions.
        ia.aImproperTorsion.clear();   // the BondTorsion term
        ia.forceImproper14.clear(); // the other three atoms
        if (x.xbonds3Atoms.isValid()) {
            for (int i2=0; i2<3; i2++) {
                for (int i3=0; i3<3; i3++) {
                    if (i3==i2) continue;
                    for (int i4=0; i4<3; i4++) {
                        if (i4==i2 || i4==i3) continue;
                        const AtomIndexTriple bx(x.xbonds3Atoms[i2],
                                                 x.xbonds3Atoms[i3],
                                                 x.xbonds3Atoms[i4]);

                        // Not heap allocated; just a reference if non-null.
                        const BondTorsion* bt = getAmberImproperTorsion(
                                        getAtomClassIndex(bx[0]),
                                        getAtomClassIndex(bx[1]),
                                        c1,
                                        getAtomClassIndex(bx[2]));
                        if (bt) {
                            ia.forceImproper14.push_back(IncludedAtomIndexTriple(
                                                atoms[bx[0]].inclAtomIndex,
                                                atoms[bx[1]].inclAtomIndex,
                                                atoms[bx[2]].inclAtomIndex));
                            ia.aImproperTorsion.push_back(bt);
                        }
                    }
                }
            }
        }
    }


    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 7.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 7. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 7. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }



// GMolModel - Same for AllAtomIndex ....this can be optimised

    for (DuMM::IncludedAtomIndex iax(0); iax < AllAtoms.size(); ++iax) {
        const DuMM::AtomIndex ax = getAtomIndexOfAllAtom(iax);
        const CrossBodyBondInfo& x = crossBodyBondInfo[ax]; // computed above
        IncludedAtom& ia = mutableThis->updAllAtom(iax);

        ia.scale12All.clear(); ia.scale13All.clear();
        ia.scale14All.clear(); ia.scale15All.clear();

        for (int j=0; j < (int)x.xshortPath12All.size(); ++j) {
            const DuMM::AtomIndex a2x = x.xshortPath12All[j];
            assert(atoms[a2x].AllnonbondAtomIndex.isValid());
            ia.scale12All.push_back(atoms[a2x].AllnonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath13All.size(); ++j) {
            const DuMM::AtomIndex a3x = x.xshortPath13All[j][1];
            assert(atoms[a3x].AllnonbondAtomIndex.isValid());
            ia.scale13All.push_back(atoms[a3x].AllnonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath14All.size(); ++j) {
            const DuMM::AtomIndex a4x = x.xshortPath14All[j][2];
            assert(atoms[a4x].AllnonbondAtomIndex.isValid());
            ia.scale14All.push_back(atoms[a4x].AllnonbondAtomIndex);
        }
        for (int j=0; j < (int)x.xshortPath15All.size(); ++j) {
            const DuMM::AtomIndex a4x = x.xshortPath15All[j][3];
            assert(atoms[a4x].AllnonbondAtomIndex.isValid());
            ia.scale15All.push_back(atoms[a4x].AllnonbondAtomIndex);
        }


        const DuMM::AtomClassIndex c1 = getAtomClassIndex(ax);

        // Save a BondStretch entry for each cross-body 1-2 bond
        ia.force12All.resize(x.xbond12All.size());
        ia.stretchAll.resize(x.xbond12All.size());

        for (int b12=0; b12 < (int)x.xbond12All.size(); ++b12) {
            const DuMM::AtomIndex bx = x.xbond12All[b12];
            ia.force12All[b12] = atoms[bx].AllAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx);
            ia.stretchAll[b12] = getBondStretch(c1, c2);
        }



        // Save a BondBend entry for each cross-body 1-3 bond

        ia.force13All.resize(x.xbond13All.size());
        ia.bendAll.resize(x.xbond13All.size());
        for (int b13=0; b13 < (int)x.xbond13All.size(); ++b13) {
            const AtomIndexPair& bx = x.xbond13All[b13];
            ia.force13All[b13][0] = atoms[bx[0]].AllAtomIndex;
            ia.force13All[b13][1] = atoms[bx[1]].AllAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx[0]);
            const DuMM::AtomClassIndex c3 = getAtomClassIndex(bx[1]);
            ia.bendAll[b13] = getBondBend(c1, c2, c3);
        }


        // Save a BondTorsion entry for each cross-body 1-4 bond
        ia.force14All.resize(x.xbond14All.size());
        ia.torsionAll.resize(x.xbond14All.size());
        for (int b14=0; b14 < (int)x.xbond14All.size(); ++b14) {
            const AtomIndexTriple& bx = x.xbond14All[b14];
            ia.force14All[b14][0] = atoms[bx[0]].AllAtomIndex;
            ia.force14All[b14][1] = atoms[bx[1]].AllAtomIndex;
            ia.force14All[b14][2] = atoms[bx[2]].AllAtomIndex;

            const DuMM::AtomClassIndex c2 = getAtomClassIndex(bx[0]);
            const DuMM::AtomClassIndex c3 = getAtomClassIndex(bx[1]);
            const DuMM::AtomClassIndex c4 = getAtomClassIndex(bx[2]);
            ia.torsionAll[b14] = getBondTorsion(c1, c2, c3, c4);
        }

        ia.aImproperTorsionAll.clear();   // the BondTorsion term
        ia.forceImproper14All.clear(); // the other three atoms
        if (x.xbonds3AtomsAll.isValid()) {
            for (int i2=0; i2<3; i2++) {
                for (int i3=0; i3<3; i3++) {
                    if (i3==i2) continue;
                    for (int i4=0; i4<3; i4++) {
                        if (i4==i2 || i4==i3) continue;
                        const AtomIndexTriple bx(x.xbonds3AtomsAll[i2],
                                                 x.xbonds3AtomsAll[i3],
                                                 x.xbonds3AtomsAll[i4]);

                        const BondTorsion* bt = getAmberImproperTorsion(
                                        getAtomClassIndex(bx[0]),
                                        getAtomClassIndex(bx[1]),
                                        c1,
                                        getAtomClassIndex(bx[2]));
                        if (bt) {
                            ia.forceImproper14All.push_back(IncludedAtomIndexTriple(
                                                atoms[bx[0]].AllAtomIndex,
                                                atoms[bx[1]].AllAtomIndex,
                                                atoms[bx[2]].AllAtomIndex));
                            ia.aImproperTorsionAll.push_back(bt);
                        }
                    }
                }
            }
        }
    }

    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 8.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }



        /////////////////////////////
        // Fill in GBSA parameters //
        /////////////////////////////

    if (gbsaGlobalScaleFactor != 0 && getNumNonbondAtoms()) {
        // These parameters will be used either by our local GBSA implemention or
        // by OpenMM's if we're using that. Note that we're working only with
        // nonbond atoms so must refer to them by nonbond index rather than
        // atom index or included atom index.
        mutableThis->gbsaAtomicPartialCharges.resize(getNumNonbondAtoms());
        mutableThis->gbsaAtomicNumbers.resize(getNumNonbondAtoms());
        mutableThis->gbsaNumberOfCovalentBondPartners.resize(getNumNonbondAtoms());
        mutableThis->atomicNumberOfHCovalentPartner.resize(getNumNonbondAtoms());

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 8.1. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8.1. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

        for (DuMM::NonbondAtomIndex nbx(0); nbx < getNumNonbondAtoms(); ++nbx) {
            const DuMMAtom& atom = getAtom(getAtomIndexOfNonbondAtom(nbx));
            const ChargedAtomType& atype      = chargedAtomTypes[atom.chargedAtomTypeIndex];
            const Real             charge     = atype.partialCharge;
            const AtomClass&       aclass     = atomClasses[atype.atomClassIx];
            const int              element    = aclass.element;

            mutableThis->gbsaAtomicPartialCharges[nbx]         = charge;
            mutableThis->gbsaAtomicNumbers[nbx]                = element;
            // Note that the number of covalent partners is topological information
            // needed by GBSA to identify the right parameters to use for this atom;
            // it does *not* depend on whether we're crossing bodies, or whether
            // covalent partners are themselves nonbonded atoms.
            mutableThis->gbsaNumberOfCovalentBondPartners[nbx] = atom.bond12.size();
            // This is used only for hydrogens which will only have one bond.
            mutableThis->atomicNumberOfHCovalentPartner[nbx] = -1;
            if (atom.bond12.size()==1) {
                const DuMMAtom& partner = getAtom(atom.bond12[0]);
                const ChargedAtomType& ptype =
                    chargedAtomTypes[partner.chargedAtomTypeIndex];
                const AtomClass& pclass = atomClasses[ptype.atomClassIx];
                mutableThis->atomicNumberOfHCovalentPartner[nbx] = pclass.element;
            }
        }

        // look up obc scale factor for each atom
        mutableThis->gbsaObcScaleFactors.resize(getNumNonbondAtoms());
        int returnValue = getObcScaleFactors(getNumNonbondAtoms(),
                                             &gbsaAtomicNumbers.front(),
                                             &mutableThis->gbsaObcScaleFactors.front());
        SimTK_ASSERT_ALWAYS(returnValue == 0, "Couldn't get GBSA scale factors.");

        mutableThis->gbsaRadii.resize(getNumNonbondAtoms());
        returnValue = getGbsaRadii(getNumNonbondAtoms(),
                                   &gbsaAtomicNumbers.front(),
                                   &gbsaNumberOfCovalentBondPartners.front(),
                                   &atomicNumberOfHCovalentPartner.front(),
                                   &mutableThis->gbsaRadii.front());
        SimTK_ASSERT_ALWAYS(returnValue == 0, "Couldn't get GBSA input radii.");

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 8.2. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8.2. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

        // Don't delete this object here.
        ObcParameters* obcParameters =
            new ObcParameters(getNumNonbondAtoms(), ObcParameters::ObcTypeII);
        obcParameters->setScaledRadiusFactors( &gbsaObcScaleFactors.front() );
        obcParameters->setAtomicRadii(&gbsaRadii.front(),
                                      SimTKOpenMMCommon::KcalAngUnits);
        obcParameters->setSolventDielectric(gbsaSolventDielectric);
        obcParameters->setSoluteDielectric(gbsaSoluteDielectric);
        mutableThis->gbsaCpuObc = new CpuObc(obcParameters); // CpuObc takes ownership of the parameters object
        gbsaCpuObc->setIncludeAceApproximation((int)gbsaIncludeAceApproximation);

        // The cpu GBSA requires vectors of pointers to locations and forces.
        // Allocate those now. Note that forces have to be zeroed each step.
        gbsaRawCoordinates.resize(3 * getNumNonbondAtoms(), NaN); // [x,y,z,x,y,z,...], Angstrom(!) units
        gbsaAtomicForces.resize  (3 * getNumNonbondAtoms(), NaN);
        mutableThis->gbsaCoordinatePointers.resize (getNumNonbondAtoms()); // [&x0,&x1,&x2...]
        mutableThis->gbsaAtomicForcePointers.resize(getNumNonbondAtoms()); // [&x0,&x1,&x2...]

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 8.3. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8.3. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

        // Load the pointers -- after this they don't change.
        for (DuMM::NonbondAtomIndex a(0); a < getNumNonbondAtoms(); ++a) {
            mutableThis->gbsaCoordinatePointers[a]  = &gbsaRawCoordinates[3*a];
            mutableThis->gbsaAtomicForcePointers[a] = &gbsaAtomicForces[3*a];
        }
    }

        ///////////////////////////////////////////
        // Initialize OpenMM if it is being used //
        ///////////////////////////////////////////

    // OpenMM has to be enabled at run time to be used. And if
    // all we can find is a slow reference implementation we still won't use it
    // unless the reference platform has been explicitly allowed.

    // Using "while" here just so we can break out; this won't ever loop. If we
    // decide to use OpenMM, the flag usingOpenMM will be set true.
    mutableThis->usingOpenMM = false;
    if (wantOpenMMAcceleration ) {
        mutableThis->openMMPlatformInUse = mutableThis->openMMPlugin.initializeOpenMM(allowOpenMMReference, this);

        if (openMMPlatformInUse.empty()) {
            if (tracing)
                std::cout << "WARNING: DuMM: failed to initialize OpenMM\n";
        }

        if (tracing)
            std::cout << "NOTE: DuMM: using OpenMM platform '" << openMMPlatformInUse << "'\n";

        mutableThis->usingOpenMM = true;
    }

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory 8.4. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 8.4. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    if (!usingOpenMM) {
        // If the caller specified how many threads to use, even if only one,
        // and says "useMultithreadedComputation" then we will use the
        // multithreaded code.
        const bool wantParallel =
            useMultithreadedComputation
            && (   numThreadsRequested > 0
                || ParallelExecutor::getNumProcessors() > 1);

        const bool anyNonbonded = getNumNonbondAtoms()
                                  && (   coulombGlobalScaleFactor !=0
                                      || vdwGlobalScaleFactor     !=0
                                      || gbsaGlobalScaleFactor    !=0);

        if (wantParallel) {
            mutableThis->usingMultithreaded = true;
            if (!anyNonbonded) {
                mutableThis->usingMultithreaded = false;
                if (tracing)
                    // EU BEGIN COMMENT
                    std::clog << "NOTE: DuMM: not using multithreading because"
                                 " there are no nonbonded or implicit solvent"
                                 " terms to calculate.\n";

            }
            // This will probably never happen.
            if (ParallelExecutor::isWorkerThread()) {
                mutableThis->usingMultithreaded = false;
                if (tracing)
                    // EU BEGIN COMMENT
                    std::clog << "NOTE: DuMM: can't use multithreading because"
                                 " the main thread is already a ParallelExecutor"
                                 " worker thread.\n";
            }
        }
    }

    if(MEMDEBUG){
        //std::cout << "DuMMRep::realizeInternalLists memory 9.\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 9. " << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory 9. " << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }


    // If we're not using openMM, but we are using multithreaded computation, create
    // Parallel2DExecutors for parallelizing expensive force calculations.
    mutableThis->numThreadsInUse = 1;
    if (usingMultithreaded) {
        const int numThreadsWanted = numThreadsRequested > 0
                            ? numThreadsRequested
                            : ParallelExecutor::getNumProcessors();
        mutableThis->numThreadsInUse = std::min(numThreadsWanted,
                                                std::max(1, getNumNonbondAtoms()/2));

        if (tracing && (numThreadsInUse < numThreadsWanted))
            std::clog << "NOTE: DuMM: reduced number of threads from "
                      << numThreadsWanted << " to " << numThreadsInUse
                      << " because there were only " << getNumNonbondAtoms()
                      << " atoms included in nonbonded force calculations.\n";

        //mutableThis->executor = new ParallelExecutor(numThreadsInUse);

        if (tracing)
            std::clog << "NOTE: DuMM: using multithreading code with "
                      << numThreadsInUse << " threads.\n";

        mutableThis->nonbondedExecutor =
            //new Parallel2DExecutor(includedBodies.size(), *executor);
            new Parallel2DExecutor(includedBodies.size(), numThreadsInUse);
        mutableThis->gbsaExecutor =
            //new Parallel2DExecutor(getNumNonbondAtoms(), *executor);
            new Parallel2DExecutor(getNumNonbondAtoms(), numThreadsInUse);

        // gmolmodel
        mutableThis->NonbondedFullExecutor =
                new Parallel2DExecutor( AllBodies.size(), numThreadsInUse) ;
    }

    if (!(usingOpenMM || usingMultithreaded)) {
        if (tracing)
            // EU BEGIN COMMENT
            std::clog << "NOTE: DuMM: using single threaded code.\n";


        // Using single threaded -- allocate global temporaries
        vdwScaleSingleThread.resize(getNumNonbondAtoms(), Real(1));
        coulombScaleSingleThread.resize(getNumNonbondAtoms(), Real(1));

    }

    //GMolModel
    vdwScaleAllSingleThread.resize(getNumAllNonbondAtoms(), Real(1));
    coulombScaleAllSingleThread.resize(getNumAllNonbondAtoms(), Real(1));


    // Create cache entries for storing position info and forces for included
    // atoms and included bodies.

    // Included atom position information is realized unconditionally when
    // realizeSubsystemPosition() is called.
    mutableThis->inclAtomStationCacheIndex  = allocateCacheEntry
       (s, Stage::Position, new Value<Vector_<Vec3> >());
    mutableThis->inclAtomPositionCacheIndex = allocateCacheEntry
       (s, Stage::Position, new Value<Vector_<Vec3> >());

    // Included atom velocity information is "lazy evaluated" because it
    // usually isn't need for anything. We'll realize it if someone
    // asks for it.
    mutableThis->inclAtomVelocityCacheIndex = allocateLazyCacheEntry
       (s, Stage::Velocity, // no earlier than this stage
        new Value<Vector_<Vec3> >()); // don't allocate yet

    // Forces and potential energy here can be calculated any time
    // after Position stage has been realized, but we won't calculate
    // them until Dynamics stage unless someone asks for them earlier.
    mutableThis->inclAtomForceCacheIndex = allocateCacheEntry
       (s, Stage::Position, Stage::Dynamics,
        new Value<Vector_<Vec3> >());

    // Allocate cache entry to hold rigid body moments and forces
    // on included bodies, resulting from the included atom forces above.
    mutableThis->inclBodyForceCacheIndex  = allocateCacheEntry
       (s, Stage::Position, Stage::Dynamics,
        new Value<Vector_<SpatialVec> >());

    mutableThis->energyCacheIndex = allocateCacheEntry
       (s, Stage::Position, Stage::Dynamics, new Value<Real>());

    mutableThis->internalListsRealized = true;

    if(MEMDEBUG){
        std::cout << "DuMMRep::realizeInternalLists memory .\n" << exec_molmodel("free") << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory .\n" << getLinuxMemoryUsageFromProc_m() << " kB" << std::endl << std::flush;
        std::cout << "DuMMRep::realizeInternalLists memory .\n" << getResourceUsage_m() << " kB" << std::endl << std::flush;
    }

    return 0;
}
//.............................REALIZE TOPOLOGY.................................


//------------------------------------------------------------------------------
//                             REALIZE POSITION
//------------------------------------------------------------------------------
// Here we calculate for every included atom: (1) its body station re-expressed
// in Ground, and (2) its absolute Cartesian position in Ground. Forces and
// energy *could* be calculated here, but we'll hold off until we get to
// Dynamics stage unless someone asks for them earlier.
// Cost is 18 flops per atom plus bookkeeping.
int DuMMForceFieldSubsystemRep::
realizeSubsystemPositionImpl(const State& s) const {
    const MultibodySystem&        mbs    = getMultibodySystem();
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    Vector_<Vec3>& inclAtomStation_G = updIncludedAtomStationCache(s);
    Vector_<Vec3>& inclAtomPos_G     = updIncludedAtomPositionCache(s);

    inclAtomStation_G.resize(getNumIncludedAtoms());
    inclAtomPos_G.resize(getNumIncludedAtoms());

    for (DuMMIncludedBodyIndex dbx(0); dbx < includedBodies.size(); ++dbx) {
        const IncludedBody&      inclBod = includedBodies[dbx];
        const MobilizedBodyIndex mbx     = inclBod.mobodIx;
        const MobilizedBody&     mobod   = matter.getMobilizedBody(mbx);

        const Transform&    X_GB  = mobod.getBodyTransform(s);
        const Rotation&     R_GB = X_GB.R();
        const Vec3&         p_GB = X_GB.p();

        for (DuMM::IncludedAtomIndex iax=inclBod.beginIncludedAtoms;
             iax != inclBod.endIncludedAtoms; ++iax)
        {
            const Vec3& station_B = getIncludedAtomStation(iax);

            // atomic coordinates with respect to Ground frame
            const Vec3 p_BS_G = R_GB * station_B;  // 15 flops
            inclAtomStation_G[iax] = p_BS_G;
            inclAtomPos_G[iax]     = p_GB + p_BS_G;  //  3 flops
        }
    }

    return 0;
}
//.............................REALIZE POSITION.................................



//------------------------------------------------------------------------------
//                        REALIZE ATOM VELOCITY CACHE
//------------------------------------------------------------------------------
// Compute included atom velocities from body velocities, if they haven't
// already been computed since velocities were last changed.
// Cost is 12 flops per atom plus bookkeeping.
void DuMMForceFieldSubsystemRep::
realizeIncludedAtomVelocityCache(const State& state) const {
    if (isIncludedAtomVelocityCacheRealized(state))
        return; // nothing to do

    const MultibodySystem&        mbs    = getMultibodySystem();
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    const Vector_<Vec3>& inclAtomBodyStation_G = getIncludedAtomStationsInG(state);

    Vector_<Vec3>& inclAtomVelocityCache = updIncludedAtomVelocityCache(state);
    inclAtomVelocityCache.resize(getNumAtoms());

    for (DuMMIncludedBodyIndex dbx(0); dbx < getNumIncludedBodies(); ++dbx) {
        const IncludedBody&      inclBod = includedBodies[dbx];
        const MobilizedBodyIndex mbx     = inclBod.mobodIx;
        const MobilizedBody&     mobod   = matter.getMobilizedBody(mbx);

        const SpatialVec&   V_GB = mobod.getBodyVelocity(state);
        const Vec3&         w    = V_GB[0]; // angular velocity
        const Vec3&         v    = V_GB[1]; // linear velocity

        for (DuMM::IncludedAtomIndex iax=inclBod.beginIncludedAtoms;
             iax != inclBod.endIncludedAtoms; ++iax)
        {
            inclAtomVelocityCache[iax] =
                v + w % inclAtomBodyStation_G[iax]; // 12 flops
        }
    }

    markIncludedAtomVelocityCacheRealized(state);

}
//...........................REALIZE ATOM VELOCITY CACHE........................



//------------------------------------------------------------------------------
//                             REALIZE VELOCITY
//------------------------------------------------------------------------------
int DuMMForceFieldSubsystemRep::realizeSubsystemVelocityImpl(const State& state) const
{
    // NOTHING TO DO
    return 0;
}
//.............................REALIZE VELOCITY.................................





//------------------------------------------------------------------------------
//                             CALC BOND STRETCH
//------------------------------------------------------------------------------
// Helper routine for realizeDynamics. This should be called only once it has
// been determined that we need to calculate bond stretch terms. We're given an
// included atom 1 and we process all the bond stretch terms involving that
// atom and one other immediately-adjacent atom. To avoid double-counting the
// 1-2 bond array force12 for this atom only contain bonds for which this
// atom's original atom index was less than the atom index of the other atom.
// (The included atom index values that are in force12 won't necessarily be
// ordered the same way.)
void DuMMForceFieldSubsystemRep::calcBondStretch
   (DuMM::IncludedAtomIndex             a1num,
    const Vector_<Vec3>&                inclAtomStation_G,
    const Vector_<Vec3>&                inclAtomPos_G,
    Real                                bondStretchScaleFactor,
    Real                                customBondStretchScaleFactor,
    Vector_<SpatialVec>&                inclBodyForces_G,
    Real&                               energy) const
{

    const IncludedAtom&           a1    = getIncludedAtom(a1num);
    const DuMMIncludedBodyIndex   b1    = a1.inclBodyIndex;
    const Vec3& a1Station_G = inclAtomStation_G[a1num];
    const Vec3& a1Pos_G     = inclAtomPos_G[a1num];

    for (DuMM::IncludedAtomIndex b12(0); b12 < a1.force12.size(); ++b12) {
        const DuMM::IncludedAtomIndex a2num = a1.force12[b12];

        const IncludedAtom&  a2 = getIncludedAtom(a2num);
        const Vec3& a2Station_G = inclAtomStation_G[a2num];
        const Vec3& a2Pos_G     = inclAtomPos_G[a2num];
        const Vec3  r           = a2Pos_G - a1Pos_G;
        const Real  d           = r.norm();

        const BondStretch& bs = *a1.stretch[b12];

        Real eStretch, fStretch;
        if (bs.hasBuiltinTerm()) {
            const Real x = d - bs.d0;
            const Real kx = bondStretchScaleFactor * bs.k * x;
            eStretch =  kx*x; // no factor of 1/2!
            fStretch = -2*kx; // sign is as would be applied to a2
        } else
            eStretch = fStretch = 0;

        for (int i=0; i < (int)bs.customTerms.size(); ++i) {
            const DuMM::CustomBondStretch& term = *bs.customTerms[i];
            eStretch += customBondStretchScaleFactor * term.calcEnergy(d);
            // expecting f = -dE/dx but can't check
            fStretch += customBondStretchScaleFactor * term.calcForce(d);
        }

        // Force is normally directed along the vector from atom 1
        // to atom 2. But if the atoms are exactly on top of one another
        // we're just going to use an arbitrary direction in the hope
        // that this is just some relaxation where the only thing that
        // matters is that the atoms do separate.
        const Vec3 f2 = (d==0 ? Vec3(fStretch,0,0) : (fStretch/d)*r);

        const DuMMIncludedBodyIndex b2 = a2.inclBodyIndex;
 
	assert(b2 != b1);
        energy += eStretch;

        //TRACE("calcBondStretch: ");
        //TRACE((std::to_string(energy)).c_str());
        //TRACE("\n");

        inclBodyForces_G[b2] += SpatialVec( a2Station_G % f2, f2); // 15 flops
        inclBodyForces_G[b1] -= SpatialVec( a1Station_G % f2, f2); // 15 flops

    } 
}
//.............................CALC BOND STRETCH................................



//------------------------------------------------------------------------------
//                              CALC BOND BEND
//------------------------------------------------------------------------------
// Helper routine for realizeDynamics. This should be called only
// once it has been determined that we need to calculate bond bending terms.
// We're given an atom 1 and we process all the bond bending forces
// involving that atom and two others. To avoid double-counting the 1-3 bond
// array force13 for this atom only contain bonds for which this atom's original
// atom index was less than the atom index of the final atom. (The included atom
// index values that are in force13 won't necessarily be ordered the same way.)
void DuMMForceFieldSubsystemRep::calcBondBend
   (DuMM::IncludedAtomIndex             a1num,
    const Vector_<Vec3>&                inclAtomStation_G,
    const Vector_<Vec3>&                inclAtomPos_G,
    Real                                bondBendScaleFactor,
    Real                                customBondBendScaleFactor,
    Vector_<SpatialVec>&                inclBodyForces_G,
    Real&                               energy) const
{
    const IncludedAtom&           a1    = getIncludedAtom(a1num);
    const DuMMIncludedBodyIndex   b1    = a1.inclBodyIndex;
    const Vec3& a1Station_G = inclAtomStation_G[a1num];
    const Vec3& a1Pos_G     = inclAtomPos_G[a1num];

    for (int b13=0; b13 < (int)a1.force13.size(); ++b13) {
        const DuMM::IncludedAtomIndex a2num = a1.force13[b13][0];
        const DuMM::IncludedAtomIndex a3num = a1.force13[b13][1];

        const IncludedAtom& a2 = getIncludedAtom(a2num);
        const IncludedAtom& a3 = getIncludedAtom(a3num);

        // TODO: These might be the same body but for now we don't care.
        const Vec3& a2Station_G = inclAtomStation_G[a2num];
        const Vec3& a3Station_G = inclAtomStation_G[a3num];
        const Vec3& a2Pos_G     = inclAtomPos_G[a2num];
        const Vec3& a3Pos_G     = inclAtomPos_G[a3num];

        Real angle, e;
        Vec3 f1, f2, f3;
        const BondBend& bb = *a1.bend[b13];

        // atom 2 is the central one
        bb.calculateAtomForces(a2Pos_G, a1Pos_G, a3Pos_G,
                               bondBendScaleFactor, customBondBendScaleFactor,
                               angle, e, f2, f1, f3);

        const DuMMIncludedBodyIndex b2 = a2.inclBodyIndex;
        const DuMMIncludedBodyIndex b3 = a3.inclBodyIndex;

	    assert(!(b2==b1 && b3==b1)); // shouldn't be on the list if all on 1 body

        energy += e;

        //TRACE("calcBondBend: ");
        //TRACE((std::to_string(e)).c_str());
        //TRACE("\n");

        inclBodyForces_G[b1] += SpatialVec( a1Station_G % f1, f1); // 15 flops
        inclBodyForces_G[b2] += SpatialVec( a2Station_G % f2, f2); // 15 flops
        inclBodyForces_G[b3] += SpatialVec( a3Station_G % f3, f3); // 15 flops
    }
}
//..............................CALC BOND BEND..................................



//------------------------------------------------------------------------------
//                              CALC BOND TORSION
//------------------------------------------------------------------------------
// Helper routine for realizeDynamics. This should be called only
// once it has been determined that we need to calculate bond torsion terms.
// We're given an atom 1 and we process all the bond torsion forces
// involving that atom and three others. To avoid double-counting the 1-4 bond
// array force14 for this atom only contain bonds for which this atom's original
// atom index was less than the atom index of the final atom. (The included atom
// index values that are in force14 won't necessarily be ordered the same way.)
void DuMMForceFieldSubsystemRep::calcBondTorsion
   (DuMM::IncludedAtomIndex             a1num,
    const Vector_<Vec3>&                inclAtomStation_G,
    const Vector_<Vec3>&                inclAtomPos_G,
    Real                                bondTorsionScaleFactor,
    Real                                customBondTorsionScaleFactor,
    Vector_<SpatialVec>&                inclBodyForces_G,
    Real&                               energy) const
{
    const IncludedAtom&           a1    = getIncludedAtom(a1num);
    const DuMMIncludedBodyIndex   b1    = a1.inclBodyIndex;
    const Vec3& a1Station_G = inclAtomStation_G[a1num];
    const Vec3& a1Pos_G     = inclAtomPos_G[a1num];

    for (int b14=0; b14 < (int)a1.force14.size(); ++b14) {
        const DuMM::IncludedAtomIndex a2num = a1.force14[b14][0];
        const DuMM::IncludedAtomIndex a3num = a1.force14[b14][1];
        const DuMM::IncludedAtomIndex a4num = a1.force14[b14][2];

        const IncludedAtom& a2 = getIncludedAtom(a2num);
        const IncludedAtom& a3 = getIncludedAtom(a3num);
        const IncludedAtom& a4 = getIncludedAtom(a4num);

        const Vec3& a2Station_G = inclAtomStation_G[a2num];
        const Vec3& a3Station_G = inclAtomStation_G[a3num];
        const Vec3& a4Station_G = inclAtomStation_G[a4num];
        const Vec3& a2Pos_G     = inclAtomPos_G[a2num];
        const Vec3& a3Pos_G     = inclAtomPos_G[a3num];
        const Vec3& a4Pos_G     = inclAtomPos_G[a4num];

        Real angle, e;
        Vec3 f1, f2, f3, f4;
        const BondTorsion& bt = *a1.torsion[b14];
        bt.calculateAtomForces
           (a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G,
            bondTorsionScaleFactor, customBondTorsionScaleFactor,
            angle, e, f1, f2, f3, f4);

        const DuMMIncludedBodyIndex b2 = a2.inclBodyIndex;
        const DuMMIncludedBodyIndex b3 = a3.inclBodyIndex;
        const DuMMIncludedBodyIndex b4 = a4.inclBodyIndex;

	    assert(!(b2==b1 && b3==b1 && b4==b1)); // shouldn't be on the list if all on 1 body
    
        energy += e;

        //TRACE("calcBondTorsion: ");
        //TRACE((std::to_string(e)).c_str());
        //TRACE("\n");

        inclBodyForces_G[b1] += SpatialVec( a1Station_G % f1, f1); // 15 flops
        inclBodyForces_G[b2] += SpatialVec( a2Station_G % f2, f2); // 15 flops
        inclBodyForces_G[b3] += SpatialVec( a3Station_G % f3, f3); // 15 flops
        inclBodyForces_G[b4] += SpatialVec( a4Station_G % f4, f4); // 15 flops
    }
}
//..............................CALC BOND TORSION...............................



//------------------------------------------------------------------------------
//                         CALC AMBER IMPROPER TORSION
//------------------------------------------------------------------------------
// Helper routine for realizeDynamics. This should be called only once it has
// been determined that we need to calculate Amber improper torsion terms.
// We're given an atom 1 and we process all the improper torsion forces
// involving that atom and two others. Unlike the other bonded terms this one
// is counted every time it comes up.
void DuMMForceFieldSubsystemRep::calcAmberImproperTorsion
   (DuMM::IncludedAtomIndex             a1num,
    const Vector_<Vec3>&                inclAtomStation_G,
    const Vector_<Vec3>&                inclAtomPos_G,
    Real                                amberImproperTorsionScaleFactor,
    Vector_<SpatialVec>&                inclBodyForces_G,
    Real&                               energy) const
{
    const IncludedAtom&           a1    = getIncludedAtom(a1num);
    const DuMMIncludedBodyIndex   b1    = a1.inclBodyIndex;
    const Vec3& a1Station_G = inclAtomStation_G[a1num];
    const Vec3& a1Pos_G     = inclAtomPos_G[a1num];

    // Note that a1 is the *third* atom in the torsion; the other three
    // are stored for each entry in forceImproper14.
    for (int b14=0; b14 < (int)a1.forceImproper14.size(); ++b14) {
        const DuMM::IncludedAtomIndex a2num = a1.forceImproper14[b14][0];
        const DuMM::IncludedAtomIndex a3num = a1.forceImproper14[b14][1];
        const DuMM::IncludedAtomIndex a4num = a1.forceImproper14[b14][2];

        const IncludedAtom& a2 = getIncludedAtom(a2num);
        const IncludedAtom& a3 = getIncludedAtom(a3num);
        const IncludedAtom& a4 = getIncludedAtom(a4num);

        const Vec3& a2Station_G = inclAtomStation_G[a2num];
        const Vec3& a3Station_G = inclAtomStation_G[a3num];
        const Vec3& a4Station_G = inclAtomStation_G[a4num];
        const Vec3& a2Pos_G     = inclAtomPos_G[a2num];
        const Vec3& a3Pos_G     = inclAtomPos_G[a3num];
        const Vec3& a4Pos_G     = inclAtomPos_G[a4num];

        Real angle, e;
        Vec3 f1, f2, f3, f4;
        const BondTorsion& bt = *a1.aImproperTorsion[b14];

        bt.calculateAtomForces
           (a2Pos_G, a3Pos_G, a1Pos_G, a4Pos_G,
            amberImproperTorsionScaleFactor, 0, // no custom improper
            angle, e, f2, f3, f1, f4);

        const DuMMIncludedBodyIndex b2 = a2.inclBodyIndex;
        const DuMMIncludedBodyIndex b3 = a3.inclBodyIndex;
        const DuMMIncludedBodyIndex b4 = a4.inclBodyIndex;

	    assert(!(b2==b1 && b3==b1 && b4==b1)); // shouldn't be on the list if all on 1 body

        energy += e;

        inclBodyForces_G[b1] += SpatialVec( a1Station_G % f1, f1); // 15 flops
        inclBodyForces_G[b2] += SpatialVec( a2Station_G % f2, f2); // 15 flops
        inclBodyForces_G[b3] += SpatialVec( a3Station_G % f3, f3); // 15 flops
        inclBodyForces_G[b4] += SpatialVec( a4Station_G % f4, f4); // 15 flops
    }
}
//.........................CALC AMBER IMPROPER TORSION..........................



//------------------------------------------------------------------------------
//                    CALC BODY SUBSET NONBONDED FORCES
//------------------------------------------------------------------------------
// Helper routine for realizeDynamics when nonbonded forces are being calculated
// on the CPU either single-threaded or in parallel. This is just van der Waals
// and Coulomb forces, not GBSA.
// TODO: there are *no* cutoffs here.
// This is *very* expensive -- code carefully!
//
// Strategy:
//   for a single included body b
//     for each nonbond atom ab on b
//          set scale factors on atoms bonded to atom ab
//          for each body c in range [first,last]
//            for each nonbond atom ac on c
//                 compute vector r=ac-ab and distance d=|r|
//                 compute vdw forces
//                 compute Coulomb forces
//                 add force contribution to atoms
//          reset scale factors on atoms bonded to atom ab
//

void DuMMForceFieldSubsystemRep::calcBodySubsetNonbondedForces
   (DuMMIncludedBodyIndex                   dummBodIx,
    DuMMIncludedBodyIndex                   firstIx,
    DuMMIncludedBodyIndex                   lastIx,
    const Vector_<Vec3>&                    inclAtomPos_G,
    Array_<Real,DuMM::NonbondAtomIndex>&    vdwScale,    // temps: all 1s
    Array_<Real,DuMM::NonbondAtomIndex>&    coulombScale,
    Vector_<Vec3>&                          inclAtomForce_G,
    Real&                                   energy) const
{
    //TRACE("DuMMForceFieldSubsystemRep::calcBodySubsetNonbondedForces BEGIN\n");

    const IncludedBody& inclBod1 = includedBodies[dummBodIx];

    // Run through every nonbond atom that is attached to this included body.
    for (DuMM::NonbondAtomIndex nax1 = inclBod1.beginNonbondAtoms;
         nax1 != inclBod1.endNonbondAtoms; ++nax1)
    {
        DuMM::IncludedAtomIndex iax1 = getIncludedAtomIndexOfNonbondAtom(nax1);
        const IncludedAtom& a1 = getIncludedAtom(iax1);
        const ChargedAtomType& a1type = chargedAtomTypes[a1.chargedAtomTypeIndex];
        const DuMM::AtomClassIndex a1cnum = a1type.atomClassIx;

        //TRACE( (std::string(" NonbondAtomIndex ") + std::to_string(nax1)).c_str() );
        //TRACE( (std::string(" AtomIndex ") + std::to_string(getAtomIndexOfIncludedAtom(iax1))).c_str() );
        //TRACE( (std::string(" IncludedAtomIndex ") + std::to_string(iax1)).c_str() );
        //TRACE( (std::string(" ChargedAtomTypeIndex ") + std::to_string(a1.chargedAtomTypeIndex)).c_str() );
        //TRACE( (std::string(" AtomClassIndex ") + std::to_string(a1cnum)).c_str() );

        const AtomClass&           a1class = atomClasses[a1cnum];
        const Vec3&                a1Pos_G = inclAtomPos_G[iax1];

        const Real q1Fac = coulombGlobalScaleFactor
                                * CoulombFac * a1type.partialCharge;

        // Set scale factors for all closely-bonded atoms to a1 that are
        // involved in nonbond calculations.
        scaleBondedAtoms(a1,vdwScale,coulombScale);

        // We'll update this element repeatedly below so get a reference
        // to it once outside the loop.
        Vec3& afrc1_G = inclAtomForce_G[iax1];

        for (DuMMIncludedBodyIndex dbx2 = firstIx; dbx2 <= lastIx; ++dbx2) {
            assert(dbx2 != dummBodIx);
            const IncludedBody& inclBod2 = includedBodies[dbx2];

            // Run through all the nonbond atoms that are attached to this body.
            for (DuMM::NonbondAtomIndex nax2 = inclBod2.beginNonbondAtoms;
                 nax2 != inclBod2.endNonbondAtoms; ++nax2)
            {
   //             TRACE(" ");
                DuMM::IncludedAtomIndex iax2 =
                    getIncludedAtomIndexOfNonbondAtom(nax2);
                assert(iax2 != iax1);
                const IncludedAtom& a2 = getIncludedAtom(iax2);
                const ChargedAtomType& a2type  = chargedAtomTypes[a2.chargedAtomTypeIndex];
                const DuMM::AtomClassIndex a2cnum  = a2type.atomClassIx;
                //TRACE( (std::string(" NonbondAtomIndex ") + std::to_string(nax2)).c_str() );
                //TRACE( (std::string(" AtomIndex ") + std::to_string(getAtomIndexOfIncludedAtom(iax1))).c_str() );
                //TRACE( (std::string(" IncludedAtomIndex ") + std::to_string(iax2)).c_str() );
                //TRACE( (std::string(" ChargedAtomTypeIndex ") + std::to_string(a2.chargedAtomTypeIndex)).c_str() );
                //TRACE( (std::string(" AtomClassIndex ") + std::to_string(a2cnum)).c_str() );
                const AtomClass& a2class = atomClasses[a2cnum];
                const Vec3&      a2Pos_G = inclAtomPos_G[iax2];

                const Vec3  r  = a2Pos_G - a1Pos_G; // from a1 to a2 (3 flops)
                const Real  d2 = r.normSqr() ;     // 5 flops
                const Real  d = std::sqrt(d2);

                if( nonbondedMethod ==0 || ( nonbondedMethod ==1 && d <= nonbondedCutoff ) ) {

                    //TRACE( (std::string(" r ") + std::to_string(std::sqrt(d2))).c_str() );
                    const Real ood = 1 / d; // approx 40 flops
                    const Real ood2 = ood * ood;        // 1 flop

                    // Coulombic electrostatic force
                    const Real qq = coulombScale[nax2] // 2 flops
                                    * q1Fac * a2type.partialCharge;
                    //TRACE( (std::string(" a1type.partialCharge ") + std::to_string(a1type.partialCharge) + std::string(" a2type.partialCharge ") + std::to_string(a2type.partialCharge)).c_str() );
                    // e = scale*(1/(4*pi*e0)) *  q1*q2/d
                    const Real eCoulomb = qq * ood;     // 1 flop
                    // f = -[scale*(1/(4*pi*e0)) * -q1*q2/d^2] * r
                    // Note: we're missing a factor of 1/d^2 here; see below.
                    const Real fCoulomb = eCoulomb;

                    // van der Waals forces

                    // Get precomputed mixed dmin and emin. Must ask the lower-numbered atom class.
                    Real dij, eij;
                    if (a1cnum <= a2cnum) {
                        dij = a1class.vdwDij[a2cnum - a1cnum];
                        eij = a1class.vdwEij[a2cnum - a1cnum];
                    } else {
                        dij = a2class.vdwDij[a1cnum - a2cnum];
                        eij = a2class.vdwEij[a1cnum - a2cnum];
                    }
                    //TRACE( (std::string(" Dij ") + std::to_string(dij) + std::string(" eij ") + std::to_string(eij)).c_str() );

                    // 5 flops
                    const Real ddij2 = dij * dij * ood2;   // (dmin_ij/d)^2
                    const Real ddij6 = ddij2 * ddij2 * ddij2;
                    const Real ddij12 = ddij6 * ddij6;

                    // 8 flops
                    const Real eijScale = vdwGlobalScaleFactor * vdwScale[nax2] * eij;
                    //TRACE( (std::string(" vdwGlobal ") + std::to_string(vdwGlobalScaleFactor) + std::string(" vdwScale[nax2] ") + std::to_string(vdwScale[nax2])).c_str() );
                    const Real eVdw = eijScale * (ddij12 - 2 * ddij6);
                    //TRACE(" eVdw: ");
                    //TRACE((std::to_string(eVdw)).c_str());
                    // Note: factor of 1/d^2 missing here; see below.
                    const Real fVdw = 12 * eijScale * (ddij12 - ddij6);

                    // Here's where we restore the missing 1/d^2. This
                    // is the force to apply to atom 2; apply equal and
                    // opposite to atom 1.
                    const Vec3 fj = ((fCoulomb + fVdw) * ood2) * r; // 5 flops
                    // kJ (Da-nm^2/ps^2)        // 2 flops
                    energy += (eCoulomb + eVdw);
                    //TRACE(" energy: ");
                    //TRACE((std::to_string(eVdw)).c_str());
                    //TRACE(" ");
                    inclAtomForce_G[iax2] += fj;   // 3 flops
                    afrc1_G -= fj;   // 3 flops
                }
            }
        }

        // This is the end of the outer atom loop. We're done with atom a1.
        unscaleBondedAtoms(a1,vdwScale,coulombScale);
    }
    //TRACE("\nDuMMForceFieldSubsystemRep::calcBodySubsetNonbondedForces END\n");
}
//....................CALC BODY SUBSET NONBONDED FORCES.........................



//------------------------------------------------------------------------------
//                          CALC NONBONDED FORCES
//------------------------------------------------------------------------------
// This is the single-threaded nonbonded force calculation (not including GBSA).
// For each included body, we calculate its interactions with all
// higher-numbered included bodies.
void DuMMForceFieldSubsystemRep::calcNonbondedForces
   (const Vector_<Vec3>&                inclAtomPos_G,
    Vector_<Vec3>&                      inclAtomForce_G,
    Real&                               energy) const
{
    //TRACE("calcNonbondedForces() BEGIN");
    for (DuMMIncludedBodyIndex inclBodyIx(0);
         inclBodyIx < getNumIncludedBodies(); ++inclBodyIx)
    {
        calcBodySubsetNonbondedForces(
            inclBodyIx,
            DuMMIncludedBodyIndex(inclBodyIx + 1),
            DuMMIncludedBodyIndex(getNumIncludedBodies()-1),
            inclAtomPos_G,
            vdwScaleSingleThread,     // these 2 temps are indexed by nonbond
            coulombScaleSingleThread, //   atom index, *not* included atom index
            inclAtomForce_G, energy);
    }
}
//............................CALC NONBONDED FORCES.............................



//------------------------------------------------------------------------------
//                        class NonbondedForceTask
//------------------------------------------------------------------------------
// This class is used by realizeDynamics for calculating nonbonded interactions
// in multiple threads. The strategy here is to let the Parallel2DExecutor
// utility execute disjoint pairs of included bodies in parallel.


// ----------------------------------------------------------------------------
// Gmolmodel - modified for simbody ParallelExecutor update
// ----------------------------------------------------------------------------
class NonbondedForceTask : public SimTK::Parallel2DExecutor::Task {
public:
    NonbondedForceTask
       (const DuMMForceFieldSubsystemRep& dumm,
        const Vector_<Vec3>& inclAtomPos_G,
        Vector_<Vec3>& inclAtomForces_G, Real& energy)
    :   dumm(dumm), inclAtomPos_G(inclAtomPos_G),
        globalAtomForces_G(inclAtomForces_G), globalEnergy(energy)
    {
    }


    void initialize() {
  //      TRACE("NonbondedForceTask::initialize BEGIN\n");
        localEnergy = 0.0;

        // Temps for nonbonded scale factors; initialize to 1
        localVdwScale.resize(getNumNonbondAtoms(), Real(1));
        localCoulombScale.resize(getNumNonbondAtoms(), Real(1));
   //     TRACE("NonbondedForceTask::initialize END\n");
    }

    // At the end of execution, each thread adds its local energy contribution
    // to the global total. See comment above regarding forces.
    void finish() {
    //    TRACE("NonbondedForceTask::finish BEGIN\n");

        globalEnergy += localEnergy;

//        TRACE("NonbondedForceTask::finish ");
//        TRACE(std::to_string(globalEnergy).c_str());
//        TRACE(" ");
//        TRACE(std::to_string(localEnergy).c_str());
//        TRACE(" END\n");
    }

    // This is the standard unit of work for a pair of body indices. No
    // simultaneous task will be using either of these two indices, so we
    // can access global data that is indexed by them without synchronization.
    void execute(int body1, int body2) {
     //   TRACE("NonbondedForceTask::execute ");
     //   TRACE( (std::to_string(body1) + std::string(" ") + std::to_string(body2) ).c_str());
      //  TRACE(" BEGIN\n");

        // // test
        // localEnergy += 1.0;    

        dumm.calcBodySubsetNonbondedForces(
            DuMMIncludedBodyIndex(body1),
            DuMMIncludedBodyIndex(body2),   // i.e, just one body
            DuMMIncludedBodyIndex(body2),
            inclAtomPos_G,
            localVdwScale, localCoulombScale,
            globalAtomForces_G, localEnergy);

//        TRACE("NonbondedForceTask::execute ");
//        TRACE(std::to_string(localEnergy).c_str());
//       TRACE(" END\n");
    }


private:
    int getNumNonbondAtoms() const {return dumm.getNumNonbondAtoms();}

    const DuMMForceFieldSubsystemRep&   dumm;
    const Vector_<Vec3>&                inclAtomPos_G;
    Vector_<Vec3>&                      globalAtomForces_G;
    Real&                               globalEnergy;

    // Thread local temporaries.
    static thread_local Real                                  localEnergy;
    static thread_local Array_<Real, DuMM::NonbondAtomIndex>  localVdwScale;
    static thread_local Array_<Real, DuMM::NonbondAtomIndex>  localCoulombScale;



};

// because static
thread_local Real                                  NonbondedForceTask::localEnergy = 0.0;
thread_local Array_<Real, DuMM::NonbondAtomIndex>  NonbondedForceTask::localVdwScale;
thread_local Array_<Real, DuMM::NonbondAtomIndex>  NonbondedForceTask::localCoulombScale;

//..........................class NonbondedForceTask............................



//------------------------------------------------------------------------------
//                              CALC GBSA FORCES
//------------------------------------------------------------------------------
// Helper function used by realizeSubsystemDynamicsImpl().
// This calculates GBSA implicit solvent forces. This is the most expensive
// calculation -- probably 3X the other nonbonded terms.
void DuMMForceFieldSubsystemRep::calcGBSAForces
   (const Vector_<Vec3>&                inclAtomStation_G,
    const Vector_<Vec3>&                inclAtomPos_G,
    bool                                useParallel,
    Real                                gbsaGlobalScaleFac,
    Vector_<SpatialVec>&                inclBodyForces_G,
    Real&                               energy) const
{
    // 1) Populate array of atom positions for gbsa IN ANGSTROMS
    // Put atomic coordinates relative to ground in gbsaRawCoordinates.
    // Note that we pass a pre-calculated array of pointers to GBSA; it
    // was set in realizeTopology() to point into gbsaRawCoordinates.

    for (DuMM::NonbondAtomIndex nax(0); nax < getNumNonbondAtoms(); ++nax) {
        const DuMM::IncludedAtomIndex iax =
            getIncludedAtomIndexOfNonbondAtom(nax);
        // atomic coordinates with respect to Ground frame
        const Vec3 aPos_G  = inclAtomPos_G[iax]*DuMM::Nm2Ang; // Angstroms

        gbsaRawCoordinates[3 * nax + 0] = aPos_G[0];
        gbsaRawCoordinates[3 * nax + 1] = aPos_G[1];
        gbsaRawCoordinates[3 * nax + 2] = aPos_G[2];
    }


    // compute GBSA forces and energy
    const int returnValue = gbsaCpuObc->computeImplicitSolventForces
       (&gbsaCoordinatePointers.front(), &gbsaAtomicPartialCharges.front(),
        &gbsaAtomicForcePointers.front(), useParallel ? gbsaExecutor : NULL );
    SimTK_ASSERT_ALWAYS(returnValue == 0,
        "GBSA CpuObc::computeImplicitSolventForces() failed.");

    RealOpenMM gbsaEnergy = gbsaCpuObc->getEnergy();

    // 4)  apply GBSA forces to bodies

    // convert force units from kcal/mol-A to to kJ/mol-nm
    const Real KcalA2KJnm = DuMM::Kcal2KJ/DuMM::Ang2Nm;

    for (DuMMIncludedBodyIndex inclBodyIx(0);
         inclBodyIx < getNumIncludedBodies(); ++inclBodyIx)
    {
        const IncludedBody& inclBod = includedBodies[inclBodyIx];

        // Run through all the nonbond atoms that are attached to this body
        // and collect up the GBSA forces from each.
        for (DuMM::NonbondAtomIndex nax = inclBod.beginNonbondAtoms;
                nax != inclBod.endNonbondAtoms; ++nax)
        {
            const DuMM::IncludedAtomIndex iax =
                getIncludedAtomIndexOfNonbondAtom(nax);
            const Vec3& aStation_G = inclAtomStation_G[iax];  // nm

            Vec3 fGbsa(gbsaAtomicForces[3 * nax + 0],      // kcal/mol-A
                       gbsaAtomicForces[3 * nax + 1],
                       gbsaAtomicForces[3 * nax + 2]);

            // convert force units from kcal/mol-A to to kJ/mol-nm
            fGbsa *= KcalA2KJnm;
            fGbsa *= gbsaGlobalScaleFac;

            inclBodyForces_G[inclBodyIx] +=
                SpatialVec( aStation_G % fGbsa, fGbsa );
        }
    }

    // update potential energy from gbsa
    // convert kcal/mol to kJ/mol
    gbsaEnergy *= DuMM::Kcal2KJ;
    gbsaEnergy *= gbsaGlobalScaleFac;
    energy += gbsaEnergy;
}
//..............................CALC GBSA FORCES................................



//------------------------------------------------------------------------------
//                          REALIZE FORCES AND ENERGY
//------------------------------------------------------------------------------
// Here's where we calculate all the forces if they haven't already been done.
// Potential energy is calculated at the same time since that comes for free.
void DuMMForceFieldSubsystemRep::realizeForcesAndEnergy(const State& s) const
{

    if (   isIncludedAtomForceCacheRealized(s)
        && isIncludedBodyForceCacheRealized(s)
        && isEnergyCacheRealized(s))
        return; // nothing to do

    // Get access to the matter subsystem so we can access the bodies.
    // const MultibodySystem&        mbs    = getMultibodySystem();
    //const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // These are the DuMM-local cache entries that we're going to fill in here.
    Vector_<Vec3>&          inclAtomForce_G  = updIncludedAtomForceCache(s);
    Vector_<SpatialVec>&    inclBodyForces_G = updIncludedBodyForceCache(s);
    Real&                   energy           = updEnergyCache(s);

    inclAtomForce_G.resize(getNumIncludedAtoms());
    inclAtomForce_G = Vec3(0);

    inclBodyForces_G.resize(getNumIncludedBodies());
    inclBodyForces_G = SpatialVec(Vec3(0), Vec3(0));

    energy = 0;
	++forceEvaluationCount;

    // Get access to already-calculated position dependent quantities.
    const Vector_<Vec3>& inclAtomStation_G = getIncludedAtomStationsInG(s);
    const Vector_<Vec3>& inclAtomPos_G     = getIncludedAtomPositionsInG(s);


        // BONDED FORCES //

    if ( ! ( usingOpenMM && ! wantOpenMMCalcOnlyNonBonded ) ) {

        const bool doStretch = bondStretchGlobalScaleFactor != 0
                               || customBondStretchGlobalScaleFactor != 0;
        const bool doBend = bondBendGlobalScaleFactor != 0
                            || customBondBendGlobalScaleFactor != 0;
        const bool doTorsion = bondTorsionGlobalScaleFactor != 0
                               || customBondTorsionGlobalScaleFactor != 0;
        const bool doImproper = amberImproperTorsionGlobalScaleFactor != 0;

        for (DuMMIncludedBodyIndex incBodyIx(0);
             incBodyIx < getNumIncludedBodies(); ++incBodyIx) {
            const IncludedBody &inclBody = includedBodies[incBodyIx];
            assert(inclBody.isValid());

            // For each atom on this body, deal with all the bonded forces for which
            // it is atom 1 of a bond. See discussion elsewhere for how we avoid
            // double counting.
            for (DuMMBondStarterIndex bsx = inclBody.beginBondStarterAtoms;
                 bsx != inclBody.endBondStarterAtoms; ++bsx) {
                const DuMM::IncludedAtomIndex atom = bondStarterAtoms[bsx];

                // Bond stretch (1-2)
                if (doStretch)
                    calcBondStretch(atom, inclAtomStation_G, inclAtomPos_G,
                                    bondStretchGlobalScaleFactor, customBondStretchGlobalScaleFactor,
                                    inclBodyForces_G, energy);

                // Bond bend (1-2-3)
                if (doBend)
                    calcBondBend(atom, inclAtomStation_G, inclAtomPos_G,
                                 bondBendGlobalScaleFactor, customBondBendGlobalScaleFactor,
                                 inclBodyForces_G, energy);

                // Bond torsion (1-2-3-4)
                if (doTorsion)
                    calcBondTorsion(atom, inclAtomStation_G, inclAtomPos_G,
                                    bondTorsionGlobalScaleFactor, customBondTorsionGlobalScaleFactor,
                                    inclBodyForces_G, energy);

                // Amber improper torsion   2-1-3
                //                             \4
                if (doImproper)
                    calcAmberImproperTorsion(atom, inclAtomStation_G,
                                             inclAtomPos_G, amberImproperTorsionGlobalScaleFactor,
                                             inclBodyForces_G, energy);
            }
        }
        TRACE_OPENMM (("CALC BONDED with DUMM: Ebonded = " + std::to_string(energy) +  " \n").c_str());
    }




                // NONBONDED FORCES //

   // TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: trying to evaluate nonbonded\n");
    if (getNumNonbondAtoms()) {
   //     TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: getNumNonbondAtoms != 0\n");
        // We'll use GPU acceleration if possible; otherwise parallel computation;
        // otherwise serial calculation here.

        if (usingOpenMM) {
    //        TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: using OpenMM\n");
     //       TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: assert passed\n");

            // Calculate forces and energy.
            // TODO: should calculate energy only when it is asked for.

            TRACE_OPENMM (("BEFORE OpenMM\t" + std::to_string(energy) +  " \n").c_str());
            openMMPlugin.calcOpenMMEnergyAndForces(
                inclAtomStation_G, inclAtomPos_G, true /*forces*/, true /*energy*/,
                inclBodyForces_G, energy);
            TRACE_OPENMM (("AFTER OpenMM\t" + std::to_string(energy) +  " \n").c_str());

            // All done!
            markIncludedAtomForceCacheRealized(s);
            markIncludedBodyForceCacheRealized(s);
            markEnergyCacheRealized(s);
            return;
        }

        // We're not using OpenMM; calculate these terms here as best we can.
        if (usingMultithreaded) {
     //       TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: using multithreading\n");
            // Parallel calculation.
            NonbondedForceTask task
               (*this, inclAtomPos_G, inclAtomForce_G, energy);
     //       TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: about to execute nonbondedExecutor with: ");
      //      TRACE(std::to_string(numThreadsInUse).c_str());
      //      TRACE(" threads\n");
            nonbondedExecutor->execute(task, Parallel2DExecutor::HalfMatrix);

            TRACE_OPENMM (("CALC NONBONDED with DUMM Parallel2DExecutor: Energy = " + std::to_string(energy) +  " \n").c_str());

        } else {
            // Serial calculation in this thread.
      //      TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: begin serial calculation\n");
            if (!(coulombGlobalScaleFactor==0 && vdwGlobalScaleFactor==0)) {
      //          TRACE("DuMMForceFieldSubsystemRep::realizeForcesAndEnergy: scale factor not 0\n");
                calcNonbondedForces(inclAtomPos_G, inclAtomForce_G, energy);

                TRACE_OPENMM (("CALC NONBONDED with DUMM SingleThread: Energy = " + std::to_string(energy) +  " \n").c_str());

            }
        }

        // GBSA - (Generalized Born/solvent accessibility implicit) solvent model
        if (gbsaGlobalScaleFactor != 0) {
            calcGBSAForces(inclAtomStation_G, inclAtomPos_G, usingMultithreaded,
                           gbsaGlobalScaleFactor, inclBodyForces_G, energy);
            TRACE_OPENMM (("CALC GBSA with DUMM : Energy = " + std::to_string(energy) +  " \n").c_str());
        }
    }

    // Compute included body spatial forces from generated atom forces.
    for (DuMMIncludedBodyIndex dbx(0); dbx < getNumIncludedBodies(); ++dbx) {
        const IncludedBody& inclBod = includedBodies[dbx];
        for (DuMM::IncludedAtomIndex iax = inclBod.beginIncludedAtoms;
                iax != inclBod.endIncludedAtoms; ++iax)
        {
            const Vec3& aStation_G = inclAtomStation_G[iax];
            const Vec3& aFrc_G     = inclAtomForce_G[iax];
            inclBodyForces_G[dbx] += SpatialVec(aStation_G % aFrc_G, aFrc_G);
        }
    }

    TRACE_OPENMM (("FINAL ENERGY = " + std::to_string(energy) +  " \n").c_str());

    // Done.
    markIncludedAtomForceCacheRealized(s);
    markIncludedBodyForceCacheRealized(s);
    markEnergyCacheRealized(s);

}
//..........................REALIZE FORCES AND ENERGY...........................



//------------------------------------------------------------------------------
//                              REALIZE DYNAMICS
//------------------------------------------------------------------------------
// Here's where forces from this Subsystem are added into the MultibodySystem's
// global array of forces. Forces and energy will be realized here if they
// haven't already been calculated.
int DuMMForceFieldSubsystemRep::realizeSubsystemDynamicsImpl(const State& s) const
{
    // Get access to the System-global cache entry into which all forces must
    // ultimately be accumulated. Units are: kJ (torque), kJ/nm (force).
    const MultibodySystem& mbs = getMultibodySystem();
    Vector_<SpatialVec>& rigidBodyForces = // indexed by MobilizedBodyIndex
        mbs.updRigidBodyForces(s, Stage::Dynamics);

    // Ensure that our force cache contains valid forces.
    realizeForcesAndEnergy(s);
    const Vector_<SpatialVec>& inclBodyForces = // indexed by IncludedBodyIndex
        getIncludedBodyForceCache(s);

    // Apply the generated forces to the mobilized bodies.
    for (DuMMIncludedBodyIndex i(0); i < getNumIncludedBodies(); ++i) {
        const IncludedBody& inclBod = includedBodies[i];
        rigidBodyForces[inclBod.mobodIx] += inclBodyForces[i];
    }

    return 0;
}
//..............................REALIZE DYNAMICS................................



//------------------------------------------------------------------------------
//                            CALC POTENTIAL ENERGY
//------------------------------------------------------------------------------
// Return the potential energy. This can be done any time after stage Position,
// however we have to make sure it has been realized first.
Real DuMMForceFieldSubsystemRep::calcPotentialEnergy(const State& state) const {
    // Currently there is no way to compute only the energy, although it
    // would be somewhat cheaper if forces aren't needed.


    realizeForcesAndEnergy(state);

    //CalcFullPotEnergyIncludingRigidBodiesRep(state);
    return getEnergyCache(state);
}
//............................CALC POTENTIAL ENERGY.............................



//------------------------------------------------------------------------------
//                             SCALE BONDED ATOMS
//------------------------------------------------------------------------------
// Prior to performing nonbonded calculations involving a particular nonbond
// atom, we set scale factors on certain closely-bonded atoms that have been
// precalculated as needing to be scaled.
void DuMMForceFieldSubsystemRep::scaleBondedAtoms
   (const IncludedAtom& a,  // this is a nonbond atom
    Array_<Real,DuMM::NonbondAtomIndex>& vdwScale,
    Array_<Real,DuMM::NonbondAtomIndex>& coulombScale) const
{
    for (int i=0; i < (int)a.scale12.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale12[i];
        vdwScale[ix]=vdwScale12; coulombScale[ix]=coulombScale12;
    }
    for (int i=0; i < (int)a.scale13.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale13[i];
        vdwScale[ix]=vdwScale13; coulombScale[ix]=coulombScale13;
    }
    for (int i=0; i < (int)a.scale14.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale14[i];
        vdwScale[ix]=vdwScale14; coulombScale[ix]=coulombScale14;
    }
    for (int i=0; i < (int)a.scale15.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale15[i];
        vdwScale[ix]=vdwScale15; coulombScale[ix]=coulombScale15;
    }
}
//.............................SCALE BONDED ATOMS...............................



//------------------------------------------------------------------------------
//                            UNSCALE BONDED ATOMS
//------------------------------------------------------------------------------
void DuMMForceFieldSubsystemRep::unscaleBondedAtoms
   (const IncludedAtom& a,
    Array_<Real,DuMM::NonbondAtomIndex>& vdwScale,
    Array_<Real,DuMM::NonbondAtomIndex>& coulombScale) const
{
    for (int i=0; i < (int)a.scale12.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale12[i];
        vdwScale[ix]=coulombScale[ix]=1;
    }
    for (int i=0; i < (int)a.scale13.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale13[i];
        vdwScale[ix]=coulombScale[ix]=1;
    }
    for (int i=0; i < (int)a.scale14.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale14[i];
        vdwScale[ix]=coulombScale[ix]=1;
    }
    for (int i=0; i < (int)a.scale15.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale15[i];
        vdwScale[ix]=coulombScale[ix]=1;
    }
}
//............................UNSCALE BONDED ATOMS..............................



// Vdw combining functions
// -----------------------
// There are several in common use. The most common
// one, Lorentz-Berthelot is also the worst one!
// The pragmatically best seems to be the Waldman-Hagler rule, which
// we will use by default. In between is the Halgren-HHG
// rule. Another good rule is Tang-Toennies but it requires
// additional empirical data (the "sixth dispersion coefficient"
// C6) which we do not have available. An alternative to
// Tang-Toennies is Kong, which uses the Tang-Toennies radius
// formula, but Waldman-Hagler's well depth formula (and Kong
// came considerably before either of them).
//
// The Lennard-Jones 12-6 potential is specified as follows:
// Each atom type i has two parameters ri and ei, resp. the
// van der Waals radius and energy well depth. The radii are
// defined so that if two atoms of type i are separated by
// a distance dmin=2*ri, then the van der Waals energy is -ei.
// For a pair of atoms of types i and j we define an effective
// separation dmin_ij and well depth e_ij. Then if the vector
// from atom i to atom j is v, and d=|v| we have
//
//    Evdw(d) = e_ij * ( (dmin_ij/d)^12 - 2*(dmin_ij/d)^6 )
//
//    Fvdw_j(d) = -grad_j(Evdw)
//              = 12 e_ij * ( (dmin_ij/d)^12 - (dmin_ij/d)^6 ) * v/d^2
//    Fvdw_i(d) = -Fvdw_j(d)
//
// Some cautions: it is common among force fields to specify
// the vdw size (1) either by radius or diameter, and (2) by
// minimum energy or zero crossing. In the latter case the
// symbol "sigma" is used instead of "r", with r=2^(1/6) * sigma
// (that is, sigma is smaller than r). We will be using the
// "radius at minimum energy" convention; note that that has to
// be doubled to produce the dmin used in the LJ formula.



static inline Real arithmeticMean(Real a, Real b) {
    return 0.5*(a+b);
}
static inline Real geometricMean(Real a, Real b) {
    return std::sqrt(a*b);
}
static inline Real harmonicMean(Real a, Real b) {
    return (2*a*b) / (a+b);
}


// cubicMean = (a^3+b^3)/(a^2+b^2)
static inline Real cubicMean(Real a, Real b) {
    return (a*a*a+b*b*b)/(a*a+b*b);
}

// Harmonic mean of harmonic & geometric means
// hhgMean = 4ab/(sqrt(a)+sqrt(b))^2
static inline Real hhgMean(Real a, Real b) {
    return harmonicMean(harmonicMean(a,b), geometricMean(a,b));
}

// Used in AMBER, CHARMM, and MM2/3 (but MMs don't use LJ)
static inline void vdwCombineLorentzBerthelot(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = arithmeticMean(ri,rj);
    e = geometricMean(ei,ej);
}

// Used in OPLS, DANG
static inline void vdwCombineJorgensen(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = geometricMean(ri,rj);
    e = geometricMean(ei,ej);
}

// Used in MMFF, AMOEBA (but with Buffered 14-7 rather than LJ)
static inline void vdwCombineHalgrenHHG(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    r = cubicMean(ri,rj);
    e = hhgMean(ei,ej);
}

static const Real oo6  = Real(1/6.L);
static const Real oo13 = Real(1/13.L);

// This doesn't seem to be used by anyone but it should be!
// Ref: Waldman, M. & Hagler, A.T. New combining rules for
// rare gas van der Waals parameters.
// J. Comput. Chem. 14(9):1077 (1993).
static inline void vdwCombineWaldmanHagler(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    const Real ri3 = ri*ri*ri, ri6 = ri3*ri3;
    const Real rj3 = rj*rj*rj, rj6 = rj3*rj3;
    const Real er6 = geometricMean(ei*ri6, ej*rj6);
    const Real r6  = arithmeticMean(ri6, rj6);

    r = std::pow(r6, oo6);
    e = er6 / r6;
}

// This is a possible alternative to Waldman-Hagler. It uses
// the same well depth combination term as WH, but with a different
// radius combination term which is the same as Tang-Toennies.
// Ref: Kong, C.L. Combining rules for intermolecular potential
// parameters. II. Rules for the Lennard-Jones (12-6) potential
// and the Morse potential. J. Chem. Phys. 59(5):2464 (1973).
// Comparison with WH: Delhommelle, J. & Millie, P. Inadequacy of
// the Lorentz-Berthelot combining rules for accurate predictions
// of equilibrium properties by molecular simulation. Molecular
// Physics 99(8):619 (2001).

static inline void vdwCombineKong(
    Real ri, Real rj, Real ei, Real ej,
    Real& r, Real& e)
{
    const Real ri3 = ri*ri*ri, ri6 = ri3*ri3, ri12 = ri6*ri6;
    const Real rj3 = rj*rj*rj, rj6 = rj3*rj3, rj12 = rj6*rj6;
    const Real er6 = geometricMean(ei*ri6, ej*rj6);

    // calculate (ei*ri^12)^(1/13), etc.
    const Real eri12_13 = std::pow(ei*ri12, oo13);
    const Real erj12_13 = std::pow(ej*rj12, oo13);
    const Real er12_13  = arithmeticMean(eri12_13, erj12_13);
    const Real r6 =  std::pow(er12_13, 13) / er6;

    r = std::pow(r6, oo6);
    e = er6 / r6;
}


// Radii and returned diameter are given in nm, energies in kJ/mol.
void DuMMForceFieldSubsystemRep::applyMixingRule(Real ri, Real rj, Real ei, Real ej, Real& dmin, Real& emin) const
{
    Real rmin = 0;

    switch(vdwMixingRule) {
    case DuMMForceFieldSubsystem::WaldmanHagler:
        vdwCombineWaldmanHagler(ri,rj,ei,ej,rmin,emin);     break;
    case DuMMForceFieldSubsystem::HalgrenHHG:
        vdwCombineHalgrenHHG(ri,rj,ei,ej,rmin,emin);        break;
    case DuMMForceFieldSubsystem::Jorgensen:
        vdwCombineJorgensen(ri,rj,ei,ej,rmin,emin);         break;
    case DuMMForceFieldSubsystem::LorentzBerthelot:
        vdwCombineLorentzBerthelot(ri,rj,ei,ej,rmin,emin);  break;
    case DuMMForceFieldSubsystem::Kong:
        vdwCombineKong(ri,rj,ei,ej,rmin,emin);              break;
    default: assert(!"unknown vdw mixing rule");
    };

    dmin = 2*rmin;
}


void DuMMForceFieldSubsystemRep::dump() const
{
    printf("===============================================================\n");
    printf("Dump of DuMMForceFieldSubsystem:\n");
    printf("  # bodies in system=%d, with atoms=%d, included in forces=%d\n",
        getMultibodySystem().getMatterSubsystem().getNumBodies(),
        duMMSubsetOfBodies.size(), getNumIncludedBodies());
    printf("  NClusters=%d NAtoms=%d NAtomClasses=%d NChargedAtomTypes=%d NBonds=%d\n",
        clusters.size(), atoms.size(),
        atomClasses.size(), chargedAtomTypes.size(), bonds.size());
    printf("  # included atoms=%d, # nonbond atoms=%d, # bondstarter atoms=%d\n",
        getNumIncludedAtoms(), getNumNonbondAtoms(), getNumBondStarterAtoms());

    printf("\n================= DuMM BODIES ================\n");
    for (DuMMBodyIndex i(0); i < duMMSubsetOfBodies.size(); ++i) {
        printf("  DuMMBody %d:\n", (int)i);
        duMMSubsetOfBodies[i].dump();
    }
    printf("\n================= INCLUDED BODIES ================\n");
    for (DuMMIncludedBodyIndex i(0); i < getNumIncludedBodies(); ++i) {
        printf("  IncludedBody %d:\n", (int)i);
        getIncludedBody(i).dump();
    }
    printf("\n================= ALL ATOMS ================\n");
    for (DuMM::AtomIndex i(0); i < (int)atoms.size(); ++i) {
        printf("  Atom %d: ", (int)i);
        atoms[i].dump();
    }
    printf("\n================= INCLUDED ATOMS ================\n");
    for (DuMM::IncludedAtomIndex i(0); i < getNumIncludedAtoms(); ++i) {
        printf("  Incl Atom %d: ", (int)i);
        getIncludedAtom(i).dump();
    }
    printf("\n================= NONBOND ATOMS ================\n");
    printf("  ");
    for (DuMM::NonbondAtomIndex i(0); i < getNumNonbondAtoms(); ++i) {
        printf("%d ", (int)getIncludedAtomIndexOfNonbondAtom(i));
        if ((i+1)%21 == 0) printf("\n  ");
    }
    printf("\n");
    printf("\n================= BOND STARTER ATOMS ================\n");
    printf("  ");
    for (DuMMBondStarterIndex i(0); i < getNumBondStarterAtoms(); ++i) {
        printf("%d ", (int)getIncludedAtomIndexOfBondStarterAtom(i));
        if ((i+1)%21 == 0) printf("\n  ");
    }
    printf("\n");
    printf("\n================= CLUSTERS ================\n");
    for (DuMM::ClusterIndex i(0); i < clusters.size(); ++i) {
        printf("  Cluster %d:\n", (int)i);
        clusters[i].dump();
    }
    printf("\n================= ATOM CLASSES ================\n");
    for (DuMM::AtomClassIndex i(0); i < (int)atomClasses.size(); ++i) {
        if (!atomClasses[i].isValid()) continue;
        printf("  AtomClass %d:\n", (int)i);
        atomClasses[i].dump();
    }
    printf("\n================= CHARGED ATOM TYPES ================\n");
    for (DuMM::ChargedAtomTypeIndex i(0); i < (int)chargedAtomTypes.size(); ++i) {
        if (!chargedAtomTypes[i].isValid()) continue;
        printf("  ChargedAtomType %d:\n", (int)i);
        chargedAtomTypes[i].dump();
    }
    printf("===============================================================\n");

}

std::ostream& DuMMForceFieldSubsystemRep::generateBiotypeChargedAtomTypeSelfCode(std::ostream& os, BiotypeIndex biotypeIx) const
{
    assert(chargedAtomTypesByBiotype.find(biotypeIx) != chargedAtomTypesByBiotype.end());
    DuMM::ChargedAtomTypeIndex typeId = chargedAtomTypesByBiotype.find(biotypeIx)->second;
    os << "    dumm.setBiotypeChargedAtomType(";
    os << "DuMM::ChargedAtomTypeIndex(" << typeId << ")";

    const Biotype& biotype = Biotype::get(biotypeIx);
    String biotypeAtomName = Biotype::get(biotypeIx).getAtomName();

    os << ", Biotype::get(";
    os << "\"" << biotype.getResidueName() << "\"";
    os << ", \"" << biotype.getAtomName() << "\"";
    switch (biotype.getOrdinality()) {
        case Ordinality::Any:
            os << ", Ordinality::Any";
            break;
        case Ordinality::Initial:
            os << ", Ordinality::Initial";
            break;
        case Ordinality::Final:
            os << ", Ordinality::Final";
            break;
        default:
            assert(false);
            break;
    }
    os << ").getIndex()";

    os << ");";
    os << std::endl;

    return os;
}

void DuMMForceFieldSubsystemRep::setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex, BiotypeIndex biotypeIx)
{
    assert(biotypeIx.isValid());
    assert(chargedAtomTypeIndex.isValid());

    // OK if value is already populated
    if (chargedAtomTypesByBiotype.find(biotypeIx) != chargedAtomTypesByBiotype.end())
    {
        [[maybe_unused]] const auto oldId = chargedAtomTypesByBiotype.find(biotypeIx)->second;
        assert(oldId == chargedAtomTypeIndex);
    }

    chargedAtomTypesByBiotype[biotypeIx] = chargedAtomTypeIndex;
    assert(chargedAtomTypesByBiotype.find(biotypeIx) != chargedAtomTypesByBiotype.end());
}

DuMM::ChargedAtomTypeIndex DuMMForceFieldSubsystemRep::getBiotypeChargedAtomType(BiotypeIndex biotypeIx) const {
    assert(biotypeIx.isValid());

//std::cout << "getBiotypeChargedAtomType  biotypeIx= " <<  biotypeIx << "  atom=" << Biotype::get(biotypeIx).getAtomName()
//<< "  residue=" << Biotype::get(biotypeIx).getResidueName() << std::endl;

  //  assert (chargedAtomTypesByBiotype.find(biotypeIx) != chargedAtomTypesByBiotype.end());
    return chargedAtomTypesByBiotype.find(biotypeIx)->second;
}

    ///////////////
    // BOND BEND //
    ///////////////

// Given a central atom location c bonded to atoms at r and s,
// calculate the angle between them, the potential energy,
// and forces on each of the three atoms.
//
// Singular cases:
//   (1) |r| is 0
//   (2) |s| is 0
//   (3) r and s are aligned or anti-aligned
// In the first two cases we can't generate a torque at all because
// two of the atoms are coincident. We will do nothing here since
// the corresponding bond stretch term will eventually separate these
// atoms. In case (3) we can't find a unique direction about which to apply
// the torque. In that case we'll just make up a direction perpendicular
// to both vectors (r and s) and use it; that is the right thing to
// do if we're just doing some kind of minimization since it will push
// the zero bend angle to something nonzero.
void BondBend::calculateAtomForces
   (const Vec3& cG, const Vec3& rG, const Vec3& sG,
    const Real& builtinScale, const Real& customScale,
    Real& theta, Real& pe, Vec3& cf, Vec3& rf, Vec3& sf) const
{
    const Vec3 r = rG - cG; //               3 flops
    const Vec3 s = sG - cG; //               3 flops
    const Real rr = ~r*r, ss = ~s*s;    // |r|^2, |s|^2 ( 10 flops)

    // Check for singular cases (1) and (2).
    if (rr==0 || ss==0) {
        theta = 0; pe = 0;
        cf = rf = sf = Vec3(0); // no forces
        return;
    }

    // If we get here we know that all three atoms are
    // separated, at least a little.

    const Real rs = ~r * s; // r dot s      (5 flops)
    const Vec3 rxs = r % s; // r cross s    (9 flops)
    const Real rxslen = rxs.norm(); //      (~35 flops)
    theta = std::atan2(rxslen, rs); //       ~50 flops

    Real torque;
    if (hasBuiltinTerm()) {
        const Real bend = theta - theta0;       //   1 flop
        const Real skb  = builtinScale*k*bend;  //   2 flops
        pe     =  skb*bend; // NOTE: no factor of 1/2 (1 flop)
        torque = -2*skb;                        // 1 flop
    } else
        pe = torque = 0;

    for (int i=0; i < (int)customTerms.size(); ++i) {
        const DuMM::CustomBondBend& term = *customTerms[i];
        pe     += customScale*term.calcEnergy(theta);
        torque += customScale*term.calcTorque(theta); // expecting -dE/dtheta
    }

    // p is unit vector perpendicular to r and s. Note handling of
    // singular case (3) here.
    const UnitVec3 p = (rxslen != 0 ? UnitVec3(rxs/rxslen,true)  // ~11 flops
                                    : UnitVec3(r).perp());

    // Already checked above that rr and ss are non zero.
    rf = (torque/rr)*(r % p);          // ~20 flops
    sf = (torque/ss)*(p % s);          // ~20 flops
    cf = -(rf+sf); // makes the net force zero (6 flops)
}

    //////////////////
    // BOND TORSION //
    //////////////////

// Given atom locations r-x-y-s in the ground frame, calculate the
// torsion angle, energy and a force on each atom so that the desired
// pure torque is produced.
// This code is modeled in part after Tinker's torsion code in
// etors1.f because I couldn't figure out how to do it myself
// (sherm 060905). Thanks, Jay!
void BondTorsion::calculateAtomForces
   (const Vec3& rG, const Vec3& xG, const Vec3& yG, const Vec3& sG,
    const Real& builtinScale, const Real& customScale,
    Real& theta, Real& pe,
    Vec3& rf, Vec3& xf, Vec3& yf, Vec3& sf) const
{
    // All vectors point along the r->x->y->s direction
    const Vec3 r  = xG - rG; //               3 flops
    const Vec3 s  = sG - yG; //               3 flops
    const Vec3 xy = yG - xG; //               3 flops

    // Create a unit vector v along the axis, using increasingly
    // desperate measures in case of overlapping atoms. If we
    // don't have a real axis (i.e., atoms x and y overlap)
    // we'll signal that with oov==0 (see below). We don't care
    // much what happens in that case, but we hope to do something
    // remotely plausible so a stuck minimization will have some
    // hope of getting unstuck.

    const Real vv = ~xy*xy;                     //   5 flops
    const Real oov = (vv==0 ? Real(0)
                            : 1/std::sqrt(vv)); // ~40 flops
    const UnitVec3 v =
        (oov != 0 ? UnitVec3(xy*oov,true)       //   4 flops
                   : ((r%s).norm() != 0 ? UnitVec3(r % s)
                                        : UnitVec3(r).perp()));

    // Calculate plane normals. Axis vector v serves as the "x"
    // axis of both planes. Vectors r (r->x) and s (y->s) are in
    // the plane in a vaguely "y axis" way, so t=rXv is the "z" axis
    // (plane normal) for the first plane and u=vXs is the plane normal
    // for the second. When those normals are aligned theta is 0.
    const Vec3 t = r % v, u = v % s; // 18 flops

    // If either r or s are aligned with the axis, we can't generate
    // a torque so we're done.
    const Real tt = ~t*t, uu = ~u*u; // 10 flops
    if (tt == 0 || uu == 0) {
        theta=0; pe=0; rf=xf=yf=sf=Vec3(0);
        return;
    }

    const Vec3 txu = t % u;                 //   9 flops
    const Real ootu = 1/std::sqrt(tt*uu);   // ~40 flops
    const Real cth = (~t*u)*ootu;           //   6 flops
    const Real sth = (~v*txu)*ootu;         //   6 flops
    theta = std::atan2(sth,cth);            // ~50 flops

    Real torque = 0;
    pe = 0;
    for (int i=0; i < (int)terms.size(); ++i) {
        pe     += terms[i].energy(theta);
        torque += terms[i].torque(theta);
    }
    pe     *= builtinScale;
    torque *= builtinScale;
    for (int i=0; i < (int)customTerms.size(); ++i) {
        const DuMM::CustomBondTorsion& term = *customTerms[i];
        pe     += customScale*term.calcEnergy(theta);
        torque += customScale*term.calcTorque(theta); // expecting -dE/dtheta
    }

    const Vec3 ry = yG-rG;    // from r->y        3 flops
    const Vec3 xs = sG-xG;    // from x->s        3 flops
    const Vec3 dedt =  (torque/tt)*(t % v);  // ~20 flops
    const Vec3 dedu = -(torque/uu)*(u % v);  // ~21 flops

    rf = dedt % v; // 9 flops
    sf = dedu % v; // 9 flops
    if (oov==0) {
        xf = -rf;   // No axis; this is just desperation.
        yf = -sf;   // At least it keeps the forces summing to 0.
    } else {
        xf = ((ry % dedt) + (dedu % s))*oov;
        yf = ((dedt % r) + (xs % dedu))*oov;
    }
}


    //////////
    // ATOM //
    //////////


void IncludedAtom::dump() const {
    printf(" includedAtomIx=%d includedBodyIx=%d (atomIx=%d)\n",
            (int)inclAtomIndex, (int)inclBodyIndex, (int)atomIndex);
    printf(" chargedAtomType=%d\n", (int)chargedAtomTypeIndex);

    printf("\n          force 1-2 (IncludedAtomIndex):");
    for (int i=0; i < (int)force12.size(); ++i)
        printf(" %d", (int)force12[i]);
    printf("\n          force 1-3 (IncludedAtomIndex):");
    for (int i=0; i < (int)force13.size(); ++i)
        printf(" %d-%d", (int)force13[i][0], (int)force13[i][1]);
    printf("\n          force 1-4 (IncludedAtomIndex):");
    for (int i=0; i < (int)force14.size(); ++i)
        printf(" %d-%d-%d", (int)force14[i][0], (int)force14[i][1],
                            (int)force14[i][2]);
    printf("\n          force 1-5 (IncludedAtomIndex):");
    for (int i=0; i < (int)force15.size(); ++i)
        printf(" %d-%d-%d-%d", (int)force15[i][0], (int)force15[i][1],
                               (int)force15[i][2], (int)force15[i][3]);
    printf("\n          forceImproper 1-4 (IncludedAtomIndex):");
    for (int i=0; i < (int)forceImproper14.size(); ++i)
        printf(" %d- %d-x-%d", (int)forceImproper14[i][0],
                               (int)forceImproper14[i][1],
                               (int)forceImproper14[i][2]);

    printf("\n          scale 1-2 (NonbondAtomIndex):");
    for (int i=0; i < (int)scale12.size(); ++i)
        printf(" %d", (int)scale12[i]);
    printf("\n          scale 1-3 (NonbondAtomIndex):");
    for (int i=0; i < (int)scale13.size(); ++i)
        printf(" %d", (int)scale13[i]);
    printf("\n          scale 1-4 (NonbondAtomIndex):");
    for (int i=0; i < (int)scale14.size(); ++i)
        printf(" %d", (int)scale14[i]);
    printf("\n          scale 1-5 (NonbondAtomIndex):");
    for (int i=0; i < (int)scale15.size(); ++i)
        printf(" %d", (int)scale15[i]);

    printf("\n");

    printf("    1-2 stretch:");
    for (int i=0; i < (int)stretch.size(); ++i)
        printf(" (%g,%g)", stretch[i]->k, stretch[i]->d0);
    printf("\n    1-3 bend:");
    for (int i=0; i < (int)bend.size(); ++i)
        printf(" (%g,%g)", bend[i]->k, bend[i]->theta0);
    printf("\n    1-4 torsion:\n");
    for (int i=0; i < (int)torsion.size(); ++i) {
        const BondTorsion& bt = *torsion[i];
        printf("     ");
        for (int j=0; j<(int)bt.terms.size(); ++j) {
            const TorsionTerm& tt = bt.terms[j];
            printf(" (%d:%g,%g)", tt.periodicity,
                                  tt.amplitude, tt.theta0);
        }
        printf("\n");
    }
    printf("\n    Amber improper torsion:\n");
    for (int i=0; i < (int)aImproperTorsion.size(); ++i) {
        const BondTorsion& bt = *aImproperTorsion[i];
        for (int j=0; j<(int)bt.terms.size(); ++j) {
            const TorsionTerm& tt = bt.terms[j];
            printf(" (%d:%g,%g)", tt.periodicity, tt.amplitude, tt.theta0);
        }
        printf("\n");
    }
    printf("\n");
}

    /////////////
    // CLUSTER //
    /////////////


void Cluster::attachToBody(MobilizedBodyIndex bnum, const Transform& X_BR, DuMMForceFieldSubsystemRep& mm) {
    assert(!isAttachedToBody());
    mobodIx = bnum;
    placement_B = X_BR;

    // Tell all the atoms directly contained in this cluster that they are
    // now attached to the body also. This will fail if any of the atoms are
    // alread attached -- no polygamy.
    AtomPlacementSet::const_iterator ap = directAtomPlacements.begin();
    while (ap != directAtomPlacements.end()) {
        DuMMAtom& a = mm.updAtom(ap->atomIndex);
        a.attachToBody(bnum, X_BR*ap->station);
        ++ap;
    }

    // Now do the same for our contained groups, who will in turn notify their
    // own atoms and subgroups.
    ClusterPlacementSet::const_iterator cp = directClusterPlacements.begin();
    while (cp != directClusterPlacements.end()) {
        Cluster& c = mm.updCluster(cp->clusterIndex);
        c.attachToBody(bnum, X_BR*cp->placement, mm);
        ++cp;
    }
}

// Return true if this cluster contains (directly or indirectly) any atom which has already
// been attached to a body. If so return one of the attached atoms and its body, which can
// be helpful in error messages.
bool Cluster::containsAnyAtomsAttachedToABody(DuMM::AtomIndex& atomIndex, MobilizedBodyIndex& bodyIx,
                                              const DuMMForceFieldSubsystemRep& mm) const
{
    const AtomPlacementSet& myAtoms   = getAllContainedAtoms();
    AtomPlacementSet::const_iterator ap = myAtoms.begin();
    while (ap != myAtoms.end()) {
        const DuMMAtom& a = mm.getAtom(ap->atomIndex);
        if (a.isAttachedToBody()) {
            atomIndex = ap->atomIndex;
            bodyIx = a.getMobodIndex();
            return true;
        }
        ++ap;
    }
    atomIndex = DuMM::InvalidAtomIndex;
    bodyIx = InvalidMobilizedBodyIndex;
    return false;
}

// Place an atom in this cluster. To be valid, the atom must not
// already be
//   (a) in any of the trees of which this group is apart, or
//   (b) attached to a body.
// TODO: (c) at the moment we don't allow placing an atom in a group unless
//           that group is a top-level group (i.e., it has no parents).
// If this group is already attached to a body, then we will update
// the atom entry to note that it is now attached to the body also.
// EU BEGIN COMMENT
void Cluster::placeAtom(DuMM::AtomIndex atomIndex, const Vec3& station,
                        DuMMForceFieldSubsystemRep& mm)
// EU BEGIN
//void Cluster::placeAtom(DuMM::AtomIndex atomIndex, Vec3 station,
//                        DuMMForceFieldSubsystemRep& mm)
// EU END
{
    assert(isTopLevelCluster()); // TODO
    assert(!mm.getAtom(atomIndex).isAttachedToBody());
    assert(!containsAtom(atomIndex));

    std::pair<AtomPlacementSet::iterator, bool> ret;
    ret = directAtomPlacements.insert(AtomPlacement(atomIndex,station));
    assert(ret.second); // must not have been there already

    ret = allAtomPlacements.insert(AtomPlacement(atomIndex,station));
    assert(ret.second); // must not have been there already

    if (isAttachedToBody())
        mm.updAtom(atomIndex).attachToBody(mobodIx, placement_B*station);
}

// Place a child cluster in this parent cluster. To be valid, the child
// must not
//   (a) already be contained in the parent group or one of the parent's subgroups, or
//   (b) contain any atoms which are already present in the parent or any
//       of the parent's subgroups, or
//   (c) already be attached to a body.
// TODO: (d) at the moment we don't allow adding a child group unless
//           the parent (this) group is a top-level group (i.e., it has no parents).
// If the parent is already attached to a body, then we will update
// the child to note that it is now attached to the body also (and it
// will update its contained atoms).
void Cluster::placeCluster(DuMM::ClusterIndex childClusterIndex,
                           const Transform& placement,
                           DuMMForceFieldSubsystemRep& mm)
{   assert(isTopLevelCluster()); // TODO

    Cluster& child = mm.updCluster(childClusterIndex);
    assert(!child.isAttachedToBody());
    assert(!containsCluster(childClusterIndex));

    // Make sure the new child cluster doesn't contain any atoms which are already in
    // any of the trees to which the parent cluster (this) is associated.
    // TODO: for now we need only look at the parent since we know it is top level.
    const AtomPlacementSet& childsAtoms  = child.getAllContainedAtoms();
    AtomPlacementSet&       parentsAtoms = updAllContainedAtoms();

    // Make sure none of the child's atoms are already in the parent.
    AtomPlacementSet::const_iterator ap = childsAtoms.begin();
    while (ap != childsAtoms.end()) {
        [[maybe_unused]] auto ret = parentsAtoms.insert(AtomPlacement(ap->atomIndex, placement*ap->station));
        assert(ret.second); // mustn't have been there already
        ++ap;
    }

    const ClusterPlacementSet& childsClusters  = child.getAllContainedClusters();
    ClusterPlacementSet&       parentsClusters = updAllContainedClusters();

    // Make sure none of the child's atoms are already in the parent.
    ClusterPlacementSet::const_iterator cp = childsClusters.begin();
    while (cp != childsClusters.end()) {
        [[maybe_unused]] auto ret = parentsClusters.insert(ClusterPlacement(cp->clusterIndex, placement*cp->placement));
        assert(ret.second); // mustn't have been there already
        ++cp;
    }

    noteNewChildCluster(childClusterIndex, placement);
    child.noteNewParentCluster(clusterIndex, placement);

    if (isAttachedToBody())
        child.attachToBody(mobodIx, placement_B*placement, mm);

    //TODO: check for loops
}



// Calculate the composite mass properties for this cluster, transformed into
// the indicated frame.
MassProperties Cluster::calcMassProperties
   (const Transform& tr, const DuMMForceFieldSubsystemRep& mm) const
{
    Real    mass = 0;
    Vec3    com(0);
    Inertia inertia(0);
	Inertia inertia_Spheres(0);

    // Calculate the mass properties in the local frame and transform last.
    AtomPlacementSet::const_iterator aap = allAtomPlacements.begin();
    while (aap != allAtomPlacements.end()) {
        const Real ma = mm.getElement(mm.getAtomElementNum(aap->atomIndex))
                                                                .getMass();
        // Get the mass of the nucleus
        SimTK::Real femto2nano = 0.000001;                                                                
        Real ra = ma;
        ra = std::pow(ra, 1.0 / 3.0);
        ra *= 1.2 * 0.01;

        mass    += ma;
        com     += ma*aap->station;
        Inertia pointMass = Inertia(aap->station, ma);

        // Accumulate point masses inertia
        inertia += pointMass;

        // Accumulate spherical inertia
        Inertia sphericalInertia = UnitInertia::sphere(ra);
        sphericalInertia *= ma;
		sphericalInertia += pointMass;
        inertia_Spheres += sphericalInertia;

        ++aap;
    }
    com /= mass;

    //std::cout << "[INERTIA Cluster::calcMassProperties] " << inertia <<" " << inertia_Spheres << std::endl;

    //return MassProperties(mass,com,inertia).calcTransformedMassProps(tr);
    return MassProperties(mass,com,inertia_Spheres).calcTransformedMassProps(tr);
}

    ///////////////
    // DUMM BODY //
    ///////////////

void DuMMBody::realizeTopologicalCache(const DuMMForceFieldSubsystemRep& mm) {
    // nothing
}





//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//		        GMOLMODEL - EXTRA FUNCTIONALITIES
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// int DuMMForceFieldSubsystemRep::realizeSubsystemTopologyImpl(State& s)
// was hijacked by adding the needed topology lists in order to calculate full 
// potential energy at any given time --> see up


//------------------------------------------------------------------------------
//                        class NonbondedEnergyTask
//------------------------------------------------------------------------------
// This class is used by CalcFullPotEnergyIncludingRigidBodies for calculating
// nonbonded energy separated in vdw and coulomb term in multiple threads, using
// Parallel2DExecutor.

class NonbondedFullEnergyTask : public SimTK::Parallel2DExecutor::Task {
public:
    NonbondedFullEnergyTask
            (const DuMMForceFieldSubsystemRep& dumm,
             const Vector_<Vec3>& AllAtomPos_G,
             Real& eVdW, Real& eCoulomb)
            :   dumm(dumm), AllAtomPos_G(AllAtomPos_G),
                eCoulomb(eCoulomb), eVdW(eVdW)
    {
    }


    void initialize() {

        localVdwEnergyAll = 0.0;
        localCoulombEnergyAll = 0.0;

        localVdwScaleAll.resize( getNumAllNonbondAtoms(), Real(1) );
        localCoulombScaleAll.resize( getNumAllNonbondAtoms(), Real(1) );

    }
    void finish() {
        eVdW += localVdwEnergyAll;
        eCoulomb += localCoulombEnergyAll;
    }

    void execute(int body1, int body2) {

        dumm.CalcFullPotEnergyNonbonded(
                DuMMIncludedBodyIndex(body1),
                DuMMIncludedBodyIndex(body2),   // i.e, just one body
                DuMMIncludedBodyIndex(body2),
                AllAtomPos_G,
                localVdwScaleAll, localCoulombScaleAll,
                localVdwEnergyAll, localCoulombEnergyAll);
    }


private:
    int getNumAllNonbondAtoms() const {return dumm.getNumAllNonbondAtoms();}

    const DuMMForceFieldSubsystemRep&   dumm;
    const Vector_<Vec3>&                AllAtomPos_G;
    Real&                               eCoulomb;
    Real&                               eVdW;

    // Thread local temporaries.
    static thread_local Real                                  localVdwEnergyAll;
    static thread_local Real                                  localCoulombEnergyAll;
    static thread_local Array_<Real, DuMM::NonbondAtomIndex>  localVdwScaleAll;
    static thread_local Array_<Real, DuMM::NonbondAtomIndex>  localCoulombScaleAll;

};

// because static
thread_local Real                                  NonbondedFullEnergyTask::localVdwEnergyAll = 0.0;
thread_local Real                                  NonbondedFullEnergyTask::localCoulombEnergyAll = 0.0;
thread_local Array_<Real, DuMM::NonbondAtomIndex>  NonbondedFullEnergyTask::localVdwScaleAll;
thread_local Array_<Real, DuMM::NonbondAtomIndex>  NonbondedFullEnergyTask::localCoulombScaleAll;

//........................class NonbondedFullEnergyTask.........................



//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - INCLUDING RIGID BODIES  
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the potential terms values that will not be 
// calculated, since they reffere to atoms inside rigid bodies and are constant 
// during dynamics. However, we need to evaluate the correct potential values 
// when comparing different constrainted regimens of the same compound.
// All variable's names were kept as in MolModel.
// TODO: Stretch custom terms are not implemented.
// TODO: GBSA potential is not calculated.
Real DuMMForceFieldSubsystemRep::
CalcFullPotEnergyIncludingRigidBodiesRep(const State& s) const {


    Real eStretch = 0;		// Bonds Stretch Potential
    Real eBend = 0;		    // Bonds Bend Potential
    Real eTorsion = 0;		// Bonds Torsion Potential
    Real eImproper = 0;		// Bond Improper Potential
    Real eVdW = 0; 		    // Van der Waals Potential
    Real eCoulomb = 0;		// Coulomb Potential

    Real eTotal = 0;

    const MultibodySystem&        mbs    = getMultibodySystem();
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    Vector_<Vec3>  AllAtomStation_G = getIncludedAtomStationCache(s);
    Vector_<Vec3>  AllAtomPos_G     = getIncludedAtomPositionCache(s);     
    AllAtomStation_G.resize( getNumAllAtoms() );
    AllAtomPos_G.resize( getNumAllAtoms() );

    // iterate all bodies
    for (DuMMIncludedBodyIndex dbx(0); dbx < AllBodies.size(); ++dbx) {

        const IncludedBody&      currentBody  = AllBodies[dbx];
        const MobilizedBodyIndex mbx     = currentBody.mobodIx;
        const MobilizedBody&     mobod   = matter.getMobilizedBody(mbx);

        const Transform&    X_GB  = mobod.getBodyTransform(s);
        const Rotation&     R_GB  = X_GB.R();
        const Vec3&         p_GB  = X_GB.p();


	    // First we make sure we have updated all atoms positions
	    int iax_count = 0;
        for (DuMM::IncludedAtomIndex iax=currentBody.beginAllAtoms;
             iax != currentBody.endAllAtoms; ++iax)
        {
            const Vec3& station_B_All = getAllAtomStation(iax);
            // atomic coordinates with respect to Ground frame
            const Vec3 p_BS_G = R_GB * station_B_All; 
            AllAtomStation_G[iax] = p_BS_G;
            AllAtomPos_G[iax]     = p_GB + p_BS_G;

            iax_count++;
        }

	    // iterate all bond forming atoms within currentBody
	    for (DuMMBondStarterIndex bsx = currentBody.beginAllBondStarterAtoms;
    		bsx != currentBody.endAllBondStarterAtoms; ++bsx)
    	{

	    const DuMM::IncludedAtomIndex a1num = AllbondStarterAtoms[bsx];

	    //Calculate Bonded Terms
	
	    if ( bondStretchGlobalScaleFactor !=0 
	    || customBondStretchGlobalScaleFactor != 0 ){
            eStretch += CalcFullPotEnergyBondStretch( a1num, AllAtomPos_G );
	    }
	    if ( bondBendGlobalScaleFactor !=0 
	    || customBondTorsionGlobalScaleFactor != 0 ){
            eBend += CalcFullPotEnergyBondBend( a1num, AllAtomPos_G );
	    }
	    if ( bondTorsionGlobalScaleFactor !=0 
	    || customBondTorsionGlobalScaleFactor != 0 ){
            eTorsion += CalcFullPotEnergyBondTorsion( a1num, AllAtomPos_G );
	    }
	    if ( amberImproperTorsionGlobalScaleFactor !=0 ){
            eImproper += CalcFullPotEnergyBondImproper( a1num, AllAtomPos_G );
	    }


	    }

    }

     //Calculate Nonbonded Terms
     if (usingMultithreaded) {

        NonbondedFullEnergyTask NonbondedFullTask
                        (*this, AllAtomPos_G, eVdW, eCoulomb);


        //NonbondedFullExecutor -> execute(NonbondedFullTask, Parallel2DExecutor::HalfMatrix); // inter body
        NonbondedFullExecutor -> execute(NonbondedFullTask, Parallel2DExecutor::HalfPlusDiagonal ); // inter and inside bodies
     }

     else {

         if (!(coulombGlobalScaleFactor == 0 && vdwGlobalScaleFactor == 0)) {
             CalcFullPotEnergyNonbondedSingleThread(AllAtomPos_G, eVdW, eCoulomb);
         }
     }

    eTotal = eStretch + eBend + eTorsion + eImproper + eVdW + eCoulomb;

/*
    TRACE( ("FUll Potential Bond Stretch: " + std::to_string( eStretch ) + "\n").c_str() );
    TRACE( ("FUll Potential Bond Bend   : " + std::to_string( eBend ) + "\n").c_str() );
    TRACE( ("FUll Potential Bond Torsion: " + std::to_string( eTorsion ) + "\n").c_str() );
    TRACE( ("FUll Potential Bond Impoper: " + std::to_string( eImproper ) + "\n").c_str() );
    TRACE( ("FUll Potential VdW         : " + std::to_string( eVdW ) + "\n").c_str() );
    TRACE( ("FUll Potential Coulomnb    : " + std::to_string( eCoulomb ) + "\n").c_str() );
    TRACE( ("FUll Potential TOTAL       : " + std::to_string( eTotal ) + "\n").c_str() );
//*/
    return eTotal;

}
//.........CALC FULL POTENTIAL ENERGY - INCLUDING RIGID BODIES .................




//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - NONBONDED TERM
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the full potential nonbonded terms
// the original function calcBodySubsetNonbondedForces was slightly modified to
// iterate Allbodies AllNonbonded lists and to calculate only energy
void DuMMForceFieldSubsystemRep::CalcFullPotEnergyNonbonded
        (DuMMIncludedBodyIndex                   dummBodIx,
         DuMMIncludedBodyIndex                   firstIx,
         DuMMIncludedBodyIndex                   lastIx,
         const Vector_<Vec3>&                    AllAtomPos_G,
         Array_<Real,DuMM::NonbondAtomIndex>&    vdwScaleAll,    // temps: all 1s
         Array_<Real,DuMM::NonbondAtomIndex>&    coulombScaleAll,
         Real&                                   eVdW,
         Real&                                   eCoulomb) const
{

    const IncludedBody& inclBod1 = AllBodies[dummBodIx];


    for (DuMM::NonbondAtomIndex nax1 = inclBod1.beginAllNonbondAtoms;
         nax1 != inclBod1.endAllNonbondAtoms; ++nax1)
    {
        DuMM::IncludedAtomIndex iax1 = getAllAtomIndexOfNonbondAtom(nax1);
        const IncludedAtom& a1 = getAllAtom(iax1);
        const ChargedAtomType& a1type = chargedAtomTypes[a1.chargedAtomTypeIndex];
        const DuMM::AtomClassIndex a1cnum = a1type.atomClassIx;


        const AtomClass&           a1class = atomClasses[a1cnum];
        const Vec3&                a1Pos_G = AllAtomPos_G[iax1];

        const Real q1Fac = coulombGlobalScaleFactor
                           * CoulombFac * a1type.partialCharge;

        scaleAllBondedAtoms(a1,vdwScaleAll,coulombScaleAll);

        for (DuMMIncludedBodyIndex dbx2 = firstIx; dbx2 <= lastIx; ++dbx2) {
            // REC BUG RESTORE We also want atoms inside the bodies: assert(dbx2 != dummBodIx);
            const IncludedBody& inclBod2 = AllBodies[dbx2];


            for (DuMM::NonbondAtomIndex nax2 = inclBod2.beginAllNonbondAtoms;
                 nax2 != inclBod2.endAllNonbondAtoms; ++nax2)
            {

                DuMM::IncludedAtomIndex iax2 =
                        getAllAtomIndexOfNonbondAtom(nax2);

                // assert(iax2 != iax1); REC BUG RESTORE
                if ((dummBodIx == dbx2) && (iax2 <= iax1)) {continue;} // don't double count pairs

                const IncludedAtom& a2 = getAllAtom(iax2);
                const ChargedAtomType& a2type  = chargedAtomTypes[a2.chargedAtomTypeIndex];
                const DuMM::AtomClassIndex a2cnum  = a2type.atomClassIx;

                const AtomClass& a2class = atomClasses[a2cnum];
                const Vec3&      a2Pos_G = AllAtomPos_G[iax2];

                const Vec3  r  = a2Pos_G - a1Pos_G;
                const Real d2 = r.normSqr();
                const Real  d = std::sqrt(d2);

                if( nonbondedMethod ==0 || ( nonbondedMethod ==1 && d <= nonbondedCutoff ) ) {

                    //TRACE( (std::string(" r ") + std::to_string(std::sqrt(d2))).c_str() );
                    const Real ood = 1 / d; // approx 40 flops
                    const Real ood2 = ood * ood;        // 1 flop


                    // Coulombic electrostatic force
                    const Real qq = coulombScaleAll[nax2] // 2 flops
                                    * q1Fac * a2type.partialCharge;
                    // e = scale*(1/(4*pi*e0)) *  q1*q2/d
                    eCoulomb += qq * ood;     // 1 flop


                    Real dij, eij;
                    if (a1cnum <= a2cnum) {
                        dij = a1class.vdwDij[a2cnum - a1cnum];
                        eij = a1class.vdwEij[a2cnum - a1cnum];
                    } else {
                        dij = a2class.vdwDij[a1cnum - a2cnum];
                        eij = a2class.vdwEij[a1cnum - a2cnum];
                    }

                    const Real ddij2 = dij * dij * ood2;   // (dmin_ij/d)^2
                    const Real ddij6 = ddij2 * ddij2 * ddij2;
                    const Real ddij12 = ddij6 * ddij6;

                    const Real eijScale = vdwGlobalScaleFactor * vdwScaleAll[nax2] * eij;

                    eVdW += eijScale * (ddij12 - 2 * ddij6);
                }
            }
        }

        unscaleAllBondedAtoms(a1,vdwScaleAll,coulombScaleAll);
    }
}
//..............CALC FULL POTENTIAL ENERGY - NONBONDED TERM ....................

//------------------------------------------------------------------------------
//           CALC FULL POTENTIAL ENERGY NONBONDED TERM SINGLE THREAD
//------------------------------------------------------------------------------
// This is the single-threaded nonbonded force calculation (not including GBSA).
// For each included body, we calculate its interactions with all
// higher-numbered included bodies.
void DuMMForceFieldSubsystemRep::CalcFullPotEnergyNonbondedSingleThread
        (const Vector_<Vec3>&                AllAtomPos_G,
         Real&                               eVdW,
         Real&                               eCoulomb) const
{

    for (DuMMIncludedBodyIndex inclBodyIx(0);
         inclBodyIx < getNumAllBodies(); ++inclBodyIx)
    {
        CalcFullPotEnergyNonbonded(
                inclBodyIx,
                DuMMIncludedBodyIndex(inclBodyIx + 1),
                DuMMIncludedBodyIndex(getNumAllBodies()-1),
                AllAtomPos_G,
                vdwScaleAllSingleThread,
                coulombScaleAllSingleThread,
                eVdW, eCoulomb);
    }
}
//...........CALC FULL POTENTIAL ENERGY NONBONDED TERM SINGLE THREAD............


//------------------------------------------------------------------------------
//                        SCALE ALL BONDED ATOMS
//------------------------------------------------------------------------------
// Some modifications to scaleBondedAtoms in order to scale All Nonbonded atoms

void DuMMForceFieldSubsystemRep::scaleAllBondedAtoms
        (const IncludedAtom& a,  // this is a nonbond atom
         Array_<Real,DuMM::NonbondAtomIndex>& vdwScaleAll,
         Array_<Real,DuMM::NonbondAtomIndex>& coulombScaleAll) const
{
    for (int i=0; i < (int)a.scale12All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale12All[i];
        vdwScaleAll[ix]=vdwScale12; coulombScaleAll[ix]=coulombScale12;
    }
    for (int i=0; i < (int)a.scale13All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale13All[i];
        vdwScaleAll[ix]=vdwScale13; coulombScaleAll[ix]=coulombScale13;
    }
    for (int i=0; i < (int)a.scale14All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale14All[i];
        vdwScaleAll[ix]=vdwScale14; coulombScaleAll[ix]=coulombScale14;
    }
    for (int i=0; i < (int)a.scale15All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale15All[i];
        vdwScaleAll[ix]=vdwScale15; coulombScaleAll[ix]=coulombScale15;
    }

}
//.........................SCALE ALL BONDED ATOMS...............................


//------------------------------------------------------------------------------
//                       UNSCALE ALL BONDED ATOMS
//------------------------------------------------------------------------------
// Some modifications to unscaleBondedAtoms in order to unscale All Nonbonded atoms
void DuMMForceFieldSubsystemRep::unscaleAllBondedAtoms
        (const IncludedAtom& a,
         Array_<Real,DuMM::NonbondAtomIndex>& vdwScaleAll,
         Array_<Real,DuMM::NonbondAtomIndex>& coulombScaleAll) const
{
    for (int i=0; i < (int)a.scale12All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale12All[i];
        vdwScaleAll[ix]=coulombScaleAll[ix]=1;
    }
    for (int i=0; i < (int)a.scale13All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale13All[i];
        vdwScaleAll[ix]=coulombScaleAll[ix]=1;
    }
    for (int i=0; i < (int)a.scale14All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale14All[i];
        vdwScaleAll[ix]=coulombScaleAll[ix]=1;
    }
    for (int i=0; i < (int)a.scale15All.size(); ++i) {
        const DuMM::NonbondAtomIndex ix = a.scale15All[i];
        vdwScaleAll[ix]=coulombScaleAll[ix]=1;
    }
}
//........................UNSCALE ALL BONDED ATOMS..............................



//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - BOND STRETCH TERM 
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the full potential bond stretch terms
// All variable's names were kept as in MolModel. 
// Warning: Be sure that atom positions were updated before function call
// TODO: Custom terms are not implemented.
Real DuMMForceFieldSubsystemRep::
CalcFullPotEnergyBondStretch(const DuMM::IncludedAtomIndex a1num, 
Vector_<Vec3>  AllAtomPos_G ) const {

    Real eCurrentStretch = 0;
    	
    const IncludedAtom&  a1    	     = getAllAtom( a1num );
    const Vec3&          a1Pos_G     = AllAtomPos_G[ a1num ] ;

    for (DuMM::IncludedAtomIndex b12(0); b12 < a1.force12All.size(); ++b12) 
    {
        const DuMM::IncludedAtomIndex a2num = a1.force12All[b12]; 
      
        const Vec3& a2Pos_G     = AllAtomPos_G[a2num];
        const Vec3  r           = a2Pos_G - a1Pos_G;
        const Real  d           = r.norm();
        const BondStretch& bs = *a1.stretchAll[b12];

	if (bs.hasBuiltinTerm()) {
            const Real x = d - bs.d0;
            const Real kx = bondStretchGlobalScaleFactor * bs.k * x;
            eCurrentStretch +=  kx*x; // no factor of 1/2!     
        } 
        else
            SimTK_ASSERT_ALWAYS(1, "CustomBondStretch not implemented");
    }

    return eCurrentStretch;
}
//.............CALC FULL POTENTIAL ENERGY - BOND STRETCH TERM...................


//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - BOND BEND TERM 
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the full potential bond bend terms 
// Warning: Be sure that atom positions were updated before function call
Real DuMMForceFieldSubsystemRep::
CalcFullPotEnergyBondBend(const DuMM::IncludedAtomIndex a1num, 
Vector_<Vec3>  AllAtomPos_G ) const {

    Real eCurrentBend = 0;
    	
    const IncludedAtom&  a1    	     = getAllAtom( a1num );
    const Vec3&          a1Pos_G     = AllAtomPos_G[ a1num ] ;

    for (DuMM::IncludedAtomIndex b13(0); b13 < a1.force13All.size(); ++b13) 
    {
        const DuMM::IncludedAtomIndex a2num = a1.force13All[b13][0];
        const DuMM::IncludedAtomIndex a3num = a1.force13All[b13][1];
      
        const Vec3& a2Pos_G     = AllAtomPos_G[a2num];
        const Vec3& a3Pos_G     = AllAtomPos_G[a3num];

	Vec3 f1, f2, f3;  // we don't need them actually
        Real  angle, e;

        const BondBend& bb = *a1.bendAll[b13];

        bb.calculateAtomForces(a2Pos_G, a1Pos_G, a3Pos_G,
                               bondBendGlobalScaleFactor, 
                               customBondBendGlobalScaleFactor, angle, 
			       e, f2, f1, f3);
	eCurrentBend += e;
	
    }

    return eCurrentBend;
}
//.............CALC FULL POTENTIAL ENERGY - BOND BEND TERM......................


//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - BOND TORSION TERM 
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the full potential bond torsion terms 
// Warning: Be sure that atom positions were updated before function call
Real DuMMForceFieldSubsystemRep::
CalcFullPotEnergyBondTorsion(const DuMM::IncludedAtomIndex a1num, 
Vector_<Vec3>  AllAtomPos_G ) const {

    Real eCurrentTorsion = 0;
    	
    const IncludedAtom&  a1    	     = getAllAtom( a1num );
    const Vec3&          a1Pos_G     = AllAtomPos_G[ a1num ] ;

    for (DuMM::IncludedAtomIndex b14(0); b14 < a1.force14All.size(); ++b14) 
    {
        const DuMM::IncludedAtomIndex a2num = a1.force14All[b14][0];
        const DuMM::IncludedAtomIndex a3num = a1.force14All[b14][1];
        const DuMM::IncludedAtomIndex a4num = a1.force14All[b14][2];
      
        const Vec3& a2Pos_G     = AllAtomPos_G[a2num];
        const Vec3& a3Pos_G     = AllAtomPos_G[a3num];
        const Vec3& a4Pos_G     = AllAtomPos_G[a4num];

	Vec3 f1, f2, f3, f4;  // we don't need them actually
        Real  torsion, e;

        const BondTorsion& bt = *a1.torsionAll[b14];

        bt.calculateAtomForces(a1Pos_G, a2Pos_G, a3Pos_G, a4Pos_G,
                               bondTorsionGlobalScaleFactor, 
                               customBondTorsionGlobalScaleFactor, torsion, 
			       e, f1, f2, f3, f4);
	eCurrentTorsion += e;
	
    }

    return eCurrentTorsion;
}
//.............CALC FULL POTENTIAL ENERGY - BOND TORSION TERM...................


//------------------------------------------------------------------------------
//          CALC FULL POTENTIAL ENERGY - BOND IMPROPER TERM 
//------------------------------------------------------------------------------
// GMolModel needs to keep track of the full potential bond torsion terms 
// Warning: Be sure that atom positions were updated before function call
Real DuMMForceFieldSubsystemRep::
CalcFullPotEnergyBondImproper(const DuMM::IncludedAtomIndex a1num, 
Vector_<Vec3>  AllAtomPos_G ) const {

    Real eCurrentImproper = 0;
    	
    const IncludedAtom&  a1    	     = getAllAtom( a1num );
    const Vec3&          a1Pos_G     = AllAtomPos_G[ a1num ] ;

    for (DuMM::IncludedAtomIndex b14(0); b14 < a1.forceImproper14All.size(); ++b14) 
    {
        const DuMM::IncludedAtomIndex a2num = a1.forceImproper14All[b14][0];
        const DuMM::IncludedAtomIndex a3num = a1.forceImproper14All[b14][1];
        const DuMM::IncludedAtomIndex a4num = a1.forceImproper14All[b14][2];
      
        const Vec3& a2Pos_G     = AllAtomPos_G[a2num];
        const Vec3& a3Pos_G     = AllAtomPos_G[a3num];
        const Vec3& a4Pos_G     = AllAtomPos_G[a4num];

	Vec3 f1, f2, f3, f4;  // we don't need them actually
        Real improper, e;

        const BondTorsion& bt = *a1.aImproperTorsionAll[b14];

        bt.calculateAtomForces(a2Pos_G, a3Pos_G, a1Pos_G, a4Pos_G,
                               amberImproperTorsionGlobalScaleFactor, 
                               0, // no custom improper in MolModel - why ?
			       improper, e, f2, f3, f1, f4);
	eCurrentImproper += e;
	
    }

    return eCurrentImproper;
}
//.............CALC FULL POTENTIAL ENERGY - BOND IMPROPER TERM..................

