#ifndef SimTK_MOLMODEL_COMPOUNDSYSTEM_H_
#define SimTK_MOLMODEL_COMPOUNDSYSTEM_H_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/MolecularMechanicsSystem.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"

#include <map>

namespace SimTK {

/**
 * \brief Derived class of MolecularMechanicsSystem that knows how to model molmodel Compounds
 *
 * \todo merge this class with MolecularMechanicsSystem
 */
class SimTK_MOLMODEL_EXPORT CompoundSystem : public MolecularMechanicsSystem {
public:

    /** @class SimTK::CompoundSystem::CompoundIndex
     * Compound::Index type is an integer index into subcompounds of a Compound.  It is NOT
     * instrinsic to the subcompound, but represents the relationship between a subcompound
     * and precisely one of its parent compounds.
     */
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(CompoundSystem,CompoundIndex);

    /// default constructor
    CompoundSystem() {}

    /// destructor
    ~CompoundSystem() {
        for (int i=0; i < (int)compounds.size(); ++i)
            delete compounds[i];
    }

    /**
     * Install a new Compound into this system. We take over ownership of the Compound's
     * representation from the given handle, leaving that handle as a reference to our
     * new Compound.
     * It is an error if the given handle wasn't the owner of the Compound.
     */
    void adoptCompound(
        Compound& child, ///< Compound to be incorporated
        const Transform& compoundTransform = Transform() ///< location and orientation of the Compound
        ) 
    {
        // const Compound::Index id((int)compounds.size());

        compounds.push_back(new Compound((CompoundRep*)0)); // grow
        Compound& c = *compounds.back(); // refer to the empty handle we just created

        child.disown(c); // transfer ownership to c

        // Now tell the Compound object its owning CompoundSystem and id within
        // that System.
        // c.setCompoundSystem(*this, id);
        c.setMultibodySystem(*this);

        // Save transform
        // March 6, 2008 -- adjust internal Transform of Compound, rather than 
        // saving the Transform in CompoundSystem
        c.setTopLevelTransform(compoundTransform * c.getTopLevelTransform());
        // assert((int) compoundTransforms.size() == (int) id);
        // compoundTransforms.push_back(compoundTransform);

        // return id;
    }

    /** Instantiate a Simbody model representing the adopted Compounds, using
    the same base body-to-ground connection type (Free or Weld mobilizer) for 
    all top-level Compounds. If you want some compounds free and others welded,
    use modelOneCompound() repeatedly instead. 
    @param[in] mobilizedBodyType    Value must be exactly "Free" or "Weld"
                                    including capitalization. The default is
                                    "Free".
    @bug "Weld" will only be applied to those compounds for which the mobilizer 
    would otherwise have been "Free" (six degrees of freedom). For compounds of 
    just a few atoms other mobilizers may be used and "Weld" is ignored. **/
    void modelCompounds(String mobilizedBodyType = "Free");

    /** Build the Simbody model one compound at a time to allow differing
    base body-to-ground connection types. If you use this method you must 
    use it for all compounds; you can't mix with modelCompounds(). For 
    example:
    @code
      CompoundSystem sys; // ... with compounds already adopted.
      for (CompoundSystem::CompoundIndex c(0); c < sys.getNumCompounds(); ++c) 
          modelOneCompound(c, "Weld");
    @endcode 

    @param[in] mobilizedBodyType    Value must be exactly "Free" or "Weld"
                                    including capitalization. The default is
                                    "Free".
    @bug "Weld" will only be applied if the mobilizer would otherwise have been
    "Free" (six degrees of freedom). For compounds of just a few atoms other
    mobilizers may be used. **/
    void modelOneCompound(CompoundIndex compoundId, 
                          String        mobilizedBodyType = "Free");

    /**
     * \return number of top-level Compounds adopted by this CompoundSystem
     */
    size_t             getNumCompounds() const {return compounds.size();}

    /// \return read-only reference to an adopted Compound
    const Compound& getCompound(CompoundIndex i ///< integer index of Compound
        ) const {return *compounds.at(i);}

    /// \return mutable reference to an adopted Compound
    Compound&       updCompound(CompoundIndex i ///< integer index of Compound
        )       {return *compounds.at(i);}

private:


    void setClusterCompound(DuMM::ClusterIndex clusterIx, const Compound& compound) 
    {
        assert(! clusterHasCompound(clusterIx));
        compoundPtrsByClusterIndex[clusterIx] = &compound;
        assert(clusterHasCompound(clusterIx));
    }

    bool clusterHasCompound(DuMM::ClusterIndex clusterIx) {
        return compoundPtrsByClusterIndex.find(clusterIx) != compoundPtrsByClusterIndex.end();
    }

    void generateTopologyFromCompounds();

    // suppress
    CompoundSystem(const CompoundSystem&);
    CompoundSystem& operator=(const CompoundSystem&);

    // std::vector<Transform> compoundTransforms;

    // retarded visual studio compiler complains about being unable to 
    // export private stl class members
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    std::vector<Compound*> compounds;
    std::map<DuMM::ClusterIndex, const Compound*> compoundPtrsByClusterIndex;
    std::map<DuMM::ClusterIndex, int> clusterAtomCounts;
    std::map<DuMM::ClusterIndex, MobilizedBodyIndex> clusterBodies;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_COMPOUNDSYSTEM_H_

