#ifndef SimTK_MOLMODEL_COMPOUNDSYSTEM_H_
#define SimTK_MOLMODEL_COMPOUNDSYSTEM_H_

#include <map>

#include "molmodel/internal/Compound.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"
#include "molmodel/internal/MolecularMechanicsSystem.h"
#include "molmodel/internal/common.h"

#include "SimTKsimbody.h"


namespace SimTK {

enum class RootMobility : std::uint8_t {
    Free = 0,
    Cartesian,
    Weld,
    FreeLine,
    Ball,
    Pin
};

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
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(CompoundSystem, CompoundIndex);

    /// default constructor
    CompoundSystem() {
    }

    /// destructor
    ~CompoundSystem() {
    }

    /**
     * Install a new Compound into this system. We take over ownership of the Compound's
     * representation from the given handle, leaving that handle as a reference to our
     * new Compound.
     * It is an error if the given handle wasn't the owner of the Compound.
     */
    void adoptCompound(
        Compound& child,                                 ///< Compound to be incorporated
        const Transform& compoundTransform = Transform() ///< location and orientation of the Compound
    ) {
        // const Compound::Index id((int)compounds.size());

        // Create a new empty Compound and get a reference to it
        compounds.push_back(new Compound((CompoundRep*)0)); // grow
        Compound& newCompound = *compounds.back();

        // Transfer ownership to the supplied new empty handle
        child.disown(newCompound);

        // Now tell the Compound object its owning CompoundSystem and id within
        // that System.
        // c.setCompoundSystem(*this, id);
        newCompound.setMultibodySystem(*this);

        // Save transform
        // March 6, 2008 -- adjust internal Transform of Compound, rather than
        // saving the Transform in CompoundSystem
        newCompound.setTopLevelTransform(compoundTransform * newCompound.getTopLevelTransform());
        // std::cout << "SP_NEW  CompoundSystem::adoptCompound Top transforms:" << std::endl;
        // std::cout << compoundTransform;
        // std::cout << newCompound.getTopLevelTransform();


        // assert((int) compoundTransforms.size() == (int) id);
        // compoundTransforms.push_back(compoundTransform);

        // return id;
    }

    /**
     * \return number of top-level Compounds adopted by this CompoundSystem
     */
    size_t getNumCompounds() const {
        return compounds.size();
    }

    /// \return read-only reference to an adopted Compound
    const Compound& getCompound(CompoundIndex i ///< integer index of Compound
    ) const {
        return *compounds.at(i);
    }

    /// \return mutable reference to an adopted Compound
    Compound& updCompound(CompoundIndex i ///< integer index of Compound
    ) {
        return *compounds.at(i);
    }

    private:
    // suppress
    CompoundSystem(const CompoundSystem&);
    CompoundSystem& operator=(const CompoundSystem&);

    std::vector<Compound*> compounds;
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_COMPOUNDSYSTEM_H_
