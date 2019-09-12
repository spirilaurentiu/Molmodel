#ifndef SIMTK_MOLMODEL_VANDERWAALS_FORCE_H_
#define SIMTK_MOLMODEL_VANDERWAALS_FORCE_H_

#include "SimTKsimbody.h"
#include <vector>
#include <iterator>

namespace SimTK {

namespace md {
    const Real angstroms = 0.10; // angstroms to nanometers
    const Real nanometers = 1.0;
    const Real square_nanometers = 1.0;

    const Real kilocalories_per_mole = 4.184; // calories to joules
    const Real kilojoules_per_mole = 1.0;

    const Real daltons = 1.0;
}

class VanDerWaalsForce : public Force::Custom::Implementation
{
public:
    typedef AtomSubsystem::AtomIndex AtomIndex;

    // Definitions to elide previous use of Boost::Units
    typedef Real length_t;
    typedef Real area_t;
    typedef Real energy_t;
    typedef Real inverse_area_t;
    typedef Vec3 location_t;
    typedef Vec3 force_t;

    class VdwAtom
    {
    public:
        VdwAtom() 
            : rMin(1.90*md::angstroms), wellDepth(0.1*md::kilocalories_per_mole)
        {}
        
        VdwAtom(AtomIndex atomIndex, length_t rMin, energy_t wellDepth) 
            : atomIndex(atomIndex), rMin(rMin), wellDepth(wellDepth)
        {}

        AtomIndex atomIndex;
        length_t rMin; // nanometers, minimum energy half-distance
        energy_t wellDepth; // kilojoules per mole
    };

    // Sphere for containing atoms with forces sort of like other atoms
    struct WallSphere
    {
        length_t radius;
        location_t center;

        length_t vdwRadius;
        energy_t wellDepth;
    };

    VanDerWaalsForce(const AtomSubsystem& atomSubsystem, length_t cutoff = 0.0 * nanometers) 
        : atomSubsystem(atomSubsystem), cutoff(cutoff), cutoffSquared(cutoff * cutoff)
    {}

    bool dependsOnlyOnPositions() const {return true;}

    void addAtom(AtomSubsystem::AtomIndex atomIndex, length_t rMin, energy_t wellDepth) 
    {
        if (AtomSubsystem::AtomIndex(vdwAtoms.size()) <= atomIndex)
            vdwAtoms.resize(atomIndex + 1);
        
        vdwAtoms[atomIndex] = VdwAtom(atomIndex, rMin, wellDepth);
    }

    void addBoundarySphere(
            const location_t& center, 
            length_t radius,
            length_t vdwRadius,
            energy_t wellDepth) 
    {
        WallSphere sphere;
        sphere.center = center;
        sphere.radius = radius;
        sphere.vdwRadius = vdwRadius;
        sphere.wellDepth = wellDepth;
        wallSpheres.push_back(sphere);
    }

    void decorateAtoms(DecorationSubsystem& decorations, Real scale = 1.0) const {
        for (AtomIndex a(0); a < vdwAtoms.size(); ++a) 
        {      
            const VdwAtom& vdwAtom = vdwAtoms[a];
            
            if (!vdwAtom.atomIndex.isValid()) continue;
                        
            const AtomSubsystem::Atom& atom = atomSubsystem.getAtom(a);

            Vec3 color(0.5, 0.5, 0.5); // default color gray

            // Color argon cyan
            // Argon mass (39.9) is very close to Calcium (40.1)
            if (atom.getMass() < 39.5) 
                ;
            else if (atom.getMass() < 40.0) 
                color = Vec3(0.5, 1, 1); // argon

            decorations.addBodyFixedDecoration(
                atom.getMobilizedBodyIndex(), 
                atom.getStationInBodyFrame(),
                DecorativeSphere(vdwAtom.rMin * scale).setColor(color)
            );
        }
    }

    void decorateBoundingSpheres(DecorationSubsystem& decorations) const {
        for (size_t s(0); s < wallSpheres.size(); ++s) 
        {      
            Vec3 color(1.0, 1.0, 0.6); // pale yellow

            const WallSphere& sphere = wallSpheres[s];

            decorations.addBodyFixedDecoration(
                atomSubsystem.getMultibodySystem().getMatterSubsystem().getGround().getMobilizedBodyIndex(),
                sphere.center,
                DecorativeSphere(sphere.radius).setColor(color).setOpacity(0.2)
            );
        }
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
    {
        // Forces on pairs of atoms
        for ( AtomSubsystem::PairIterator pair = 
                    atomSubsystem.pairBegin(state, cutoff);
              pair != atomSubsystem.pairEnd();
              ++pair) 
        {
            const VdwAtom& vdwAtom1 = vdwAtoms[pair->first];
            const VdwAtom& vdwAtom2 = vdwAtoms[pair->second];

            if (!vdwAtom1.atomIndex.isValid()) continue;
            if (!vdwAtom2.atomIndex.isValid()) continue;
            
            calcForce(vdwAtom1, vdwAtom2, state, bodyForces);
        }

        // Apply any bounding spheres
        for ( std::vector<WallSphere>::const_iterator sphereI = wallSpheres.begin();
              sphereI != wallSpheres.end(); ++sphereI) {
                for ( std::vector<VdwAtom>::const_iterator atomI = vdwAtoms.begin();
                    atomI != vdwAtoms.end();  ++atomI) {
                        if ( ! atomI->atomIndex.isValid() ) continue;
                        calcForce(*sphereI, *atomI, state, bodyForces);
                }
        }
    }

    Real calcPotentialEnergy(const State& state) const {
        energy_t energy = 0.0 * md::kilojoules_per_mole;

        for (   AtomSubsystem::PairIterator pair = 
                        atomSubsystem.pairBegin(state, cutoff);
                pair != atomSubsystem.pairEnd();
                ++pair) 
        {
            const VdwAtom& vdwAtom1 = vdwAtoms[pair->first];
            const VdwAtom& vdwAtom2 = vdwAtoms[pair->second];

            if (!vdwAtom1.atomIndex.isValid()) continue;
            if (!vdwAtom2.atomIndex.isValid()) continue;
            
            energy += calcPotentialEnergy(vdwAtom1, vdwAtom2, state);
        }

        // Apply any bounding spheres
        for ( std::vector<WallSphere>::const_iterator sphereI = wallSpheres.begin();
              sphereI != wallSpheres.end(); ++sphereI) {
                for ( std::vector<VdwAtom>::const_iterator atomI = vdwAtoms.begin();
                    atomI != vdwAtoms.end();  ++atomI) {
                        if ( ! atomI->atomIndex.isValid() ) continue;
                        energy += calcPotentialEnergy(*sphereI, *atomI, state);
                }
        }

        return energy;
    }

protected:

    struct ForceAndEnergy {
        force_t force;
        energy_t energy;
    };

    // calcForceAndEnergy.  
    // Agnostic to DuMMForceFieldSubsystem vs AtomSubsystem.  Location in ground are arguments.
    // Agnostic to combining rules: precombined values are arguments to calcForceAndEnergy
    // Bonded atom scaling should be performed outside of this routine
    // Return force and energy on second location
    // TODO softer potential with non-infinity center
    ForceAndEnergy calcForceAndEnergy( 
        const location_t& r, // vector from atom1 to atom2
        Real dij, Real eij) const // pre-combined Dij and eij
    {
        ForceAndEnergy forceAndEnergy;

        const area_t d2 = r.normSqr(); // 5 flops

        // If cutoff is present, force drops to zero abruptly at and beyond the cutoff radius.
        // But force is not otherwise affected, since within the cutoff radius
        // the truncated potential is a constant offset to the potential energy.
        if ( (cutoff > 0.0) && (d2 > cutoffSquared) )
            return forceAndEnergy;

        // Don't need to calculate 1/d unless we compute Coulomb at the same time.
        inverse_area_t ood2(1.0 / d2);        // 1 flop

        // should be dimensionless
        const Real ddij2  = dij*dij*ood2;        // (dmin_ij/d)^2 (2 flops)
        const Real ddij6  = ddij2*ddij2*ddij2;   // 2 flops
        const Real ddij12 = ddij6*ddij6;         // 1 flop

        // Force
        const energy_t fVdw = 12.0 * eij * (ddij12 - ddij6);   // factor of 1/d^2 missing (3 flops)
        forceAndEnergy.force = fVdw * ood2 * r;      // to apply to atom j on b2 (5 flops)

        // Energy
        forceAndEnergy.energy = eij * (ddij12 - 2*ddij6); // 3 flops

        // offset potential to reach zero at cutoff
        if (cutoff > 0.0) 
        {
            energy_t truncationOffset = 0.0;
            assert(d2 < cutoffSquared); // should have short circuited above
            const Real off2 = dij * dij / cutoffSquared; // off2 = (rmin/cutoff)^2
            const Real off6 = off2 * off2 * off2;
            truncationOffset = eij * ( (off6 - 2.0) * off6);
            assert(truncationOffset < 0.0); // unreasonable to put cutoff within repulsive region

            forceAndEnergy.energy -= truncationOffset;
        }

        return forceAndEnergy;
    }

    // Calculate force between two atoms
    void calcForce(const VdwAtom& vdwAtom1, const VdwAtom& vdwAtom2, const State& state, Vector_<SpatialVec>& bodyForces) const
    {
        Real vdwScale = 1.0;  // TODO scale bonded atoms

        const AtomSubsystem::Atom& atom1 = atomSubsystem.getAtom(vdwAtom1.atomIndex);
        const AtomSubsystem::Atom& atom2 = atomSubsystem.getAtom(vdwAtom2.atomIndex);

        const Vec3& a1Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom1.atomIndex);
        const Vec3& a2Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom2.atomIndex);

        // TODO better combining rules
        const length_t dij = vdwAtom1.rMin + vdwAtom2.rMin;
        const energy_t eij = std::sqrt(vdwAtom1.wellDepth * vdwAtom2.wellDepth);

        const force_t fj = calcForceAndEnergy(a2Pos_G - a1Pos_G, dij, eij).force * vdwScale;

        const Vec3& a1Station_G = atomSubsystem.getAtomBodyStationInGround(state, vdwAtom1.atomIndex);
        const Vec3& a2Station_G = atomSubsystem.getAtomBodyStationInGround(state, vdwAtom2.atomIndex);

        bodyForces[atom2.getMobilizedBodyIndex()] += SpatialVec( a2Station_G % fj, fj);   // 15 flops
        bodyForces[atom1.getMobilizedBodyIndex()] -= SpatialVec( a1Station_G % fj, fj);   // 15 flops
    }

    energy_t calcPotentialEnergy(const VdwAtom& vdwAtom1, const VdwAtom& vdwAtom2, const State& state) const 
    {
        Real vdwScale = 1.0;  // TODO scale bonded atoms

        const AtomSubsystem::Atom& atom1 = atomSubsystem.getAtom(vdwAtom1.atomIndex);
        const AtomSubsystem::Atom& atom2 = atomSubsystem.getAtom(vdwAtom2.atomIndex);

        const Vec3& a1Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom1.atomIndex);
        const Vec3& a2Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom2.atomIndex);

        // TODO better combining rules
        const length_t dij = vdwAtom1.rMin + vdwAtom2.rMin;
        const energy_t eij = std::sqrt(vdwAtom1.wellDepth * vdwAtom2.wellDepth);

        return calcForceAndEnergy(a2Pos_G - a1Pos_G, dij, eij).energy * vdwScale;
    }

    // Calculate force between one atom and a bounding sphere
    void calcForce(const WallSphere& sphere, const VdwAtom& vdwAtom2, const State& state, Vector_<SpatialVec>& bodyForces) const
    {
        const AtomSubsystem::Atom& atom2 = atomSubsystem.getAtom(vdwAtom2.atomIndex);

        const Vec3& a2Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom2.atomIndex);

        // TODO better combining rules
        const length_t dij = sphere.vdwRadius + vdwAtom2.rMin;
        const energy_t eij = std::sqrt(sphere.wellDepth * vdwAtom2.wellDepth);

        Vec3 dv = a2Pos_G - sphere.center; // Atom relative to center of bounding sphere

        // Distance from atom to closest edge of bounding sphere
        Real centerDist = dv.norm();
        UnitVec3 centerDir(dv);
        
        // Vector to closest edge of sphere
        Vec3 rA = (centerDist - sphere.radius) * centerDir;
        // And to furthest edge of sphere, to avoid singularity at center
        Vec3 rB = (sphere.radius + centerDist) * centerDir;

        force_t fj = calcForceAndEnergy(rA, dij, eij).force;
        fj += calcForceAndEnergy(rB, dij, eij).force;

        const Vec3& a2Station_G = atomSubsystem.getAtomBodyStationInGround(state, vdwAtom2.atomIndex);

        // No force on the sphere - it's just a forcefield
        bodyForces[atom2.getMobilizedBodyIndex()] += SpatialVec( a2Station_G % fj, fj);   // 15 flops
    }

    energy_t calcPotentialEnergy(const WallSphere& sphere, const VdwAtom& vdwAtom2, const State& state) const
    {
        const AtomSubsystem::Atom& atom2 = atomSubsystem.getAtom(vdwAtom2.atomIndex);

        const Vec3& a2Pos_G = atomSubsystem.getAtomLocationInGround(state, vdwAtom2.atomIndex);

        // TODO better combining rules
        const length_t dij = sphere.vdwRadius + vdwAtom2.rMin;
        const energy_t eij = std::sqrt(sphere.wellDepth * vdwAtom2.wellDepth);

        Vec3 dv = a2Pos_G - sphere.center; // Atom relative to center of bounding sphere

        // Distance from atom to closest edge of bounding sphere
        Real centerDist = dv.norm();
        UnitVec3 centerDir(dv);
        
        // Vector to closest edge of sphere
        Vec3 rA = (sphere.radius - centerDist) * centerDir;
        // And to furthest edge of sphere, to avoid singularity at center
        Vec3 rB = -(sphere.radius + centerDist) * centerDir;

        energy_t energy = calcForceAndEnergy(rA, dij, eij).energy;
        energy += calcForceAndEnergy(rB, dij, eij).energy;

        return energy;
    }

    const AtomSubsystem& atomSubsystem;

    // TODO - place these data into State
    std::vector<VdwAtom> vdwAtoms;
    length_t cutoff;
    area_t cutoffSquared;
    std::vector<WallSphere> wallSpheres;
};

} // namespace SimTK

#endif /* SIMTK_MOLMODEL_VANDERWAALS_FORCE_H_ */
