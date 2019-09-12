#include "molmodel/internal/VanderWallSphere.h"

namespace SimTK {

struct ForceAndEnergy {
    Vec3 force;
    Real energy;
};

class VanderWallSphereImpl : public Force::Custom::Implementation {
public:
	VanderWallSphereImpl(DuMMForceFieldSubsystem& dmm, Vec3 c, Real r, Real v, Real dpth = 0.0)
	  : Force::Custom::Implementation(),
	    dumm(dmm),
	    center(c),
	    radius(r),
	    vdwRadius(v),
	    vdwWellDepth(dpth)
	{}
	
	bool dependsOnlyOnPositions() const {return true;}
	
	void calcForce(
			const State& state, 
			Vector_< SpatialVec >& bodyForces, 
			Vector_< Vec3 >& particleForces, 
			Vector& mobilityForces) const
	{
	    const MultibodySystem& mbs    = static_cast<const MultibodySystem&>(dumm.getSystem()); // my owner
	    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

		for (DuMM::AtomIndex a(0); a < dumm.getNumAtoms(); ++a) 
		{
			MobilizedBodyIndex body = dumm.getAtomBody(a);
			
			const Transform&  G_X_B  = matter.getMobilizedBody(body).getBodyTransform(state);
			Vec3 B_v_atom = dumm.getAtomStationOnBody(a);
			Vec3 G_v_atom = G_X_B * B_v_atom; // Atom location in ground frame
			
			Vec3 dv = G_v_atom - center; // Atom relative to center of bounding sphere
			
            DuMM::AtomClassIndex atomClassIndex = dumm.getAtomClassIndex(a);
            
            // Force from near side of sphere
            const Vec3 fj = compAtomForceAndEnergy(dv, atomClassIndex).force;
            bodyForces[body] -= SpatialVec( G_v_atom % fj, fj);
		}
	}
	
	Real calcPotentialEnergy(const State & state) const
	{
	    const MultibodySystem& mbs    = static_cast<const MultibodySystem&>(dumm.getSystem()); // my owner
	    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

		Real energy = 0.0;
		
		for (DuMM::AtomIndex a(0); a < dumm.getNumAtoms(); ++a) 
		{
			MobilizedBodyIndex body = dumm.getAtomBody(a);
			
			const Transform&  G_X_B  = matter.getMobilizedBody(body).getBodyTransform(state);
			Vec3 B_v_atom = dumm.getAtomStationOnBody(a);
			Vec3 G_v_atom = G_X_B * B_v_atom; // Atom location in ground frame
			
			Vec3 dv = G_v_atom - center; // Atom relative to center of bounding sphere
			
            DuMM::AtomClassIndex atomClassIndex = dumm.getAtomClassIndex(a);
            
            energy += compAtomForceAndEnergy(dv, atomClassIndex).energy;
            
		}
		
		return energy;
	}
	
private:
    
    ForceAndEnergy compAtomForceAndEnergy(const Vec3& atomPosRelCenter, DuMM::AtomClassIndex atomClassIndex) const
    {
        ForceAndEnergy answer;
        
        // Distance from atom to closest edge of bounding sphere
        Real centerDist = atomPosRelCenter.norm();
        UnitVec3 centerDir(atomPosRelCenter);
        
        // Vector to closest edge of sphere
        Vec3 rA = (radius - centerDist) * centerDir;
        // And to furthest edge of sphere, to avoid singularity at center
        Vec3 rB = -(radius + centerDist) * centerDir;

        // Treat atom as particle i, wall as particle j
        Real rj = vdwRadius;
        Real ej = vdwWellDepth;

        // Tget atom class parameters di, ei
        Real ri = dumm.getVdwRadius(atomClassIndex);
        Real ei = dumm.getVdwWellDepth(atomClassIndex);

        // Use Lorentz Bethelot mixing rule
        Real dij = 0.5 * (ri + rj); // arithmetic mean
        Real eij = std::sqrt(ei * ej); // geometric mean

        const Real d2A = rA.normSqr(); // 5 flops
        const Real d2B = rB.normSqr(); // 5 flops

        Real oodA = 1.0/std::sqrt(d2A);
        Real ood2A = oodA*oodA;
        Real oodB = 1.0/std::sqrt(d2B);
        Real ood2B = oodB*oodB;

        const Real ddij2A  = dij*dij*ood2A;        // (dmin_ij/d)^2 (2 flops)
        const Real ddij6A  = ddij2A*ddij2A*ddij2A;   // 2 flops
        const Real ddij12A = ddij6A*ddij6A;         // 1 flop
        const Real ddij2B  = dij*dij*ood2B;        // (dmin_ij/d)^2 (2 flops)
        const Real ddij6B  = ddij2B*ddij2B*ddij2B;   // 2 flops
        const Real ddij12B = ddij6B*ddij6B;         // 1 flop

        const Real eijScale = dumm.getVdwGlobalScaleFactor() * eij; // 1 flops

        answer.energy =  eijScale * (ddij12A - 2*ddij6A); // 3 flops
        answer.energy +=  eijScale * (ddij12B - 2*ddij6B); // 3 flops
        
        const Real fVdwA = 12 * eijScale * (ddij12A - ddij6A);   // factor of 1/d^2 missing (3 flops)
        const Real fVdwB = 12 * eijScale * (ddij12B - ddij6B);   // factor of 1/d^2 missing (3 flops)
        answer.force = fVdwA*ood2A * rA;      // to apply to atom j on b2 (5 flops)
        answer.force += fVdwB*ood2B * rB;
        
        return answer;
    }
	
	DuMMForceFieldSubsystem& dumm;
	Vec3 center;
	Real radius;
	Real vdwRadius;
	Real vdwWellDepth;
	
};

VanderWallSphere::VanderWallSphere(GeneralForceSubsystem& forces, DuMMForceFieldSubsystem& dumm, Vec3 center, Real radius, Real vdwRadius, Real wellDepth)
	  : Force::Custom( forces, new VanderWallSphereImpl(dumm, center, radius, vdwRadius, wellDepth) )
	{}

} // namespace SimTK

