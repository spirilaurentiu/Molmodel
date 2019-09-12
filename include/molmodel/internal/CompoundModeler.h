#ifndef SimTK_MOLMODEL_COMPOUND_MODELER_H_
#define SimTK_MOLMODEL_COMPOUND_MODELER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
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


#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/AtomSubsystem.h"

namespace SimTK {

/// Turns a compound into a multibody system
class SimTK_MOLMODEL_EXPORT CompoundModeler 
{
public:

    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(CompoundModeler, RigidUnitIndex);
    
    // RigidUnit data structure is for use in modelCompounds() method
    class RigidUnit {
    public:
        RigidUnit(RigidUnitIndex ix, CompoundModeler& parent) 
            : myIndex(ix), compoundModeler(&parent)
        {}
        
        MassProperties calcMassProperties() const {
            Real mass = 0;
            Vec3 com(0);
            Inertia inertia(0);
            
            std::set<AtomSubsystem::AtomIndex>::const_iterator a;
            for (a = clusterAtoms.begin(); a != clusterAtoms.end(); ++a) 
            {
                const AtomSubsystem::Atom& ssAtom = 
                        compoundModeler->getAtomSubsystem().getAtom(*a);
                
                Real ma = ssAtom.getMass();
                mass += ma;
                com += ma * ssAtom.getStationInBodyFrame();
                inertia += Inertia(ssAtom.getStationInBodyFrame(), ma);
            }
            if (mass > 0) com /= mass;
            
            return MassProperties(mass, com, inertia);
        }

        RigidUnit& setMobilizedBodyIndex(MobilizedBodyIndex ix) 
        {
            assert(!bodyIx.isValid());
            assert(ix.isValid());
            bodyIx = ix;
            std::set<AtomSubsystem::AtomIndex>::iterator a;
            for (a = clusterAtoms.begin(); a != clusterAtoms.end(); ++a) 
                compoundModeler->updAtomSubsystem().setAtomMobilizedBodyIndex(*a, ix);

            return *this;
        }
        
        void buildUp(
                AtomSubsystem::AtomIndex seedAtomIx, 
                CompoundRep& compoundRep);
        
        // TODO: these members should be private with accessors
        RigidUnitIndex parentId; // InvalidId implies parented to Ground
        Compound::BondCenterIndex inboardBondCenterIndex;
        Angle inboardBondDihedralAngle;
        Transform frameInTopCompoundFrame; // useful intermediate computation
        Transform frameInParentFrame; // what we ultimately want
        
        CompoundModeler* compoundModeler;
        
        MobilizedBodyIndex bodyIx; // populated toward the end
        RigidUnitIndex myIndex;

        Compound::BondIndex inboardBondIndex;

        std::set<AtomSubsystem::AtomIndex> clusterAtoms;
    };
    
    
    // AtomBonding data structure represents one Atom in the modelCompounds() method
    class AtomBonding {
    public:
        AtomBonding() {}
        AtomBonding(Compound::AtomIndex compoundIx)
            : compoundAtomIndex(compoundIx) {}

        AtomSubsystem::AtomIndex parentAtomIndex; // for finding rootiest atom in cluster
        RigidUnitIndex clusterIx;

        // Store only those bonds that are part of the multibody tree structure
        std::set<AtomSubsystem::AtomIndex> treeBonds;
        std::set<AtomSubsystem::AtomIndex> freeTreeBonds;

        Vec3 locationInBodyFrame;
        
        Compound::AtomIndex getCompoundAtomIndex() const {
            return compoundAtomIndex;
        }
        
        AtomBonding& setCompoundAtomIndex(Compound::AtomIndex ix) {
            compoundAtomIndex = ix;
            return (*this);
        }
        
    private:
        Compound::AtomIndex compoundAtomIndex;        
    };
    
    
    
    CompoundModeler(AtomSubsystem& atomSubsystem)
            : atomSubsystem(atomSubsystem)
    {}
    
    const AtomSubsystem& getAtomSubsystem() const {
        return atomSubsystem;
    }
    AtomSubsystem& updAtomSubsystem() {
        return atomSubsystem;
    }
    const AtomBonding& getAtomBonding(AtomSubsystem::AtomIndex a) const {
        return atomBondings[a];
    }
    AtomBonding& updAtomBonding(AtomSubsystem::AtomIndex a) {
        return atomBondings[a];
    }
    const RigidUnit& getRigidUnit(RigidUnitIndex a) const {
        return rigidUnits[a];
    }
    RigidUnit& updRigidUnit(RigidUnitIndex a) {
        return rigidUnits[a];
    }
    
    // This is intended to replace CompoundSystem.modelCompounds()
    void model(Compound& compound);
        
private:
    AtomSubsystem& atomSubsystem;
    std::vector<RigidUnit> rigidUnits;
    std::vector<AtomBonding> atomBondings;
    
};
    
} // namespace SimTK

#endif // SimTK_MOLMODEL_COMPOUND_MODELER_H_
