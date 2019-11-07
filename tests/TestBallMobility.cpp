/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

#include "SimTKmolmodel.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;

using namespace SimTK;
using namespace std;

// define CREATE_VIZ_WINDOW to see animated window of simulation
// undefine for automated nightly builds
#define CREATE_VIZ_WINDOW

#ifndef ANG_360_TO_180
#define ANG_360_TO_180(x) (((x)>180) ? ((x)-360) : (x))
#endif

class HPeroxide : public Compound {
public:
    HPeroxide(DuMMForceFieldSubsystem& dumm)
    {

        // THe way the molecule is built here replicates the way is built
        // in Gmolmodel but at a smaller scale.

        // Set Gmol atoms Molmodel types
        Compound::SingleAtom h1comp = UnivalentAtom("AAAA", Element(1, "Hydrogen", "H", 1));
        h1comp.setDefaultInboardBondLength(0.1112);

        Compound::SingleAtom o3comp = BivalentAtom("AAAB", Element(8, "Oxygen", "O", 14));
        o3comp.setDefaultInboardBondLength(0.19);

        Compound::SingleAtom o2comp = BivalentAtom("AAAC", Element(8, "Oxygen", "O", 14));
        o2comp.setDefaultInboardBondLength(0.19);

        Compound::SingleAtom h4comp = UnivalentAtom("AAAD", Element(1, "Hydrogen", "H", 1));
        h4comp.setDefaultInboardBondLength(0.1112);

        // Add Biotypes
        BiotypeIndex h1BiotypeIndex = Biotype::defineBiotype(Element(1, "H", "H", 1), 1, "TD0", "AAAA", Ordinality::Any);
        BiotypeIndex o3BiotypeIndex = Biotype::defineBiotype(Element(8, "O", "O", 14), 2, "TD0", "AAAB", Ordinality::Any);
        BiotypeIndex o2BiotypeIndex = Biotype::defineBiotype(Element(8, "O", "O", 14), 2, "TD0", "AAAC", Ordinality::Any);
        BiotypeIndex h4BiotypeIndex = Biotype::defineBiotype(Element(1, "H", "H", 1), 1, "TD0", "AAAD", Ordinality::Any);

        // Add AtomClasses
        DuMM::AtomClassIndex h1DAIx = dumm.getNextUnusedAtomClassIndex();
        const char *h1AtomClassName = "topTD0H_0";
        dumm.defineAtomClass(h1DAIx, h1AtomClassName, 1, 1
                , 0.05000000074505806 / 10.0, 0.0);

        DuMM::AtomClassIndex o3DAIx = dumm.getNextUnusedAtomClassIndex();
        const char *o3AtomClassName = "topTD0O_1";
        dumm.defineAtomClass(o3DAIx, o3AtomClassName, 8, 2
                , 0.05000000074505806 / 10.0, 0.0);

        DuMM::AtomClassIndex o2DAIx = dumm.getNextUnusedAtomClassIndex();
        const char *o2AtomClassName = "topTD0O_2";
        dumm.defineAtomClass(o2DAIx, o2AtomClassName, 8, 2
                , 0.05000000074505806 / 10.0, 0.0);

        DuMM::AtomClassIndex h4DAIx = dumm.getNextUnusedAtomClassIndex();
        const char *h4AtomClassName = "topTD0H_3";
        dumm.defineAtomClass(h4DAIx, h4AtomClassName, 1, 1
                , 0.05000000074505806 / 10.0, 0.0);

        // Define ChargedAtomTypes and BiotypeChargedAtomType
        DuMM::ChargedAtomTypeIndex h1ChargedDAIx = dumm.getNextUnusedChargedAtomTypeIndex();
        std::string h1ChargedATypeName = "TD0AAAAH";
        dumm.defineChargedAtomType(h1ChargedDAIx, h1ChargedATypeName.c_str(), h1DAIx, 0);
        dumm.setBiotypeChargedAtomType(h1ChargedDAIx, h1BiotypeIndex);

        DuMM::ChargedAtomTypeIndex o3ChargedDAIx = dumm.getNextUnusedChargedAtomTypeIndex();
        std::string o3ChargedATypeName = "TD0AAABO";
        dumm.defineChargedAtomType(o3ChargedDAIx, o3ChargedATypeName.c_str(), o3DAIx, 0);
        dumm.setBiotypeChargedAtomType(o3ChargedDAIx, o3BiotypeIndex);

        DuMM::ChargedAtomTypeIndex o2ChargedDAIx = dumm.getNextUnusedChargedAtomTypeIndex();
        std::string o2ChargedATypeName = "TD0AAACO";
        dumm.defineChargedAtomType(o2ChargedDAIx, o2ChargedATypeName.c_str(), o2DAIx, 0);
        dumm.setBiotypeChargedAtomType(o2ChargedDAIx, o2BiotypeIndex);

        DuMM::ChargedAtomTypeIndex h4ChargedDAIx = dumm.getNextUnusedChargedAtomTypeIndex();
        std::string h4ChargedATypeName = "TD0AAADH";
        dumm.defineChargedAtomType(h4ChargedDAIx, h4ChargedATypeName.c_str(), h4DAIx, 0);
        dumm.setBiotypeChargedAtomType(h4ChargedDAIx, h4BiotypeIndex);

        // Add bond params
        dumm.defineBondStretch_KA(h1DAIx, o3DAIx, 83.659999999999997, 1.54);
        dumm.defineBondStretch_KA(o3DAIx, o2DAIx, 83.659999999999997, 1.54);
        dumm.defineBondStretch_KA(o2DAIx, h4DAIx, 83.659999999999997, 1.54);

        // Add angle params
        dumm.defineBondBend_KA(h1DAIx, o3DAIx, o2DAIx, 43.460000000000001,
                ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * 1.570797));
        dumm.defineBondBend_KA(o3DAIx, o2DAIx, h4DAIx, 43.460000000000001,
                ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * 1.570797));

        // Add dihedral parameters
        dumm.defineBondTorsion_KA(h1DAIx, o3DAIx, o2DAIx, h4DAIx,
                1, 0.0, ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * 0));

        // Add this atoms to this Compound
        setBaseAtom(o3comp);
        setAtomBiotype("AAAB", "TD0", "AAAB");
        convertInboardBondCenterToOutboard();

        bondAtom(h1comp, "AAAB/bond2", 0.149, 0);
        setAtomBiotype("AAAA", "TD0", "AAAA");

        bondAtom(o2comp, "AAAB/bond1", 0.149, 0);
        setAtomBiotype("AAAC", "TD0", "AAAC");

        bondAtom(h4comp, "AAAC/bond2", 0.149, 0);
        setAtomBiotype("AAAD", "TD0", "AAAD");

        // Set a configuration
        std::map<AtomIndex, Vec3> atomTargets;
        atomTargets.insert(pair<AtomIndex, Vec3> (Compound::AtomIndex(1), Vec3( 0, -1, 0)));
        atomTargets.insert(pair<AtomIndex, Vec3> (Compound::AtomIndex(0), Vec3( 0,  0, 0)));
        atomTargets.insert(pair<AtomIndex, Vec3> (Compound::AtomIndex(2), Vec3( 1,  0, 0)));
        atomTargets.insert(pair<AtomIndex, Vec3> (Compound::AtomIndex(3), Vec3( 1,  1, 0)));

        matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0);

        // Ultimaately what we want to try: Bond Mobility
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            //setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
            setBondMobility(BondMobility::Ball, Compound::BondIndex(r));
        }

    }
};

int main() {

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);

    DuMMForceFieldSubsystem dumm(system);

    HPeroxide hper(dumm);

    //dumm.loadAmber99Parameters();

    system.adoptCompound(hper);
    
    system.modelCompounds();

    State state = system.realizeTopology();

    dumm.dump();

#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system);
#endif

    //RungeKuttaMersonIntegrator integrator(system);
    VerletIntegrator integrator(system);
    integrator.initialize(state);

#ifdef CREATE_VIZ_WINDOW
    display.report(integrator.getState());
#endif

    Real RT = 300 * SimTK_BOLTZMANN_CONSTANT_MD;
    Random::Gaussian random;
    Real timeInterval = 0.02;
    Real time = 0.0;
    for (int i = 0; i < 2; i++){
        time = i * timeInterval;

        // Get detM
        int nu = state.getNU();
        SimTK::Vector V(nu);
        system.realize(state, SimTK::Stage::Position);
        matter.realizeArticulatedBodyInertias(state);

        SimTK::Vector DetV(nu);
        SimTK::Real D0 = 1.0;
        matter.calcDetM(state, V, DetV, &D0);

        // Initialize velocities
        if( !(i % 1)) {
            /*
            SimTK::Vector U = state.updU();
            for (int j = 0; j < state.getNU(); ++j) {
                U[j] = random.getValue();
            }
            */

            double sqrtRT = std::sqrt(RT);
            SimTK::Vector V(nu);
            SimTK::Vector SqrtMInvV(nu);
            for (int i=0; i < nu; ++i){
                V[i] = random.getValue();
            }
            matter.multiplyBySqrtMInv(state, V, SqrtMInvV);

            SqrtMInvV *= (sqrtRT); // Set stddev according to temperature
            system.realize(state, SimTK::Stage::Velocity);

        }

        // Integrate trajectory
        integrator.stepTo(time);

#ifdef CREATE_VIZ_WINDOW
        display.report(integrator.getState());
#endif

    }

}

