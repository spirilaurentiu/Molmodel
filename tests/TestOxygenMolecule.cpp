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

#define ASSERT(x) SimTK_ASSERT_ALWAYS(x, "Assertion failed")

// define CREATE_VIZ_WINDOW to see animated window of simulation
// undefine for automated nightly builds
// #define CREATE_VIZ_WINDOW

class Oxygen2 : public Compound {
public:
    Oxygen2() {
        if (! Biotype::exists("Oxygen2", "O") )
            Biotype::defineBiotype(Element::Oxygen(), 1, "Oxygen2", "O");

        BiotypeIndex biotypeIx = Biotype::get("Oxygen2", "O").getIndex();

        setBaseAtom( UnivalentAtom("O1", Element::Oxygen()) );
        bondAtom( UnivalentAtom("O2", Element::Oxygen()), "O1/bond", 0.13);

        setBiotypeIndex("O1", biotypeIx);
        setBiotypeIndex("O2", biotypeIx);
    }
};

void testWater() {
    Compound water;
    Real angle = 105 * SimTK::Deg2Rad;
    Real length = 0.09;
    water.setBaseAtom(BivalentAtom("O", Element::Oxygen(), angle));
    water.bondAtom(UnivalentAtom("H1", Element::Hydrogen()), "O/bond1", length);
    water.bondAtom(UnivalentAtom("H2", Element::Hydrogen()), "O/bond2", length);
    ASSERT((water.calcDefaultAtomLocationInGroundFrame("O") - Vec3(0,0,0)).norm() < 0.01);
    ASSERT((water.calcDefaultAtomLocationInGroundFrame("H1") - Vec3(length,0,0)).norm() < 0.01);
    ASSERT((water.calcDefaultAtomLocationInGroundFrame("H2") - Vec3(std::cos(angle)*length,std::sin(angle)*length,0)).norm() < 0.01);
    water.writeDefaultPdb(std::cout);
}

int main() {

    testWater();

    Oxygen2 oxygen1;
    oxygen1.writeDefaultPdb(std::cout);

    Oxygen2 oxygen2;

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);

    // ifstream tinkerStream("C:/cygwin/home/cmbruns/svn/molmodel/resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();

    DuMM::ChargedAtomTypeIndex atomTypeId(4000);
    DuMM::AtomClassIndex atomClassId(24); // carbonyl oxygen
    dumm.defineChargedAtomType(atomTypeId, "Oxygen2 O",   atomClassId,  0.00);
    BiotypeIndex biotypeIx = Biotype::get("Oxygen2", "O").getIndex();
    dumm.setBiotypeChargedAtomType(atomTypeId, biotypeIx);
    dumm.defineBondStretch_KA(atomClassId, atomClassId, 570, 1.3);

    system.adoptCompound(oxygen1);
    system.adoptCompound(oxygen2, Vec3(0.5, 0, 0));

    Methane methane;
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), "Methane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), "Ethane C",   DuMM::AtomClassIndex(1),  0.03);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), "Ethane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), Biotype::MethaneH().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), Biotype::EthaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), Biotype::EthaneH().getIndex());
    system.adoptCompound(methane, Vec3(0, 0.5, 0));

    system.modelCompounds();        

    State state = system.realizeTopology();

    system.realize(state, Stage::Position);
    oxygen1.writePdb(state, std::cout);

#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system, 0.1);
#endif

    RungeKuttaMersonIntegrator study(system);
    study.initialize(state);


#ifdef CREATE_VIZ_WINDOW
    display.report(study.getState());
#endif

    Real timeInterval = 0.05;
    for (Real time=0.0; time < (10 * timeInterval); time += timeInterval) // picoseconds
    {
        study.stepTo(time);

#ifdef CREATE_VIZ_WINDOW
        display.report(study.getState());
#endif

    }

}

