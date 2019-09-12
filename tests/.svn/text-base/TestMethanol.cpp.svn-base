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
// #define CREATE_VIZ_WINDOW

class Methanol : public Compound {
public:
    Methanol() 
    {
        setBaseCompound("methyl", MethylGroup() );
        bondCompound("hydroxyl", AlcoholOHGroup(), "methyl/bond", 0.141);

        nameAtom("C", "methyl/C");
        nameAtom("1HC", "methyl/H1");
        nameAtom("2HC", "methyl/H2");
        nameAtom("3HC", "methyl/H3");
        nameAtom("O", "hydroxyl/O");
        nameAtom("HO", "hydroxyl/H");


        if (! Biotype::exists("Methanol", "C") )
            Biotype::defineBiotype(Element::Carbon(), 4, "Methanol", "C");
        if (! Biotype::exists("Methanol", "HC") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "Methanol", "HC");
        if (! Biotype::exists("Methanol", "O") )
            Biotype::defineBiotype(Element::Oxygen(), 2, "Methanol", "O");
        if (! Biotype::exists("Methanol", "HO") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "Methanol", "HO");

        setBiotypeIndex( "C", Biotype::get("Methanol", "C").getIndex() );
        setBiotypeIndex( "1HC", Biotype::get("Methanol", "HC").getIndex() );
        setBiotypeIndex( "2HC", Biotype::get("Methanol", "HC").getIndex() );
        setBiotypeIndex( "3HC", Biotype::get("Methanol", "HC").getIndex() );
        setBiotypeIndex( "O", Biotype::get("Methanol", "O").getIndex() );
        setBiotypeIndex( "HO", Biotype::get("Methanol", "HO").getIndex() );
    }
};

// Rubber
static const Real rubber_density = 1100.;  // kg/m^3
static const Real rubber_young   = 0.01e9; // pascals (N/m)
static const Real rubber_poisson = 0.5;    // ratio
static const Real rubber_planestrain =
                    rubber_young/(1.-rubber_poisson*rubber_poisson);
static const Real rubber_dissipation = 0.005;

int main() {

    Methanol methanol1, methanol2, methanol3;

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);

    DuMMForceFieldSubsystem dumm(system);
    HuntCrossleyContact     contact(system);

    contact.addHalfSpace(matter.Ground(), UnitVec3(0,1,0), -3, rubber_planestrain, rubber_dissipation);
    decorations.addBodyFixedDecoration(GroundIndex, Transform(Vec3(0, -3, 0)), DecorativeBrick(Vec3(7,.02,7)).setColor(Yellow).setOpacity(0.25));

    contact.addHalfSpace(matter.Ground(), UnitVec3(1,0,0), -3, rubber_planestrain, rubber_dissipation);
    decorations.addBodyFixedDecoration(GroundIndex, Transform(Vec3(-3, 0, 0)), DecorativeBrick(Vec3(.02, 7 ,7)).setColor(Yellow).setOpacity(0.25));

    contact.addHalfSpace(matter.Ground(), UnitVec3(1,0,0), 3, rubber_planestrain, rubber_dissipation);
    decorations.addBodyFixedDecoration(GroundIndex, Transform(Vec3(3, 0, 0)), DecorativeBrick(Vec3(.02,7,7)).setColor(Yellow).setOpacity(0.25));

    // ifstream tinkerStream("C:/cygwin/home/cmbruns/svn/molmodel/resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();

    // hydroxyl hydrogen
    DuMM::AtomClassIndex amberHOAtomClassIndex(31);
    DuMM::ChargedAtomTypeIndex methanolHOAtomTypeIndex(4003);
    dumm.defineChargedAtomType(methanolHOAtomTypeIndex, "Methanol OH", amberHOAtomClassIndex, 0.4); // from serine OH
    dumm.setBiotypeChargedAtomType( methanolHOAtomTypeIndex, Biotype::get("Methanol", "HO").getIndex() );

    // oxygen
    DuMM::AtomClassIndex amberOAtomClassIndex(22);
    DuMM::ChargedAtomTypeIndex methanolOAtomTypeIndex(4002);
    dumm.defineChargedAtomType(methanolOAtomTypeIndex, "Methanol O", amberOAtomClassIndex, -0.6); // from serine OH
    dumm.setBiotypeChargedAtomType( methanolOAtomTypeIndex, Biotype::get("Methanol", "O").getIndex() );

    // methyl hydrogen atoms
    DuMM::AtomClassIndex amberHCAtomClassIndex(35);
    DuMM::ChargedAtomTypeIndex methanolHCAtomTypeIndex(4001);
    dumm.defineChargedAtomType(methanolHCAtomTypeIndex, "Methanol HC", amberHCAtomClassIndex, 0.025);
    dumm.setBiotypeChargedAtomType( methanolHCAtomTypeIndex, Biotype::get("Methanol", "HC").getIndex() );

    // carbon atom
    DuMM::AtomClassIndex amberCTAtomClassIndex(1);
    DuMM::ChargedAtomTypeIndex methanolCAtomTypeIndex(4000);
    dumm.defineChargedAtomType(methanolCAtomTypeIndex, "Methanol C", amberCTAtomClassIndex, 0.125);
    dumm.setBiotypeChargedAtomType( methanolCAtomTypeIndex, Biotype::get("Methanol", "C").getIndex() );

    system.adoptCompound(methanol1, Vec3(-0.5, 0, 0));
    system.adoptCompound(methanol2, Vec3(   0, 0, 0));
    system.adoptCompound(methanol3, Vec3( 0.5, 0, 0));

    system.modelCompounds();        

    State state = system.realizeTopology();

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

