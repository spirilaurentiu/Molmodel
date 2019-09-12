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

int main() {
    Biotype::initializePopularBiotypes();

    AminoAcidResidue::Alanine alanine1, alanine2;
    alanine1.assignBiotypes();
    alanine2.assignBiotypes();

    AminoAcidResidue::Serine serine1, serine2;
    serine1.assignBiotypes();
    serine2.assignBiotypes();

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    // DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem dumm(system);
    GeneralForceSubsystem    forces(system);

    // ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();


    // system.adoptCompound(alanine1);
    // system.adoptCompound(alanine2);

    // system.adoptCompound(serine1);
    // system.adoptCompound(serine2);


    DuMM::ChargedAtomTypeIndex ethaneCType(4000);
    dumm.defineChargedAtomType(ethaneCType, "ethane C", DuMM::AtomClassIndex(1), -0.075);
    dumm.setBiotypeChargedAtomType(ethaneCType, Biotype::get("ethane", "C").getIndex());

    DuMM::ChargedAtomTypeIndex ethaneHType(4001);
    dumm.defineChargedAtomType(ethaneHType, "ethane C", DuMM::AtomClassIndex(35), 0.025);
    dumm.setBiotypeChargedAtomType(ethaneHType, Biotype::get("ethane", "H").getIndex());

    Ethane ethane1, ethane2;
    // system.adoptCompound(ethane1);
    // system.adoptCompound(ethane2);

    // Protein protein("ACDEFGHIKLMNPQRSTVWY");
    Protein protein("AAAAAAAAAA");
    // Protein protein("APF");
    protein.assignBiotypes();
    system.adoptCompound( protein, Vec3(1.5, 0, 0) );

    DuMM::AtomClassIndex argonClass(200);
    dumm.defineAtomClass_KA(argonClass, "argon", 18, 0, 1.88,  0.2);
    DuMM::ChargedAtomTypeIndex argonType(4010);
    dumm.defineChargedAtomType(argonType, "argon", argonClass, 0.0);
    dumm.setBiotypeChargedAtomType(argonType, Biotype::get("argon", "argon").getIndex());
    Argon argon1, argon2;
    // system.adoptCompound( argon1 );
    // system.adoptCompound( argon2, Vec3(0.7, 0, 0) );

    system.modelCompounds();  

    Force::GlobalDamper(forces, matter, 0.1);

    State state = system.realizeTopology();

#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system, 0.1);
#endif

    RungeKuttaMersonIntegrator study(system);
    // CPodesIntegrator study(system, CPodes::BDF, CPodes::Newton);
    study.setAccuracy(1e-2);
    state.updU() = 1.;
    state.updU()(0,6) = 0; // no rigid body motion
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

