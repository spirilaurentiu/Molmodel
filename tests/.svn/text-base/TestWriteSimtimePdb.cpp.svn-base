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

    Methane methane;
    Ethane ethane;
    Protein protein("ECG");

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);

    // ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();

    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), "Methane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), "Ethane C",   DuMM::AtomClassIndex(1),  0.03);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), "Ethane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), Biotype::MethaneH().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), Biotype::EthaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), Biotype::EthaneH().getIndex());
    protein.assignBiotypes();

    system.adoptCompound(methane, Vec3(-0.5, 0, 0));
    system.adoptCompound(ethane, Vec3(   0, 0, 0));
    system.adoptCompound(protein, Vec3( 0.5, 0, 0));

    system.modelCompounds();        

    State state = system.realizeTopology();
#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system);
#endif


    RungeKuttaMersonIntegrator study(system);
    study.initialize(state);

    methane.writeDefaultPdb(cout);

#ifdef CREATE_VIZ_WINDOW
    display.report(study.getState());
#endif

    for (Real time=0.0; time < 0.2; time += 0.0500) { // picoseconds
        study.stepTo(time);

        // writePdb() is what this test case means to exercise
        methane.writePdb(study.getState(), cout);

#ifdef CREATE_VIZ_WINDOW
        display.report(study.getState());
#endif

    }

}

