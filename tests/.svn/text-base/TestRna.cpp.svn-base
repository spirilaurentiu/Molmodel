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

int main() 
{
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem dumm(system);

    dumm.loadAmber99Parameters();

    RNA rna("GCAU");
    rna.assignBiotypes();
    system.adoptCompound(rna);
    rna.writeDefaultPdb(cout);

    RibonucleotideResidue::Cytidylate cytidylate;
    cytidylate.assignBiotypes();
    // system.adoptCompound(cytidylate);

    Real desiredAngle = 200*Deg2Rad;
    cout << "desired chi = " << DuMM::Rad2Deg * desiredAngle << endl;
    cytidylate.setDefaultDihedralAngle("chi", desiredAngle);
    cout << "default chi = " << DuMM::Rad2Deg * cytidylate.calcDefaultDihedralAngle("chi") << endl;

    system.modelCompounds();  

    State state = system.realizeTopology();

#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system);
#endif

    RungeKuttaMersonIntegrator study(system);

    study.initialize(state);

#ifdef CREATE_VIZ_WINDOW
    display.report(study.getState());
#endif

    Real timeInterval = 0.0020; // picoseconds
    for (Real time=0.0; time < 10 * timeInterval; time += timeInterval) { // picoseconds
        study.stepTo(time);

        rna.writePdb(study.getState(), cout);

#ifdef CREATE_VIZ_WINDOW
        display.report(study.getState());
#endif

        // cout << "state chi1 = " << DuMM::Rad2Deg * cytidylate.getDihedralAngle(study.getState(), "chi1") << endl;
        // cout << "default chi1 = " << DuMM::Rad2Deg * cytidylate.getDefaultDihedralAngle("chi1") << endl;
    }

}

