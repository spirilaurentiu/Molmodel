/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

using namespace SimTK;
using namespace std;

int main() {
try {

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem dumm(system);
    dumm.loadAmber99Parameters();

    dumm.setGbsaGlobalScaleFactor(1.0);
    dumm.setAmberImproperTorsionGlobalScaleFactor(0.0);
    dumm.setBondBendGlobalScaleFactor(0.0);
    dumm.setBondStretchGlobalScaleFactor(0.0);
    dumm.setBondTorsionGlobalScaleFactor(0.0);
    dumm.setCoulombGlobalScaleFactor(0.0);
    dumm.setVdwGlobalScaleFactor(0.0);

    dumm.setGbsaIncludeAceApproximationOff();

    // A lone serine with end caps is used for comparison with Mark Friedrichs' Tinker run
    Protein serine("S");
    serine.writeDefaultPdb(cout);
    serine.assignBiotypes();
    system.adoptCompound(serine);

    system.modelCompounds();
    State state = system.realizeTopology();

    RungeKuttaMersonIntegrator study(system);
    study.initialize(state);

    Real timeInterval = 0.0020; // picoseconds
    for (Real time=0.0; time < 10 * timeInterval; time += timeInterval) { // picoseconds
        study.stepTo(time);

        // serine.writePdb(study.getState(), cout);

        system.realize(study.getState(), Stage::Dynamics);
        Real energy = system.calcPotentialEnergy(study.getState());
        cout << "Potential energy = " << energy << "KJ/mol (" << energy * DuMM::KJ2Kcal << " kcal/mol)" << endl;
    }

    return 0;
}

catch (const std::exception& e)
{
    printf("EXCEPTION THROWN: %s\n", e.what());
    return 1;
}
catch (...)
{
    printf("UNKNOWN EXCEPTION THROWN\n");
}    return 1;
}
