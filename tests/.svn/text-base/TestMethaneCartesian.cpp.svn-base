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

    Methane methane1, methane2;

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);

    dumm.loadAmber99Parameters();
    dumm.loadTestMoleculeParameters();

    system.adoptCompound(methane1, Vec3(   0, 0, 0));
    system.adoptCompound(methane2, Vec3( 0.4, 0, 0));

    methane1.setBondMobility( BondMobility::Free, Compound::BondIndex(0) );
    methane1.setBondMobility( BondMobility::Free, Compound::BondIndex(1) );
    methane1.setBondMobility( BondMobility::Free, Compound::BondIndex(2) );
    methane1.setBondMobility( BondMobility::Free, Compound::BondIndex(3) );

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

