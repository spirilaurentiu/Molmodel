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

int main() {

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    // DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem dumm(system);

    // ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();

    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), "Ethane C",   DuMM::AtomClassIndex(1),  0.03);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), "Ethane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), Biotype::EthaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), Biotype::EthaneH().getIndex());

    cout << "Default coordinates:" << endl;
    std::vector<Ethane> ethanes;
    for (int a = 0; a <= 13; a += 13) {
        Ethane ethane;

        cout << a << " degrees:" << endl;
        ethane.setDefaultTorsionAngle(a * Deg2Rad);
        ethane.writeDefaultPdb(cout);
        cout << "END" << endl;

        ethanes.push_back(ethane);
    }

    std::vector<Ethane>::iterator ethaneI;
    for (ethaneI = ethanes.begin(); ethaneI != ethanes.end();  ++ethaneI)
    {
        system.adoptCompound(*ethaneI);
    }

    system.modelCompounds();        

    State state = system.realizeTopology();

    RungeKuttaMersonIntegrator study(system);
    study.initialize(state);

    cout << "Simbody coordinates:" << endl;
    for (ethaneI = ethanes.begin(); ethaneI != ethanes.end();  ++ethaneI)
    {
        MobilizedBodyIndex bodyId = ethaneI->getAtomMobilizedBodyIndex(ethaneI->getAtomIndex("C2"));
        const MobilizedBody::Pin& pin = (MobilizedBody::Pin&) matter.getMobilizedBody(bodyId);
        Real bodyQ = pin.getAngle(state) / Deg2Rad;

        cout << "Measured torsion angle = " << ethaneI->calcDihedralAngle(study.getState(), "torsion") / Deg2Rad << " degrees" << endl;
        cout << "Mobilizer angle = " << bodyQ << " degrees" << endl;
        ethaneI->writePdb(study.getState(), cout);
        cout << "END" << endl;
    }

}

