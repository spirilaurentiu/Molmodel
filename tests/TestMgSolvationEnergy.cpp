/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
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
#include <fstream>

using namespace SimTK;
using namespace std;

namespace SimTK {

}

int main() {
try
  { 
		
	    CompoundSystem system; // molecule-specialized simbody System
        SimbodyMatterSubsystem matter(system);
	    DuMMForceFieldSubsystem dumm(system); // molecular force field
	    
		// MagnesiumIon mg1, mg2;
		MagnesiumIon::setAmberLikeParameters(dumm);
		MagnesiumIon mg1;

		system.adoptCompound(mg1, Vec3(-0.5, 1, 0));
				
		system.modelCompounds();
		const State& state = system.realizeTopology();
		
		system.realize(state, Stage::Dynamics);
		Vec3 mg1Force = system.getRigidBodyForces(state, Stage::Dynamics)
			[mg1.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))][1];
		cout << "Force on first magnesium = " << mg1Force << endl;
		
		cout << "Mobilized body Id " << mg1.getAtomMobilizedBodyIndex(Compound::AtomIndex(0)) << endl;

		Real JouleToCalorie = 1.0 / 4.184;
		Real potentialEnergy = system.calcPotentialEnergy(state) * JouleToCalorie;
		cout << "Energy = " << potentialEnergy << " kcal/mole" << endl;
		
		// Experimental value is -455.5 kcal/mole, from Rosseinsky (1965) Chem. Rev. 65: 467-490
		// accuracy is greatly affected by estimate of H+ solvation free energy
		Real expectedEnergy = -455.5;
		Real deltaV = potentialEnergy - expectedEnergy;
		
		// ensure potential is within 20 kcals of experimental value, or 400 kcal^2 of squared value
		SimTK_ASSERT_ALWAYS(deltaV*deltaV < 400, (String("Expected Mg+2 solvation energy of ") + String(expectedEnergy) + ", found energy of " + String(potentialEnergy)).c_str() );
  }
catch (const std::exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
  }
catch (...)
  {
    printf("UNKNOWN EXCEPTION THROWN\n");
  }    return 0;
}

