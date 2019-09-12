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
	    
	    dumm.loadAmber99Parameters();
	    
	    // dumm.setGbsaGlobalScaleFactor(0.97);
	    
	    std::ifstream is("/home/cmbruns/svn/SimTK/molmodel/resources/SantaLuciaParams/2MG/test.param");
	    dumm.populateFromTinkerParameterFile(is);
		
	    RNA rna ("AA");
	    // rna.writeDefaultPdb(cout);
	    NaPhosphodiesterLinkage po2("", "", 'X');
	    
		// MagnesiumIon mg1, mg2;
		MagnesiumIon::setAmberLikeParameters(dumm);
		MagnesiumIon mg1, mg2;
		ZincIon::setAmberLikeParameters(dumm);

		// system.adoptCompound(p2, Vec3( 0.5, 0, 0));
		system.adoptCompound(mg1, Vec3(-0.5, 0.6, 0.6));
		// system.adoptCompound(mg2, Vec3(-0.5, -1, 0));
		system.adoptCompound(rna, Vec3(-0.5, 0, 0));
		// system.adoptCompound(po2, Vec3(-0.5, 0, 0));
		
		std::ofstream of("test.pdb");
		system.addEventReporter(new PeriodicPdbWriter(system, of, 0.020));
	    system.addEventHandler(
	    		new VelocityRescalingThermostat(system, 293.15, 0.050));
		
	    // rna.setBondMobility(BondMobility::Free , "0/O3*" ,"1/P");
	    
	    // rna.writeDefaultPdb(cout);
	    
		system.modelCompounds();
		
		cerr << "O3* mobilized body index = " << rna.getAtomMobilizedBodyIndex(rna.getAtomIndex("0/O3*")) << endl;
		
		system.realizeTopology();
		
		system.realize(system.updDefaultState(), Stage::Dynamics);
		Vec3 mg1Force = system.getRigidBodyForces(system.updDefaultState(), Stage::Dynamics)
			[mg1.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))][1];
		cout << "Force on first magnesium = " << mg1Force << endl;
		
		cout << "Mobilized body Id " << mg1.getAtomMobilizedBodyIndex(Compound::AtomIndex(0)) << endl;
		
		// LocalEnergyMinimizer::minimizeEnergy(system, system.updDefaultState(), 15.0);
		
		VerletIntegrator integrator(system);
		TimeStepper timeStepper(system, integrator);
		timeStepper.initialize(system.updDefaultState());
		
		double endTime = 0.0100;
		
		// while ( timeStepper.getTime() < endTime )
		timeStepper.stepTo(endTime);
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

