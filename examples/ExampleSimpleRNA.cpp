/* -------------------------------------------------------------------------- *
 *                   SimTK Molmodel Example: Simple RNA                       *
 * -------------------------------------------------------------------------- *
 * This is the second example from the Molmodel User's Guide. It creates a    *
 * small RNA molecule, simulates it, generates a live animation while it is   *
 * running and generates a multi-frame pdb file to look at later.             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"

#include <iostream>
#include <fstream>
#include <exception>

using namespace SimTK;

// We'll create a file of this name in the current working directory.
const char* PdbFileName = "simpleRNA.pdb";

int main() {
try {
	// molecule-specialized Simbody System
    CompoundSystem         system; 
	SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem  forces(system); // for Nose-Hoover
    DecorationSubsystem    decorations(system);

	// molecular force field
    DuMMForceFieldSubsystem forceField(system); 
    forceField.loadAmber99Parameters();

    // This is adenine-uracil-guanine.
    RNA rna("AUG");
    rna.assignBiotypes();
    system.adoptCompound(rna);

    // Use a Nose-Hoover thermostat, which is a continuous force element
    // rather than a periodic discontinuous "disturbance" like the 
    // VelocityRescalingThermostat. We're using the default relaxation time.
    NoseHooverThermostat thermo(forces, matter, 293.15);

    // Show me a movie with a frame every 20fs.
    system.addEventReporter(new Visualizer::Reporter(system, 0.020) );

    // Generate a pdb file with a frame every 100fs.
    std::ofstream pdbFile; pdbFile.open(PdbFileName);
    system.addEventReporter(new PeriodicPdbWriter(system, pdbFile, 0.100));

	// finalize multibody system
    system.modelCompounds(); 

	// Instantiate simbody model
	system.realizeTopology();
	State state = system.getDefaultState();

	// Relax the structure before dynamics run
	LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

    // Simulate it.
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0); // 10 ps

    pdbFile.close();
    std::cout << "\nWrote pdb file "
              << Pathname::getCurrentWorkingDirectory() << PdbFileName << "\n";
    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" << std::endl;
    return 1;
}

}

