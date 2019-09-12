/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Simple Protein                     *
 * -------------------------------------------------------------------------- *
 * This is the first example from the Molmodel User's Guide. It creates a     *
 * small protein (a five-residue peptide), simulates it and generates a live  *
 * animation while it is running.                                             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>
using namespace SimTK;

int main() {
try {
    // molecule-specialized simbody System
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);

    // molecular force field
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    // Create a five-residue protein. Neutral end caps are automatically
    // added unless you suppress them. This is 
    // (Ace-) Ser-Ile-Met-Thr-Lys (-Nac).
    Protein protein("SIMTK");

    protein.assignBiotypes();
    system.adoptCompound(protein);

    // finalize mapping of atoms to bodies
    system.modelCompounds(); 

    // show me a movie
    system.addEventReporter(new Visualizer::Reporter(system, 0.020));

    // Maintain a constant temperature. (This isn't a very good
    // thermostat -- consider NoseHooverThermostat instead.)
    system.addEventHandler(new VelocityRescalingThermostat(
		   system,  293.15, 0.1));

    // Instantiate simbody model and get default state
    State state = system.realizeTopology();

    // Relax the structure before dynamics run
    LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

    // Simulate it.
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(20.0); // 20ps

    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" 
              << std::endl;
    return 1;
}

}


