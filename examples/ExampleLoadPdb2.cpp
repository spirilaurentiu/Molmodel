#include "Molmodel.h"

#include <iostream>
#include <exception>
#include <fstream>

using namespace SimTK;

int main() {
try {
    // Load the PDB file and construct the system.
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

	std::ifstream pdbFile("1AKGtrim.pdb");
	Protein protein(pdbFile);

    protein.assignBiotypes();
	system.adoptCompound(protein);
    system.modelCompounds();

    system.addEventHandler(new VelocityRescalingThermostat(
        system, 293.15, 0.1));

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.025) );

    system.realizeTopology();
    
    // Create an initial state for the simulation.
    
    State state = system.getDefaultState();
    LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);
    
    // Simulate it.
    
    VerletIntegrator integ(system);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);

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
