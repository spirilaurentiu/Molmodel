#include "Molmodel.h"

#include <iostream>
#include <exception>

using namespace SimTK;

int main() { 
try {
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);

    forceField.loadAmber99Parameters();

    Protein protein("SIMTK");
    protein.assignBiotypes();
    system.adoptCompound(protein);

    for ( Compound::BondIndex bondIx(0); 
          bondIx < protein.getNumBonds(); 
          ++bondIx)
    {
        // set all bonds rigid
        protein.setBondMobility(
            BondMobility::Rigid, 
            bondIx);
    }

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.020) );

    // finalize multibody system
    system.modelCompounds(); 

	// Maintain a constant temperature
	system.addEventHandler(new VelocityRescalingThermostat(
		   system,  293.15, 0.1));

	// Instantiate simbody model
	system.realizeTopology();
	State state = system.getDefaultState();

	// Relax the structure before dynamics run
	LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

    // Simulate it.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(20.0);

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

