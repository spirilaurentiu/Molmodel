/* -------------------------------------------------------------------------- *
 *       SimTK Molmodel installation test program (without Visualizer)        *
 * -------------------------------------------------------------------------- *
 * Run this program to verify that Molmodel has been installed properly. This *
 * requires that Simbody is also installed. However, we do not use the        *
 * Visualizer here in case you don't have that available. To test for         *
 * complete functionality run MolmodelInstallTest instead.                    *
 *                                                                            *
 * Authors: Michael Sherman                                                   *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>
using namespace SimTK;

int main() {
try {
    std::clog << 
        "MolmodelInstallTestNoViz: you should see a pdb file of a small protein\n"
        " sent to stdout.\n";

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

    // Generate a pdb frame every 100fs.
    system.addEventReporter(new PeriodicPdbWriter(system, std::cout, 0.100));

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
    ts.stepTo(1.0); // 1ps (10 frame)

    std::clog << "MolmodelInstallTestNoViz: done.\n";

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

