#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h" // for vtk visualization

#include <iostream>
#include <fstream>
#include <exception>

using namespace SimTK;
using namespace std;

int main() { try 
{

    CompoundSystem system; // molecule-specialized simbody System
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    TinkerDuMMForceFieldSubsystem dumm(system); // molecular force field

    dumm.loadAmber99Parameters();

	ifstream rnaFile("1L2X.pdb");
	PdbStructure rnaStructure(rnaFile);
	RNA rna( rnaStructure );
    rna.assignBiotypes();
    system.adoptCompound(rna);

    system.updDefaultSubsystem().addEventReporter(
        new VTKEventReporter(system, 0.020));

    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology(); 

    // Simulate it.
    
    VerletIntegrator integ(system);
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

