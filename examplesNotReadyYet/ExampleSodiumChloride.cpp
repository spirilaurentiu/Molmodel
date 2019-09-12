#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h" // for vtk visualization

#include <iostream>
#include <exception>

using namespace SimTK;
using namespace std;

int main() { try 
{
    // molecule-specialized simbody System
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
 
    // molecular force field
    TinkerDuMMForceFieldSubsystem dumm(system); 

	dumm.loadAmber99Parameters();

	SodiumIon sodium1, sodium2, sodium3;
	ChlorideIon chloride1, chloride2, chloride3;
	MagnesiumIon magnesium1, magnesium2, magnesium3;

	// dumm.setGbsaGlobalScaleFactor(0);
	// dumm.setCoulombGlobalScaleFactor(0);

	LithiumIon::setAmberLikeParameters(dumm);
	SodiumIon::setAmberLikeParameters(dumm);
	PotassiumIon::setAmberLikeParameters(dumm);
	MagnesiumIon::setAmberLikeParameters(dumm);
	ChlorideIon::setAmberLikeParameters(dumm);

    // system.adoptCompound(magnesium1, Vec3(-0.3, 0, 0)); 
    // system.adoptCompound(magnesium2, Vec3( 0.3, 0, 0)); 
    system.adoptCompound(sodium1, Vec3(-0.3, 0, 0.3)); 
    // system.adoptCompound(sodium2, Vec3( 0.3, 0, 0.3)); 
    system.adoptCompound(chloride1, Vec3(0, -0.3, 0)); 
    // system.adoptCompound(chloride2, Vec3(0,  0.3, 0)); 
    // system.adoptCompound(chloride3, Vec3(0,  0.3, 0.3)); 

    system.updDefaultSubsystem().addEventReporter(new VTKEventReporter(system,
        0.010));

    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology(); 

    // Simulate it.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(500.0);

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
