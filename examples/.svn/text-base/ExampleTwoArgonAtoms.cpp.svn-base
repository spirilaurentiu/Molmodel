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

    // GeneralForceSubsystem forces(system);
    // VanderWallSphere boundary(forces, forceField, Vec3(0,0,0), 0.50, 0.2, 0.001);
    
    // Define an atom class for argon
    forceField.defineAtomClass_KA( 
        forceField.getNextUnusedAtomClassIndex(), 
        "argon", 
        18, 
        0, 
        1.88, 
        0.0003832
        );
    forceField.defineChargedAtomType(
        forceField.getNextUnusedChargedAtomTypeIndex(), 
        "argon", 
        forceField.getAtomClassIndex("argon"), 
        0.0 // neutral charge
        );

    if (! Biotype::exists("argon", "argon"))
        Biotype::defineBiotype(Element::Argon(), 0, "argon", "argon");

    forceField.setBiotypeChargedAtomType( forceField.getChargedAtomTypeIndex("argon"), Biotype::get("argon", "argon").getIndex() );
    forceField.setGbsaGlobalScaleFactor(0);

    Argon argonAtom1, argonAtom2; // two argon atoms

    // place first argon atom, units are nanometers
    system.adoptCompound(argonAtom1, Vec3(-0.3, 0, 0)); 

    // place second argon atom, units are nanometers
    system.adoptCompound(argonAtom2, Vec3( 0.3, 0, 0)); 

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.50) );
    
    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology(); 

    // Simulate it.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(2000.0);

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
