#include "Molmodel.h"

#include <iostream>
#include <exception>
#include <fstream>

using namespace SimTK;

int main() { 
try {
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system); 

    // Atom classes are available, but not charged atom types for ethane
    // in standard Amber force field
    forceField.loadAmber99Parameters();

    if (! Biotype::exists("ethane", "C"))
        Biotype::defineBiotype(Element::Carbon(), 4, "ethane", "C");
    if (! Biotype::exists("ethane", "H"))
        Biotype::defineBiotype(Element::Hydrogen(), 1, "ethane", "H");

    forceField.defineChargedAtomType(
        forceField.getNextUnusedChargedAtomTypeIndex(), 
        "ethane C", 
		forceField.getAtomClassIndex("CT"), 
        -0.060 // made up
        );
    forceField.setBiotypeChargedAtomType(
		forceField.getChargedAtomTypeIndex("ethane C"),
		Biotype::get("ethane", "C").getIndex() );

    forceField.defineChargedAtomType(
        forceField.getNextUnusedChargedAtomTypeIndex(), 
        "ethane H", 
		forceField.getAtomClassIndex("HC"), 
        0.020 // made up, use net neutral charge
        );
    forceField.setBiotypeChargedAtomType( 
		forceField.getChargedAtomTypeIndex("ethane H"),
		Biotype::get("ethane", "H").getIndex() );

    Ethane ethane1, ethane2;
    ethane1.writeDefaultPdb(std::cout);

    // place first ethane, units are nanometers
    // skew it a little to break strict symmetry
    system.adoptCompound(ethane1,   Transform(Vec3(-0.5, 0, 0)) 
                                  * Transform(Rotation(0.1, YAxis)) ); 

    // place second ethane, units are nanometers
    system.adoptCompound(ethane2, Vec3( 0.5, 0, 0)); 

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.050) );

    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology();

    // Simulate it.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(200.0);

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

