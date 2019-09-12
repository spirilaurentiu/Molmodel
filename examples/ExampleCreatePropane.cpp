#include "Molmodel.h"

#include <iostream>
#include <fstream>

using namespace SimTK;

// Propane is a three carbon linear alkane
// C(3)H(8), or CH3-CH2-CH3
class Propane : public Molecule 
{
public:
    // constructor
    Propane()
    {
        setCompoundName("Propane");

        // First atom
        setBaseAtom( AliphaticCarbon("C1") );
		setAtomBiotype("C1", "propane", "C1_or_C3");
        convertInboardBondCenterToOutboard(); // this is the root of the parent compound

        // Second atom
        bondAtom( AliphaticCarbon("C2"), "C1/bond1" );
		setAtomBiotype("C2", "propane", "C2");

        // Third atom
        bondAtom( AliphaticCarbon("C3"), "C2/bond2" );
		setAtomBiotype("C3", "propane", "C1_or_C3");

        // First methyl hydrogens
        bondAtom( AliphaticHydrogen("H11"), "C1/bond2" );
        bondAtom( AliphaticHydrogen("H12"), "C1/bond3" );
        bondAtom( AliphaticHydrogen("H13"), "C1/bond4" );
		setAtomBiotype("H11", "propane", "H1_or_H3");
		setAtomBiotype("H12", "propane", "H1_or_H3");
		setAtomBiotype("H13", "propane", "H1_or_H3");

        // Second methylene hydrogens
        bondAtom( AliphaticHydrogen("H21"), "C2/bond3" );
        bondAtom( AliphaticHydrogen("H22"), "C2/bond4" );
		setAtomBiotype("H21", "propane", "H2");
		setAtomBiotype("H22", "propane", "H2");

        // Third methyl hydrogens
        bondAtom( AliphaticHydrogen("H31"), "C3/bond2" );
        bondAtom( AliphaticHydrogen("H32"), "C3/bond3" );
        bondAtom( AliphaticHydrogen("H33"), "C3/bond4" );
		setAtomBiotype("H31", "propane", "H1_or_H3");
		setAtomBiotype("H32", "propane", "H1_or_H3");
		setAtomBiotype("H33", "propane", "H1_or_H3");
    }

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
	// Must be called AFTER first propane is declared, so Biotypes and atom classes will be defined
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& forceField)
    {
        forceField.defineChargedAtomType(
            forceField.getNextUnusedChargedAtomTypeIndex(),
            "propane C1_or_C3", 
			forceField.getAtomClassIndex("CT"), // amber tetrahedral carbon
            -0.060 // made up
            );
        forceField.setBiotypeChargedAtomType( 
			forceField.getChargedAtomTypeIndex("propane C1_or_C3"),
			Biotype::get("propane", "C1_or_C3").getIndex() );

        forceField.defineChargedAtomType(
            forceField.getNextUnusedChargedAtomTypeIndex(),
            "propane C2", 
			forceField.getAtomClassIndex("CT"), // amber tetrahedral carbon
            -0.040 // made up
            );
        forceField.setBiotypeChargedAtomType( 
			forceField.getChargedAtomTypeIndex("propane C2"),
			Biotype::get("propane", "C2").getIndex() );

        forceField.defineChargedAtomType(
            forceField.getNextUnusedChargedAtomTypeIndex(),
            "propane H1_or_H3", 
			forceField.getAtomClassIndex("HC"), // amber tetrahedral carbon
            0.020 // made up, use net neutral charge
            );
        forceField.setBiotypeChargedAtomType( 
			forceField.getChargedAtomTypeIndex("propane H1_or_H3"),
			Biotype::get("propane", "H1_or_H3").getIndex() );

        forceField.defineChargedAtomType(
            forceField.getNextUnusedChargedAtomTypeIndex(),
            "propane H2", 
            forceField.getAtomClassIndex("HC"), // amber tetrahedral carbon
            0.020 // made up, use net neutral charge
            );
        forceField.setBiotypeChargedAtomType( 
			forceField.getChargedAtomTypeIndex("propane H2"),
			Biotype::get("propane", "H2").getIndex() );
   }
};

int main() { 
try {
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system); 

    // Atom classes are available, but not charged atom types for propane
    // in standard Amber force field
    forceField.loadAmber99Parameters();

    Propane propane1, propane2;

    Propane::setAmberLikeParameters(forceField);

    // place first propane, units are nanometers
    // skew it a little to break strict symmetry
    system.adoptCompound(propane1,   Transform(Vec3(-0.5, 0, 0)) 
                                   * Transform(Rotation(0.1, YAxis)) ); 

    // place second propane, units are nanometers
    system.adoptCompound(propane2, Vec3( 0.5, 0, 0)); 

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.100) );

    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology();

    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(100.0);

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

