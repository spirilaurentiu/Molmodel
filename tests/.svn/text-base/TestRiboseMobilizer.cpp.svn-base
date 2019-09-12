#include "SimTKmolmodel.h"

#include "molmodel/internal/VanderWallSphere.h"
#include "molmodel/internal/PeriodicVmdReporter.h"

#include <iostream>
#include <exception>
#include <fstream>

using namespace SimTK;
using namespace std;

class PeriodicTorsionMeasurer : public PeriodicEventReporter {
public:
	PeriodicTorsionMeasurer(
	    std::vector<Angle>& torsions,
		const Compound& compound,
		const char*  atom1Name,
		const char*  atom2Name,
		const char*  atom3Name,
		const char*  atom4Name,
        Real interval) 
        : PeriodicEventReporter(interval), 
          torsions(torsions),
          compound(compound),
          atom1Index( compound.getAtomIndex(atom1Name) ),
          atom2Index( compound.getAtomIndex(atom2Name) ),
          atom3Index( compound.getAtomIndex(atom3Name) ),
          atom4Index( compound.getAtomIndex(atom4Name) )
    {	
    }

    void handleEvent(const State& state) const 
    {
        Vec3 p1 = compound.calcAtomLocationInGroundFrame(state, atom1Index);
        Vec3 p2 = compound.calcAtomLocationInGroundFrame(state, atom2Index);
        Vec3 p3 = compound.calcAtomLocationInGroundFrame(state, atom3Index);
        Vec3 p4 = compound.calcAtomLocationInGroundFrame(state, atom4Index);
        
        torsions.push_back( calcDihedralAngle(p1, p2, p3, p4) );
    }
    

private:
    const Compound& compound;
    Compound::AtomIndex atom1Index;
    Compound::AtomIndex atom2Index;
    Compound::AtomIndex atom3Index;
    Compound::AtomIndex atom4Index;
    
    std::vector<Angle>& torsions;
};


std::ostream& printHistogram(const std::vector<Angle>& torsions, std::ostream& os) 
{
    int nBins = 30;
    std::vector<int> bin(nBins);
    
    // initialize counts to zero
    for (int b = 0; b < nBins; ++b) bin[b] = 0;
    
    // accumulate counts
    Angle binWidth = 360.0 / nBins;
    int maxBinSize = 0;
    for (size_t t = 0; t < torsions.size(); ++t)
    {
        Angle torsion = torsions[t] * Rad2Deg;
        while (torsion < 0) torsion += 360.0;
        while (torsion >= 360.0) torsion -= 360.0;
        int binNumber = static_cast<int>((torsion / 360.0) * nBins);
        ++bin[binNumber];
        if (bin[binNumber] > maxBinSize) maxBinSize = bin[binNumber];
    }
    
    // print histogram
    
    // counts per hash mark
    float poundSize = maxBinSize / 40.0f;
    
    for (int b = 0; b < nBins; ++b) 
    {
        os << setw(6) << fixed << setprecision(1) << b * binWidth;
        os << " -- ";
        os << setw(6) << fixed << setprecision(1) << (b + 1) * binWidth;
        os << "   ";
        os << setw(5) << bin[b];
        os << "   ";
        
        for (int h = 0; h < (int)(bin[b] / poundSize); ++h )
            os << "#";
        
        os << std::endl;
    }

	return os;
}


int main() { try 
{
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system); 

    // Atom classes are available, but not charged atom types for ethane
    // in standard Amber force field
    dumm.loadAmber99Parameters();

    GeneralForceSubsystem forces(system);
    // Sulfur vdw params are 0.2 nm, 1.046 
    VanderWallSphere boundary(forces, dumm, Vec3(0,0,0), 1.00, 0.2, 1.046);
    
    if (! Biotype::exists("ethane", "C"))
        Biotype::defineBiotype(Element::Carbon(), 4, "ethane", "C");
    if (! Biotype::exists("ethane", "H"))
        Biotype::defineBiotype(Element::Hydrogen(), 1, "ethane", "H");

    dumm.defineChargedAtomType(
        DuMM::ChargedAtomTypeIndex(5000), 
        "ethane C", 
        DuMM::AtomClassIndex(1), // "CT" type in amber 
        -0.060 // made up
        );
    dumm.setBiotypeChargedAtomType( DuMM::ChargedAtomTypeIndex(5000), Biotype::get("ethane", "C").getIndex() );

    dumm.defineChargedAtomType(
        DuMM::ChargedAtomTypeIndex(5001), 
        "ethane H", 
        DuMM::AtomClassIndex(34), // "HC" type in amber 
        0.020 // made up, use net neutral charge
        );
    dumm.setBiotypeChargedAtomType( DuMM::ChargedAtomTypeIndex(5001), Biotype::get("ethane", "H").getIndex() );

    Ethane ethane1, ethane2;
    ethane2.setPdbResidueNumber(2);

    // place first ethane, units are nanometers
    // skew it a little to break strict symmetry
    system.adoptCompound(ethane1, Transform(Vec3(-0.5, 0, 0)) * Transform(Rotation(0.1, YAxis)) ); 

    // place second ethane, units are nanometers
    system.adoptCompound(ethane2, Vec3( 0.5, 0, 0)); 

	// Write PDB coordinates to the screen
	/*
    system.addEventReporter(new PeriodicPdbWriter(system, std::cout, 0.050));
	/* */
    
    system.addEventHandler(
    	new VelocityRescalingThermostat(system, 5000.0, 0.100));

    std::vector<Angle> torsions;
    system.addEventReporter(
            new PeriodicTorsionMeasurer(torsions, ethane1, "H1", "C1", "C2", "H4", 0.030) );
        
	// Send coordinates to VMD visualization program
    /* 
	system.addEventReporter(
             new PeriodicVmdReporter(system, 0.008, 3001, true) ); 
    /* */

    system.modelCompounds(); // finalize multibody system

    State state = system.realizeTopology();

    // Simulate it.
    
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(0.200);

    printHistogram(torsions, cerr);
    
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

