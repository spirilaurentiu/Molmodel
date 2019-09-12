#include "Molmodel.h"

#include <iostream>
#include <fstream>

using namespace SimTK;

int main()
{
try {
   // Initialize simbody objects
   CompoundSystem system;
   SimbodyMatterSubsystem matter(system);
   DecorationSubsystem decoration(system);
   DuMMForceFieldSubsystem forces(system);
   forces.loadAmber99Parameters();

   // Initialize molmodel objects
   RNA rna("AAA");

   // Set bond mobilities
   //Compound& residue1 = rna.updResidue(Compound::Index(0));
   //Compound& residue2 = rna.updResidue(Compound::Index(1));
   //Compound& residue3 = rna.updResidue(Compound::Index(2));

   // Set first residue to Euclidean mobilities
   rna.setResidueBondMobility(ResidueInfo::Index(0), BondMobility::Free);

   // for (Compound::BondIndex bond(0); bond < residue1.getNumBonds(); ++bond)
   //    residue1.setBondMobility(BondMobility::Free, bond);

   // Leave second residue at default combination of Torsion and Rigid mobilities

   // Set third residue to Rigid
   rna.setResidueBondMobility(ResidueInfo::Index(2), BondMobility::Rigid);
   //for (Compound::BondIndex bond(0); bond < residue3.getNumBonds(); ++bond)
   //    residue3.setBondMobility(BondMobility::Rigid, bond);

   system.adoptCompound(rna);
   system.modelCompounds();

   // Write pdb file to use for initial topology in VMD
   std::ofstream pdbOut("ThreeAdenylates.pdb");
   rna.writeDefaultPdb(pdbOut);

   // Maintain temperature
   system.addEventHandler(new VelocityRescalingThermostat(
           system,  293.15, 0.1));

   // Send coordinates to VMD
   system.addEventReporter(new PeriodicVmdReporter(system, 0.015, 3000, true));

   // system.addEventReporter(new PeriodicPdbWriter(
   //       system, cout, 0.015));

   // Visualizer viz(system);
   // system.addEventReporter( new Visualizer::Reporter(viz, 0.015) );

   // Prepare for molecular dynamics
   system.realizeTopology();
   State state = system.getDefaultState();

   // Relax the structure before dynamics run
   LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

   VerletIntegrator integrator(system);
   TimeStepper timeStepper(system, integrator);
   timeStepper.initialize(state);

   // Start simulation
   timeStepper.stepTo(500.0);

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

