/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Fold Polyalanine                   *
 * -------------------------------------------------------------------------- *
 * This example attempts to fold polyalanine from an extended configuration.  *
 * We expect a sequence of alanines to form a helix.                          *
 *                                                                            *
 * The strategy we'll follow here is a crude annealing protocol in which      *
 * we'll start at a high temperature and progressively lower the temperature  *
 * as we simulate. We hope to see the polyalanine settle into a helix         *
 * when the temperature gets low enough. As we go along we'll track the       *
 * potential energy and save the lowest-energy state we've seen so far. In    *
 * the end we'll report that state as the result.                             *
 *                                                                            *
 * This example demonstrates                                                  *
 *   - modeling a protein directly from sequence                              *
 *   - the Nose-Hoover thermostat                                             *
 *   - use of an event handler to interrupt the simulation to make            *
 *     discontinuous changes (here lowering the requested temperature)        *
 *   - how to use OpenMM for GPU acceleration if it is available, falling     *
 *     back to multithreading otherwise                                       *
 *   - how to remove overall rigid body momentum to keep the molecule from    *
 *     drifting                                                               *
 *   - how to generate a live movie using the Simbody Visualizer and how to   *
 *     write to a pdb file for later visualization in VMD or another advanced *
 *     molecule visualizer                                                    *
 *   - use of Simbody's CPU and real timers, and Pathname class to find out   *
 *     what directory we're running in (that's where the pdb file will be).   *
 *                                                                            *
 * Author: Michael Sherman                                                    *
 * -------------------------------------------------------------------------- */


#include "Molmodel.h"

#include <iostream>
#include <fstream>
#include <exception>
using namespace SimTK;

// Comment this out if it causes trouble for you.
#define TRY_TO_USE_OPENMM


// See below.
static Real startCPU, startRT; // timers
static State minState;
static Real minPE = Infinity;


// This is a reporter so we can get some output during the simulation and
// capture the lowest-energy state we've seen so far.
class SaySomething : public PeriodicEventReporter {
public:
	SaySomething(const CompoundSystem& system, 
                 const NoseHooverThermostat& thermo,
                 Real reportInterval) 
        :   PeriodicEventReporter(reportInterval), system(system), thermo(thermo) {}

    void handleEvent(const State& state) const {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        system.realize(state, Stage::Dynamics);
        const Real pe = system.calcPotentialEnergy(state);

		std::cout << "TIME = " << state.getTime() 
              << " PE=" << pe
              << " T=" << thermo.getCurrentTemperature(state)
              << " RT=" << realTime()-startRT << "s"
              << " CPU=" << cpuTime()-startCPU << "s\n"; 
 
        if (pe < minPE) {
            std::cout << "***** best so far *****\n";
            minState = state;
            minPE = pe;
        }
    }
private:
    const CompoundSystem&               system;
    const NoseHooverThermostat&         thermo;
};

// This event handler implements the annealing protocol by changing
// the temperature of the thermal bath that interacts with our peptide.
class ChangeTemperature : public PeriodicEventHandler {
public:
    ChangeTemperature(const CompoundSystem& system,                  
                      const NoseHooverThermostat& thermo,
                      Real interval,
                      Real frac=.75, // fraction by which we reduce T
                      Real minT=10)  // temperature at which we'll quit
    :   PeriodicEventHandler(interval), system(system), thermo(thermo),
        frac(frac), minT(minT) {}

    // The Simbody TimeStepper invokes this method periodically, based on
    // the time interval that was provided in the constructor.
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const 
    {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

        const Real bathT = thermo.getBathTemperature(state);
        if (state.getTime() == 0) {
            std::cout << "STARTING T=" << bathT << "K. Reduce to " << frac << " every " 
                      << getEventInterval() << "ps until " << minT << "K." << std::endl;
            return;
        }
        
        const Real newBathT = std::max(minT, bathT*frac);
        if (bathT == newBathT) {
            if (bathT <= minT) shouldTerminate=true;
            return;
        }

        thermo.setBathTemperature(state, newBathT);
        std::cout << "Changed bath T to " << newBathT << std::endl;
    }
private:
    const CompoundSystem&       system;
    const NoseHooverThermostat& thermo;
    const Real frac;
    const Real minT;
};

int main() {
try {

	// molecule-specialized Simbody System
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    DecorationSubsystem decorations(system);

	// Molecular force field. GBSA solvation is used by default.
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    // optional -- try OpenMM accelerations
    #ifdef TRY_TO_USE_OPENMM
        forceField.setUseOpenMMAcceleration(true);
    #endif
    // This causes tracing even if we aren't trying to use OpenMM.
    forceField.setTracing(true); // log OpenMM info to console

    // Multithreading is used by default unless you disable it.
    //forceField.setUseMultithreadedComputation(false);
    forceField.setNumThreadsRequested(0); // default

    // You can change to vacuum by uncommenting this line.
    //forceField.setGbsaGlobalScaleFactor(0);

    std::cout << "num threads before realizeTopology (should be zero): " 
              << forceField.getNumThreadsInUse() << std::endl;

    // The Protein constructor always adds neutral end caps unless
    // you suppress them.

    //Protein protein("SIMTK");
    //Protein protein("AAAAAA"); // 6
    //Protein protein("AAAAAAAAAA"); // 10
    //Protein protein("AAAAAAAAAAAA"); // 12
    Protein protein("AAAAAAAAAAAAAAA"); // 15


    // This is the Tryptophan cage 1L2Y.pdb.
    // (ACE) ASN LEU TYR ILE GLN TRP LEU LYS ASP GLY GLY PRO SER          
    // SER GLY ARG PRO PRO PRO SER (NAC)
    //Protein protein("NLYIQWLKDGGPSSGRPPPS");

    protein.assignBiotypes();
    system.adoptCompound(protein);

	// finalize multibody system
    system.modelCompounds(); 

    // Create a PeriodicPdbWriter (a Molmodel-provided Reporter) for 
    // generating pdb frames. We'll keep a pointer to it here so we can generate 
    // a few frames before and after the simulation. We'll ask the TimeStepper
    // to call this every 100fs.
    std::cout << "\nNOTE: Writing pdb file 'polyalanine.pdb' in " 
              << Pathname::getCurrentWorkingDirectory() << "\n\n";
    std::ofstream pdbFile; pdbFile.open("polyalanine.pdb");
    PeriodicPdbWriter* pdbWriter = new PeriodicPdbWriter(system, pdbFile, 0.1);
    system.addEventReporter(pdbWriter); // System takes ownership

    // Create a Visualizer for live animation while simulating and ask
    // the TimeStepper to invoke it every 100fs.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.1) );


    const Real startT = 1000;
    NoseHooverThermostat thermo(forces, matter, startT, .5);

    // We'll stay at a given temperature for a time proportional to
    // the number of residues we're trying to fold.
    const Real timePerResidue = 1.; // ps
    system.addEventHandler(
        new ChangeTemperature(system, thermo, 
                timePerResidue*protein.getNumResidues(), 
                0.9,    // reduce by setting newT = 90% * oldT
                150));  // stop when we get to this temperature

    // We'll remove overall rigid body momentum every 10ps so the
    // molecule doesn't run off the screen.
    system.addEventHandler(new MassCenterMotionRemover(system, 10));

    // Output some stats every 1ps, and remember the lowest energy
    // state we've seen.
    system.addEventReporter(new SaySomething(system, thermo, 1));


	// Instantiate Simbody model
	system.realizeTopology();
	State state = system.getDefaultState();

    std::cout << "num threads in use after realizeTopology: " 
              << forceField.getNumThreadsInUse() << std::endl;

    system.realize(state, Stage::Position);
    viz.report(state);             // first animation frame
    pdbWriter->handleEvent(state); // pdb frame

	// Relax the structure before dynamics run
	LocalEnergyMinimizer::minimizeEnergy(system, state, 10);

    system.realize(state, Stage::Position);
    viz.report(state);             // post-relaxation frame
    pdbWriter->handleEvent(state); // pdb frame

    startCPU = cpuTime();   // start the cpu timer
    startRT  = realTime(); // and the real timer


    // Simulate it.
    //VerletIntegrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1000.0);


    std::cout << "Done. " << integ.getTime() << "ps in elapsed(s)=" 
              << realTime()-startRT << " CPU=" << cpuTime()-startCPU
              << std::endl;

    State finalState = minState;
    system.realize(finalState, Stage::Dynamics);
    std::cout << "Best PE=" << system.calcPotentialEnergy(finalState) << std::endl;
    viz.report(finalState);             // show best conformation
    pdbWriter->handleEvent(finalState); // pdb frame

	LocalEnergyMinimizer::minimizeEnergy(system, finalState, 5);

    system.realize(finalState, Stage::Dynamics);
    std::cout << "Final PE=" << system.calcPotentialEnergy(finalState) << std::endl;
    viz.report(finalState);             // post-minimization animation frame
    pdbWriter->handleEvent(finalState); // final pdb frame

    pdbFile.close();
    std::cout << "\nWrote pdb file "
              << Pathname::getCurrentWorkingDirectory() << "polyalanine.pdb\n";
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

