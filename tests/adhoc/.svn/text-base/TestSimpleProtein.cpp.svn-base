#include "SimTKmolmodel.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <ctime>
#include <algorithm>

using namespace SimTK;
using namespace std;

std::vector<State> allStates;

class NoiseMaker : public PeriodicEventReporter 
{
public:
	NoiseMaker(CompoundSystem& system, const Compound& compound, const Force::Thermostat& thermo, 
		std::ostream& pdbFile, Real reportInterval) 
        : PeriodicEventReporter(reportInterval), system(system), compound(compound), 
		  thermostat(thermo), pdbFile(pdbFile), tempOut("tempOut.txt")
    {
	}

    void handleEvent(const State& state) const 
    {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
		static Real lastTime = -Infinity;
		Real t = state.getTime();
		system.realize(state, Stage::Dynamics);
		const Real T = thermostat.getCurrentTemperature(state);
		const Real Eb = thermostat.calcBathEnergy(state);
		const Real Es = system.calcEnergy(state);

		static Real Tavg;
		const Real movingAvgT = 0.9;
		if (t==0) Tavg = T;
		else Tavg = movingAvgT*Tavg + (1-movingAvgT)*T;
		std::cout << "TIME = " << t << " TEMP=" << T << std::endl;
		std::cout << "   TAVG=" << Tavg << " Es+Eb=" << Es+Eb << std::endl;
		std::cout << " z=" << thermostat.getChainState(state) << std::endl;

		std::ostream& o = tempOut;
		o << t << " " << T;
		const Vector z = thermostat.getChainState(state);
		for (int i=0; i < z.size(); ++i)
			o << " " << z[i];
		o << std::endl;

		if (t-lastTime >= 0.0999) {
		    //allStates.push_back(state);
			compound.writePdb(state, pdbFile);
			pdbFile << "END\n";
			lastTime = t;
		}
    }
private:
    CompoundSystem& system;
	const Compound& compound;
	const Force::Thermostat& thermostat;
	std::ostream& pdbFile;
	mutable std::ofstream tempOut;
};

class ChangeTemperature : public PeriodicEventHandler 
{
public:
    ChangeTemperature(CompoundSystem& system, 
						const Force::Thermostat& nht,
						Real reportInterval) 
        : PeriodicEventHandler(reportInterval), system(system), thermostat(nht)
    {}


    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const
    {
		const Real t = state.getTime();
		const Real Tbath = thermostat.getBathTemperature(state);
		system.realize(state, Stage::Velocity);
		const Real Tactual = thermostat.getCurrentTemperature(state);
		std::cout << "BATH TEMP IS " << Tbath << " ACTUAL IS " << Tactual << std::endl;
		if (t < 50) return;

		Real newT = Tbath*.8;
		std::cout << "--> CHANGE TEMP TO " << newT << std::endl;
		thermostat.setBathTemperature(state, newT);
    }
private:
    CompoundSystem& system;
	//const NoseHooverThermostat* thermostat;
	const Force::Thermostat& thermostat;
};

int main() {
try {
	// molecule-specialized simbody System
    CompoundSystem system;

	SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
	//Force::GlobalDamper(forces, matter, 0.1);

	// molecular force field
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

	GeneralForceSubsystem forces(system);
	const Real temp = 300;
	const Real tRelax = 1/(2*Pi); // ps
	const Real thermalMass = 100;
	Force::Thermostat thermostat(forces, matter, SimTK_BOLTZMANN_CONSTANT_MD,
								 temp, tRelax);
	//thermostat.setDefaultNumChains(3);
    system.addEventHandler(new ChangeTemperature(system, thermostat, 20));

	//bool useNoseHoover = true;
	//if (useNoseHoover) {
	//	//NoseHooverThermostat* thermostat = 
	//	//	new NoseHooverThermostat(matter, forces, forceField, temp, thermalMass, 2);
	//	//Force::Custom(forces, thermostat);  

	//} else {
	//	// Use velocity rescaling thermostat.
	//	system.addEventHandler(new VelocityRescalingThermostat(
	//		   system,  500, /*0.1*/1));
	//}

	std::ofstream pdbOut("pdbOut.pdb");

	//std::ifstream pdbIn("villin.pdb");
	//Protein protein(pdbIn, .02);

	//Protein protein("MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF");
	//  123456789012345678901234567890123456
	// Villin headpiece: 1vii.pdb is 36 residues

    //Protein protein("SIMTK");
    Protein protein("AAAAAAAAAA");
	protein.setCompoundBondMobility(BondMobility::Torsion);
    protein.assignBiotypes();

    system.adoptCompound(protein);

	// finalize multibody system
    system.modelCompounds(); 

    system.addEventReporter(
		new NoiseMaker(system, protein, thermostat, pdbOut, .05));

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, /*0.020*/.5) );

	//forceField.setUseOpenMMAcceleration(true);
    forceField.setTraceOpenMM(true);

	//forceField.setAllGlobalScaleFactors(0);
	//forceField.setVdwGlobalScaleFactor(1);
	//forceField.setCoulombGlobalScaleFactor(1);
	//forceField.setGbsaGlobalScaleFactor(1);

	//


	// Instantiate simbody model
	system.realizeTopology();
	State state = system.getDefaultState();

    if (forceField.isUsingOpenMM())
        std::cout << "**** USING OpenMM -- Platform " << forceField.getOpenMMPlatformInUse() << std::endl;
    else
        std::cout << "**** NOT USING OpenMM\n";

	system.realizeModel(state);
	Random::Uniform rand;
	for (int i=6; i < state.getNU(); ++i)
		state.updU()[i] = 0.001*(rand.getValue()-.5);
	system.realize(state, Stage::Velocity);
	Vec3 comVel = matter.calcSystemMassCenterVelocityInGround(state);
	std::cout << "COM vel=" << comVel << std::endl;
	Vec3::updAs(&state.updU()[3]) = -comVel;
	system.realize(state, Stage::Velocity);
	comVel = matter.calcSystemMassCenterVelocityInGround(state);
	std::cout << "COM vel=" << comVel << std::endl;

	//std::cout << "INITIAL U's " << state.getU() << std::endl;

	std::cout << forceField.getNumAtoms() << " ATOMS on " << matter.getNumBodies() << " BODIES.\n";

	// Relax the structure before dynamics run
	//LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

	//allStates.clear();
	//allStates.reserve(10000);

    // Simulate it.
	const Real runTime = 200; // ps
	std::cout << "START " << forceField.getNumAtoms() << " ATOMS.\n";
	const clock_t start = std::clock();

	//const int nSteps = 1000;
	//for (int i=0; i<nSteps; ++i) {
	//	system.realize(state);
	//	state.updQ()[0] = state.getQ()[0];
	//}
	//std::cout << nSteps << "evals took " << (double)(std::clock()-start)/CLOCKS_PER_SEC << "s.\n";

    //CPodesIntegrator integ(system);
    //VerletIntegrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
	integ.setAccuracy(.3e-1);
	//integ.setAccuracy(1e-3);
    TimeStepper ts(system, integ);
	try {
		ts.initialize(state);
		ts.stepTo(runTime);

	} catch(const std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
	}
	catch(...) {
		std::cerr << "ERROR: An unknown exception was raised" << std::endl;
	}
	std::cout << runTime << "ps took " << (double)(std::clock()-start)/CLOCKS_PER_SEC << "s.\n";
	std::cout << "N time steps taken (attempted)=" << integ.getNumStepsTaken() << "(" << integ.getNumStepsAttempted()
		<< ")\n";
    std::cout << integ.getNumIterations() << " iterations (conv=" << integ.getNumConvergentIterations() 
              << " div=" << integ.getNumDivergentIterations() << ")\n";
	while (true) {
		std::cout << "Watch animation? " << std::flush;
		std::string s; std::cin >> s;
		//for (unsigned i=0; i < allStates.size(); ++i)
			//viz->report(allStates[i]);
	}

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

