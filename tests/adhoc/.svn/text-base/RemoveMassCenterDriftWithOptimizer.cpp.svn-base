
/* This example removes overall system momentum using an optimizer -- that works,
but it is a very bad way to do it. Instead, use the CompoundSystem method
removeSystemRigidBodyMomentum() which does it analytically. (That is actually
inherited from CompoundSystem's parent class.)
*/

#include "Molmodel.h"

#include <iostream>
#include <exception>
#include <fstream>
#include <ctime>

#define SIMDURATION 1000        // ps

using namespace SimTK;
using std::cout; using std::endl;

// See below.
static void removeSystemMomentum(const CompoundSystem& system, State& state); 
static std::clock_t start = 0; 

// This is a reporter so we can get some output during the simulation.
// Watch the COM velocity.
class SaySomething : public PeriodicEventReporter {
public:
	SaySomething(CompoundSystem& system, Real reportInterval) 
        :   PeriodicEventReporter(reportInterval), system(system) {}

    void handleEvent(const State& state) const {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

        // This realize() is here to make sure we can reference 
        // velocity-related quantities.
        system.realize(state, Stage::Velocity);
        const Real KE = system.calcKineticEnergy(state);
        const int  N  = matter.getNumMobilities(); // total # of dofs
        const Real T  = (2*KE)/(N*SimTK_BOLTZMANN_CONSTANT_MD);
        const Vec3 Vcom = matter.calcSystemMassCenterVelocityInGround(state);

		std::cout << "TIME = " << state.getTime() << " TEMP=" << T
                  << " Elapsed(s)=" << double(std::clock()-start)/CLOCKS_PER_SEC
                  << std::endl;
        std::cout << " COM vel=" << Vcom.norm() << endl;
    }
private:
    CompoundSystem& system;
};

// This is an event handler that can be called periodically to remove
// most of the rigid body linear and angular momentum. If you allow this
// to build up you might also want a handler that shifts the COM back
// to (0,0,0) once in a while.
class KillMomentum : public PeriodicEventHandler {
public:
    KillMomentum(CompoundSystem& system, Real interval)
    :   PeriodicEventHandler(interval), system(system) {}

    // For convenience.
    void killMomentum(State& state) {
        bool dummy;
        handleEvent(state, 0, dummy);
    }

    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const 
    {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

        system.realize(state, Stage::Velocity);
        // A SpatialVec has both an angular and linear subcomponent.
        SpatialVec sysMom = matter.calcSystemMomentumAboutGroundOrigin(state);
            cout << "REMOVING MOMENTUM before=" << sysMom.norm() << endl;
        removeSystemMomentum(system, state);
        system.realize(state, Stage::Velocity);
        sysMom = matter.calcSystemMomentumAboutGroundOrigin(state);
            cout << "... after=" << sysMom.norm() << endl;
    }
private:
    CompoundSystem& system;
};

int main() { 
	try {
		// Load the PDB file and construct the system ... and specify the forcefield
		CompoundSystem system;
		SimbodyMatterSubsystem matter(system);
		DecorationSubsystem decorations(system);
		DuMMForceFieldSubsystem forceField(system);
		forceField.loadAmber99Parameters();
        //forceField.setAllGlobalScaleFactors(0);
		//PDBReader pdb("20ala.min.pdb");
		//pdb.createCompounds(system);
        Protein pdb("AAAAA");
        system.adoptCompound(pdb);
		system.modelCompounds();

	    const Real temp = 300;

        KillMomentum* killer = new KillMomentum(system, 10); // once per 10 ps
 		system.addEventHandler(killer);
   	
		system.addEventHandler(new VelocityRescalingThermostat(system, temp, 0.1));

        // Show me a movie
        Visualizer viz(system);
        system.addEventReporter( new Visualizer::Reporter(viz, 0.010) );

        system.addEventReporter(new SaySomething(system,1.));
		
		// write output to pdb
		std::ofstream pdbfile;
        pdbfile.open("c:/temp/output.pdb");
		system.addEventReporter(new PeriodicPdbWriter(system, pdbfile, 0.1)); 
		system.realizeTopology();
		
		// Create an initial state for the simulation.
		State state = system.getDefaultState();
		//pdb.createState(system, state);
		LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

        // Set all generalized speeds to random values to test
        // the momentum killer.
        Random::Uniform random;
        for (int i=0; i < state.getNU(); ++i)
            state.updU() = random.getValue();

        // Kill any rigid body linear or angular momentum.
        killer->killMomentum(state);
		
        start = std::clock();   // start the wallclock timer

		// Choose an integrator.

        // Note: Verlet at low accuracy is a very poor integrator.
        // You get about the same amount of drift with RK4 at 1e-2
        // accuracy as Verlet at 1e-4.

		//VerletIntegrator integ(system); integ.setAccuracy(1e-4);
        RungeKuttaMersonIntegrator integ(system); integ.setAccuracy(1e-2);

        // Simulate.
		TimeStepper ts(system, integ);
		ts.initialize(state);
		ts.stepTo(SIMDURATION);

        std::cout << "Done. " << SIMDURATION << " ps in elapsed(s)=" 
                  << double(std::clock()-start)/CLOCKS_PER_SEC
                  << endl;
		pdbfile.close();
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

// This satisfies the interface required by SimTK::Optimizer for an
// objective function to optimize. This one takes six scalar parameters
// (the base body angular and linear velocities) and uses them to drive
// the overall system momentum to zero.
class RemoveSystemMomentum : public OptimizerSystem {
public:
    RemoveSystemMomentum(const CompoundSystem& system, const State& initState) :
        OptimizerSystem(6), system(system), state(initState) {
    }
    int objectiveFunc(const Vector& parameters, bool new_parameters, Real& f) const {
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        const MobilizedBody& baseBody = matter.getMobilizedBody(MobilizedBodyIndex(1)); 
        baseBody.setUFromVector(state, parameters);
        system.realize(state, Stage::Velocity);
        const SpatialVec sysMom = matter.calcSystemMomentumAboutGroundOrigin(state);
        f = sysMom.norm();
        return 0;
    }

private:
    const CompoundSystem& system;
    mutable State state; // a temporary state for use while optimizing
};

// This static function invokes the Optimizer on the above objective.
static void removeSystemMomentum(const CompoundSystem& system, State& state) {
    RemoveSystemMomentum objective(system,state);
    Optimizer opt(objective);
    opt.useNumericalGradient(true);

    // Get initial base body generalized speeds.
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    const MobilizedBody& baseBody = matter.getMobilizedBody(MobilizedBodyIndex(1)); 
    Vector vel = baseBody.getUAsVector(state);

    // Optimize.
    try {opt.optimize(vel);} catch (...) {}

    // Copy optimized parameters back into base body generalized speeds.
    baseBody.setUFromVector(state, vel);
    system.realize(state, Stage::Velocity);
}


