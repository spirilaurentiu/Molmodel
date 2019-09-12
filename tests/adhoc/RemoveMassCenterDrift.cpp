#include "Molmodel.h"

#include <iostream>
#include <exception>
#include <fstream>
#include <ctime>
#include <cmath>

#define SIMDURATION 1000        // ps

using namespace SimTK;
using std::cout; using std::endl;

// See below.
static std::clock_t start = 0; 

// This is a reporter so we can get some output during the simulation.
// Watch the COM velocity.
class SaySomething : public PeriodicEventReporter {
public:
	SaySomething(const CompoundSystem& system, const VelocityRescalingThermostat& rescaler, Real reportInterval) 
        :   PeriodicEventReporter(reportInterval), system(system), rescaler(rescaler) {}

    void handleEvent(const State& state) const {
		const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

        // This realize() is here to make sure we can reference 
        // velocity-related quantities.
        system.realize(state, Stage::Velocity);
        const Real T  = rescaler.calcCurrentTemperature(state);
        const Vec3 com = system.calcSystemMassCenterLocation(state);
        const Vec3 Vcom = matter.calcSystemMassCenterVelocityInGround(state);

		std::cout << "TIME = " << state.getTime() << " TEMP=" << T
                  << " Elapsed(s)=" << double(std::clock()-start)/CLOCKS_PER_SEC
                  << std::endl;
        std::cout << " COM pos=" << com << " |vel|=" << Vcom.norm() << endl;
    }
private:
    const CompoundSystem&               system;
    const VelocityRescalingThermostat&  rescaler;
};

class SaveState : public PeriodicEventReporter {
public:
	SaveState(Real reportInterval) 
    :   PeriodicEventReporter(reportInterval) 
    {saved.reserve(5000);}

    void clear() {saved.clear(); saved.reserve(5000);}

    void handleEvent(const State& state) const {
		saved.push_back(state);
    }

    const std::vector<State>& getStates() const {return saved;}
private:
    mutable std::vector<State> saved;
};

// Where we want the molecule's COM to stay.
static Vec3 WantCOM(1,2,3);
int main() { 
	try {
		// Load the PDB file and construct the system ... and specify the forcefield
		CompoundSystem system;
		SimbodyMatterSubsystem matter(system);
        GeneralForceSubsystem forces(system);
		DecorationSubsystem decorations(system);
		DuMMForceFieldSubsystem forceField(system);
		forceField.loadAmber99Parameters();
        //forceField.setAllGlobalScaleFactors(0);
		//PDBReader pdb("20ala.min.pdb");
		//pdb.createCompounds(system);
        //Protein pdb("AEEEA");
        //pdb.setTopLevelTransform (Transform(Rotation(Pi/3, Vec3(1,2,3)), Vec3(10,5,6)));
        //system.adoptCompound(pdb);
        //Protein pdb2("AAAAA");
        //pdb2.setTopLevelTransform(Transform(Rotation(Pi/6, Vec3(.1,2,.3)), Vec3(8,5,6)));
        //system.adoptCompound(pdb2);
        //Protein pdb3("AAAAA");
        //pdb3.setTopLevelTransform(Transform(Rotation(Pi/4, Vec3(1,.2,3)), Vec3(8,6,6)));
        //system.adoptCompound(pdb3);
        Protein pdb4("KKK");
        pdb4.setTopLevelTransform(Transform(Rotation(Pi/5, Vec3(1,2,.03)), Vec3(10,4,6)));
        system.adoptCompound(pdb4);
		system.modelCompounds();

        matter.updGround().addBodyDecoration(WantCOM, DecorativeFrame().setColor(Red));

        //forceField.setUseOpenMMAcceleration(true);
        //forceField.setTraceOpenMM(true);

	    const Real temp = 300;

        MassCenterMotionRemover& remover = 
            *new MassCenterMotionRemover(system, 10);
 		system.addEventHandler(&remover);
        remover.setDesiredMassCenterLocation(WantCOM);
        remover.setMassCenterLocationTolerance(1e-3);
        //remover.disableAngularMomentumRemoval(true);
        remover.enableMassCenterCorrection(true);

   	
        VelocityRescalingThermostat& rescaler = 
            *new VelocityRescalingThermostat(system, temp, 0.1, 6);
		system.addEventHandler(&rescaler);

        //NoseHooverThermostat noseHoover(forces,matter, temp, 1, 6);

        // Show me a movie
        Visualizer viz(system);
        system.addEventReporter( new Visualizer::Reporter(viz, 0.050) );

        SaveState* saveState = new SaveState(.05);
        //system.addEventReporter(saveState);

        system.addEventReporter(new SaySomething(system,rescaler,1.));
		
		// write output to pdb
		std::ofstream pdbfile;
        pdbfile.open("c:/temp/output.pdb");
		//system.addEventReporter(new PeriodicPdbWriter(system, pdbfile, 0.1)); 
		system.realizeTopology();
		
		// Create an initial state for the simulation.
		State state = system.getDefaultState();
		//pdb.createState(system, state);
		LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);

        system.realize(state, Stage::Velocity);

        cout << "Num thermal dofs=" << rescaler.calcNumThermalDofs(state)
             << " after excluding " << rescaler.getNumExcludedDofs() 
             << " rigid body dofs\n";

        cout << "BEFORE u=" << state.getU()
             << "\nTemp=" << rescaler.calcCurrentTemperature(state) 
             << " |mom|=" << matter.calcSystemCentralMomentum(state).norm() << endl;

        rescaler.rescale(state); // randomizes initially
        system.realize(state, Stage::Velocity);

        cout << "RESCALED u=" << state.getU() 
             << "\nTemp=" << rescaler.calcCurrentTemperature(state) 
             << " |mom|=" << matter.calcSystemCentralMomentum(state).norm() << endl;

        SpatialVec sysCMOM = matter.calcSystemCentralMomentum(state);
        cout << "Central mom=" << sysCMOM << endl;

        // Kill any rigid body linear or angular momentum.
        remover.removeSystemMomentum(state);
        system.realize(state, Stage::Velocity);

        sysCMOM = matter.calcSystemCentralMomentum(state);
        cout << "Now central mom=" << sysCMOM << endl;

        cout << "NO MOMENTUM u=" << state.getU()
             << "\nTemp=" << rescaler.calcCurrentTemperature(state) 
             << " |mom|=" << matter.calcSystemCentralMomentum(state).norm() << endl;

        rescaler.rescale(state);
        system.realize(state, Stage::Velocity);

        cout << "FINAL RESCALED u=" << state.getU()
             << "\nTemp=" << rescaler.calcCurrentTemperature(state) 
             << " |mom|=" << matter.calcSystemCentralMomentum(state).norm() << endl;


		
        start = std::clock();   // start the wallclock timer

		// Choose an integrator.

        // Note: Verlet at low accuracy is a very poor integrator.
        // You get about the same amount of drift with RK4 at 1e-2
        // accuracy as Verlet at 1e-4. But we want lots of drift here to
        // show the COM mover working!

		VerletIntegrator integ(system); integ.setAccuracy(1e-1);
        //RungeKuttaMersonIntegrator integ(system); integ.setAccuracy(1e-2);

        // Simulate.
		TimeStepper ts(system, integ);
		ts.initialize(state);
		ts.stepTo(SIMDURATION);

        std::cout << "Done. " << SIMDURATION << " ps in elapsed(s)=" 
                  << double(std::clock()-start)/CLOCKS_PER_SEC
                  << endl;
		pdbfile.close();

        for (unsigned frame=0; frame < saveState->getStates().size(); ++frame)
            viz.report(saveState->getStates()[frame]);
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


