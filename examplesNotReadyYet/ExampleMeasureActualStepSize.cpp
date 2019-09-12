#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"

#include <fstream>

#include <iostream>
#include <exception>

using namespace SimTK;
using namespace std;


class WritePdbReporter : public PeriodicEventReporter {
public:
    WritePdbReporter(
        const MultibodySystem& system, 
        std::ostream& outputStream,
        Real interval) 
        : PeriodicEventReporter(interval), 
          system(system), 
          outputStream(outputStream),
		  currentModelNumber(1)
    {}

	void addCompound(const Compound& compound) {
		compoundPointers.push_back(&compound);
	}

    void handleEvent(const State& state) const {
		outputStream << "MODEL " << currentModelNumber << endl;
        system.realize(state, Stage::Position);
		for (int c = 0; c < (int)compoundPointers.size(); ++c)
			compoundPointers[c]->writePdb(state, outputStream);
		outputStream << "ENDMDL" << endl;

		++currentModelNumber;
    }

private:
	mutable int currentModelNumber;
    const MultibodySystem& system;
	std::vector<const Compound*> compoundPointers;
    std::ostream& outputStream;
};


void integrate(Real accuracy, BondMobility::Mobility mobility, int integratorIndex, std::ostream& os)
{
    cout << "accuracy = " << accuracy << endl;
	cout << "integrator index = " << integratorIndex << endl;
	cout << "mobility index = " << mobility << endl;
    
    cout << "constructing system..." << endl;
    
    CompoundSystem system; // molecule-specialized simbody System
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    TinkerDuMMForceFieldSubsystem dumm(system); // molecular force field
    
    // Two very small proteins
    Protein protein1("SPFAK");
    Protein protein2("SPFAK");
    
	protein1.setPdbChainId('A');
	protein2.setPdbChainId('B');

    // link proteins to forcefields
    dumm.loadAmber99Parameters();
    protein1.assignBiotypes();
    protein2.assignBiotypes();
    
    // Set mobility of all bonds in both structures
	// For "Torsion" mobility, leave all bonds at default values
	if (mobility != BondMobility::Torsion) {
		for (Compound::BondIndex b(0); b < protein1.getNBonds(); ++b) {
			protein1.setBondMobility(mobility, b);
		}
		for (Compound::BondIndex b(0); b < protein2.getNBonds(); ++b) {
			protein2.setBondMobility(mobility, b);
		}
	}
    
    system.adoptCompound(protein1);
    system.adoptCompound(protein2, Vec3(0.5, 0.5, 0) );
    
    // system.updDefaultSubsystem().addEventReporter(new VTKEventReporter(system, 0.010));

	ofstream pdbStream("trajectory.pdb");
	WritePdbReporter writePdbReporter(system, pdbStream, 0.010);
	writePdbReporter.addCompound(protein1);
	writePdbReporter.addCompound(protein2);
	system.updDefaultSubsystem().addEventReporter(&writePdbReporter);
    
    system.updDefaultSubsystem().addEventHandler(new VelocityRescalingThermostat(system, 293.15, 0.020));
    
    system.modelCompounds(); // finalize multibody system
    
    system.realizeTopology();
    
    State& state = system.updDefaultState();
    
    
    cout << "minimizing energy..." << endl;
        
    LocalEnergyMinimizer::minimizeEnergy(system, state, 100.0);
    
    // Choice of integrators
    Integrator* integratorPtr;
    switch (integratorIndex) {
            case 0:
                integratorPtr = new RungeKuttaMersonIntegrator(system);
                break;
            case 1:
                integratorPtr = new VerletIntegrator(system);
                break;
            case 2:
                integratorPtr = new CPodesIntegrator(system);
                break;
            case 3:
                integratorPtr = new ExplicitEulerIntegrator(system);
                break;
    }
 
    Integrator& study = *integratorPtr;

    study.setAccuracy(accuracy);
    TimeStepper timeStepper(system, study);

    timeStepper.initialize(state);

    // cause timeStepper to pause after each integration step
    timeStepper.setReportAllSignificantStates(true);
    study.setReturnEveryInternalStep(true);
    
    cout << "integrating..." << endl;
    
	long int previousForceEvaluationCount = dumm.getForceEvaluationCount();

    std::vector<Real> integratorStepSizes;
    Real endTime = 5.000; // picoseconds
    // Real endTime = 0.050; // picoseconds - fast time for debugging
    while(timeStepper.getTime() < endTime)
    {
        timeStepper.stepTo(endTime);
    
        // store integrator step size
        integratorStepSizes.push_back(study.getPreviousStepSizeTaken());
    }
       
    cout << "analyzing statistics..." << endl;

	// Dump time steps
	//ofstream stepsOut ("timeSteps.dat");
 //   for (int t = 0; t < (int)integratorStepSizes.size(); ++t) {
	//	stepsOut << integratorStepSizes[t] << endl;
	//}	
	//stepsOut.close();

    // Mean step time
    Real meanTimeStep = 0.0;
    for (int t = 0; t < (int)integratorStepSizes.size(); ++t) {
        meanTimeStep += integratorStepSizes[t];
    }
    meanTimeStep /= integratorStepSizes.size();
    
    // Standard deviation of step time
    Real timeStepDev = 0.0;
    for (int t = 0; t < (int)integratorStepSizes.size(); ++t) {
        Real diff = meanTimeStep - integratorStepSizes[t];
        timeStepDev += diff * diff;
    }
    timeStepDev /= (integratorStepSizes.size() - 1);
    timeStepDev = sqrt(timeStepDev);
    
    // Median
    std::sort(integratorStepSizes.begin(), integratorStepSizes.end());
    Real medianTimeStep = integratorStepSizes[(int)(integratorStepSizes.size() / 2)];
    
    os << study.getMethodName();
    os << "\t";
    
    switch (mobility) {
        case BondMobility::Free: os << "All atom"; break;
        case BondMobility::Torsion: os << "Torsion"; break;
        case BondMobility::Rigid: os << "Rigid"; break;
    }
    os << "\t";
    
    os << accuracy;
    os << "\t";
    os << dumm.getForceEvaluationCount() - previousForceEvaluationCount;
    os << "\t";
    os << integratorStepSizes.size();
    os << "\t";
    os << 1000.0 * meanTimeStep;
    os << "\t";
    os << "" << 1000.0 * timeStepDev;
    os << "\t";
    os << 1000.0 * medianTimeStep;
    os << std::endl;
    
    delete integratorPtr;
    
    cout << "Done." << endl;
}

int main() { try 
{
    ofstream os("integ.dat");

    // Vary accuracy to measure effect on integrator step size
    Real accuracyList[] = {1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6};
    // Real accuracyList[] = {3e-3, 1e-3, 3e-4}; // short list
    
    // Compare different modeling strategies
	BondMobility::Mobility mobilityList[] = {BondMobility::Free, BondMobility::Torsion, BondMobility::Rigid};

    // Loop over integrators, accuracies, and mobilities
    for (int integratorIndex = 0; integratorIndex <= 1; ++integratorIndex) 
    {
		// skip CPODES, it's way too slow for molecules, taking >60 force evaluations per step
		if (integratorIndex == 2) continue;

        for (int mIx = 1; mIx <= 1; ++mIx)
        {
            BondMobility::Mobility mobility = mobilityList[mIx];
            for (int aIx = 0; aIx <= 8; ++aIx)
            {
                integrate(accuracyList[aIx], mobility, integratorIndex, os);
                // integrate(3e-3, mobility, integratorIndex, os);
            }
        }
    }

    os.close();

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

