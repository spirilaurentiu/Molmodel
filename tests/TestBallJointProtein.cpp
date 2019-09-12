/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKmolmodel.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;

using namespace SimTK;
using namespace std;

// define CREATE_VIZ_WINDOW to see animated window of simulation
// undefine for automated nightly builds
 #define CREATE_VIZ_WINDOW

void testIssue907()
{
    CompoundSystem system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem forces(system);

    DuMMForceFieldSubsystem dumm(system);

    dumm.loadAmber99Parameters();
    dumm.setCoulombGlobalScaleFactor(0);//meterReader.globalCoulombScaleFactor);
    dumm.setBondTorsionGlobalScaleFactor(1);//meterReader.globalBondTorsionScaleFactor);
    dumm.setGbsaGlobalScaleFactor(0);//meterReader.globalGbsaScaleFactor);
    dumm.setVdwGlobalScaleFactor(0);//meterReader.globalVdwScaleFactor);
    dumm.setBondStretchGlobalScaleFactor(0);//meterReader.globalBondStretchScaleFactor);
    dumm.setBondBendGlobalScaleFactor(0);//meterReader.globalBondBendScaleFactor);
    dumm.setAmberImproperTorsionGlobalScaleFactor(0);//meterReader.globalAmberImproperTorsionScaleFactor);
     
    Biopolymer myMolecule;
    // myMolecule = Protein(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //50 //ok
    myMolecule = Protein("AAA");
    myMolecule.assignBiotypes();

    system.adoptCompound(myMolecule,Vec3(  0,0,0));
    cout<<"about to modelCompounds"<<endl;
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<" time:"<< asctime (timeinfo)   <<endl;
    system.modelCompounds();
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<" time:"<< asctime (timeinfo)   <<endl;
    cout<<"done with modelCompounds"<<endl;
    State state = system.realizeTopology();
    system.realize(state,Stage::Position);
    VerletIntegrator integrator(system,.002);

    TimeStepper ts(system,integrator);
    ts.initialize(state);
    //# putting writeDefaultPdb here was fine to 2200 residues
    cout <<"1 about to write default pdb coords to  test.pdb"<<endl;
    //myMolecule.writeDefaultPdb("test.pdb",Vec3(0));
    cout <<"1 just wrote default pdb coords to  test.pdb"<<endl;
    cout<<"[Repel.h:ConstrainedDynamics] Starting dynamics now."<<endl;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<" time:"<< asctime (timeinfo)   <<endl;
    //cout <<"2 about to write default pdb coords to  test.pdb"<<endl;
    //myMolecule.writeDefaultPdb("test.pdb",Vec3(0));
    //cout <<"2 just wrote default pdb coords to  test.pdb"<<endl;
    ts.stepTo(.1);

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    cout<<" time:"<< asctime (timeinfo)   <<endl;


    cout<<"done."<<endl;
}

void testShortProtein() 
{
    // Biotype::initializePopularBiotypes();

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);
    DecorationSubsystem decorations(system);

    // ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    // dumm.populateFromTinkerParameterFile(tinkerStream);
    // tinkerStream.close();
    dumm.loadAmber99Parameters();

    // Protein protein("ACDEFGHIKLMNPQRSTVWY");
    ///Protein protein("SPFGR");
    Protein protein("AAAAA");
    protein.assignBiotypes();

    // This is what we want to test
    for (unsigned int r = 0 ; r < protein.getNumBonds(); r++){
        protein.setBondMobility(BondMobility::Ball, Compound::BondIndex(r));
    }

    //protein.writeDefaultPdb(cout);

    system.adoptCompound(protein);

    system.modelCompounds();  

    State state = system.realizeTopology();
 
#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system);
    //system.addEventReporter( new Visualizer::Reporter(display, 0.005) );
#endif

    //RungeKuttaMersonIntegrator integrator(system);
    VerletIntegrator integrator(system);

    integrator.initialize(state);

    int numberOfSteps = 3;

    /*
#ifdef CREATE_VIZ_WINDOW
    numberOfSteps = 10;
#endif

    TimeStepper ts(system, integrator);
    ts.initialize(state);
    ts.stepTo(numberOfSteps * 0.0050);

#ifdef CREATE_VIZ_WINDOW
    display.report(integrator.getState());
#endif
*/
    Real RT = 300 * SimTK_BOLTZMANN_CONSTANT_MD;
    Random::Gaussian random;
    Real timeInterval = 0.02;
    Real time = 0.0;
    for (int i = 0; i < 100; i++){
        time = i * timeInterval;

        // Get detM
        int nu = state.getNU();
        SimTK::Vector V(nu);
        system.realize(state, SimTK::Stage::Position);
        matter.realizeArticulatedBodyInertias(state);

        SimTK::Real detM = 1.0;
        SimTK::Vector DetV(nu);
        SimTK::Real* D0 = new SimTK::Real(1.0);
        matter.calcDetM(state, V, DetV, D0);

        // Initialize velocities
        if( !(i % 10)) {
            /*
            SimTK::Vector U = state.updU();
            for (int j = 0; j < state.getNU(); ++j) {
                U[j] = random.getValue();
            }
            */

            double sqrtRT = std::sqrt(RT);
            SimTK::Vector V(nu);
            SimTK::Vector SqrtMInvV(nu);
            for (int i=0; i < nu; ++i){
                V[i] = random.getValue();
            }
            matter.multiplyBySqrtMInv(state, V, SqrtMInvV);

            SqrtMInvV *= (sqrtRT); // Set stddev according to temperature
            system.realize(state, SimTK::Stage::Velocity);

        }

        // Integrate trajectory
        integrator.stepTo(time);

#ifdef CREATE_VIZ_WINDOW
        display.report(integrator.getState());
#endif

    }

}

int main() {
    testShortProtein();
    //testIssue907();
}

