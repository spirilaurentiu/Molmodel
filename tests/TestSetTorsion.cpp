#include "SimTKmolmodel.h"
#include <iostream>

using namespace std;
using namespace SimTK;

void verifyAngle(Real standardAngle, Real testAngle) 
{
    // Normalize to (-Pi:Pi] range
    while (standardAngle <= -Pi) standardAngle += 2*Pi;
    while (standardAngle > Pi) standardAngle -= 2*Pi;

    // Normalize to (-Pi:Pi] range
    while (testAngle <= -Pi) testAngle += 2*Pi;
    while (testAngle > Pi) testAngle -= 2*Pi;

    cout << "Desired angle = " << standardAngle*DuMM::Rad2Deg << "\t";
    cout << "Test angle = " << testAngle*DuMM::Rad2Deg << endl;

    // Make sure set and get values are "pretty much" the same
    Real error = testAngle - standardAngle;
    error *= error;
    assert(error < 1e-6); // within one milliradian
}

int main() 
{
// try {

    int granularity = 11; // number of steps around the torsion
    Real angleStep = 2 * Pi / granularity;

    // Use integer as loop index to ensure that top value is hit
    for (int i = -granularity; i <= granularity; ++i)
    {
        // Dec 10, 2007 - works with ethane, not with cytdine
        RibonucleotideResidue::Cytidylate molecule;
        String torsionName("chi");

        // RiboseCore molecule("1", "1", '1');
        // String torsionName("delta");
        // molecule.assignBiotypes();

        // Ethane molecule;
        // String torsionName("torsion");

        Real desiredAngle = i * angleStep;

        // 1) Set default angle
        molecule.setDefaultDihedralAngle(torsionName.c_str(), desiredAngle);
        Real measuredAngle = molecule.calcDefaultDihedralAngle(torsionName.c_str());

        // molecule.writeDefaultPdb(cout);

        verifyAngle(desiredAngle, measuredAngle);

        // 2) Set angle in simbody model
        CompoundSystem system;
        SimbodyMatterSubsystem matter(system);
        DuMMForceFieldSubsystem dumm(system);
        dumm.loadAmber99Parameters();
        dumm.loadTestMoleculeParameters();

        molecule.assignBiotypes();

        system.adoptCompound( molecule );

        system.modelCompounds();

        State state = system.realizeTopology();
        system.realizeModel(state);
 
        Real initialAngle = molecule.calcDihedralAngle(state, torsionName.c_str());

        verifyAngle(desiredAngle, initialAngle);

        molecule.setDihedralAngle(state, torsionName.c_str(), desiredAngle);

        Real angleAfterSet = molecule.calcDihedralAngle(state, torsionName.c_str());

        verifyAngle(desiredAngle, angleAfterSet);
    }

    cout << "PASSED" << endl;
    return 0;

}

//catch (const std::exception& e)
//{
//printf("EXCEPTION THROWN: %s\n", e.what());
//}
//catch (...)
//{
//printf("UNKNOWN EXCEPTION THROWN\n");
//}    return 0;
//}
