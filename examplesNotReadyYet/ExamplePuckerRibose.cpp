#include "SimTKmolmodel.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace SimTK;

// Coefficients for ribose pucker torsions
Real a[] = {39.2, 37.6, 35.8, 36.4, 36.8};
Real b[] = {-161.8, 52.8, -91.3, 124.8, -17.9};
String nuName[] = {"nu0", "nu1", "nu2", "nu3", "nu4"};

int main() 
{
    RibonucleotideResidue::Adenylate adenylate;
    adenylate.setPdbResidueNumber(1);

    // Cap 3' end with a phosphate
	adenylate.bondCompound(
			"3PrimePhosphate", 
			ThreePrimeRnaPhosphateGroup("  A", "  A", 'A'), 
			"bondNext", 
			0.1610);

    // Cap 5' end with a phosphate
    adenylate.convertInboardBondCenterToOutboard();
	adenylate.bondCompound(
		"5PrimePhosphate", 
		FivePrimeRnaPhosphateGroup("  A", "  A", 'A'), 
		"bondPrevious");

    // Place base so it doesn't crash during flexing
    adenylate.setDefaultDihedralAngle("chi", 150.0*Deg2Rad);

    // Place 2' hydroxyl so it points away
    adenylate.defineDihedralAngle("HO2", "C3*", "C2*", "O2*", "HO2'");
    adenylate.setDefaultDihedralAngle("HO2", 180.0*Deg2Rad);

    // TODO - transform coordinates to match a fixed ribose location
    String referenceRiboseString = ""
"ATOM      8  C4* A       1      -0.700   3.825   0.740  1.00  0.00           C\n"
"ATOM     10  O4* A       1       0.086   4.267  -0.344  1.00  0.00           O\n"
"ATOM     11  C3* A       1      -2.137   3.837   0.126  1.00  0.00           C\n"
"ATOM     14  C2* A       1      -1.967   4.817  -1.078  1.00  0.00           C\n"
"ATOM     16  C1* A       1      -0.501   4.525  -1.531  1.00  0.00           C\n";
    istringstream referenceRiboseStream(referenceRiboseString);
    PdbStructure referenceRiboseStructure = PdbStructure(referenceRiboseStream);
    Compound::AtomTargetLocations referenceRiboseTargets = adenylate.createAtomTargets(referenceRiboseStructure);
    adenylate.matchDefaultTopLevelTransform(referenceRiboseTargets);
    cout << " Found " << referenceRiboseTargets.size() << " atom targets" << endl;

    ofstream os("ribose_trajectory.pdb");

    int modelNumber = 1;
    for (Angle p = -20; p < 230.0; p += 1.0)
    {
        for (int nu = 0; nu <= 4; ++nu)
        {
            Angle torsionAngle = a[nu]*Deg2Rad * sin( (p - b[nu])*Deg2Rad );
            adenylate.setDefaultDihedralAngle(nuName[nu], torsionAngle);
        }
        adenylate.matchDefaultTopLevelTransform(referenceRiboseTargets);
	// cout << adenylate.getTopLevelTransform() <<  endl;

        os << "MODEL " << modelNumber << endl;
        adenylate.writeDefaultPdb(os);
        os << "ENDMDL" << endl;
        ++modelNumber;
    }

    return 0;
}

