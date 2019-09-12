
#include "SimTKmolmodel.h"
#include "molmodel/internal/Pdb.h"

#include <sstream>

using namespace SimTK;
using namespace std;

int main() 
{      
    CompoundSystem system;
    DuMMForceFieldSubsystem dumm(system);
    dumm.loadAmber99Parameters();

    const char* uudPdbString = ""
"ATOM      1  P     G     1     -10.655  -7.494   9.276  1.00  0.00           P  \n"
"ATOM      2  OP1   G     1     -11.450  -7.642   8.036  1.00  0.00           O  \n"
"ATOM      3  OP2   G     1      -9.248  -7.952   9.308  1.00  0.00           O  \n"
"ATOM      4  O5*   G     1     -10.689  -5.944   9.713  1.00  0.00           O  \n"
"ATOM      5  C5*   G     1     -11.850  -5.146   9.440  1.00  0.00           C  \n"
"ATOM      6 H5*1   G     1     -12.691  -5.536   9.997  1.00  0.00           H  \n"
"ATOM      7 H5*2   G     1     -12.064  -5.200   8.363  1.00  0.00           H  \n"
"ATOM      8  C4*   G     1     -11.678  -3.694   9.863  1.00  0.00           C  \n"
"ATOM      9  H4*   G     1     -12.581  -3.135   9.623  1.00  0.00           H  \n"
"ATOM     10  O4*   G     1     -11.419  -3.590  11.266  1.00  0.00           O  \n"
"ATOM     11  C3*   G     1     -10.481  -3.050   9.181  1.00  0.00           C  \n"
"ATOM     12  H3*   G     1      -9.675  -3.764   9.008  1.00  0.00           H  \n"
"ATOM     13  O3*   G     1     -10.924  -2.423   7.974  1.00  0.00           O  \n"
"ATOM     14  C2*   G     1     -10.085  -1.969  10.167  1.00  0.00           C  \n"
"ATOM     15  H2*   G     1      -9.021  -1.742  10.071  1.00  0.00           H  \n"
"ATOM     16  C1*   G     1     -10.357  -2.648  11.503  1.00  0.00           C  \n"
"ATOM     17  H1*   G     1     -10.681  -1.907  12.234  1.00  0.00           H  \n"
"ATOM     18  O2*   G     1     -10.882  -0.790  10.004  1.00  0.00           O  \n"
"ATOM     19 HO2'   G     1     -10.518  -0.117  10.584  1.00  0.00           H  \n"
"ATOM     20  N9    G     1      -9.145  -3.330  11.991  1.00  0.00           N  \n"
"ATOM     21  C8    G     1      -8.821  -4.645  11.961  1.00  0.00           C  \n"
"ATOM     22  H8    G     1      -9.493  -5.380  11.521  1.00  0.00           H  \n"
"ATOM     23  N7    G     1      -7.681  -4.993  12.454  1.00  0.00           N  \n"
"ATOM     24  C5    G     1      -7.164  -3.763  12.870  1.00  0.00           C  \n"
"ATOM     25  C4    G     1      -8.055  -2.738  12.589  1.00  0.00           C  \n"
"ATOM     26  C6    G     1      -5.926  -3.463  13.497  1.00  0.00           C  \n"
"ATOM     27  N3    G     1      -7.907  -1.418  12.839  1.00  0.00           N  \n"
"ATOM     28  C2    G     1      -6.740  -1.143  13.429  1.00  0.00           C  \n"
"ATOM     29  N1    G     1      -5.797  -2.102  13.744  1.00  0.00           N  \n"
"ATOM     30  H1    G     1      -4.939  -1.812  14.187  1.00  0.00           H  \n"
"ATOM     31  N2    G     1      -6.432   0.114  13.748  1.00  0.00           N  \n"
"ATOM     32  H21   G     1      -5.549   0.320  14.193  1.00  0.00           H  \n"
"ATOM     33  H22   G     1      -7.081   0.860  13.544  1.00  0.00           H  \n"
"ATOM     34  O6    G     1      -5.036  -4.227  13.814  1.00  0.00           O  \n"
"END";

    const char* uudPdbString2 = ""
"ATOM      5  C5*   G     1       0.000   1.000   0.000  1.00  0.00           C  \n"
"ATOM      8  C4*   G     1       0.000   0.000   0.000  1.00  0.00           C  \n"
"ATOM     10  O4*   G     1       1.000   0.000   0.000  1.00  0.00           O  \n"
"ATOM     11  C3*   G     1       0.000   0.000   1.000  1.00  0.00           C  \n"
"END";

    // std::ifstream inFileStream( "1UUD.pdb"); 
    std::istringstream inFileStream(uudPdbString);

    PdbStructure pdbStructure(inFileStream);
    RNA myRNA("G");// GGC");

    Compound::AtomTargetLocations atomTargets = myRNA.createAtomTargets(pdbStructure);

    //for (Compound::BondIndex bondIx(0); bondIx < myRNA.getNumBonds(); ++bondIx) 
    //{
    //    myRNA.setBondMobility(BondMobility::Rigid, bondIx);
    //}

    // Four steps to a perfect match
    myRNA.matchDefaultAtomChirality(atomTargets);
    myRNA.matchDefaultBondLengths(atomTargets);
    myRNA.matchDefaultBondAngles(atomTargets);
    myRNA.matchDefaultDihedralAngles(atomTargets);
    myRNA.matchDefaultTopLevelTransform(atomTargets);

    Real residual = myRNA.getTransformAndResidual(atomTargets).residual;
    Transform transform = myRNA.getTransformAndResidual(atomTargets).transform;

    cout << "Transform = " << transform << endl;

    cout << "Residual = " << residual << endl;

    myRNA.writeDefaultPdb(std::cout);

    // Make sure H3* atom is within 0.1 Angstrom of target location
    Vec3 finalLocH3 = myRNA.getTopLevelTransform() * myRNA.calcDefaultAtomLocationInCompoundFrame("0/H3*");
    Vec3 initialLocH3 = 0.10 * Vec3(-9.675,  -3.764,   9.008);
    Vec3 errVec = finalLocH3 - initialLocH3;
    Real error = std::sqrt( dot(errVec, errVec) );

    if (error < 0.01) 
    {
        cout << "PASSED" << endl;

        return 0;
    }

    else 
    {
        cout << "FAILED" << endl;    
        cerr << "Initial H3* location was " << initialLocH3 << endl;
        cerr << "Final H3* location was " << finalLocH3 << endl;

        return 1;
    }
}

