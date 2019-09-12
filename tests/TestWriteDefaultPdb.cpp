
#include "SimTKmolmodel.h"
#include "molmodel/internal/Pdb.h"

#include <sstream>

using namespace SimTK;
using namespace std;

int main() 
{      
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
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

    // std::ifstream inFileStream( "1UUD.pdb"); 
    std::istringstream inFileStream(uudPdbString);

    PdbStructure pdbStructure(inFileStream);
    RNA myRNA("G");// GGC");

    Compound::AtomTargetLocations atomTargets = myRNA.createAtomTargets(pdbStructure);

    // Four steps to a perfect match
    myRNA.matchDefaultBondLengths(atomTargets);
    myRNA.matchDefaultBondAngles(atomTargets);
    myRNA.matchDefaultDihedralAngles(atomTargets);
    myRNA.matchDefaultTopLevelTransform(atomTargets);

    Real residual = myRNA.getTransformAndResidual(atomTargets).residual;
    Transform transform = myRNA.getTransformAndResidual(atomTargets).transform;

    cout << "Transform = " << transform << endl;

    cout << "Residual = " << residual << endl;

    myRNA.writeDefaultPdb(std::cout);

    std::ostringstream oldWayOutString;
    myRNA.writeDefaultPdb(oldWayOutString);
    oldWayOutString << "END" << endl;

    std::ostringstream newWayOutString;
    PdbStructure pdbOut(myRNA);
    pdbOut.write(newWayOutString);

    // Dynamic configuration
    myRNA.assignBiotypes();
    system.adoptCompound(myRNA);
    system.modelCompounds(); // finalize multibody system
    State state = system.realizeTopology(); 
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    std::ostringstream statelyOldPdbOutString;
    myRNA.writePdb(ts.getState(), statelyOldPdbOutString);
    statelyOldPdbOutString << "END" << endl;

    std::ostringstream statelyNewPdbOutString;
    PdbStructure statelyPdbOut(ts.getState(), myRNA);
    statelyPdbOut.write(statelyNewPdbOutString);


    if ( oldWayOutString.str() != newWayOutString.str() ) 
    {
        cout << "FAILED" << endl;    

        cout << "OLD WAY ********************" << endl;
        cout << oldWayOutString.str() << endl;

        cout << "NEW WAY ********************" << endl;
        cout << newWayOutString.str() << endl;

        return 1;
    }

    else if ( statelyOldPdbOutString.str() != statelyNewPdbOutString.str() ) 
    {
        cout << "FAILED" << endl;    

        cout << "OLD WAY ********************" << endl;
        cout << statelyOldPdbOutString.str() << endl;

        cout << "NEW WAY ********************" << endl;
        cout << statelyNewPdbOutString.str() << endl;

        return 1;
    }

    else
    {
        cout << "PASSED" << endl; 

        return 0;
    }
}
