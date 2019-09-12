#include "SimTKmolmodel.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace std; 

void testFitMagnesiums() 
{
    const char* inPdbString =
    "ATOM      1 MG+2 MG  X   0     304.790-310.727 -20.189  1.00  0.00          MG\n"
    "ATOM      1 MG+2 MG  X   1     180.134 396.232  25.745  1.00  0.00          MG\n";

    CompoundSystem system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem forces(system);
    DuMMForceFieldSubsystem dumm(system);

    MagnesiumIon myMagnesiumIonVec[2];
    for (int i = 0; i < (2); i++) {
        myMagnesiumIonVec[i].setPdbChainId("X");
        myMagnesiumIonVec[i].setPdbResidueNumber(i);

        system.adoptCompound(myMagnesiumIonVec[i],Vec3(0,0,i));
    }

    istringstream previousFrameFile (inPdbString);
    assert(previousFrameFile.good());
    PdbStructure pdbStructure(previousFrameFile);
    pdbStructure.write(std::cout);
    for (int i = 0; i < (2); i++) {
        Compound::AtomTargetLocations magnesiumAtomTargets = myMagnesiumIonVec[i].createAtomTargets(pdbStructure); 
        cout << "[Repel.h] magnesiumAtomTargets: " << magnesiumAtomTargets[Compound::AtomIndex(0)]<<endl;
        cout << "[Repel.h] magnesiumAtomTargets: " << magnesiumAtomTargets.size()       <<endl;
        myMagnesiumIonVec[i].fitDefaultConfiguration(magnesiumAtomTargets,.3);
        //myMagnesiumIonVec[i].matchDefaultTopLevelTransform(magnesiumAtomTargets);

        myMagnesiumIonVec[i].writeDefaultPdb(std::cout);
    }
}

void testFitWonkyNucleotideChirality() 
{
    const char* input_pdb_string = "\
ATOM   2250  P   G   A  71       1.419  15.496 135.643  1.00  0.00           P\n\
ATOM   2251  OP1 G   A  71       0.096  14.848 135.790  1.00  0.00           O\n\
ATOM   2252  OP2 G   A  71       1.885  15.911 134.300  1.00  0.00           O\n\
ATOM   2253  O5* G   A  71       1.456  16.786 136.608  1.00  0.00           O\n\
ATOM   2254  C5* G   A  71       2.351  16.848 137.710  1.00  0.00           C\n\
ATOM   2255 H5*1 G   A  71       3.397  16.761 137.343  1.00  0.00           H\n\
ATOM   2256 H5*2 G   A  71       2.133  16.013 138.412  1.00  0.00           H\n\
ATOM   2257  C4* G   A  71       2.177  18.165 138.426  1.00  0.00           C\n\
ATOM   2258  H4* G   A  71       1.647  18.882 137.762  1.00  0.00           H\n\
ATOM   2259  O4* G   A  71       1.399  17.960 139.635  1.00  0.00           O\n\
ATOM   2260  C3* G   A  71       3.535  18.671 138.898  1.00  0.00           C\n\
ATOM   2261  H3* G   A  71       4.292  17.863 138.783  1.00  0.00           H\n\
ATOM   2262  O3* G   A  71       3.921  19.795 138.116  1.00  0.00           O\n\
ATOM   2263  C2* G   A  71       3.250  19.118 140.327  1.00  0.00           C\n\
ATOM   2264  H2* G   A  71       4.178  19.033 140.934  1.00  0.00           H\n\
ATOM   2265  C1* G   A  71       2.246  18.060 140.781  1.00  0.00           C\n\
ATOM   2266  H1* G   A  71       1.669  17.694 139.905  1.00  0.00           H\n\
ATOM   2267  O2* G   A  71       2.624  20.384 140.308  1.00  0.00           O\n\
ATOM   2268 HO2' G   A  71       2.623  20.710 139.405  1.00  0.00           H\n\
ATOM   2269  N9  G   A  71       2.814  17.296 141.904  1.00  0.00           N\n\
ATOM   2270  C8  G   A  71       2.264  17.040 143.133  1.00  0.00           C\n\
ATOM   2271  H8  G   A  71       1.288  17.395 143.429  1.00  0.00           H\n\
ATOM   2272  N7  G   A  71       3.026  16.327 143.914  1.00  0.00           N\n\
ATOM   2273  C5  G   A  71       4.163  16.094 143.148  1.00  0.00           C\n\
ATOM   2274  C4  G   A  71       4.043  16.683 141.920  1.00  0.00           C\n\
ATOM   2275  C6  G   A  71       5.354  15.397 143.405  1.00  0.00           C\n\
ATOM   2276  N3  G   A  71       4.944  16.672 140.908  1.00  0.00           N\n\
ATOM   2277  C2  G   A  71       6.089  16.030 141.085  1.00  0.00           C\n\
ATOM   2278  N1  G   A  71       6.261  15.404 142.354  1.00  0.00           N\n\
ATOM   2279  H1  G   A  71       7.132  14.913 142.501  1.00  0.00           H\n\
ATOM   2280  N2  G   A  71       7.000  16.002 140.104  1.00  0.00           N\n\
ATOM   2281  H21 G   A  71       6.814  16.469 139.229  1.00  0.00           H\n\
ATOM   2282  H22 G   A  71       7.873  15.513 140.239  1.00  0.00           H\n\
ATOM   2283  O6  G   A  71       5.574  14.827 144.472  1.00  0.00           O\n";

    RNA rna("G");
    rna.renumberPdbResidues(71);
    rna.setPdbChainId("A");

    istringstream inString(input_pdb_string);
    Compound::AtomTargetLocations atomTargets = 
        rna.createAtomTargets(PdbStructure(inString));

    for (Compound::AtomIndex t(0); t < atomTargets.size(); ++t) {
        // cout << "atom " << t << ": " << rna.getAtomName(t) << endl;
    }

    cout << atomTargets.size() << " atom targets found" << endl;

    // Match with idealized geometry
    rna.matchDefaultConfiguration(atomTargets, Compound::Match_Idealized);
    cout << "Match idealized residual = ";
    cout << rna.getTransformAndResidual(atomTargets).residual;
    cout << " nanometers" << endl;

    // Match with exact geometry
    rna.matchDefaultConfiguration(atomTargets, Compound::Match_Exact);
    cout << "Match exact residual = ";
    cout << rna.getTransformAndResidual(atomTargets).residual;
    cout << " nanometers" << endl;

    cout << "top level transform = " << rna.getTopLevelTransform() << endl;

    // Write input and output coordinates for comparison
    ofstream pdb_output_stream("test.pdb");
    pdb_output_stream << input_pdb_string << endl;
    rna.renumberPdbResidues(72);
    rna.writeDefaultPdb(pdb_output_stream);
    pdb_output_stream.close();

    // Here is an automated test; "exact" residual should be less than 0.05 angstroms (0.005 nanometers)
    SimTK_ASSERT_ALWAYS(rna.getTransformAndResidual(atomTargets).residual < 0.005, "Residual too large");
}


int main() {
    try {
        testFitMagnesiums();
        testFitWonkyNucleotideChirality();
        cout << "PASSED" << endl;
        return 0;
    } catch (exception exc) {
        cout << "FAILED" << endl;
        return 1;
    }
}
