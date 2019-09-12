/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
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
#include "molmodel/internal/Pdb.h"

#include <sstream>

using namespace SimTK;
using namespace std;

int main() {

    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);
    dumm.loadAmber99Parameters();

try
 {
       String inputPdb=""
"ATOM      1  O5'   G B   1      -9.159 -10.595  12.585  1.00  0.00           O\n"
"ATOM      2  C5'   G B   1     -10.504 -10.828  12.144  1.00  0.00           C\n"
"ATOM      3  C4'   G B   1     -11.413  -9.637  12.402  1.00  0.00           C\n"
"ATOM      4  O4'   G B   1     -11.529  -9.358  13.803  1.00  0.00           O\n"
"ATOM      5  C3'   G B   1     -10.853  -8.372  11.768  1.00  0.00           C\n"
"ATOM      6  O3'   G B   1     -11.450  -8.222  10.474  1.00  0.00           O\n"
"ATOM      7  C2'   G B   1     -11.366  -7.276  12.681  1.00  0.00           C\n"
"ATOM      8  O2'   G B   1     -12.692  -6.870  12.323  1.00  0.00           O\n"
"ATOM      9  C1'   G B   1     -11.337  -7.955  14.044  1.00  0.00           C\n"
"ATOM     10  N9    G B   1     -10.055  -7.715  14.737  1.00  0.00           N\n"
"ATOM     11  C8    G B   1      -9.097  -8.600  15.109  1.00  0.00           C\n"
"ATOM     12  N7    G B   1      -8.062  -8.138  15.727  1.00  0.00           N\n"
"ATOM     13  C5    G B   1      -8.349  -6.771  15.783  1.00  0.00           C\n"
"ATOM     14  C6    G B   1      -7.594  -5.706  16.345  1.00  0.00           C\n"
"ATOM     15  O6    G B   1      -6.520  -5.759  16.912  1.00  0.00           O\n"
"ATOM     16  N1    G B   1      -8.239  -4.486  16.187  1.00  0.00           N\n"
"ATOM     17  C2    G B   1      -9.460  -4.306  15.567  1.00  0.00           C\n"
"ATOM     18  N2    G B   1      -9.911  -3.053  15.508  1.00  0.00           N\n"
"ATOM     19  N3    G B   1     -10.177  -5.303  15.036  1.00  0.00           N\n"
"ATOM     20  C4    G B   1      -9.568  -6.502  15.177  1.00  0.00           C\n"
"ATOM     21  H5'   G B   1     -10.920 -11.665  12.687  1.00  0.00           H\n"
"ATOM     22 H5''   G B   1     -10.482 -11.059  11.069  1.00  0.00           H\n"
"ATOM     23  H4'   G B   1     -12.402  -9.841  11.995  1.00  0.00           H\n"
"ATOM     24  H3'   G B   1      -9.764  -8.388  11.720  1.00  0.00           H\n"
"ATOM     25 HO2'   G B   1     -12.766  -6.945  11.369  1.00  0.00           H\n"
"ATOM     26  H1'   G B   1     -12.156  -7.577  14.656  1.00  0.00           H\n"
"ATOM     27  H8    G B   1      -9.198  -9.663  14.894  1.00  0.00           H\n"
"ATOM     28  H1    G B   1      -7.762  -3.680  16.564  1.00  0.00           H\n"
"ATOM     29  H21   G B   1      -9.367  -2.301  15.908  1.00  0.00           H\n"
"ATOM     30  H22   G B   1     -10.795  -2.855  15.063  1.00  0.00           H\n"
"ATOM     31 HO5'   G B   1      -8.573 -10.999  11.941  1.00  0.00           H\n"
"ATOM     32  P     G B   2     -10.655  -7.494   9.276  1.00  0.00           P\n"
"ATOM     33  OP1   G B   2     -11.450  -7.642   8.036  1.00  0.00           O\n"
"ATOM     34  OP2   G B   2      -9.248  -7.952   9.308  1.00  0.00           O\n"
"ATOM     35  O5'   G B   2     -10.689  -5.944   9.713  1.00  0.00           O\n"
"ATOM     36  C5'   G B   2     -11.850  -5.146   9.440  1.00  0.00           C\n"
"ATOM     37  C4'   G B   2     -11.678  -3.694   9.863  1.00  0.00           C\n"
"ATOM     38  O4'   G B   2     -11.419  -3.590  11.266  1.00  0.00           O\n"
"ATOM     39  C3'   G B   2     -10.481  -3.050   9.181  1.00  0.00           C\n"
"ATOM     40  O3'   G B   2     -10.924  -2.423   7.974  1.00  0.00           O\n"
"ATOM     41  C2'   G B   2     -10.085  -1.969  10.167  1.00  0.00           C\n"
"ATOM     42  O2'   G B   2     -10.882  -0.790  10.004  1.00  0.00           O\n"
"ATOM     43  C1'   G B   2     -10.357  -2.648  11.503  1.00  0.00           C\n"
"ATOM     44  N9    G B   2      -9.145  -3.330  11.991  1.00  0.00           N\n"
"ATOM     45  C8    G B   2      -8.821  -4.645  11.961  1.00  0.00           C\n"
"ATOM     46  N7    G B   2      -7.681  -4.993  12.454  1.00  0.00           N\n"
"ATOM     47  C5    G B   2      -7.164  -3.763  12.870  1.00  0.00           C\n"
"ATOM     48  C6    G B   2      -5.926  -3.463  13.497  1.00  0.00           C\n"
"ATOM     49  O6    G B   2      -5.036  -4.227  13.814  1.00  0.00           O\n"
"ATOM     50  N1    G B   2      -5.797  -2.102  13.744  1.00  0.00           N\n"
"ATOM     51  C2    G B   2      -6.740  -1.143  13.429  1.00  0.00           C\n"
"ATOM     52  N2    G B   2      -6.432   0.114  13.748  1.00  0.00           N\n"
"ATOM     53  N3    G B   2      -7.907  -1.418  12.839  1.00  0.00           N\n"
"ATOM     54  C4    G B   2      -8.055  -2.738  12.589  1.00  0.00           C\n"
"ATOM     55  H5'   G B   2     -12.691  -5.536   9.997  1.00  0.00           H\n"
"ATOM     56 H5''   G B   2     -12.064  -5.200   8.363  1.00  0.00           H\n"
"ATOM     57  H4'   G B   2     -12.581  -3.135   9.623  1.00  0.00           H\n"
"ATOM     58  H3'   G B   2      -9.675  -3.764   9.008  1.00  0.00           H\n"
"ATOM     59  H2'   G B   2      -9.021  -1.742  10.071  1.00  0.00           H\n"
"ATOM     60 HO2'   G B   2     -10.518  -0.117  10.584  1.00  0.00           H\n"
"ATOM     61  H1'   G B   2     -10.681  -1.907  12.234  1.00  0.00           H\n"
"ATOM     62  H8    G B   2      -9.493  -5.380  11.521  1.00  0.00           H\n"
"ATOM     63  H1    G B   2      -4.939  -1.812  14.187  1.00  0.00           H\n"
"ATOM     64  H21   G B   2      -5.549   0.320  14.193  1.00  0.00           H\n"
"ATOM     65  H22   G B   2      -7.081   0.860  13.544  1.00  0.00           H\n"
"ATOM     66  P     C B   3      -9.944  -2.320   6.701  1.00  0.00           P\n"
"ATOM     67  OP1   C B   3     -10.708  -1.742   5.573  1.00  0.00           O\n"
"ATOM     68  OP2   C B   3      -9.258  -3.622   6.537  1.00  0.00           O\n"
"ATOM     69  O5'   C B   3      -8.854  -1.235   7.180  1.00  0.00           O\n"
"ATOM     70  C5'   C B   3      -9.116   0.164   7.039  1.00  0.00           C\n"
"ATOM     71  C4'   C B   3      -7.921   1.010   7.472  1.00  0.00           C\n"
"ATOM     72  O4'   C B   3      -7.587   0.777   8.845  1.00  0.00           O\n"
"ATOM     73  C3'   C B   3      -6.673   0.649   6.681  1.00  0.00           C\n"
"ATOM     74  O3'   C B   3      -6.580   1.521   5.552  1.00  0.00           O\n"
"ATOM     75  C2'   C B   3      -5.552   0.993   7.643  1.00  0.00           C\n"
"ATOM     76  O2'   C B   3      -5.209   2.382   7.575  1.00  0.00           O\n"
"ATOM     77  C1'   C B   3      -6.164   0.634   8.992  1.00  0.00           C\n"
"ATOM     78  N1    C B   3      -5.806  -0.744   9.381  1.00  0.00           N\n"
"ATOM     79  C2    C B   3      -4.572  -0.943   9.986  1.00  0.00           C\n"
"ATOM     80  O2    C B   3      -3.819  -0.005  10.168  1.00  0.00           O\n"
"ATOM     81  N3    C B   3      -4.233  -2.208  10.357  1.00  0.00           N\n"
"ATOM     82  C4    C B   3      -5.067  -3.235  10.145  1.00  0.00           C\n"
"ATOM     83  N4    C B   3      -4.704  -4.460  10.526  1.00  0.00           N\n"
"ATOM     84  C5    C B   3      -6.338  -3.035   9.520  1.00  0.00           C\n"
"ATOM     85  C6    C B   3      -6.665  -1.780   9.156  1.00  0.00           C\n"
"ATOM     86  H5'   C B   3      -9.978   0.427   7.652  1.00  0.00           H\n"
"ATOM     87 H5''   C B   3      -9.343   0.379   5.995  1.00  0.00           H\n"
"ATOM     88  H4'   C B   3      -8.150   2.065   7.332  1.00  0.00           H\n"
"ATOM     89  H3'   C B   3      -6.659  -0.403   6.395  1.00  0.00           H\n"
"ATOM     90  H2'   C B   3      -4.681   0.364   7.446  1.00  0.00           H\n"
"ATOM     91 HO2'   C B   3      -5.283   2.738   8.463  1.00  0.00           H\n"
"ATOM     92  H1'   C B   3      -5.805   1.330   9.750  1.00  0.00           H\n"
"ATOM     93  H41   C B   3      -3.808  -4.607  10.967  1.00  0.00           H\n"
"ATOM     94  H42   C B   3      -5.325  -5.242  10.374  1.00  0.00           H\n"
"ATOM     95  H5    C B   3      -7.018  -3.868   9.344  1.00  0.00           H\n"
"ATOM     96  H6    C B   3      -7.623  -1.590   8.672  1.00  0.00           H\n"
"END\n";


     std::istringstream inStream(inputPdb);
     PdbStructure pdbStructure(inStream);
     RNA myRNA("GGC");
     myRNA.setPdbChainId("B");
     myRNA.writeDefaultPdb(cout); // OK

     Compound::AtomTargetLocations atomTargets = myRNA.createAtomTargets(pdbStructure);

     cout << atomTargets.size() << endl;
     SimTK_ASSERT_ALWAYS( (atomTargets.size() == 95), "Wrong number of atom matches" ); // pdb lacks H2', has HO5'
     cout << atomTargets.size() << endl;
     // assert(atomTargets.size() == 16);
     myRNA.writeDefaultPdb(cout);

     // Four steps to a perfect match
     myRNA.matchDefaultBondLengths(atomTargets);
     myRNA.matchDefaultBondAngles(atomTargets);
     myRNA.matchDefaultDihedralAngles(atomTargets);
     myRNA.matchDefaultTopLevelTransform(atomTargets);

     Real residual = myRNA.getTransformAndResidual(atomTargets).residual;
     cout << "residual = " << residual << " nanometers" << endl;

     myRNA.assignBiotypes();

     myRNA.writeDefaultPdb(cout);
     system.adoptCompound(myRNA);
     system.modelCompounds();
     system.realizeTopology();

     const State& state = system.getDefaultState();
     VerletIntegrator study(system);

     study.initialize(state);
     study.setAccuracy(0.01);

     TimeStepper myTimeStepper(system,study);

     Real mytime=0;
     Real interval = 0.02; //picoseconds
     for (int k=0; k<3; k++) {

       cout<< "step "<<k<<endl;
       mytime = k*interval;
       while (myTimeStepper.getTime() < mytime)
           myTimeStepper.stepTo(mytime);
       cout << "check 9"<<endl;

       myRNA.writePdb(myTimeStepper.getState(), cout);

     }

     cout << "PASSED" << endl;
     return 0;
  }
catch (const std::exception& e)
 {
   printf("EXCEPTION THROWN: %s\n", e.what());
 }
catch (...)
 {
   printf("UNKNOWN EXCEPTION THROWN\n");
 }    return 0;
}
