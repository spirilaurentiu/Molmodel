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

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

void testInputMatchesOutput() {
    // A very small protein
    String inputPdb = ""\
    "ATOM      1  N   HIS    19      28.165  29.227  23.618  1.00 91.78           N\n"
    "ATOM      2  CA  HIS    19      27.004  29.173  22.731  1.00 91.74           C\n"
    "ATOM      3  C   HIS    19      26.321  27.818  22.666  1.00 81.05           C\n"
    "ATOM      4  O   HIS    19      25.105  27.739  22.805  1.00 89.56           O\n"
    "ATOM      5  CB  HIS    19      27.248  29.629  21.270  1.00 98.21           C\n"
    "ATOM      6  CG  HIS    19      27.954  30.936  21.082  1.00100.00           C\n"
    "ATOM      7  ND1 HIS    19      28.852  31.112  20.015  1.00100.00           N\n"
    "ATOM      8  CD2 HIS    19      27.882  32.111  21.798  1.00100.00           C\n"
    "ATOM      9  CE1 HIS    19      29.310  32.368  20.116  1.00100.00           C\n"
    "ATOM     10  NE2 HIS    19      28.753  32.997  21.176  1.00100.00           N\n"
    "ATOM     11  N   SER    20      27.057  26.778  22.303  1.00 58.48           N\n"
    "ATOM     12  CA  SER    20      26.412  25.501  22.190  1.00 45.35           C\n"
    "ATOM     13  C   SER    20      26.089  24.909  23.558  1.00 30.40           C\n"
    "ATOM     14  O   SER    20      26.808  25.067  24.542  1.00 29.60           O\n"
    "ATOM     15  CB  SER    20      27.206  24.513  21.344  1.00 49.69           C\n"
    "ATOM     16  OG  SER    20      26.466  23.310  21.049  1.00 22.13           O\n"
    "TER\n"
    "END\n";

    istringstream inStream(inputPdb);
    PdbStructure pdbStructure(inStream);

    ostringstream outStream;
    pdbStructure.write(outStream);

    ASSERT(inputPdb == outStream.str());
}


// SimTK core bug 653
// Avoid distorting input geometry
void testMatchDefaultBreaksPlanarity() {
    // A very small RNA
    String pdbString = ""\
    "ATOM      1  O5' A       1       0.000   0.000   0.000  1.00  0.00           O\n"
    "ATOM      2  C5' A       1      -0.731   1.221   0.000  1.00  0.00           C\n"
    "ATOM      3  H5' A       1      -0.380   1.863  -0.837  1.00  0.00           H\n"
    "ATOM      4 H5'' A       1      -0.570   1.749   0.965  1.00  0.00           H\n"
    "ATOM      5  C4' A       1      -2.202   0.927  -0.173  1.00  0.00           C\n"
    "ATOM      6  H4' A       1      -2.775   1.880  -0.152  1.00  0.00           H\n"
    "ATOM      7  O4' A       1      -2.416   0.262  -1.447  1.00  0.00           O\n"
    "ATOM      8  C3' A       1      -2.649  -0.071   0.888  1.00  0.00           C\n"
    "ATOM      9  H3' A       1      -1.861  -0.843   1.026  1.00  0.00           H\n"
    "ATOM     10  O3' A       1      -2.867   0.610   2.118  1.00  0.00           O\n"
    "ATOM     11  C2' A       1      -3.982  -0.569   0.339  1.00  0.00           C\n"
    "ATOM     12  H2' A       1      -4.164  -1.607   0.693  1.00  0.00           H\n"
    "ATOM     13  C1' A       1      -3.707  -0.583  -1.164  1.00  0.00           C\n"
    "ATOM     14  H1' A       1      -4.620  -0.263  -1.711  1.00  0.00           H\n"
    "ATOM     15  O2' A       1      -4.991   0.378   0.624  1.00  0.00           O\n"
    "ATOM     16 HO2' A       1      -5.049   0.485   1.576  1.00  0.00           H\n"
    "ATOM     17  N9  A       1      -3.416  -1.965  -1.639  1.00  0.00           N\n"
    "ATOM     18  C8  A       1      -2.272  -2.678  -1.387  1.00  0.00           C\n"
    "ATOM     19  H8  A       1      -1.653  -2.510  -0.518  1.00  0.00           H\n"
    "ATOM     20  N7  A       1      -2.001  -3.578  -2.293  1.00  0.00           N\n"
    "ATOM     21  C5  A       1      -3.039  -3.450  -3.209  1.00  0.00           C\n"
    "ATOM     22  C4  A       1      -3.909  -2.469  -2.818  1.00  0.00           C\n"
    "ATOM     23  C6  A       1      -3.333  -4.123  -4.406  1.00  0.00           C\n"
    "ATOM     24  N3  A       1      -5.033  -2.063  -3.454  1.00  0.00           N\n"
    "ATOM     25  C2  A       1      -5.370  -2.657  -4.589  1.00  0.00           C\n"
    "ATOM     26  N1  A       1      -4.448  -3.699  -5.015  1.00  0.00           N\n"
    "ATOM      1  H2  A       1      -6.262  -2.352  -5.114  1.00  0.00           H\n"
    "ATOM     28  N6  A       1      -2.564  -5.104  -4.896  1.00  0.00           N\n"
    "ATOM     29  H61 A       1      -1.732  -5.389  -4.399  1.00  0.00           H\n"
    "ATOM     30  H62 A       1      -2.817  -5.558  -5.761  1.00  0.00           H\n"
    "ATOM     31  P   A       1       1.610   0.000   0.000  1.00  0.00           P\n"
    "ATOM     32  OP1 A       1       2.103   0.698   1.208  1.00  0.00           O\n"
    "ATOM     33  OP2 A       1       2.103   0.698  -1.208  1.00  0.00           O\n"
    "ATOM     34  OP3 A       1       2.103  -1.395   0.000  1.00  0.00           O\n";
    istringstream pdbStream(pdbString);
	PdbStructure pdbStructure(pdbStream);

    // 1) First reproduce undesired behavior
    RNA mol1("A");
	Compound::AtomTargetLocations atomTargets = mol1.createAtomTargets(pdbStructure); 

	mol1.matchDefaultAtomChirality(atomTargets);
    mol1.matchDefaultBondLengths(atomTargets);
    mol1.matchDefaultBondAngles(atomTargets);
    mol1.matchDefaultDihedralAngles(atomTargets, Compound::DistortPlanarBonds);
    mol1.matchDefaultTopLevelTransform(atomTargets);

    // Most distorted part in problem report is distance between atoms N1 and C2
    const ResidueInfo& res1 = mol1.getResidue( ResidueInfo::Index(0) );
    Vec3 atomN1Pos = mol1.calcDefaultAtomFrameInCompoundFrame(res1.getAtomIndex("N1")).p();
    Vec3 atomC2Pos = mol1.calcDefaultAtomFrameInCompoundFrame(res1.getAtomIndex("C2")).p();
    Real bondLength = (atomN1Pos - atomC2Pos).norm();

    ASSERT(bondLength > 0.20); // distorted


    // 2) Repair with extra parameter on matchDefaultAtomChirality

    RNA mol2("A");
	atomTargets = mol2.createAtomTargets(pdbStructure); 

	mol2.matchDefaultAtomChirality(atomTargets, 0.20);
    mol2.matchDefaultBondLengths(atomTargets);
    mol2.matchDefaultBondAngles(atomTargets);
    mol2.matchDefaultDihedralAngles(atomTargets, Compound::DistortPlanarBonds);
    mol2.matchDefaultTopLevelTransform(atomTargets);

    // Most distorted part in problem report is distance between atoms N1 and C2
    const ResidueInfo& res2 = mol2.getResidue( ResidueInfo::Index(0) );
    atomN1Pos = mol2.calcDefaultAtomFrameInCompoundFrame(res2.getAtomIndex("N1")).p();
    atomC2Pos = mol2.calcDefaultAtomFrameInCompoundFrame(res2.getAtomIndex("C2")).p();
    bondLength = (atomN1Pos - atomC2Pos).norm();

    cout << "bondLength = " << bondLength << endl;
    mol2.writeDefaultPdb(cout);

    ASSERT(bondLength < 0.20); // not distorted


    // 3) Repair by not setting torsion of bonded planar atoms

    RNA mol3("A");
	atomTargets = mol3.createAtomTargets(pdbStructure); 

	mol3.matchDefaultAtomChirality(atomTargets);
    mol3.matchDefaultBondLengths(atomTargets);
    mol3.matchDefaultBondAngles(atomTargets);
    mol3.matchDefaultDihedralAngles(atomTargets, Compound::KeepPlanarBonds);
    mol3.matchDefaultTopLevelTransform(atomTargets);

    // Most distorted part in problem report is distance between atoms N1 and C2
    const ResidueInfo& res3 = mol3.getResidue( ResidueInfo::Index(0) );
    atomN1Pos = mol3.calcDefaultAtomFrameInCompoundFrame(res3.getAtomIndex("N1")).p();
    atomC2Pos = mol3.calcDefaultAtomFrameInCompoundFrame(res3.getAtomIndex("C2")).p();
    bondLength = (atomN1Pos - atomC2Pos).norm();

    ASSERT(bondLength < 0.20); // not distorted
}

int main() {
    try {
        testInputMatchesOutput();
        testMatchDefaultBreaksPlanarity();
    }
    catch (const std::exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        return 1;
    }
    catch (...) {
        printf("UNKNOWN EXCEPTION THROWN\n");
        return 1;
    }    
    
    cout << "Done" << endl;
    return 0;
}

