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
try
  { 
      // TODO - this should be a very small protein
      // Use slight modifications of default configuration to test the import of state
      String inputPdb = ""
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
"END\n";


      std::istringstream inStream(inputPdb);
      PdbStructure pdbStructure(inStream);

      Protein protein("HS");
      protein.writeDefaultPdb(cout); // OK

      protein.updResidue(ResidueInfo::Index(1)).setPdbResidueNumber(19);
      protein.updResidue(ResidueInfo::Index(2)).setPdbResidueNumber(20);

      Compound::AtomTargetLocations atomTargets = protein.createAtomTargets(pdbStructure);

      SimTK_ASSERT_ALWAYS(atomTargets.size() == 16, "Wrong number of atoms used for matching");

      protein.writeDefaultPdb(cout);

      // Four steps to a perfect match
      protein.matchDefaultBondLengths(atomTargets);
      protein.matchDefaultBondAngles(atomTargets);
      protein.matchDefaultDihedralAngles(atomTargets);
      protein.matchDefaultTopLevelTransform(atomTargets);

      Real residual = protein.getTransformAndResidual(atomTargets).residual;
      cout << "residual = " << residual << " nanometers" << endl;

      SimTK_ASSERT_ALWAYS(residual < 0.02, "Structure matching was too inaccurate");

      protein.writeDefaultPdb(cout);

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

