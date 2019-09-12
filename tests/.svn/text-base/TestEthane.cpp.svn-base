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

#include <iostream>
#include <fstream>

using namespace SimTK;

void testConstructArgon() {
    Argon argon;

    // std::ofstream os("C:/test_argon.pdb");
    std::ostream& os = std::cout;

    argon.writeDefaultPdb(os);
}

void testConstructMethane() {
     Methane m;

     // std::ofstream os("C:/test_methane.pdb");
     std::ostream& os = std::cout;

     m.writeDefaultPdb(os);
}

void testConstructEthane() {
     Ethane e;

     // std::ofstream os("C:/test_ethane.pdb");
     std::ostream& os = std::cout;

     e.writeDefaultPdb(os);

     e.setDefaultTorsionAngle(10*Deg2Rad);
     e.writeDefaultPdb(os);

     e.setDefaultTorsionAngle(20*Deg2Rad);
     e.writeDefaultPdb(os);

     e.setDefaultTorsionAngle(30*Deg2Rad);
     e.writeDefaultPdb(os);

     e.setDefaultTorsionAngle(60*Deg2Rad);
     e.writeDefaultPdb(os);
}

void testConstructSerine() {
    AminoAcidResidue::Serine ser;

    // std::ofstream os("C:/test_serine.pdb");
    std::ostream& os = std::cout;

    ser.writeDefaultPdb(os);
}

void testConstructProtein() {
    Protein p("ACDEFGHIKLMNPQRSTVWY");
    // Protein p("AAPAA");

    // std::ofstream os("C:/test_protein.pdb");
    std::ostream& os = std::cout;

    p.writeDefaultPdb(os);
}

void testConstructRna() {
    RibonucleotideResidue::Adenylate a;

    // std::ofstream os("C:/test_rna.pdb");
    std::ostream& os = std::cout;

    a.writeDefaultPdb(os);
}

int main() {
try
  { 
      testConstructArgon();
      testConstructMethane();
      testConstructEthane();
      testConstructSerine();
      testConstructProtein();
      testConstructRna();
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

