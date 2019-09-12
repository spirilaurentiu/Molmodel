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

/// Ensure that pdb coordinates of single atoms match those defined when the atom is placed

int main() 
{
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);

    dumm.loadAmber99Parameters();

	MagnesiumIon m1;
	MagnesiumIon::setAmberLikeParameters(dumm);
	ChlorideIon c1;
	ChlorideIon::setAmberLikeParameters(dumm);

	Vec3 v1(0.73, -4.62, 0.001);
	system.adoptCompound(m1, v1);

	Vec3 v2(2, 0.00, -3.5);
	system.adoptCompound(c1, Vec3(2, 0.00, -3.5));


    system.modelCompounds();  

    State state = system.realizeTopology();
    RungeKuttaMersonIntegrator study(system);

    study.initialize(state);

	double x, y, z;

	ostringstream s4;
	c1.writePdb(study.getState(), s4);
    std::istringstream(s4.str().substr(30, 8)) >> x;
    std::istringstream(s4.str().substr(38, 8)) >> y;
    std::istringstream(s4.str().substr(46, 8)) >> z;
	SimTK_ASSERT_ALWAYS( (Vec3(x, y, z) - 10*v2).scalarNormSqr() < 0.0001 ,
		"Output PDB coordinates did not match input");

	ostringstream s3;
	c1.writeDefaultPdb(s3);
    std::istringstream(s3.str().substr(30, 8)) >> x;
    std::istringstream(s3.str().substr(38, 8)) >> y;
    std::istringstream(s3.str().substr(46, 8)) >> z;
	SimTK_ASSERT_ALWAYS( (Vec3(x, y, z) - 10*v2).scalarNormSqr() < 0.0001 ,
		"Output PDB coordinates did not match input");

	ostringstream s2;
	m1.writePdb(study.getState(), s2);
    std::istringstream(s2.str().substr(30, 8)) >> x;
    std::istringstream(s2.str().substr(38, 8)) >> y;
    std::istringstream(s2.str().substr(46, 8)) >> z;
	SimTK_ASSERT_ALWAYS( (Vec3(x, y, z) - 10*v1).scalarNormSqr() < 0.0001 ,
		"Output PDB coordinates did not match input");

	ostringstream s1;
	m1.writeDefaultPdb(s1);
    std::istringstream(s1.str().substr(30, 8)) >> x;
    std::istringstream(s1.str().substr(38, 8)) >> y;
    std::istringstream(s1.str().substr(46, 8)) >> z;
	SimTK_ASSERT_ALWAYS( (Vec3(x, y, z) - 10*v1).scalarNormSqr() < 0.0001 ,
		"Output PDB coordinates did not match input");

}

