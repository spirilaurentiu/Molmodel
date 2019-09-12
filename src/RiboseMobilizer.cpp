/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors: Ajay Seth                                                    *
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

#include "molmodel/internal/RiboseMobilizer.h"
#include <cmath>

namespace SimTK {

/* static */ std::vector< const Function* > PseudorotationMobilizer::createFunctions(Angle amplitude, Angle phase) 
{
	std::vector< const Function* > functions;
	functions.push_back(new ZeroFunction ); // x rotation
	functions.push_back(new ZeroFunction ); // y rotation
	
	// Z-rotation is the only degree of freedom for this mobilizer
	functions.push_back(new SinusoidFunction(amplitude, phase) ); // z rotation
	
	functions.push_back(new ZeroFunction ); // x translation
	functions.push_back(new ZeroFunction ); // y translation
	functions.push_back(new ZeroFunction ); // z translation
	
	return functions;
}

/* static */ std::vector< std::vector<int> > PseudorotationMobilizer::createCoordIndices() 
{
	std::vector<int> oneCoordinateVec;
	oneCoordinateVec.push_back(0); // first and only generalized coordinate index
	
	std::vector< std::vector<int> > coordIndices;
	// All six functions take the one generalized coordinate
	coordIndices.push_back(std::vector<int>()); // 1, empty
	coordIndices.push_back(std::vector<int>()); // 2, empty
	coordIndices.push_back(oneCoordinateVec);   // 3, one coordinate, zero
	coordIndices.push_back(std::vector<int>()); // 4, empty
	coordIndices.push_back(std::vector<int>()); // 5, empty
	coordIndices.push_back(std::vector<int>()); // 6, empty
	
	return coordIndices;
}

} // namespace SimTK
