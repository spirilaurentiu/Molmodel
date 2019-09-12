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

#ifndef SimTK_RIBOSE_MOBILIZER_H_
#define SimTK_RIBOSE_MOBILIZER_H_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/units.h"
#include <cmath>

namespace SimTK {

/// Function f(x) = 0 for all x
class ZeroFunction : public Function::Constant {
public:
	ZeroFunction() : Function::Constant(0, 0) {}
};

/// Implements a simple functional relationship, y = amplitude * sin(x - phase)
class SimTK_MOLMODEL_EXPORT SinusoidFunction : public Function {
private:
	Angle amplitude;
	Angle phase;
public:

	//Default constructor
	SinusoidFunction()
		: amplitude(1.0), phase(0.0*Deg2Rad) {}
	
	//Convenience constructor to specify the slope and Y-intercept of the linear r
	SinusoidFunction(Angle amp, Angle phi)
	: amplitude(amp), phase(phi) {}
	
	Real calcValue(const Vector& x) const{
		return amplitude*std::sin(x[0] - phase);
	}
	
	Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const{
		Real deriv = 0;
	
		assert(1 == x.size());
		
		// Derivatives 1, 5, 9, 13, ... are cos()
		if      ( 1 == (derivComponents.size() % 4) ){
			deriv = amplitude*std::cos(x[0] - phase);
		}
		// Derivatives 2, 6, 10, 14, ... are -sin()
		else if ( 2 == (derivComponents.size() % 4) ){
			deriv = -amplitude*std::sin(x[0] - phase);
		}
		// Derivatives 3, 7, 11, 15, ... are -cos()
		else if ( 3 == (derivComponents.size() % 4) ){
			deriv = -amplitude*std::cos(x[0] - phase);
		}
		// Derivatives 0, 4, 8, 12, ... are sin()
		else if ( 0 == (derivComponents.size() % 4) ){
			deriv = amplitude*std::sin(x[0] - phase);
		}
		else assert(false);
		
		return deriv;
	}
	
	int getArgumentSize() const{
		return 1;
	}
	
	int getMaxDerivativeOrder() const{
		return 1000;
	}
};
	
class SimTK_MOLMODEL_EXPORT PseudorotationMobilizer : public MobilizedBody::FunctionBased
{
public:
	PseudorotationMobilizer(
			MobilizedBody& parent,
			const Transform& inbFrame,
			const Body& body,
			const Transform& outbFrame,
			Angle amplitude,
			Angle phase
			)
	: MobilizedBody::FunctionBased(
			parent,
			inbFrame,
			body,
			outbFrame,
			1,
			createFunctions(amplitude, phase),
			createCoordIndices() // TODO ask Ajay
			)
	{}	
	
private:
	
    static std::vector< const Function* > createFunctions(Angle amplitude, Angle phase);
    static std::vector< std::vector<int> > createCoordIndices();	
	
};


class SimTK_MOLMODEL_EXPORT RiboseNu3Mobilizer : public PseudorotationMobilizer
{
public:
	RiboseNu3Mobilizer(
			MobilizedBody& parent,
			const Transform& inbFrame,
			const Body& body,
			const Transform& outbFrame)
	: PseudorotationMobilizer(
			parent,
			inbFrame,
			body,
			outbFrame,
			36.4*Deg2Rad,
			124.8*Deg2Rad)
	{}
			
};

} // namespace SimTK

#endif
