
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef __RealSimTk_H_
#define __RealSimTk_H__

// Set RealOpenMMType to 2 for double precision, 1 for float

#ifndef RealOpenMMType
#define RealOpenMMType 2
#endif 

#if RealOpenMMType == 1 

#define RealOpenMM     float
#define SQRT           std::sqrtf
#define POW            std::powf
#define SIN            std::sinf
#define COS            std::cosf
#define TAN            std::tanf

// LOG is used in Vishal's gpu code; modifying LOG -> LN 
#define LN             std::logf

#define EXP            std::expf
#define FABS           std::fabsf
#define ACOS           std::acosf
#define ASIN           std::asinf
#define ATAN           std::atanf
#define TANH           std::tanhf

#define ATOF           std::atoff

#define PI_M           3.141592653589f
#define TWO_SIX        1.122462048309372981f
#define RADIAN        57.29577951308f
#define LOG_TEN        2.302585092994045684f
#define SQRT_TWO       1.41421356237309504f

#else

#define RealOpenMM     double
#define SQRT           std::sqrt
#define POW            std::pow
#define SIN            std::sin
#define COS            std::cos
#define TAN            std::tan

// LOG is used in Vishal's gpu code; modifying LOG -> LN 
#define LN             std::log

#define EXP            std::exp
#define FABS           std::fabs
#define ACOS           std::acos
#define ASIN           std::asin
#define ATAN           std::atan
#define TANH           std::tanh

#define ATOF           std::atof

#define PI_M           3.141592653589
#define TWO_SIX        1.122462048309372981
#define RADIAN        57.29577951308
#define LOG_TEN        2.302585092994045684
#define SQRT_TWO       1.41421356237309504
#define RADIAN_INVERSE 0.01745329252

#endif

#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define MATRIXDOT3(u,v) u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + \
                        u[3]*v[3] + u[4]*v[4] + u[5]*v[5] + \
                        u[6]*v[6] + u[7]*v[7] + u[8]*v[8]


#endif // __RealSimTk_H__
