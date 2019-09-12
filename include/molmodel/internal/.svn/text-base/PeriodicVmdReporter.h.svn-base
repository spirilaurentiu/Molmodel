/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#ifndef SimTK_MOLMODEL_PERIODICVMDREPORTER_H_
#define SimTK_MOLMODEL_PERIODICVMDREPORTER_H_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/VmdConnection.h"

namespace SimTK {

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.
class SimTK_MOLMODEL_EXPORT PeriodicVmdReporter : public PeriodicEventReporter {
public:
    PeriodicVmdReporter(
        const CompoundSystem& system, 
        Real interval,
        int localSocketNumber,
        bool blockWaitingForVmdConnection = false
        ) 
        : PeriodicEventReporter(interval), 
          system(system), 
          vmdConnection(localSocketNumber),
          blockWaitingForVmdConnection(blockWaitingForVmdConnection)
    {}


    void handleEvent(const State& state) const;

private:
    const CompoundSystem& system;
    mutable VmdConnection vmdConnection;
    bool blockWaitingForVmdConnection;
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_PERIODICVMDREPORTER_H_
