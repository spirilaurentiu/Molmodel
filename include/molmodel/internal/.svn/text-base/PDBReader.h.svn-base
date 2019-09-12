#ifndef SimTK_MOLMODEL_PDB_READER_H_
#define SimTK_MOLMODEL_PDB_READER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "molmodel/internal/common.h"
#include "molmodel/internal/CompoundSystem.h"

namespace SimTK {

/**
 * This class parses PDB files, then constructs Systems and States based on them for use
 * in Simbody.
 */

class SimTK_MOLMODEL_EXPORT PDBReader {
public:
    /**
     * Create a PDBReader object that encapsulates the information in a PDB file.
     */
    explicit PDBReader(std::string filename);
    ~PDBReader();
    /**
     * Create one or more Compounds representing the protein or nucleic acid described in the PDB file.  You
     * pass to this method a CompoundSystem which has all the necessary subsystems.  It then creates
     * Compounds based on the sequence in the PDB file and adds them to the CompoundSystem.
     */
    void createCompounds(CompoundSystem& system);
    /**
     * Create a State reflecting the structure loaded from the PDB file.  Before calling this method,
     * you must first call createCompounds(), then call modelCompounds() on the CompoundSystem.  You
     * pass it the System that was created by createCompounds() and a State in which to store the result.
     * It performs a nonlinear optimization to find the set of internal coordinates for the System which
     * most closely match the structure specified in the PDB file.
     * 
     * Because this method actually has to solve an optimization problem, it may take a long time
     * for large proteins.
     */
    Real createState(const CompoundSystem& system, State& state) const;
private:
    class PDBReaderImpl;
    PDBReaderImpl* impl;
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_PDB_READER_H_
