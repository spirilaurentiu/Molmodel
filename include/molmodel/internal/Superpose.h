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


/**@file
 *
 * Kabsch superposition algorithm implementation.
 */

#ifndef MOLMODEL_SUPERPOSE_H_
#define MOLMODEL_SUPERPOSE_H_

#include <array>
#include <vector>

#include "SimTKmath.h"


namespace SimTK {

/// Matched pair of 3D vectors to be used in least-squares superposition
class Vec3Pair {
    public:
    Vec3Pair() = default;

    Vec3Pair(const Vec3& source, const Vec3& target, Real weight = 1.0) {
        this->source = source;
        this->target = target;
        this->weight = weight;
    }

    [[nodiscard]] auto getSource() const -> const Vec3& {
        return source;
    }
    [[nodiscard]] auto getTarget() const -> const Vec3& {
        return target;
    }
    [[nodiscard]] auto getWeight() const -> Real {
        return weight;
    }

    private:
    Vec3 source;
    Vec3 target;
    Real weight;
};

class TransformAndResidual {
    public:
    TransformAndResidual() = default;

    TransformAndResidual(const Transform& transform, Real residual) {
        this->transform = transform;
        this->residual = residual;
    }

    Transform transform;
    Real residual;
};

class Kabsch78 {
    public:
    using VectorSet = std::vector<Vec3Pair>;

    struct AtomSet {
        std::vector<SimTK::Real> sourceX;
        std::vector<SimTK::Real> sourceY;
        std::vector<SimTK::Real> sourceZ;

        std::vector<SimTK::Real> targetX;
        std::vector<SimTK::Real> targetY;
        std::vector<SimTK::Real> targetZ;

        int n;
    };

    /**
     * Compute the transformation that orients the first (source) set of vectors
     * to match as closely as possible the second (target) set.
     * Using a weighted least-squares criterion.
     *
     * \return Transform that, when applied to source vectors, minimizes weighted
     * least-squares residual with respect to the target vectors.
     */
    static auto superpose(const VectorSet& vectors) -> TransformAndResidual;

    static auto superpose_unweighted(const AtomSet& vectors) -> TransformAndResidual;

    private:
    struct Centroids {
        Vec3 sourceCentroid;
        Vec3 targetCentroid;
    };

    static void accumulate(const AtomSet& atoms, Real& E0_out, Mat<3, 3>& R_out, Centroids& centroids);

    // Drop-in replacement for the Eigen call + sort block.
    // RtR must be symmetric PSD (which R^T*R always is).
    // Returns mu[0] >= mu[1] >= mu[2] >= 0 and unit eigenvectors a[0..2].
    // Already sorted - remove the sort block that follows in the original.
    static void eigenSymm3(const Mat<3, 3>& RtR, std::array<Real, 3>& mu, std::array<Vec3, 3>& a);
};

} // namespace SimTK

#endif // MOLMODEL_SUPERPOSE_H_