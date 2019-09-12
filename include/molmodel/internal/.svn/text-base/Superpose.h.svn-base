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

#include "SimTKmath.h"
#include <vector>

namespace SimTK {

/// Matched pair of 3D vectors to be used in least-squares superposition
class Vec3Pair 
{
public:
    Vec3Pair(const Vec3& s, const Vec3& t, Real w = 1.0)
        : source(s), target(t), weight(w) {}

    const Vec3& getSource() const {return source;}
    const Vec3& getTarget() const {return target;}
    Real getWeight() const {return weight;}

private:
    Vec3 source;
    Vec3 target;
    Real weight;
};

class TransformAndResidual {
public:
    TransformAndResidual(const Transform& t, Real r)
        : transform(t), residual(r) {}

    Transform transform;
    Real residual;
};

class SimTK_MOLMODEL_EXPORT Kabsch78 {
public:
    typedef std::vector<Vec3Pair> VectorSet ;

    /**
     * Compute the transformation that orients the first (source) set of vectors
     * to match as closely as possible the second (target) set.
     * Using a weighted least-squares criterion.
     *
     * \return Transform that, when applied to source vectors, minimizes weighted least-squares residual with respect to the target vectors.
     */
    static TransformAndResidual superpose(const VectorSet& vectors) 
    {
        // a) Remove any translation between the two given
        // vector sets x(n) and y(n), and determine
        // E0 = 1/2 SUM[ w(n)*(x(n)^2 + y(n)^2) ]
        // and R, r(ij) = SUM(n)[ w(n)*y(ni)*x(nj) ]

        // 1) Identify center of mass of each set of vectors
        Real totalMass = 0.0;
        Vec3 sourceCentroid(0.0, 0.0, 0.0);
        Vec3 targetCentroid(0.0, 0.0, 0.0);
        VectorSet::const_iterator vI;
        for (vI = vectors.begin(); vI != vectors.end(); ++vI)
        {
            totalMass += vI->getWeight();
            sourceCentroid += vI->getWeight() * vI->getSource();
            targetCentroid += vI->getWeight() * vI->getTarget();
        }
        if (totalMass != 0.0) {
            sourceCentroid /= totalMass;
            targetCentroid /= totalMass;
        }

        // Form R matrix from Kabsch paper
        Mat<3,3> R(0.0);
        Real E0 = 0.0; // Initial residual, see Kabsch
        for (vI = vectors.begin(); vI != vectors.end(); ++vI)
        {
            Vec3 x = vI->getSource() - sourceCentroid;
            Vec3 y = vI->getTarget() - targetCentroid;
            E0 += 0.5 * vI->getWeight() * ( dot(x, x) + dot(y, y) );
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                      R[i][j] += vI->getWeight() * y[i] * x[j];
        }

        // b) Form ~RR, determine eigenvalues mu(k) and the
        // mutually orthogonal eigenvectors a(k) and
        //    sort so that mu1 >= mu2 >= mu3.
        // Set
        //    a3 == a1 cross a2
        // to be sure to have a right handed system
        Vector_<std::complex<double> > mu0; // eigenvalues, complex, unsorted
        Matrix_<std::complex<double> > a0; // eigenvectors, complex, unsorted
        Eigen eigen( Matrix( R.transpose() * R ) );
        eigen.getAllEigenValuesAndVectors(mu0, a0);

        // use only real component of results
        Vec3 a1[3]; // eigenvectors, real, unsorted
        Real mu1[3]; // eigenvalues, real, unsorted
        for (int i = 0; i < 3; ++i) {
            mu1[i] = mu0[i].real();
            for (int j = 0; j < 3; ++j) 
                // Note swapping of indices: It appears that the columns of a0 are eigenvectors
                a1[j][i] = a0[i][j].real();
        }

        // sort indices of eigenvalues, from largest eigenvalue to smallest
        int indMax(0);
        int indMin(0);
        for (int i = 0; i < 3; ++i) {
            if ( mu1[i] > mu1[indMax] ) {
                indMax = i;
            }
            if ( mu1[i] <= mu1[indMin] ) {
                indMin = i;
            }
        }
        assert (indMin != indMax);
        int indMid = 3 - indMax - indMin; // too clever...
        assert(indMid >= 0);
        assert(indMid <= 2);
        assert(indMid != indMax);
        assert(indMid != indMin);
        int indSort[3] = {indMax, indMid, indMin};

        Real mu[3]; // eigenvalues, sorted
        Vec3 a[3]; // eigenvectors, sorted
        for (int i = 0; i < 3; ++i) {
            mu[i] = mu1[indSort[i]];
            a[i] = Vec3(UnitVec3(a1[indSort[i]]));
        }

        a[2] = cross(a[0], a[1]); // force right handed system
        
        // c) Determine Ra(k) (k = 1,2,3), normalize the first
        // two vectors to obtain b1, b2, and set b3 == b1 cross b2.
        // This will also take care of the case mu2 > mu3 = 0.
        Real sigma[] = {1.0, 1.0, 1.0};
        Vec3 b[3];
        b[0] = Vec3(UnitVec3(R * a[0]));
        b[1] = Vec3(UnitVec3(R * a[1]));
        b[2] = cross(b[0], b[1]);

        if ( dot(b[2], R * a[2]) < 0 )
            sigma[2] = -1.0;

        // d) Form U according to eq. 7:
        // u(ij) = SUM(k)[ b(ki)a(kj) ]
        // where b(k) = R*a(k)/(sigma(k)sqrt(mu(k)))
        Mat<3,3> U(0.0);   // initialize 2-d array
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    U[i][j] += b[k][i] * a[k][j];
        
        // Compute residual error
        Real residualError = E0 - sigma[0] * std::sqrt(mu[0]) - sigma[1] * std::sqrt(mu[1]) - sigma[2] * std::sqrt(std::abs(mu[2]));
        // E = 1/2 SUM(over n)[w(n)*(Ux(n) - y(n))^2]
        // variance would be sqrt( 1/SUM(w(n)) * SUM(w(n)*(Ux(n) - y(n))^2) )
        Real variance = 0.0;
        if (totalMass > 0)
            variance = std::sqrt(2 * residualError / totalMass);

        Rotation rotation(U);

        // No rotation for single point overlay
        if (vectors.size() < 2)
            rotation = Rotation();

        Transform transform1(-sourceCentroid);
        Transform transform2(rotation);
        Transform transform3(targetCentroid);

        return TransformAndResidual(transform3 * transform2 * transform1, variance);
    }
};

} // namespace SimTK

#endif // MOLMODEL_SUPERPOSE_H_
