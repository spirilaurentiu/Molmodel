#include "Superpose.h"

using namespace SimTK;

auto Kabsch78::superpose(const VectorSet& vectors) -> TransformAndResidual {
    // a) Remove any translation between the two given
    // vector sets x(n) and y(n), and determine
    // E0 = 1/2 SUM[ w(n)*(x(n)^2 + y(n)^2) ]
    // and R, r(ij) = SUM(n)[ w(n)*y(ni)*x(nj) ]

    // 1) Identify center of mass of each set of vectors
    Real totalMass = 0.0;
    Vec3 sourceCentroid(0.0, 0.0, 0.0);
    Vec3 targetCentroid(0.0, 0.0, 0.0);
    VectorSet::const_iterator vI;
    for (vI = vectors.begin(); vI != vectors.end(); ++vI) {
        totalMass += vI->getWeight();
        sourceCentroid += vI->getWeight() * vI->getSource();
        targetCentroid += vI->getWeight() * vI->getTarget();
    }
    if (totalMass != 0.0) {
        sourceCentroid /= totalMass;
        targetCentroid /= totalMass;
    }

    // Form R matrix from Kabsch paper
    Mat<3, 3> R(0.0);
    Real E0 = 0.0; // Initial residual, see Kabsch
    for (vI = vectors.begin(); vI != vectors.end(); ++vI) {
        Vec3 x = vI->getSource() - sourceCentroid;
        Vec3 y = vI->getTarget() - targetCentroid;
        E0 += 0.5 * vI->getWeight() * (dot(x, x) + dot(y, y));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R[i][j] += vI->getWeight() * y[i] * x[j];
            }
        }
    }

    // b) Form ~RR, determine eigenvalues mu(k) and the
    // mutually orthogonal eigenvectors a(k) and
    //    sort so that mu1 >= mu2 >= mu3.
    // Set
    //    a3 == a1 cross a2
    // to be sure to have a right handed system
    Vector_<std::complex<SimTK::Real>> mu0; // eigenvalues, complex, unsorted
    Matrix_<std::complex<SimTK::Real>> a0;  // eigenvectors, complex, unsorted
    Eigen eigen(Matrix(R.transpose() * R));
    eigen.getAllEigenValuesAndVectors(mu0, a0);

    // use only real component of results
    Vec3 a1[3];  // eigenvectors, real, unsorted
    Real mu1[3]; // eigenvalues, real, unsorted
    for (int i = 0; i < 3; ++i) {
        mu1[i] = mu0[i].real();
        for (int j = 0; j < 3; ++j) {
            // Note swapping of indices: It appears that the columns of a0 are eigenvectors
            a1[j][i] = a0[i][j].real();
        }
    }

    // sort indices of eigenvalues, from largest eigenvalue to smallest
    int indMax(0);
    int indMin(0);
    for (int i = 0; i < 3; ++i) {
        if (mu1[i] > mu1[indMax]) {
            indMax = i;
        }
        if (mu1[i] <= mu1[indMin]) {
            indMin = i;
        }
    }
    assert(indMin != indMax);
    int indMid = 3 - indMax - indMin; // too clever...
    assert(indMid >= 0);
    assert(indMid <= 2);
    assert(indMid != indMax);
    assert(indMid != indMin);
    int indSort[3] = {indMax, indMid, indMin};

    Real mu[3]; // eigenvalues, sorted
    Vec3 a[3];  // eigenvectors, sorted
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

    if (dot(b[2], R * a[2]) < 0) {
        sigma[2] = -1.0;
    }

    // d) Form U according to eq. 7:
    // u(ij) = SUM(k)[ b(ki)a(kj) ]
    // where b(k) = R*a(k)/(sigma(k)sqrt(mu(k)))
    Mat<3, 3> U(0.0); // initialize 2-d array
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                U[i][j] += b[k][i] * a[k][j];
            }
        }
    }

    // Compute residual error
    const auto residualError = E0 - (sigma[0] * std::sqrt(mu[0])) - (sigma[1] * std::sqrt(mu[1]))
                               - (sigma[2] * std::sqrt(std::abs(mu[2])));

    // E = 1/2 SUM(over n)[w(n)*(Ux(n) - y(n))^2]
    // variance would be sqrt( 1/SUM(w(n)) * SUM(w(n)*(Ux(n) - y(n))^2) )
    Real variance = 0.0;
    if (totalMass > 0) {
        variance = std::sqrt(2 * residualError / totalMass);
    }

    Rotation rotation(U);

    // No rotation for single point overlay
    if (vectors.size() < 2) {
        rotation = Rotation();
    }

    Transform transform1(-sourceCentroid);
    Transform transform2(rotation);
    Transform transform3(targetCentroid);

    return {transform3 * transform2 * transform1, variance};
}

auto Kabsch78::superpose_unweighted(const AtomSet& vectors) -> TransformAndResidual {
    Real E0 = 0.0;    // Initial residual, see Kabsch
    Mat<3, 3> R(0.0); // Form R matrix from Kabsch paper
    Centroids centroids;

    // a) Remove any translation between the two given
    // vector sets x(n) and y(n), and determine
    // E0 = 1/2 SUM[ w(n)*(x(n)^2 + y(n)^2) ]
    // and R, r(ij) = SUM(n)[ w(n)*y(ni)*x(nj) ]
    accumulate(vectors, E0, R, centroids);

    // b) Form ~RR, determine eigenvalues mu(k) and the
    // mutually orthogonal eigenvectors a(k) and
    //    sort so that mu1 >= mu2 >= mu3.
    // Set
    //    a3 == a1 cross a2
    // to be sure to have a right handed system
    std::array<Real, 3> mu;
    std::array<Vec3, 3> a;
    eigenSymm3(R.transpose() * R, mu, a); // sorted descending, real, no heap

    a[2] = cross(a[0], a[1]); // force right handed system

    // c) Determine Ra(k) (k = 1,2,3), normalize the first
    // two vectors to obtain b1, b2, and set b3 == b1 cross b2.
    // This will also take care of the case mu2 > mu3 = 0.
    std::array<Real, 3> sigma = {1.0, 1.0, 1.0};
    std::array<Vec3, 3> b;
    b[0] = Vec3(UnitVec3(R * a[0]));
    b[1] = Vec3(UnitVec3(R * a[1]));
    b[2] = cross(b[0], b[1]);

    if (dot(b[2], R * a[2]) < 0) {
        sigma[2] = -1.0;
    }

    // d) Form U according to eq. 7:
    // u(ij) = SUM(k)[ b(ki)a(kj) ]
    // where b(k) = R*a(k)/(sigma(k)sqrt(mu(k)))
    Mat<3, 3> U(0.0); // initialize 2-d array
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                U[i][j] += b[k][i] * a[k][j];
            }
        }
    }

    const auto residualError =
        E0 - std::sqrt(mu[0]) - std::sqrt(mu[1]) - (sigma[2] * std::sqrt(std::abs(mu[2])));
    const Real variance = (vectors.n > 0) ? std::sqrt(2.0 * residualError / vectors.n) : 0.0;

    Rotation rotation(U);

    // No rotation for single point overlay
    if (vectors.n < 2) {
        rotation = Rotation();
    }

    Transform transform1(-centroids.sourceCentroid);
    Transform transform2(rotation);
    Transform transform3(centroids.targetCentroid);

    return {transform3 * transform2 * transform1, variance};
}

void Kabsch78::accumulate(const AtomSet& atoms, Real& E0_out, Mat<3, 3>& R_out, Centroids& centroids) {
    const int n = atoms.n;

    // Raw restricted pointers - compiler can assume no aliasing
    const Real* __restrict__ srcX = atoms.sourceX.data();
    const Real* __restrict__ srcY = atoms.sourceY.data();
    const Real* __restrict__ srcZ = atoms.sourceZ.data();
    const Real* __restrict__ targetX = atoms.targetX.data();
    const Real* __restrict__ targetY = atoms.targetY.data();
    const Real* __restrict__ targetZ = atoms.targetZ.data();

    // All 18 accumulators as plain locals - compiler keeps these in registers.
    // Deliberately NOT in a struct/array: named scalars are easier for the
    // register allocator and prevent false dependency chains.
    Real sxs = 0;
    Real sys = 0;
    Real szs = 0;
    Real txs = 0;
    Real tys = 0;
    Real tzs = 0;
    Real srcSq = 0;
    Real tgtSq = 0;
    Real Rxx = 0;
    Real Rxy = 0;
    Real Rxz = 0;
    Real Ryx = 0;
    Real Ryy = 0;
    Real Ryz = 0;
    Real Rzx = 0;
    Real Rzy = 0;
    Real Rzz = 0;

    for (int i = 0; i < n; ++i) {
        const Real xi = srcX[i];
        const Real yi = srcY[i];
        const Real zi = srcZ[i];
        const Real pi = targetX[i];
        const Real qi = targetY[i];
        const Real ri = targetZ[i];

        sxs += xi;
        sys += yi;
        szs += zi;
        txs += pi;
        tys += qi;
        tzs += ri;

        srcSq += (xi * xi) + (yi * yi) + (zi * zi);
        tgtSq += (pi * pi) + (qi * qi) + (ri * ri);

        Rxx += pi * xi;
        Rxy += pi * yi;
        Rxz += pi * zi;
        Ryx += qi * xi;
        Ryy += qi * yi;
        Ryz += qi * zi;
        Rzx += ri * xi;
        Rzy += ri * yi;
        Rzz += ri * zi;
    }

    // --- Everything below is O(1), outside the hot loop ---
    const auto invN = 1.0 / static_cast<Real>(n);
    const Real x_mean = sxs * invN;
    const Real y_mean = sys * invN;
    const Real z_mean = szs * invN;
    const Real p_mean = txs * invN;
    const Real q_mean = tys * invN;
    const Real r_mean = tzs * invN;
    const auto nf = static_cast<Real>(n);

    // Parallel-axis correction
    R_out[0][0] = Rxx - (nf * p_mean * x_mean);
    R_out[0][1] = Rxy - (nf * p_mean * y_mean);
    R_out[0][2] = Rxz - (nf * p_mean * z_mean);
    R_out[1][0] = Ryx - (nf * q_mean * x_mean);
    R_out[1][1] = Ryy - (nf * q_mean * y_mean);
    R_out[1][2] = Ryz - (nf * q_mean * z_mean);
    R_out[2][0] = Rzx - (nf * r_mean * x_mean);
    R_out[2][1] = Rzy - (nf * r_mean * y_mean);
    R_out[2][2] = Rzz - (nf * r_mean * z_mean);

    E0_out = 0.5
             * (srcSq - nf * (x_mean * x_mean + y_mean * y_mean + z_mean * z_mean) + tgtSq
                - nf * (p_mean * p_mean + q_mean * q_mean + r_mean * r_mean));

    centroids.sourceCentroid = Vec3(x_mean, y_mean, z_mean);
    centroids.targetCentroid = Vec3(p_mean, q_mean, r_mean);
}

void Kabsch78::eigenSymm3(const Mat<3, 3>& RtR, std::array<Real, 3>& mu, std::array<Vec3, 3>& a) {
    const Real m00 = RtR[0][0];
    const Real m11 = RtR[1][1];
    const Real m22 = RtR[2][2];
    const Real m01 = RtR[0][1];
    const Real m02 = RtR[0][2];
    const Real m12 = RtR[1][2];

    // --- Eigenvalues: trigonometric method (Smith 1961) ---
    // Shift by mean eigenvalue to improve numerical conditioning
    const Real tr = m00 + m11 + m22;
    const Real q = tr * (1.0 / 3.0);
    const Real b00 = m00 - q;
    const Real b11 = m11 - q;
    const Real b22 = m22 - q;

    // ||M - qI||_F^2
    const Real p2 =
        (b00 * b00) + (b11 * b11) + (b22 * b22) + (2.0 * ((m01 * m01) + (m02 * m02) + (m12 * m12)));

    if (p2 < 1e-30) {
        // All eigenvalues equal - pathological case (perfect sphere)
        mu[0] = mu[1] = mu[2] = q;
        a[0] = Vec3(1, 0, 0);
        a[1] = Vec3(0, 1, 0);
        a[2] = Vec3(0, 0, 1);
        return;
    }

    const Real p = std::sqrt(p2 * (1.0 / 6.0));
    const Real invP3 = 1.0 / (p * p * p);

    // det((M - qI)/p) / 2  -- lies in [-1, 1] for real symmetric matrices
    const Real r = std::clamp(
        0.5 * invP3
            * (b00 * (b11 * b22 - m12 * m12) - m01 * (m01 * b22 - m12 * m02) + m02 * (m01 * m12 - b11 * m02)),
        -1.0,
        1.0);

    const Real phi = std::acos(r) * (1.0 / 3.0); // in [0, pi/3]

    // Eigenvalues sorted descending (guaranteed by phi in [0, pi/3])
    mu[0] = q + (2.0 * p * std::cos(phi));
    mu[2] = q + (2.0 * p * std::cos(phi + (2.0 * M_PI / 3.0)));
    mu[1] = tr - mu[0] - mu[2]; // exact: avoids third trig call

    // --- Eigenvectors: nullspace via row cross-products ---
    // For eigenvalue lam, rows of (M - lam*I) all lie perpendicular
    // to the eigenvector. Take all three pairwise cross products,
    // use the longest (best-conditioned).
    //
    // Rows: r0=(d0,m01,m02), r1=(m01,d1,m12), r2=(m02,m12,d2)
    for (int k = 0; k < 3; ++k) {
        const Real d0 = m00 - mu[k];
        const Real d1 = m11 - mu[k];
        const Real d2 = m22 - mu[k];

        const Vec3 v01((m01 * m12) - (m02 * d1), // r0 x r1
                       (m02 * m01) - (d0 * m12),
                       (d0 * d1) - (m01 * m01));

        const Vec3 v02((m01 * d2) - (m02 * m12), // r0 x r2
                       (m02 * m02) - (d0 * d2),
                       (d0 * m12) - (m01 * m02));

        const Vec3 v12((d1 * d2) - (m12 * m12), // r1 x r2
                       (m12 * m02) - (m01 * d2),
                       (m01 * m12) - (d1 * m02));

        const Real n01 = dot(v01, v01);
        const Real n02 = dot(v02, v02);
        const Real n12 = dot(v12, v12);

        const Vec3& best = (n01 >= n02 && n01 >= n12) ? v01 : (n02 >= n12) ? v02 : v12;
        a[k] = Vec3(UnitVec3(best));
    }
}
