#include "molmodel/internal/bondGeometry.h"

namespace SimTK {

// Re-using the fast horizontal sum for 3 elements to avoid memory spills
inline double fast_hsum3(__m256d v) {
    __m128d low = _mm256_castpd256_pd128(v);
    __m128d high = _mm256_extractf128_pd(v, 1);
    __m128d combined = _mm_add_pd(low, high);
    double res[2];
    _mm_storeu_pd(res, combined);
    return res[0] + res[1];
}

// Helper for 3D Cross Product using AVX
// Logic: (a1b2 - a2b1, a2b0 - a0b2, a0b1 - a1b0)
inline __m256d avx_cross(const __m256d a, const __m256d b) {
    __m256d a_yzx = _mm256_permute4x64_pd(a, _MM_SHUFFLE(3, 0, 2, 1));
    __m256d b_zxy = _mm256_permute4x64_pd(b, _MM_SHUFFLE(3, 1, 0, 2));
    __m256d a_zxy = _mm256_permute4x64_pd(a, _MM_SHUFFLE(3, 1, 0, 2));
    __m256d b_yzx = _mm256_permute4x64_pd(b, _MM_SHUFFLE(3, 0, 2, 1));
    return _mm256_sub_pd(_mm256_mul_pd(a_yzx, b_zxy), _mm256_mul_pd(a_zxy, b_yzx));
}

// Inline dot product for first 3 elements
inline double avx_dot3(const __m256d a, const __m256d b) {
    __m256d mul = _mm256_mul_pd(a, b);
    // Sum the first three doubles
    return ((double*)&mul)[0] + ((double*)&mul)[1] + ((double*)&mul)[2];
}


Angle calcAngle(const Vec3& p1, const Vec3& p2, const Vec3& p3) {
    // 1. Load points into registers
    // Note: If Vec3 is already 4-double aligned, use _mm256_load_pd
    __m256d v1 = _mm256_set_pd(0.0, p1[2], p1[1], p1[0]);
    __m256d v2 = _mm256_set_pd(0.0, p2[2], p2[1], p2[0]);
    __m256d v3 = _mm256_set_pd(0.0, p3[2], p3[1], p3[0]);

    // 2. Compute vectors relative to the vertex (p2)
    __m256d v21 = _mm256_sub_pd(v1, v2);
    __m256d v23 = _mm256_sub_pd(v3, v2);

    // 3. Compute dot product and squared magnitudes
    // Use the same multiplication result for dot(v21, v21) and dot(v23, v23)
    double d21_23 = fast_hsum3(_mm256_mul_pd(v21, v23));
    double d21_21 = fast_hsum3(_mm256_mul_pd(v21, v21));
    double d23_23 = fast_hsum3(_mm256_mul_pd(v23, v23));

    // 4. Calculate cosine using 1 sqrt and 1 div
    // cos(theta) = (v21 . v23) / sqrt(|v21|^2 * |v23|^2)
    double denomSq = d21_21 * d23_23;
    
    // Guard against division by zero for overlapping atoms
    if (denomSq < 1e-18) return 0.0; 

    double cosAngle = d21_23 / std::sqrt(denomSq);

    // 5. Clamp and result
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    return std::acos(cosAngle);
}

/*
    Dihedral angles are define in terms of four atomic positions

    1
    \
        \
        2------3
                \
                \
                4

    calcDihedralAngle returns a dihedral angle in radians, 
    in the range (-Pi, Pi]
    given three unit vectors
    pointing in the direction of the 1->2 axis, the 2->3 axis, and the 3->4 axis, respectively.
    The answer is unchanged if the order of atoms is reversed from (1,2,3,4) to (4,3,2,1)
    This method can be used as a helper method for dihedral angles express using either 4 atomic
    location, or two bond center orientations.
*/
Angle calcDihedralAngle(const UnitVec3& b12, const UnitVec3& b23, const UnitVec3& b34) {
    __m256d v12 = _mm256_set_pd(0.0, b12[2], b12[1], b12[0]);
    __m256d v23 = _mm256_set_pd(0.0, b23[2], b23[1], b23[0]);
    __m256d v34 = _mm256_set_pd(0.0, b34[2], b34[1], b34[0]);

    // Normal vectors: n1 = b12 x b23, n2 = b23 x b34
    __m256d vn1_raw = avx_cross(v12, v23);
    __m256d vn2_raw = avx_cross(v23, v34);

    // Normalize n1
    double dot_n1 = avx_dot3(vn1_raw, vn1_raw);
    __m256d vn1 = _mm256_div_pd(vn1_raw, _mm256_set1_pd(std::sqrt(dot_n1)));

    // Normalize n2
    double dot_n2 = avx_dot3(vn2_raw, vn2_raw);
    __m256d vn2 = _mm256_div_pd(vn2_raw, _mm256_set1_pd(std::sqrt(dot_n2)));

    double cosAngle = avx_dot3(vn1, vn2);

    // Clamp for acos
    if (cosAngle > 1.0) cosAngle = 1.0;
    else if (cosAngle < -1.0) cosAngle = -1.0;

    double angle = std::acos(cosAngle);

    // Sign bit logic: dot(n1, b34)
    if (avx_dot3(vn1, v34) < 0) angle = -angle;

    return angle;
}

Angle calcDihedralAngle(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4) {
    // 1. Load and Subtrack (Cheap)
    __m256d v1 = _mm256_set_pd(0.0, p1[2], p1[1], p1[0]);
    __m256d v2 = _mm256_set_pd(0.0, p2[2], p2[1], p2[0]);
    __m256d v3 = _mm256_set_pd(0.0, p3[2], p3[1], p3[0]);
    __m256d v4 = _mm256_set_pd(0.0, p4[2], p4[1], p4[0]);

    __m256d b12 = _mm256_sub_pd(v2, v1);
    __m256d b23 = _mm256_sub_pd(v3, v2);
    __m256d b34 = _mm256_sub_pd(v4, v3);

    // 2. Cross Products (Raw vectors, skip bond normalization!)
    __m256d vn1 = avx_cross(b12, b23);
    __m256d vn2 = avx_cross(b23, b34);

    // 3. One-shot normalization for the final cosine
    double dot_n1n1 = avx_dot3(vn1, vn1);
    double dot_n2n2 = avx_dot3(vn2, vn2);
    double dot_n1n2 = avx_dot3(vn1, vn2);

    // Combined: only 1 sqrt and 1 div total
    double cosAngle = dot_n1n2 / std::sqrt(dot_n1n1 * dot_n2n2);

    // 4. Cleanup and Sign
    if (cosAngle > 1.0) cosAngle = 1.0;
    else if (cosAngle < -1.0) cosAngle = -1.0;
    
    double angle = std::acos(cosAngle);
    if (avx_dot3(vn1, b34) < 0) angle = -angle;

    return angle;
}

} // namespace SimTK

