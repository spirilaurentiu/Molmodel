#include "molmodel/internal/bondGeometry.h"

namespace SimTK {

// Dihedral angles are define in terms of four atomic positions
// 
//   1
//    \
//     \
//      2------3
//              \
//               \
//                4
//
// calcDihedralAngle returns a dihedral angle in radians, 
// in the range (-Pi, Pi]
// given three unit vectors
// pointing in the direction of the 1->2 axis, the 2->3 axis, and the 3->4 axis, respectively.
// The answer is unchanged if the order of atoms is reversed from (1,2,3,4) to (4,3,2,1)
// This method can be used as a helper method for dihedral angles express using either 4 atomic
// location, or two bond center orientations.
Angle calcDihedralAngle(const UnitVec3& bond12, const UnitVec3& bond23, const UnitVec3& bond34)
{
    // 12 and 34 vectors must not be colinear with 23
    assert( dot (bond12, bond23) < 0.999 );
    assert( dot (bond12, bond23) > -0.999 );
    assert( dot (bond23, bond34) < 0.999 );
    assert( dot (bond23, bond34) > -0.999 );

    // Normal vectors to the junctions
    UnitVec3 n1(bond12 % bond23);
    UnitVec3 n2(bond23 % bond34);

    Real cosAngle = dot(n1, n2);

    // Sometimes roundoff gives an proposed cosine a smidgen over 1.0
    assert(cosAngle < 1.1);
    assert(cosAngle > -1.1);
    if (cosAngle > 1.0) cosAngle = 1.0;
    if (cosAngle < -1.0) cosAngle = -1.0;

    Angle nominalDihedralAngle = std::acos( cosAngle );

    // Determine sign
    if (dot(n1, bond34) < 0) nominalDihedralAngle *= -1;

    // avoid NaN-ish values
    assert( nominalDihedralAngle > -4);
    assert( nominalDihedralAngle < 4);
    assert( ! (nominalDihedralAngle > 4) );
    assert( ! (nominalDihedralAngle < -4) );

    return nominalDihedralAngle;
}

Angle calcDihedralAngle(const Vec3& atomPos1, const Vec3& atomPos2, const Vec3& atomPos3, const Vec3& atomPos4)
{
    UnitVec3 bond12(atomPos2 - atomPos1);
    UnitVec3 bond23(atomPos3 - atomPos2);
    UnitVec3 bond34(atomPos4 - atomPos3);
    
    return calcDihedralAngle(bond12, bond23, bond34);
}

} // namespace SimTK

