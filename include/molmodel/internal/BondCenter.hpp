#pragma once

#include "molmodel/internal/common.h"
#include "molmodel/internal/units.h"
#include "SimTKcommon/internal/UnitVec.h"

namespace SimTK {

#ifndef EXPLORER
    #define EXPLORER 0
#endif

//==============================================================================
//                           CLASS BondCenter
//==============================================================================
/** 
 * BondCenter is a class representing one half of a covalent
 * bond between two atoms.  BondCenter belongs to one atom and
 * represents one possible direction for a covalent bond.
**/

/*! <!-- BondCenter is a class representing one half of a covalent
 * bond between two atoms.  BondCenter belongs to one atom and
 * represents one possible direction for a covalent bond. -->
*/
class BondCenter {
public:
    enum Chirality {RightHanded, LeftHanded, Planar};

    // Use Paul's method
    // Given two other bond directions, and two bond angles,
    // determine direction of third bond
    static UnitVec3 getBondDirection(
        const UnitVec3& bondDir_1, Angle theta1,
        const UnitVec3& bondDir_2, Angle theta2,
        Chirality chirality
        ) 
    {
        if(EXPLORER){
            std::cout<<"BondCenter::getBondDirection() a1 "<<bondDir_1<<" a2 "<<bondDir_2<<std::endl;
            std::cout<<__FILE__<<":"<<__LINE__<<" a1 theta1 a2 theta2 chirality "<< bondDir_1<<", "<< theta1<<", "<<bondDir_2<<", "<<theta2<<", "<<chirality<<std::endl;
        }

        // Complete basis with third vector
        const UnitVec3 bondDir_12_normal(bondDir_1 % bondDir_2);

        if (chirality == Planar) {

            if(EXPLORER){
                std::cout<<"  BondCenter::getBondDirection() chirality == Planar"<<std::endl;
            }
            assert( theta1 >= 0 );
            // assert( theta1 >= -180*Deg2Rad );
            assert( theta1 <= 180*Deg2Rad );
            assert( theta2 >= 0 );
            assert( theta2 <= 180*Deg2Rad );

            // TODO - decide direction of theta1 based on estimate of theta2

            // estimate angle between first two bond centers
            Angle theta12 = dot(bondDir_1, bondDir_2);
            if (theta12 > 1.0) theta12 = 1.0;
            if (theta12 < -1.0) theta12 = -1.0;
            theta12 = std::acos(theta12);

            // is this bond on opposite side from bond2 with respect to bond1?
            // two cases, 1) 1-3 + 1-2 == 2-3
            Angle oppositeError1 = std::abs(theta2 - theta1 - theta12);
            // 2) 1-3 + 1-2 == 360 - 2-3
            Angle oppositeError2 = std::abs(360*Deg2Rad - theta12 - theta1 - theta2);

            // or is this bond on the same side as bond2? => 2-3 angle == difference between 1-3 and 1-2
            Angle adjacentError = std::abs(std::abs(theta12 - theta1) - theta2);

            // if this is on the opposite side from bond2, flip theta1
            if ( (adjacentError > oppositeError1) || (adjacentError > oppositeError2) )
                theta1 = -theta1; // move to opposite hemisphere from a2 w.r.t. a1

            UnitVec3 direction( Rotation( theta1, bondDir_12_normal ) * bondDir_1 );
            UnitVec3 sanityCheck( Rotation( theta12, bondDir_12_normal) * bondDir_1 );

            return direction;
        }

        else { // non-planar chirality

            // Bond angles are strictly positive
            assert(theta1 >= 0);
            assert(theta1 <= 180*Deg2Rad);
            assert(theta2 >= 0);
            assert(theta2 <= 180*Deg2Rad);


            // Compute coefficients of new bond direction in basis a1,a2,a3
            Real cosTheta12 = dot(bondDir_1, bondDir_2);

            assert(cosTheta12 != 0); // non colinear
            assert(cosTheta12 >= -1);
            assert(cosTheta12 <= 1);

            Real cosTheta1 = std::cos(theta1);
            Real cosTheta2 = std::cos(theta2);

            if(EXPLORER){
                std::cout<<"BondCenter::getBondDirection() cosTheta1 "
                    <<cosTheta1<<" cosTheta2 "<<cosTheta2<<" cosTheta "<<cosTheta12<<std::endl;
            }            
            Real sinSquaredTheta12 = 1.0 - (cosTheta12*cosTheta12);

            // a1 and a2 must not be parallel
            assert(sinSquaredTheta12 > 0);

            Real v1 = (cosTheta1 - cosTheta12*cosTheta2) / sinSquaredTheta12;
            Real v2 = (cosTheta2 - cosTheta12*cosTheta1) / sinSquaredTheta12;
            Real v3Squared = 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta12);

            if(EXPLORER){
                std::cout<<__FILE__<<":"<<__LINE__<<" v1 v2 = "<<v1<<", "<<v2<<std::endl;
                std::cout<<"  BondCenter::getBondDirection() v1 "<<v1<<" v2 "<<v2<<" v3Sq  "<<v3Squared<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" v1 v2 cosTheta = "<<v1<<", "<<v2<<", "<<cosTheta12<<std::endl; 
                std::cout<<__FILE__<<":"<<__LINE__<<" 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta) = "<<1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta12)<<std::endl;
            }

            // no solutions for certain sets of angles
            //assert(v3Squared >= 0);
            Real v3; // EU
            if (!(v3Squared >= 0)){
                //RESTORE std::cout<<__FILE__<<":"<<__LINE__<<" No solution found .. v1 v2 cosTheta = "<<v1<<", "<<v2<<", "<<cosTheta<<". Consider a larger planarity threshold."<<std::endl;
                //RESTORE std::cout<<__FILE__<<":"<<__LINE__<<" 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta) = "<<1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta)<<std::endl;
                //RESTORE RESTORE exit(1);
                v3 = -1*(std::sqrt(-1*v3Squared)); // EU
            }
            else{ v3 = std::sqrt(v3Squared); } // EU

            //Real v3 = std::sqrt(v3Squared); // ReSTORE

            if ( chirality == LeftHanded ) v3 = -v3;

            return UnitVec3(v1*bondDir_1 + v2*bondDir_2 + v3*bondDir_12_normal);
        }
 
    }

    // Default constructor to be used only for storage in std types
    BondCenter();

    BondCenter(Angle angle1, Angle angle2, int yCenter, Chirality c);

    BondCenter(Angle angle1, Angle angle2, int yCenter, Chirality c, UnitVec3 dir); // NEWMOB

    bool isBonded() const;

    BondCenter& setDefaultBondLength(mdunits::Length l) {
        assert (! isBonded() );

        defaultBondLength = l;

        return *this;
    }

    BondCenter& setDefaultDihedralAngle(Angle a) {
        assert (! isBonded() );

        defaultDihedralAngle = a;

        return *this;
    }

    mdunits::Length getDefaultBondLength() const {
        return defaultBondLength;
    }

    Angle getDefaultDihedralAngle() const {
        return defaultDihedralAngle;
    }

    int getDefaultDihedralReferenceCenter() const {
        return defaultDihedralReferenceCenter;
    }

    Angle getDefaultBond1Angle() const {return defaultBond1Angle;}

    BondCenter& setDefaultBond1Angle(Angle angle) 
    {
        assert(angle > 0);
        assert(angle <= SimTK::Pi);
        defaultBond1Angle = angle;
        return *this;
    }

    Angle getDefaultBond2Angle() const {return defaultBond2Angle;}
    
    BondCenter& setDefaultBond2Angle(Angle angle) 
    {
        assert(angle > 0);
        assert(angle <= SimTK::Pi);
        defaultBond2Angle = angle;
        return *this;
    }

    Chirality getChirality() const {return chirality;}

    BondCenter& setChirality(Chirality c) {
        chirality = c;
        return *this;
    }

    BondCenter& setInboard(bool b) {
        assert(!bonded);
        inboard = b;
        return *this;
    }

    bool isInboard() const {
        return inboard;
    }

    BondCenter& setBonded(bool b);

    // NEWMOB BEGIN
    const UnitVec3& getDirection(void) const{
        return direction;
    }

    UnitVec3& updDirection(void){
        return direction;
    }

    BondCenter& setDirection(UnitVec3 dir) {
        direction = dir;
        return *this;
    }
    //NEWMOB END

protected:

private:

    bool inboard;
    bool bonded;

    Angle defaultBond1Angle; // bond angle with first bond center (on x-axis)
    Angle defaultBond2Angle; // bond angle with second bond center (in +y half of XY plane)
    int defaultDihedralReferenceCenter; // index of other bond center to use as y-axis
    Chirality chirality;
    UnitVec3 direction; // NEWMOB

    // Compound::BondIndex bondIndex; // if bonded, the bond object containing a subcompound

    mdunits::Length defaultBondLength;
    Angle defaultDihedralAngle;
};

inline BondCenter::Chirality chiralityFromPlaneDeviation(Real deviation)
{
    return deviation < 0 ?
        BondCenter::LeftHanded :
        BondCenter::RightHanded;
}

inline BondCenter::Chirality flippedChirality(
    BondCenter::Chirality c)
{
    if (c == BondCenter::RightHanded)
        return BondCenter::LeftHanded;
    if (c == BondCenter::LeftHanded)
        return BondCenter::RightHanded;
    return c; // planar unchanged
}

} // namespace SimTK
