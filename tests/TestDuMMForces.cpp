/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-9 Stanford University and the Authors.         *
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

#include "SimTKmolmodel.h"
#include "SimTKsimbody.h"

#include "SimTKcommon/Testing.h"

#include <iostream>
#include <vector>

using namespace SimTK;
using namespace std;

static const Real angstroms = 0.10;
static const Real kilocalories_per_mole = 4.184;
static const Real degrees = Pi / 180.0;

// Base class for molecule systems used for testing individual forces
class TestSystem {
public:

    TestSystem() : system(), matter(system), dumm(system)
    {
        //dumm.setUseOpenMMAcceleration(false);
    }

    Real calcEnergy() 
    {
	    system.realize(system.updDefaultState(), Stage::Dynamics);
        return system.calcPotentialEnergy(system.getDefaultState());
    }

    SpatialVec calcForce(MobilizedBodyIndex body) 
    {
	   system.realize(system.updDefaultState(), Stage::Dynamics);
	   return system.getRigidBodyForces(system.getDefaultState(), Stage::Dynamics)[body];
    }

    const CompoundSystem& getSystem() const {return system;}
    CompoundSystem& updSystem() {return system;}

    const DuMMForceFieldSubsystem& getDuMM() const {return dumm;}
    DuMMForceFieldSubsystem& updDuMM() {return dumm;}

    const SimbodyMatterSubsystem& getMatter() const {return matter;}
    SimbodyMatterSubsystem& updMatter() {return matter;}

protected:
    CompoundSystem system;
	SimbodyMatterSubsystem matter;
    DuMMForceFieldSubsystem dumm;
};

/////////////
// Coulomb //
/////////////

// Simple test system for Coulomb force
class SodiumChlorideSystem : public TestSystem {
public:
    SodiumChlorideSystem(Real distance) 
    {
        // Turn off all but Coulomb force
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setCoulombGlobalScaleFactor(1.0);

	    SodiumIon::setAmberLikeParameters(dumm);
	    ChlorideIon::setAmberLikeParameters(dumm);

	    SodiumIon na;
	    ChlorideIon cl;
	    system.adoptCompound(na, Vec3(distance, 0, 0));
	    system.adoptCompound(cl);
	    system.modelCompounds();

	    sodiumBody = na.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));

	    system.updDefaultState() = system.realizeTopology();
    }

    Real calcForce() {
        Vec3 forceOnSodium = TestSystem::calcForce(sodiumBody)[1];
        return forceOnSodium[0]; // X-component is signed force
    }

protected:
    MobilizedBodyIndex sodiumBody;
};


void testCoulombEnergy(Real distance) 
{
    Real observedEnergy = SodiumChlorideSystem(distance).calcEnergy();

    // Real expectedEnergy0 = -SimTK_COULOMB_CONSTANT_IN_MD / distance;
    Real q1 = -1.0; // charge on chloride
    Real q2 = 1.0; // charge on sodium
    Real expectedEnergy = 138.935456 * q1 * q2 / distance;

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}


// Finite difference to prove Force = -dEnergy/dLength
void testCoulombEnergyVsForce(Real distance) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;
    // this assert is really a debugging assert, not a unit test condition
    assert( (delta/distance) < 1e-2 );

    // Estimate force from -dEnergy/dDistance
    Real energy1 = SodiumChlorideSystem(distance).calcEnergy();
    Real energy2 = SodiumChlorideSystem(distance + delta).calcEnergy();

    Real expectedForce = - (energy2 - energy1)/delta;
    Real observedForce = SodiumChlorideSystem(distance + delta/2.0).calcForce();

    SimTK_TEST(observedForce < 0.0);

    SimTK_TEST_EQ_TOL(expectedForce, observedForce, 0.1);
}

void testCoulombForce() 
{
    testCoulombEnergy(1.0);
    testCoulombEnergy(0.12345);
    testCoulombEnergy(17.623);

    testCoulombEnergyVsForce(1.0);
    testCoulombEnergyVsForce(0.12345);
    testCoulombEnergyVsForce(17.623);
}


//////////////////////////
// van der Waals forces //
//////////////////////////

class ArgonSystem : public TestSystem
{
public:
    ArgonSystem(Real distance) {
        // Turn off all but van der Waals force
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setVdwGlobalScaleFactor(1.0);

        if (! dumm.hasAtomClass("Argon") ) {
			dumm.defineAtomClass(
				dumm.getNextUnusedAtomClassIndex(),
				"Argon",
				Element::Oxygen().getAtomicNumber(),
				0, // no bonds
				1.88 * angstroms, // radius
				0.0003832 * kilocalories_per_mole // well depth
				);
		}

		if (! dumm.hasChargedAtomType("Argon") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"Argon", 
				dumm.getAtomClassIndex("Argon"),
				0.00 // no charge on symmetric molecule
				);
		}

        if (! Biotype::exists("Argon", "Ar") )
            Biotype::defineBiotype(Element::Argon(), 0, "Argon", "Ar");

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("Argon"), Biotype::get("Argon", "Ar").getIndex() );

        Argon argon1, argon2;
        argon1.setAtomBiotype("Ar", "Argon", "Ar");
        argon2.setAtomBiotype("Ar", "Argon", "Ar");

        system.adoptCompound(argon1);
        system.adoptCompound(argon2, Vec3(distance, 0, 0));
	    system.modelCompounds();

	    argon2Body = argon2.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));

	    system.updDefaultState() = system.realizeTopology();
    }

    Real calcForce() {
        Vec3 forceOnArgon2 = TestSystem::calcForce(argon2Body)[1];
        return forceOnArgon2[0]; // X-component is signed force
    }

protected:
    MobilizedBodyIndex argon2Body;
};

void testVanDerWaalsEnergy(Real distance) 
{
    Real observedEnergy = ArgonSystem(distance).calcEnergy();

    Real r = 0.5 * distance;
    Real rMin = 1.88 * angstroms;
    Real epsilon = 0.0003832 * kilocalories_per_mole;

    Real expectedEnergy = epsilon * ( std::pow((rMin/r), 12.0) - 2.0 * std::pow((rMin/r), 6.0) );

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

// Finite difference to prove Force = -dEnergy/dLength
void testVanDerWaalsEnergyVsForce(Real distance) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;
    // this assert is really a debugging assert, not a unit test condition
    assert( (delta/distance) < 1e-2 );

    // Estimate force from -dEnergy/dDistance
    Real energy1 = ArgonSystem(distance).calcEnergy();
    Real energy2 = ArgonSystem(distance + delta).calcEnergy();

    Real expectedForce = - (energy2 - energy1)/delta;
    Real observedForce = ArgonSystem(distance + delta/2.0).calcForce();

    SimTK_TEST_EQ_TOL(expectedForce, observedForce, 0.1);
}


void testVanDerWaalsForce() {
    testVanDerWaalsEnergy(5.0 * angstroms);
    testVanDerWaalsEnergy(3.7 * angstroms);
    testVanDerWaalsEnergy(3.0 * angstroms);

    testVanDerWaalsEnergyVsForce(5.0 * angstroms);
    testVanDerWaalsEnergyVsForce(3.7 * angstroms);
    testVanDerWaalsEnergyVsForce(3.0 * angstroms);
}


//////////////////
// Bond Stretch //
//////////////////

class OxygenSystem : public TestSystem {
public:
    OxygenSystem(Real bondLength) 
    {
        // Turn off all but bond stretch force
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setBondStretchGlobalScaleFactor(1.0);

        if (! dumm.hasAtomClass("O2Mol") ) {
			dumm.defineAtomClass(
				dumm.getNextUnusedAtomClassIndex(),
				"O2Mol",
				Element::Oxygen().getAtomicNumber(),
				1, // one bond
				1.70 * angstroms, // radius
				0.20 * kilocalories_per_mole // well depth
				);
		}

		if (! dumm.hasChargedAtomType("O2Mol") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"O2Mol", 
				dumm.getAtomClassIndex("O2Mol"),
				0.00 // no charge on symmetric molecule
				);
		}

        dumm.defineBondStretch(
            dumm.getAtomClassIndex("O2Mol"),
            dumm.getAtomClassIndex("O2Mol"),
            500.0 * kilocalories_per_mole / (angstroms * angstroms),
            1.21 * angstroms
        );

        if (! Biotype::exists("Oxygen Molecule", "O") )
            Biotype::defineBiotype(Element::Oxygen(), 1, "Oxygen Molecule", "O");

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("O2Mol"), Biotype::get("Oxygen Molecule", "O").getIndex() );

        Compound o2Molecule;
        o2Molecule.setBaseAtom( UnivalentAtom("O1", Element::Oxygen()) );
        o2Molecule.bondAtom( UnivalentAtom("O2", Element::Oxygen()), "O1/bond", bondLength );
        o2Molecule.setBondMobility(BondMobility::Free, "O1", "O2");
        o2Molecule.setAtomBiotype("O1", "Oxygen Molecule", "O");
        o2Molecule.setAtomBiotype("O2", "Oxygen Molecule", "O");

	    system.adoptCompound(o2Molecule);
	    system.modelCompounds();

	    o2Body = o2Molecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(1));

	    system.updDefaultState() = system.realizeTopology();
    }

    Real calcForce() {
        Vec3 forceOnO2 = TestSystem::calcForce(o2Body)[1];
        return forceOnO2[0]; // X-component is signed force
    }

protected:
    MobilizedBodyIndex o2Body;
};

void testBondStretchEnergy(Real bondLength) 
{
    Real observedEnergy = OxygenSystem(bondLength).calcEnergy();

    Real r0 = 1.21 * angstroms;
    Real kR = 500.0 * kilocalories_per_mole / (angstroms * angstroms);

    // Stiffness constant kR in MD codes like molmodel and tinker is 2.0 times what engineers would expect
    // Real expectedEnergy = 0.5 * kR * (bondLength - r0) * (bondLength - r0);
    Real expectedEnergy = kR * (bondLength - r0) * (bondLength - r0);

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

// Finite difference to prove Force = -dEnergy/dLength
void testBondStretchEnergyVsForce(Real bondLength) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;
    // this assert is really a debugging assert, not a unit test condition
    assert( (delta/bondLength) < 1e-2 );

    // Estimate force from -dEnergy/dDistance
    Real energy1 = OxygenSystem(bondLength).calcEnergy();
    Real energy2 = OxygenSystem(bondLength + delta).calcEnergy();

    Real expectedForce = - (energy2 - energy1)/delta;
    Real observedForce = OxygenSystem(bondLength + delta/2.0).calcForce();

    SimTK_TEST_EQ_TOL(expectedForce, observedForce, 0.01);
}

void testBondStretchForce()
{
    testBondStretchEnergy(1.10*angstroms);
    testBondStretchEnergy(1.35*angstroms);
    testBondStretchEnergy(1.71*angstroms);

    testBondStretchEnergyVsForce(1.10*angstroms);
    testBondStretchEnergyVsForce(1.35*angstroms);
    testBondStretchEnergyVsForce(1.71*angstroms);
}


/////////////////////////
// Custom Bond Stretch //
/////////////////////////

class HarmonicBondStretch
     : public DuMM::CustomBondStretch
{
public:
    HarmonicBondStretch(Real stiffness, Real idealLength)
        : stiffness(stiffness), idealLength(idealLength)
    {}

    Real calcEnergy(Real length) const {
        Real dR = length - idealLength;
        return stiffness * dR * dR;
    }

    Real calcForce(Real length) const {
        Real dR = length - idealLength;
        return -2.0 * stiffness * dR;
    }

private:
    Real idealLength;
    Real stiffness;
};

// Oxygen modified to use custom bond stretch force
class CustomOxygenSystem : public OxygenSystem {
public:
    CustomOxygenSystem(Real bondLength) : OxygenSystem(bondLength) 
    {
        dumm.setAllGlobalScaleFactors(0.0);
        dumm.setCustomBondStretchGlobalScaleFactor(1.0);

        dumm.defineCustomBondStretch(
            dumm.getAtomClassIndex("O2Mol"),
            dumm.getAtomClassIndex("O2Mol"),
            new HarmonicBondStretch(
                600.0 * kilocalories_per_mole / (angstroms * angstroms),
                1.21 * angstroms
            )
        );

        system.realizeTopology();
    }
};

void testCustomBondStretchEnergy(Real bondLength) 
{
    Real observedEnergy = CustomOxygenSystem(bondLength).calcEnergy();

    Real r0 = 1.21 * angstroms;
    Real kR = 600.0 * kilocalories_per_mole / (angstroms * angstroms);

    // Stiffness constant kR in MD codes like molmodel and tinker is 2.0 times what engineers would expect
    // Real expectedEnergy = 0.5 * kR * (bondLength - r0) * (bondLength - r0);
    Real expectedEnergy = kR * (bondLength - r0) * (bondLength - r0);

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

// Finite difference to prove Force = -dEnergy/dLength
void testCustomBondStretchEnergyVsForce(Real bondLength) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;
    // this assert is really a debugging assert, not a unit test condition
    assert( (delta/bondLength) < 1e-2 );

    // Estimate force from -dEnergy/dDistance
    Real energy1 = CustomOxygenSystem(bondLength).calcEnergy();
    Real energy2 = CustomOxygenSystem(bondLength + delta).calcEnergy();

    Real expectedForce = - (energy2 - energy1)/delta;
    Real observedForce = CustomOxygenSystem(bondLength + delta/2.0).calcForce();

    SimTK_TEST_EQ_TOL(expectedForce, observedForce, 0.01);
}

void testCustomBondStretchForce()
{
    testCustomBondStretchEnergy(1.10 * angstroms);
    testCustomBondStretchEnergy(1.35 * angstroms);
    testCustomBondStretchEnergy(1.70 * angstroms);

    testCustomBondStretchEnergyVsForce(1.10*angstroms);
    testCustomBondStretchEnergyVsForce(1.35*angstroms);
    testCustomBondStretchEnergyVsForce(1.71*angstroms);
}


///////////////
// Bond Bend //
///////////////

class WaterSystem : public TestSystem {
public:
    WaterSystem(Real bondAngle) 
    {
        // Turn off all but bond bend force
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setBondBendGlobalScaleFactor(1.0);

        if (! dumm.hasAtomClass("TIP3P Water O") ) {
			dumm.defineAtomClass(
				dumm.getNextUnusedAtomClassIndex(),
				"TIP3P Water O",
				Element::Oxygen().getAtomicNumber(),
				2, // one bond
				1.7683 * angstroms, // radius
				0.1520 * kilocalories_per_mole // well depth
				);
		}

        if (! dumm.hasAtomClass("TIP3P Water H") ) {
			dumm.defineAtomClass(
				dumm.getNextUnusedAtomClassIndex(),
				"TIP3P Water H",
				Element::Hydrogen().getAtomicNumber(),
				1, // one bond
				0.0001 * angstroms, // radius
				0.0000 * kilocalories_per_mole // well depth
				);
		}

		if (! dumm.hasChargedAtomType("TIP3P Water O") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"TIP3P Water O", 
				dumm.getAtomClassIndex("TIP3P Water O"),
				-0.834
				);
		}

		if (! dumm.hasChargedAtomType("TIP3P Water H") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"TIP3P Water H", 
				dumm.getAtomClassIndex("TIP3P Water H"),
				0.417
				);
		}

        dumm.defineBondStretch(
            dumm.getAtomClassIndex("TIP3P Water O"),
            dumm.getAtomClassIndex("TIP3P Water H"),
            553.0 * kilocalories_per_mole / (angstroms * angstroms),
            0.9572 * angstroms
        );

        dumm.defineBondBend(
            dumm.getAtomClassIndex("TIP3P Water H"),
            dumm.getAtomClassIndex("TIP3P Water O"),
            dumm.getAtomClassIndex("TIP3P Water H"),
            100.00 * kilocalories_per_mole / (degrees * degrees),
            104.52); // define bond bend takes angle in degrees!?!

        if (! Biotype::exists("TIP3P Water", "O") )
            Biotype::defineBiotype(Element::Oxygen(), 2, "TIP3P Water", "O");

        if (! Biotype::exists("TIP3P Water", "H") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TIP3P Water", "H");

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("TIP3P Water O"), Biotype::get("TIP3P Water", "O").getIndex() );
        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("TIP3P Water H"), Biotype::get("TIP3P Water", "H").getIndex() );

        Compound water;
        water.setBaseAtom( BivalentAtom("O", Element::Oxygen()) );
        water.bondAtom( UnivalentAtom("H1", Element::Hydrogen()), "O/bond1", 0.9572 * angstroms );
        water.bondAtom( UnivalentAtom("H2", Element::Hydrogen()), "O/bond2", 0.9572 * angstroms );
        water.setDefaultBondAngle(bondAngle, "H1", "O", "H2");

        // Unfortunately, there is not yet a BondMobility::Angle concept, so pure Pin is not available here
        water.setBondMobility(BondMobility::Rigid, "O", "H1");
        water.setBondMobility(BondMobility::Free, "O", "H2");

        water.setAtomBiotype("O", "TIP3P Water", "O");
        water.setAtomBiotype("H1", "TIP3P Water", "H");
        water.setAtomBiotype("H2", "TIP3P Water", "H");

	    system.adoptCompound(water);
	    system.modelCompounds();

	    system.updDefaultState() = system.realizeTopology();

        // debugging...
        // system.realize(system.getDefaultState(), Stage::Position);
        // cout << endl << endl;
        // water.writePdb(system.getDefaultState(), cout);
    }
};

void testBondBendEnergy(Real bondAngle) 
{
    Real observedEnergy = WaterSystem(bondAngle).calcEnergy();

    Real r0 = 104.52 * degrees;
    Real kR = 100.0 * kilocalories_per_mole / (degrees * degrees);

    // Stiffness constant kR in MD codes like molmodel and tinker is 2.0 times what engineers would expect
    // Real expectedEnergy = 0.5 * kR * (bondLength - r0) * (bondLength - r0);
    Real expectedEnergy = kR * (bondAngle - r0) * (bondAngle - r0);

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

void testBondBendForce() {
    testBondBendEnergy(104.52*degrees);
    testBondBendEnergy(90.0*degrees);
    testBondBendEnergy(85.3*degrees);
}


//////////////////////
// Custom Bond Bend //
//////////////////////

class HarmonicBondBend
     : public DuMM::CustomBondBend
{
public:
    HarmonicBondBend(Real stiffnessInKJPerNmSquared, Real idealAngleInRadians)
        : stiffness(stiffnessInKJPerNmSquared), idealAngle(idealAngleInRadians)
    {}

    Real calcEnergy(Real angle) const {
        Real dR = angle - idealAngle;
        return stiffness * dR * dR;
    }

    Real calcTorque(Real angle) const {
        Real dR = angle - idealAngle;
        return -2.0 * stiffness * dR;
    }

private:
    Real idealAngle;
    Real stiffness;
};

class CustomWaterSystem : public WaterSystem
{
public:
    CustomWaterSystem(Real bondAngle) 
        : WaterSystem(bondAngle)
    {
        dumm.defineCustomBondBend(
            dumm.getAtomClassIndex("TIP3P Water H"),
            dumm.getAtomClassIndex("TIP3P Water O"),
            dumm.getAtomClassIndex("TIP3P Water H"),
            new HarmonicBondBend(
                100.00 * kilocalories_per_mole / (degrees * degrees),
                104.52 * degrees)
        );

	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setCustomBondBendGlobalScaleFactor(1.0);
	    system.updDefaultState() = system.realizeTopology();
    }
};


// Show that HarmonicCustomBondBend forces are the same as those of internal DuMM bond bends
// It doesn't prove the forces are correct, but at least they agree.
void testCustomBondBendToBondBend(Real bondAngle)
{
    WaterSystem internalSystem(bondAngle);
    CustomWaterSystem customSystem(bondAngle);
    for (MobilizedBodyIndex b(0); b < internalSystem.getMatter().getNumBodies(); ++b) {
        SimTK_TEST_EQ(internalSystem.calcForce(b), customSystem.calcForce(b));
    }
}

void testCustomBondBendEnergy(Real bondAngle) 
{
    Real observedEnergy = CustomWaterSystem(bondAngle).calcEnergy();

    Real r0 = 104.52 * degrees;
    Real kR = 100.0 * kilocalories_per_mole / (degrees * degrees);

    // Stiffness constant kR in MD codes like molmodel and tinker is 2.0 times what engineers would expect
    // Real expectedEnergy = 0.5 * kR * (bondLength - r0) * (bondLength - r0);
    Real expectedEnergy = kR * (bondAngle - r0) * (bondAngle - r0);

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

void testCustomBondBendForce() 
{
    testCustomBondBendEnergy(90.0 * degrees);
    testCustomBondBendEnergy(104.5 * degrees);
    testCustomBondBendEnergy(110.0 * degrees);

    testCustomBondBendToBondBend(90.0 * degrees);
    testCustomBondBendToBondBend(104.5 * degrees);
    testCustomBondBendToBondBend(110.0 * degrees);
}


////////////////////////////
// Dihedral Torsion Angle //
////////////////////////////

class HydrogenPeroxideSystem : public TestSystem
{
public:
    HydrogenPeroxideSystem(Real torsion) {
        // Turn off all but bond bend force
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setBondTorsionGlobalScaleFactor(1.0);

        dumm.loadAmber99Parameters();

		if (! dumm.hasChargedAtomType("H2O2 O") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"H2O2 O", 
				dumm.getAtomClassIndex("OH"),
				-0.834
				);
		}

		if (! dumm.hasChargedAtomType("H2O2 H") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"H2O2 H", 
				dumm.getAtomClassIndex("HO"),
				0.417
				);
		}

        if (! Biotype::exists("H2O2", "O") )
            Biotype::defineBiotype(Element::Oxygen(), 2, "H2O2", "O");

        if (! Biotype::exists("H2O2", "H") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "H2O2", "H");

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("H2O2 O"), Biotype::get("H2O2", "O").getIndex() );
        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("H2O2 H"), Biotype::get("H2O2", "H").getIndex() );

        dumm.defineBondStretch(
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("OH"),
            500.0 * kilocalories_per_mole / (angstroms * angstroms),
            1.4464 * angstroms
        );

        dumm.defineBondBend(
            dumm.getAtomClassIndex("HO"),
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("OH"),
            120.0 * kilocalories_per_mole / (degrees * degrees),
            100.8982 // * degrees // defineBondBend() takes actual degrees
        );

        dumm.defineBondTorsion(
            dumm.getAtomClassIndex("HO"),
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("HO"),
            2, // periodicity
            0.53 * kilocalories_per_mole,
            111.9324 // * degrees // defineTorsion() takes actual degrees
            // 0.0 // * degrees // defineTorsion() takes actual degrees
        );

        Compound h2O2;
        h2O2.setBaseAtom( BivalentAtom("O1", Element::Oxygen()) );
        h2O2.bondAtom( BivalentAtom("O2", Element::Oxygen()), "O1/bond1", 1.4464 * angstroms );
        h2O2.bondAtom( UnivalentAtom("H1", Element::Hydrogen()), "O1/bond2", 0.9659 * angstroms );
        h2O2.bondAtom( UnivalentAtom("H2", Element::Hydrogen()), "O2/bond2", 0.9659 * angstroms );
        h2O2.setDefaultBondAngle(100.8982 * degrees, "H1", "O1", "O2");
        h2O2.setDefaultBondAngle(100.8982 * degrees, "O1", "O2", "H2");
        h2O2.setDefaultDihedralAngle(torsion, "H1", "O1", "O2", "H2");

        // Unfortunately, there is not yet a BondMobility::Angle concept, so pure Pin is not available here
        h2O2.setBondMobility(BondMobility::Rigid, "O1", "H1");
        h2O2.setBondMobility(BondMobility::Rigid, "O2", "H2");
        h2O2.setBondMobility(BondMobility::Torsion, "O1", "O2");

        h2O2.setAtomBiotype("O1", "H2O2", "O");
        h2O2.setAtomBiotype("O2", "H2O2", "O");
        h2O2.setAtomBiotype("H1", "H2O2", "H");
        h2O2.setAtomBiotype("H2", "H2O2", "H");

	    system.adoptCompound(h2O2);
	    system.modelCompounds();

        o2BodyIndex = h2O2.getAtomMobilizedBodyIndex(h2O2.getAtomIndex("O2"));

	    system.updDefaultState() = system.realizeTopology();
    }

    Real calcTorque() const 
    {
        system.realize(system.getDefaultState(), Stage::Dynamics);

        // Express rigid body forces in terms of mobilities
        Vector mobilityForces;

        matter.calcTreeEquivalentMobilityForces(
            system.getDefaultState(), 
            system.getRigidBodyForces(system.getDefaultState(), Stage::Dynamics),
            mobilityForces);

        Real torque = matter.getMobilizedBody(o2BodyIndex).getOneFromUPartition(
            system.getDefaultState(),
            0,
            mobilityForces);

        return torque;
    }

protected:
    MobilizedBodyIndex o2BodyIndex;
};

void testTorsionEnergy(Real torsion) 
{
    Real observedEnergy = HydrogenPeroxideSystem(torsion).calcEnergy();

    int periodicity = 2;
    Real halfAmplitude = 0.53 * kilocalories_per_mole;
    Real phase = 111.9324 * degrees;
    // Real phase = 0.0 * degrees;

    // Torsion is not doubled
    Real expectedEnergy = halfAmplitude * ( 1.0 + std::cos(periodicity * torsion - phase) );

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

// Finite difference to prove Force = -dEnergy/dLength
void testTorsionEnergyVsTorque(Real torsion) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;

    // Estimate force from -dEnergy/dDistance
    Real energy1 = HydrogenPeroxideSystem(torsion).calcEnergy();
    Real energy2 = HydrogenPeroxideSystem(torsion + delta).calcEnergy();

    Real expectedTorque = - (energy2 - energy1)/delta;
    Real observedTorque = HydrogenPeroxideSystem(torsion + delta/2.0).calcTorque();

    SimTK_TEST_EQ_TOL(expectedTorque, observedTorque, 0.01);
}

void testTorsionForce() {
    testTorsionEnergy(90.0 * degrees);
    testTorsionEnergy(111.932 * degrees);
    testTorsionEnergy(120.0 * degrees);
    testTorsionEnergy(191.932 * degrees);

    testTorsionEnergyVsTorque(90.0 * degrees);
    testTorsionEnergyVsTorque(111.932 * degrees);
    testTorsionEnergyVsTorque(120.0 * degrees);
    testTorsionEnergyVsTorque(191.932 * degrees);
}



///////////////////////////
// Custom Torsion forces //
///////////////////////////

class SinusoidalTorsion : public DuMM::CustomBondTorsion {
public:
    SinusoidalTorsion(int periodicity, Real halfAmplitudeInKJPerMole, Real phaseInRadians) 
        : periodicity(periodicity), halfAmplitude(halfAmplitudeInKJPerMole), phase(phaseInRadians)
    {}

    Real calcEnergy(Real torsionInRadians) const {
        return halfAmplitude * ( 1.0 + std::cos(periodicity * torsionInRadians - phase) );
    }

    Real calcTorque(Real torsionInRadians) const {
        return periodicity * halfAmplitude * ( std::sin(periodicity * torsionInRadians - phase) );
    }

private:
    int periodicity;
    Real halfAmplitude;
    Real phase;
};

class CustomHydrogenPeroxideSystem : public HydrogenPeroxideSystem
{
public:
    CustomHydrogenPeroxideSystem(Real torsion) : HydrogenPeroxideSystem(torsion) {
	    dumm.setAllGlobalScaleFactors(0.0);
        dumm.setCustomBondTorsionGlobalScaleFactor(1.0);

        dumm.defineCustomBondTorsion(
            dumm.getAtomClassIndex("HO"),
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("OH"),
            dumm.getAtomClassIndex("HO"),
            new SinusoidalTorsion(
                2, // periodicity
                1.42 * kilocalories_per_mole, // half amplitude
                73.4 * degrees) // phase
        );

	    system.updDefaultState() = system.realizeTopology();
    }
};

void testCustomTorsionEnergy(Real torsion) 
{
    Real observedEnergy = CustomHydrogenPeroxideSystem(torsion).calcEnergy();

    int periodicity = 2;
    Real halfAmplitude = 1.42 * kilocalories_per_mole;
    Real phase = 73.4 * degrees;
    // Real phase = 0.0 * degrees;

    // Torsion is not doubled
    Real expectedEnergy = halfAmplitude * ( 1.0 + std::cos(periodicity * torsion - phase) );

    SimTK_TEST_EQ(expectedEnergy, observedEnergy);
}

// Finite difference to prove Force = -dEnergy/dLength
void testCustomTorsionEnergyVsTorque(Real torsion) 
{
    // Finite difference must be small relative to total distance
    Real delta = 1e-4;

    // Estimate force from -dEnergy/dDistance
    Real energy1 = CustomHydrogenPeroxideSystem(torsion).calcEnergy();
    Real energy2 = CustomHydrogenPeroxideSystem(torsion + delta).calcEnergy();

    Real expectedTorque = - (energy2 - energy1)/delta;
    Real observedTorque = CustomHydrogenPeroxideSystem(torsion + delta/2.0).calcTorque();

    SimTK_TEST_EQ_TOL(expectedTorque, observedTorque, 0.01);
}

void testCustomTorsionForce() 
{
    testCustomTorsionEnergy(90.0 * degrees);
    testCustomTorsionEnergy(111.932 * degrees);
    testCustomTorsionEnergy(120.0 * degrees);
    testCustomTorsionEnergy(191.932 * degrees);

    testCustomTorsionEnergyVsTorque(90.0 * degrees);
    testCustomTorsionEnergyVsTorque(111.932 * degrees);
    testCustomTorsionEnergyVsTorque(120.0 * degrees);
    testCustomTorsionEnergyVsTorque(191.932 * degrees);
}




int main() 
{
    SimTK_START_TEST("TestDuMMForces");

    SimTK_SUBTEST(testCoulombForce);
    SimTK_SUBTEST(testVanDerWaalsForce);
    SimTK_SUBTEST(testBondStretchForce);
    SimTK_SUBTEST(testCustomBondStretchForce);
    SimTK_SUBTEST(testBondBendForce);
    SimTK_SUBTEST(testCustomBondBendForce);
    SimTK_SUBTEST(testTorsionForce);
    SimTK_SUBTEST(testCustomTorsionForce);

    SimTK_END_TEST();
}
