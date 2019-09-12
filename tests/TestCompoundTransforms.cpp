#include "SimTKmolmodel.h"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;
using namespace SimTK;

void testBaseAtomTransform() 
{
    Compound compound;

    Vec3 baseTransform(17, 0, 0);
    compound.setBaseAtom( QuadrivalentAtom("A", Element::Carbon()), baseTransform );

    // test initial default location
    Vec3 atomLocation = compound.calcDefaultAtomLocationInGroundFrame("A");
    Vec3 errorVec( baseTransform - atomLocation );
    SimTK_ASSERT_ALWAYS(dot(errorVec, errorVec) < 0.1, "setBaseAtom() transform not applied correctly in default location");
    SimTK_ASSERT_ALWAYS(dot(atomLocation, atomLocation) > 0.1, "setBaseAtom() transform not applied correctly in default location");

    // test dynamic location
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    Amber99ForceSubsystem amber99(system);
    compound.setAtomBiotype("A", "Alanine", " CA ");

    system.adoptCompound(compound);
    system.modelCompounds();
    State state = system.realizeTopology();
    system.realize(state, Stage::Position);
    atomLocation = compound.calcAtomLocationInGroundFrame( state, Compound::AtomIndex(0) );
    errorVec = baseTransform - atomLocation;
    SimTK_ASSERT_ALWAYS(dot(errorVec, errorVec) < 0.1, "setBaseAtom() transform not applied correctly in dynamic location");
}

void testAdoptCompoundTransform()
{
    Compound compound;

    compound.setBaseAtom( QuadrivalentAtom("A", Element::Carbon()) );

    Vec3 atomLocation = compound.calcDefaultAtomLocationInGroundFrame("A");
    assert(dot(atomLocation, atomLocation) < 0.1);

    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    Amber99ForceSubsystem amber99(system);
    compound.setAtomBiotype("A", "Alanine", " CA ");

    Vec3 adoptTransform(0, -6.8, 0);
    assert(dot(adoptTransform, adoptTransform) > 0.1);

    system.adoptCompound(compound, adoptTransform);

    // adoptCompound(Transform) should affect default configuration too!!!
    atomLocation = compound.calcDefaultAtomLocationInGroundFrame("A");
    Vec3 errorVec = atomLocation - adoptTransform;
    SimTK_ASSERT_ALWAYS(dot(atomLocation, atomLocation) > 0.1, "adoptCompound() transform not applied correctly in default location");
    SimTK_ASSERT_ALWAYS(dot(errorVec, errorVec) < 0.1, "adoptCompound() transform not applied correctly in default location");

    // test dynamic location
    system.modelCompounds();
    State state = system.realizeTopology();
    system.realize(state, Stage::Position);
    atomLocation = compound.calcAtomLocationInGroundFrame( state, Compound::AtomIndex(0) );
    errorVec= adoptTransform - atomLocation;
    SimTK_ASSERT_ALWAYS(dot(errorVec, errorVec) < 0.1, "adoptCompound() transform not applied correctly in dynamic location");
}

void testMatchTopLevelTransform()
{
    Compound compound;

    compound.setBaseAtom( QuadrivalentAtom("A", Element::Carbon()) );

    Vec3 atomLocation = compound.calcDefaultAtomLocationInGroundFrame("A");
    assert(dot(atomLocation, atomLocation) < 0.1);

    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    Amber99ForceSubsystem amber99(system);
    compound.setAtomBiotype("A", "Alanine", " CA ");

    Vec3 adoptTransform(0, 0, 3.14159);
    assert(dot(adoptTransform, adoptTransform) > 0.1);

    system.adoptCompound(compound, adoptTransform);

    Vec3 atomTarget(1.7, -3.6, 0.002);
    Compound::AtomTargetLocations atomTargetLocations;
    atomTargetLocations[Compound::AtomIndex(0)] = atomTarget;

    // test getTransformAndResidual()
    Transform transform1 = compound.getTransformAndResidual(atomTargetLocations).transform;
    Vec3 errorVec = transform1 * adoptTransform - atomTarget;
    SimTK_ASSERT_ALWAYS( dot(errorVec, errorVec) < 0.1, "getTransformAndResidual() got wrong answer");

    // test matchDefaultTopLevelTransform()
    compound.matchDefaultTopLevelTransform(atomTargetLocations);
    errorVec = atomTarget - compound.calcDefaultAtomLocationInGroundFrame("A");
    SimTK_ASSERT_ALWAYS( dot(errorVec, errorVec) < 0.1, "matchDefaultTopLevelTransform got wrong answer" );
}


int main() 
{
try {
    testBaseAtomTransform();
    testAdoptCompoundTransform();
    testMatchTopLevelTransform();

    cout << "PASSED" << endl;
    return 0;
}

catch (const std::exception& e)
{
    cerr << "EXCEPTION THROWN: " << e.what() << endl;

    cerr << "FAILED" << endl;
    return 1;
}

catch (...)
{
    cerr << "UNKNOWN EXCEPTION THROWN" << endl;

    cerr << "FAILED" << endl;
    return 1;
}

}

