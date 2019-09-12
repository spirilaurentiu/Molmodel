#include "SimTKmolmodel.h"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;
using namespace SimTK;

void testCompletePdbLineShouldBeOK() 
{
    // Truncate at character 57; has coordinates but nothing after
    istringstream completeLine("ATOM    866  CA  LEU A 112      30.133  37.313  21.927  0.90 34.39           C  ");
    PdbStructure pdbStructure(completeLine);

    const PdbAtom& atom = pdbStructure.getAtom(" CA ", PdbResidueId(112, ' '), "A");
    SimTK_ASSERT_ALWAYS(atom.getOccupancy() == 0.90, "Atom occupancy incorrectly parsed");
    SimTK_ASSERT_ALWAYS(atom.getTemperatureFactor() == 34.39, "Atom temperature factor inorrectly parsed");
}

void testCoordinatesButNoOccupancyShouldBeOK() 
{
    // Truncate at character 57; has coordinates but nothing after
    istringstream trunc57("ATOM      1  N   ALA A   1     -52.630  -1.437  23.003");
    PdbStructure pdbStructure(trunc57);

    const PdbAtom& atom = pdbStructure.getAtom(" N  ", PdbResidueId(1, ' '), "A");
    SimTK_ASSERT_ALWAYS(atom.getOccupancy() == 1.00, "Atom occupancy incorrectly parsed");
    SimTK_ASSERT_ALWAYS(atom.getTemperatureFactor() == 0.00, "Atom temperature factor inorrectly parsed");
}

// Lack of coordinates should raise exception
void testMissingZCoordinateShouldRaiseException() 
{
    try {
        istringstream trunc("ATOM      1  N   ALA A   1     -52.630  -1.437");
        PdbStructure pdbStructure(trunc);
    }
    catch (...) {
        return;
    }

    throw std::logic_error("Missing z-coordinate failed to raise an exception");
}

int main() 
{
try {
    testCompletePdbLineShouldBeOK();
    testCoordinatesButNoOccupancyShouldBeOK();
    testMissingZCoordinateShouldRaiseException();

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

