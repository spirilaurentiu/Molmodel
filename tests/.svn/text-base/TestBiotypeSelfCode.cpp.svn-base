#include "SimTKmolmodel.h"
#include <iostream>
#include <fstream>

using namespace SimTK;
using namespace std;

int main() 
{
    Biotype argonBiotype = Biotype::Argon();

    DuMMForceFieldSubsystem dumm;
    ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    dumm.populateFromTinkerParameterFile(tinkerStream);

    // ofstream codeOut("biotype_test.cpp");

    Biotype::generateAllBiotypeCode(cout);
}
