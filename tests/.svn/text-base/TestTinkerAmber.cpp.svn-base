#include "SimTKmolmodel.h"
#include <iostream>
#include <fstream>
#include "TinkerAmber99.h"

using namespace SimTK;
using namespace std;

int main() {
    cout << "TESTING..." << endl;
    try {
        // For nightly builds, don't count on presence of external files
        // DuMMForceFieldSubsystem dumm("../../resources/tinker_amber99_clean.prm");
        DuMMForceFieldSubsystem dumm;
        dumm.loadAmber99Parameters();
        dumm.dumpCForceFieldParameters(cout, "populateAmber99Params");

        // std::ofstream testOutFile("C:/testAmberParams.cpp");
        // dumm.dumpCForceFieldParameters(testOutFile, "populateAmber99Params");
        // dumm.generateBiotypeChargedAtomTypeSelfCode(testOutFile);
        dumm.generateBiotypeChargedAtomTypeSelfCode(cout);

        // Create a second file, which should be identical
        DuMMForceFieldSubsystem dumm2;
        // populateAmber99Params(dumm2);
        dumm2.loadAmber99Parameters();
        // std::ofstream testOutFile2("C:/testAmberParams2.cpp");
        dumm2.dumpCForceFieldParameters(cout, "populateAmber99Params");

        cout << "PASSED" << endl;
        return 0;
    }
    catch (...) {
        cout << "FAILED" << endl;
        return 1;
    }
}

