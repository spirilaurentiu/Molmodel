
#include "molmodel/internal/Compound.h"

#include "SimTKmolmodel.h"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using namespace SimTK;
using namespace std;

void testExactMatch(const std::string& pdbFileName) 
{
    time_t startTime = time(NULL);
    
    // In order for Biotypes to work correctly, 
    // amber99 parameters must be loaded first?
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem forces(system);
    forces.loadAmber99Parameters();

    cout << "Loading Structure from File..." << endl;
    std::ifstream pdbInputStream(pdbFileName.c_str());
    assert(pdbInputStream.is_open());
    
    // Check speed of various pdb file matching operations
    PdbStructure inputPdbStructure(pdbInputStream);
    pdbInputStream.close();
    
    time_t endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating internal coordinate model..." << endl;
    RNA rna(inputPdbStructure);
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Noting target atom locations..." << endl;
    Compound::AtomTargetLocations atomTargets = 
            rna.createAtomTargets(inputPdbStructure);
    cout << atomTargets.size() << " matches found" << endl;
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching chirality..." << endl;
    rna.matchDefaultAtomChirality(atomTargets, 0.01);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching bond angles..." << endl;
    rna.matchDefaultBondAngles(atomTargets);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching bond lengths..." << endl;
    rna.matchDefaultBondLengths(atomTargets);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching dihedral angles..." << endl;
    rna.matchDefaultDihedralAngles(atomTargets, Compound::DistortPlanarBonds);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching top level transform..." << endl;
    rna.matchDefaultTopLevelTransform(atomTargets);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating simbody system..." << endl;
    system.adoptCompound(rna);
    rna.assignBiotypes();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating multibody model..." << endl;
    system.modelCompounds();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Writing PDB file" << endl;
    PdbStructure outputPdbStructure(rna);
    std::ofstream pdbOutputStream("test1.pdb");
    outputPdbStructure.write(pdbOutputStream);
    pdbOutputStream.close();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Done." << endl;
}

void testIdealizedMatch(const std::string& pdbFileName) 
{
    time_t startTime = time(NULL);
    
    // In order for Biotypes to work correctly, 
    // amber99 parameters must be loaded first?
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem forces(system);
    forces.loadAmber99Parameters();

    cout << "Loading Structure from File..." << endl;
    std::ifstream pdbInputStream(pdbFileName.c_str());
    assert(pdbInputStream.is_open());
    
    // Check speed of various pdb file matching operations
    PdbStructure inputPdbStructure(pdbInputStream);
    pdbInputStream.close();
    
    time_t endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating internal coordinate model..." << endl;
    RNA rna(inputPdbStructure);
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Noting target atom locations..." << endl;
    Compound::AtomTargetLocations atomTargets = 
            rna.createAtomTargets(inputPdbStructure);
    cout << atomTargets.size() << " matches found" << endl;
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching chirality..." << endl;
    rna.matchDefaultAtomChirality(atomTargets, 90*Deg2Rad);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching dihedral angles..." << endl;
    rna.matchDefaultDihedralAngles(atomTargets, Compound::KeepPlanarBonds);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Matching top level transform..." << endl;
    rna.matchDefaultTopLevelTransform(atomTargets);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Final fitting..." << endl;
    rna.fitDefaultConfiguration(atomTargets, 0.005);

    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating simbody system..." << endl;
    system.adoptCompound(rna);
    rna.assignBiotypes();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Creating multibody model..." << endl;
    system.modelCompounds();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Writing PDB file" << endl;
    PdbStructure outputPdbStructure(rna);
    std::ofstream pdbOutputStream("test2.pdb");
    outputPdbStructure.write(pdbOutputStream);
    pdbOutputStream.close();
    
    endTime = time(NULL);
    cout << endTime - startTime << " s elapsed time" << endl;
    startTime = endTime;
    
    cout << "Done." << endl;
}

int main(int argc, char *argv[]) 
{
    // PDB file name is required as argument
    assert(argc > 1);
    
    const std::string pdbFileName(argv[1]);
    
    testExactMatch(pdbFileName);
    testIdealizedMatch(pdbFileName);
}

