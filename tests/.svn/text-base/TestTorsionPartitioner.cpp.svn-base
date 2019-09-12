/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

#include <iostream>
#include <vector>

using std::cout;
using std::endl;

using namespace SimTK;
using namespace std;

// define CREATE_VIZ_WINDOW to see animated window of simulation
// undefine for automated nightly builds
// #define CREATE_VIZ_WINDOW

int main() {
//try
//  {
    Element e;
    e = Element::Argon();
    // cout << "Element=" << e << endl;

    AliphaticCarbon ac;
    // cout << "AliphaticCarbon " << ac << endl;

    MethylGroup methyl;
    // cout << "Methyl " << methyl << endl;

    Methane methane1;
    // cout << "Methane=" << methane1 << endl;
    Methane methane2;

    Ethane ethane1;
    Ethane ethane2;
    Ethane ethane3;
    ethane1.setDefaultTorsionAngle(55*Deg2Rad);
    // cout << "Ethane=" << ethane << endl;
    // cout << "  default torsion angle=" << ethane.getDefaultTorsionAngle() << endl;

    AminoAcidResidue::Serine serine;


    AminoAcidResidue::Alanine alanine1;
    // AminoAcidResidue::Alanine alanine2;
    // cout << "Alanine = " << alanine1 << endl;

    //cout << "PROTEIN=" << protein << endl;

    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);
    DecorationSubsystem decorations(system);

    // system.adoptCompound(methane1);
    // system.adoptCompound(ethane1);
    // system.adoptCompound(alanine1);

    // system.adoptCompound(alanine1);

    // system.adoptCompound(ethane2);
    // system.adoptCompound(ethane3);
    // system.adoptCompound(methane2);

    // system.adoptCompound(protein);

    // system.adoptCompound(alanine1);
    // system.adoptCompound(alanine2);

    dumm.loadAmber99Parameters();
    //ifstream tinkerStream("../../resources/tinker_amber99_clean.prm");
    //dumm.populateFromTinkerParameterFile(tinkerStream);
    //tinkerStream.close();

    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), "Methane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), "Ethane C",   DuMM::AtomClassIndex(1),  0.03);
    dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), "Ethane H",  DuMM::AtomClassIndex(34),  -0.01);
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5001), Biotype::MethaneH().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5002), Biotype::EthaneC().getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5003), Biotype::EthaneH().getIndex());

    // Complete C-terminal group
    // alanine1.bondAtom(UnivalentAtom("2OXT", Element::Oxygen()), "bondC", 0.1250);
    // alanine1.nameAtom("1OXT", "O");

    // Add one more hydrogen to nitrogen to create unnatural 3-body cluster
    // alanine1.convertInboardBondCenterToOutboard();
    // alanine1.bondAtom(UnivalentAtom("2HN", Element::Hydrogen()), "bondN", 0.1010);

    // build acetyl alanine
    //AminoAcidResidue::Alanine alanine1;
    alanine1.assignBiotypes();
    //AcetylResidue compound1;
    //compound1.assignBiotypes(Ordinality::Initial);
    //compound1.bondCompound("ala1", alanine1, "bondC", 0.13350);
    //system.adoptCompound(compound1);

    AcetylResidue acetyl1;
    acetyl1.assignBiotypes(Ordinality::Initial);
    acetyl1.bondCompound("ala1", alanine1, "bondC", 0.135);
    acetyl1.setBondMobility(BondMobility::Rigid, "C", "ala1/N");

    // Add some more alanines to make a peptide
    AminoAcidResidue::Alanine alanine2;
    AminoAcidResidue::Alanine alanine3;
    AminoAcidResidue::Alanine alanine4;
    alanine2.assignBiotypes();
    alanine3.assignBiotypes();
    alanine4.assignBiotypes(Ordinality::Final);
    acetyl1.bondCompound("ala2", alanine2, "ala1/bondC", 0.135);

    // system.adoptCompound(acetyl1);

    // serine.assignBiotypes();
    // system.adoptCompound();
    AminoAcidResidue::Threonine thr;
    thr.assignBiotypes();
    // system.adoptCompound(thr);

//    Protein protein("AAAAAAAAA");
    Protein protein("GPRFS");
    protein.writeDefaultPdb(cout);
    protein.assignBiotypes();
    system.adoptCompound(protein);

    // bool biotypesComplete = alanine1.assignBiotypes();
    // bool biotypesComplete = alanine1.assignBiotypes(Ordinality::ANY);
    // assert(biotypesComplete);
    // alanine2.assignBiotypes();

    cout << "Methane = " << methane1 << endl;
    cout << "Ethane = " << ethane1 << endl;
    cout << "Alanine = " << alanine1 << endl;

    // alanine1.writeDefaultPdb(cout);

    system.modelCompounds();
                

    State state = system.realizeTopology();

#ifdef CREATE_VIZ_WINDOW
    Visualizer display(system);
#endif

    RungeKuttaMersonIntegrator study(system);
    //CPodesIntegrator study(system);
    study.initialize(state);

#ifdef CREATE_VIZ_WINDOW
    display.report(study.getState());
#endif

    // getchar(); // for synchronizing movie capture

    // ethane.setDihedralAngle(state, "torsion", 60*DuMM::Deg2Rad);

    //cout << protein.getAtomLocation(state, "0/CA") << endl;

    Real timeInterval = 0.005;
    for (Real time=0.0; time < (10 * timeInterval); time += timeInterval) // picoseconds
    {

        // std::cout << "Ethane torsion angle = " << DuMM::Rad2Deg * ethane.getDihedralAngle(state, "torsion") << std::endl;

        // study.stepTo(time, Infinity);
        study.stepTo(time);

#ifdef CREATE_VIZ_WINDOW
        display.report(study.getState());
#endif

        //protein.writePdbCoordinates(state, fileName);
    }

//  }
//catch (const std::exception& e)
//  {
//    printf("EXCEPTION THROWN: %s\n", e.what());
//  }
//catch (...)
//  {
//    printf("UNKNOWN EXCEPTION THROWN\n");
//  }    return 0;
}

