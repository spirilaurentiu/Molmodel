#include "SimTKmolmodel.h"
#include <iostream>

using namespace std;
using namespace SimTK;

int main() 
{
    // Start with a simple protein
    // NOTE that the first residue, number 0, will actually be an acetyl end cap.
    // The first real amino acid, an alanine, will be residue number 1
    // Protein protein("ACDEFGHIKLMNPQRSTVWY");

    // In order to extract body ids, we need to model the protein in a system
    CompoundSystem system;
	SimbodyMatterSubsystem matter(system);
    // DecorationSubsystem decorations(system);

    DuMMForceFieldSubsystem dumm(system);

    // Use hard-coded force field parameters
    dumm.loadAmber99Parameters();

    std::vector<Protein> proteins;
    {
        Protein protein1("QISQAIKYLQNNIKGFIIRQRVNDEMKVNCATLLQAAYRGHSIRANVF");
        // Protein protein1("ACDEFGHIKLMNPQRSTVWY");
        protein1.assignBiotypes();
        proteins.push_back(protein1);


        Protein protein2("AAAA");
        protein2.assignBiotypes();
        proteins.push_back(protein2);

        // adoptCompound must be done after proteins vector has stabilized
        system.adoptCompound(proteins[0]);
        system.adoptCompound( proteins[1], Vec3(-2, 0, 0) );
    }

    // Now we can model
    system.modelCompounds();

    const Protein& protein = proteins[0];

    // Next extract the body ids and locations
    int numResidues = protein.getNumResidues() - 1;
    for (ResidueInfo::Index r(0); r < numResidues; ++r)
    {
        const ResidueInfo& residue = protein.getResidue(ResidueInfo::Index(r+1));
        for (ResidueInfo::AtomIndex ra(0); ra < residue.getNumAtoms(); ++ra)
        {
            Compound::AtomIndex a = residue.getAtomIndex(ra);

            cout << residue.getPdbResidueName();
            cout << "\t";
            cout << residue.getPdbResidueNumber();
            cout << "\t";
            cout << residue.getAtomName(ra);
            cout << "\t";
            cout << "Body Index = " << protein.getAtomMobilizedBodyIndex(a);
            cout << "\t";
            cout << protein.getAtomLocationInMobilizedBodyFrame(a);

            const MobilizedBody& body = matter.getMobilizedBody(protein.getAtomMobilizedBodyIndex(a));
            const MobilizedBody& parentBody = body.getParentMobilizedBody();

            cout << "\t";
            cout << "parent body index = " << parentBody.getMobilizedBodyIndex();

            cout << endl;
        }
    }

}
