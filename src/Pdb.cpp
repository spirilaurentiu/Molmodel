//
//
//
#include "molmodel/internal/Pdb.h"
#include "molmodel/internal/Compound.h"
#include <stdexcept>

namespace SimTK {

// Operator < is required for use as a key in a hash table.
bool PdbResidueId::operator<(const PdbResidueId& other) const {
    if (residueNumber < other.residueNumber) return true;
    else if (residueNumber == other.residueNumber) 
        return insertionCode < other.insertionCode;
    else return false;
}

// Write columns 31 to 66 of a PDB ATOM/HETATM record at this location.
std::ostream& PdbAtomLocation::writePdb(std::ostream& os, const Transform& transform) const 
{
    Vec3 modCoords = transform * coordinates;

    // X, Y, Z coordinates - 
    // convert internal nanometers to Angstroms by multiplying by ten

    // S Flores modified so that precision is reduced to accomodate larger numbers
    
    if (((modCoords[0] * 10.0) > -1000) && ((modCoords[0] * 10.0) < 10000))  
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    else if (((modCoords[0] * 10.0) > -10000) && ((modCoords[0] * 10.0) < 100000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    else if (((modCoords[0] * 10.0) > -100000) && ((modCoords[0] * 10.0) < 1000000)) {
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(1);
    }
    else {
        std::cout<<__FILE__<<":"<<__LINE__<<" : Trouble writing coordinates = "<<modCoords[0]<<", "<<modCoords[1]<<", "<<modCoords[2]<<std::endl;
        SimTK_ASSERT_ALWAYS(0,"An X-coordinate is too large to print in PDB format");}
    os << modCoords[0] * 10.0; // X

    if (((modCoords[1] * 10.0) > -1000) && ((modCoords[1] * 10.0) < 10000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    else if (((modCoords[1] * 10.0) > -10000) && ((modCoords[1] * 10.0) < 100000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    else if (((modCoords[1] * 10.0) > -100000) && ((modCoords[1] * 10.0) < 1000000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(1);
    else SimTK_ASSERT_ALWAYS(0,"A Y-coordinate is too large to print in PDB format");
    os << modCoords[1] * 10.0; // Y

    if (((modCoords[2] * 10.0) > -1000) && ((modCoords[2] * 10.0) < 10000))  
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(3);
    else if (((modCoords[2] * 10.0) > -10000) && ((modCoords[2] * 10.0) < 100000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    else if (((modCoords[2] * 10.0) > -100000) && ((modCoords[2] * 10.0) < 1000000))
        os << std::right << std::setw(8) << std::setiosflags(std::ios::fixed) << std::setprecision(1);
    else SimTK_ASSERT_ALWAYS(0,"A Z-coordinate is too large to print in PDB format");
    os << modCoords[2] * 10.0; // Z

    // occupancy
    os << std::right << std::setw(6) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << occupancy;

    // thermal parameter
    os << std::right << std::setw(6) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << temperatureFactor;

    return os;
}

/*static*/ bool PdbAtom::writeExtraPrecision = false;

/*static*/void PdbAtom::setWriteFullPrecisionLocation(bool saveFullPrecision) 
{   writeExtraPrecision = saveFullPrecision; }
/*static*/bool PdbAtom::getWriteFullPrecisionLocation() 
{   return writeExtraPrecision; }

PdbAtom::PdbAtom(const SimTK::String& name, const Element& e) : element(e) 
{
    atomName = canonicalizeAtomName(name);
}

PdbAtom::PdbAtom( 
                 const class Compound& compound, 
                 const String& n, 
                 const Transform& transform )
    : element(compound.getAtomElement(n))
{
    atomName = canonicalizeAtomName(n);
    locations.push_back( PdbAtomLocation(transform * compound.calcDefaultAtomLocationInCompoundFrame(n)) );
}

PdbAtom::PdbAtom( 
                 const State& state, 
                 const class Compound& compound, 
                 const String& n, 
                 const Transform& transform )
    : element(compound.getAtomElement(n))
{
    atomName = canonicalizeAtomName(n);
    Compound::AtomIndex atomIndex = compound.getAtomIndex(n);
  
    locations.push_back( PdbAtomLocation(transform * compound.calcAtomLocationInGroundFrame(state, atomIndex)) );
}

std::ostream& PdbAtom::write(
    std::ostream& os, 
    int& nextAtomSerialNumber, 
    const char residueName[4], 
    PdbResidueId residueId, 
    String chainId, 
    const Transform& transform) const 
{
    // remember initial stream state
    std::ios_base::fmtflags originalFlags = os.flags();

    // record name
    os << "ATOM  ";

    // serial number
    if (nextAtomSerialNumber < 100000) 
    // added a std::dec just to be sure we're back to decimal.
        {os << std::right << std::dec << std::setw(5) << nextAtomSerialNumber;}
    else 
	// Following Mike Kuiper's suggestion, if atom numbers are too large to fit in the 5-character alloted column, write in hex.
	// making sure to set os back to std::dec so as to not muck up all the other fields in the ATOM record
        {os << std::right << std::hex << std::setw(5) << nextAtomSerialNumber << std::dec;}
    ++nextAtomSerialNumber;

    os << " "; // blank at column 12

    os << std::setw(4) << atomName;

    // alternate location indicator
    os << " ";

    // residue name
    os << std::setw(3) << residueName;

    // blank at column 22
    os << " "; 

    // chain id
    //os << chainId;
    if (chainId.length() == 1) {
	 os << chainId;
    } else if (chainId.length() > 1){
	    os << String (" "); // write out single whitespace.  We have a long chain ID so we will override the PDB chain ID column.
    } else {
	    std::cout<<__FILE__":"<<__LINE__<<" chainId is less than 1 character long!"<<std::endl;
	    exit(1);
    }

    // residue number
    os << std::right << std::setw(4) << residueId.residueNumber;

    // residue insertion code at position 27
    os << residueId.insertionCode;

    // blank columns 28-30
    os << "   ";

    const PdbAtomLocation& location = locations[0]; // TODO - select particular alternate location
    location.writePdb(os, transform);

    // blank columns 67-76
    os << "          ";

    // element
    // convert to upper case
    std::string e = element.getSymbol();
    std::transform(e.begin(), e.end(), e.begin(), (int(*)(int)) toupper);
    os << std::right << std::setw(2) << e;

    os << std::endl; 

    // Samuel Flores added this ignorable high-precision coordinate output.
    // TODO: this should be optional.
    if (writeExtraPrecision) {
        os << "REMARK-SIMTK-COORDS" << std::setprecision(17);
        for (int i=0; i<3; ++i) os << " " << 10*location.getCoordinates()[i];
        os << std::endl;
    }

    // restore stream state
    os.flags(originalFlags);

    return os;
}

bool PdbAtom::hasLocation(char altLoc) const {
    return locationIndicesById.find(altLoc) != locationIndicesById.end();
}

bool PdbAtom::hasLocation() const {
    return locations.size() > 0;
}

Vec3 PdbAtom::getLocation(char altLoc) const {
    return locations[locationIndicesById.find(altLoc)->second].getCoordinates();
}

Vec3 PdbAtom::getLocation() const {
    return locations[0].getCoordinates();
}

// Process an ATOM/HETATM line to get location and other information about
// this atom, or revise the location from a REMARK-SIMTK-COORDS line that
// provides higher precision.
void PdbAtom::parsePdbLine(const String& line)
{
    if ( (line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM") ) 
    {
        if (line.length() < 53) {
            throw std::runtime_error("Pdb atom line too short");
        }

        char altLoc = line[16];

        if (! hasLocation(altLoc)) {
            double x, y, z;
            SimTK::Real tempFac = 0.0; // default value in case not present
            SimTK::Real occ = 1.0; // default value in case not present

            // Is is an error to lack coordinates
            std::istringstream zeroes(String("00000000000000000"));
            std::istringstream myx ((line.substr(30, 8)) + String("00000000000000000") ); 
            myx >>x; 
            std::istringstream(line.substr(38, 8)) >> y;
            std::istringstream myy ((line.substr(38, 8)) + String("00000000000000000") ); 
            myy >> y;
            std::istringstream myz ((line.substr(46, 8)) + String("00000000000000000") ); 
            std::istringstream(line.substr(46, 8)) >> z;
            myz >> z;
            myz >> z;
            // 0.0000001 had no effect 
            // 0.00001 had an effect
            // If line is short, set temperature factor and occupancy to default values
            if (line.length() > 66) {
                std::istringstream(line.substr(54, 6)) >> occ;
                std::istringstream(line.substr(60, 6)) >> tempFac;
            }

            locationIndicesById[altLoc] = locations.size();

            // Convert from Angstroms to nanometers
            locations.push_back(PdbAtomLocation(0.1 * SimTK::Vec3(x,y,z), altLoc, tempFac, occ));
        }
        else {
            std::cout<<__FILE__<<":"<<__LINE__<<" read line : "<<line<<std::endl;
            std::cout<<__FILE__<<":"<<__LINE__<<" duplicate location! Please check your residue numbers, atom names, etc."<<std::endl; exit(1); 
            assert(!"duplicate location");
        }
    } else if ((line.substr(0,19) == "REMARK-SIMTK-COORDS")) {
        // We require that a "REMARK-SIMTK-COORDS" follow its corresponding 
        // ATOM or HETATM record, with no other intervening ATOM or HETATM records.       
        //      syntax: REMARK-SIMTK-COORDS X Y Z
        // where REMARK-SIMTK-COORDS should be left-justified, and records should be 
        // separated by one or more whitespaces.

        Real x, y, z;
        std::stringstream(line.substr(19)) >>  x>>y>>z;
        //std::cout<< std::setprecision(17) <<__FILE__<<":"<<__LINE__<<" Parsing PDB line : >"<<line<<"< .. This is a double precision record. Coordinates are interpreted as : "<<x<<","<<y<<","<<z<<std::endl;
        const Vec3 highResNm(x/10,y/10,z/10); // convert Å to nm
        //std::cout<< std::setprecision(17) <<__FILE__<<":"<<__LINE__<<" put into nm vector : "<<highResNm<<std::endl;
        // Sanity-check the high precision coordinates. Standard PDB coordinates
        // have 6 digits of precision with no more than 3 to the right of the
        // decimal (i.e., .001Å precision), less precision if the numbers are large.
        // We'll check that the high- and low-res numbers agree to .01% or .002Å, 
        // whichever is less stringent.
        const Vec3& lowResNm = locations.back().getCoordinates();
        Vec3 tolNm;
        const Real minTolNm = .0002; // .002Å converted to nm
        for (int i=0; i<3; ++i) tolNm[i] = std::max(std::abs(lowResNm[i]/1000), minTolNm);
        const Vec3 absDiffNm = (highResNm - lowResNm).abs();

	    if (absDiffNm[0] > tolNm[0] || absDiffNm[1] > tolNm[1] || absDiffNm[2] > tolNm[2]) {
            throw std::runtime_error("The coordinates given in the REMARK-SIMTK-COORDS record "
                                     "disagree with those on the preceding ATOM or HETATM record "
                                     "by too much. Please correct this issue or remove the "
                                     "REMARK-SIMTK-COORDS record."); 
        }

        // Copy altLoc, tempFac and occ from preceding PdbAtomLocation.
        PdbAtomLocation tempPdbAtomLocation(highResNm, locations.back().getAlternateLocationIndicator(), 
                                            locations.back().getTemperatureFactor(), locations.back().getOccupancy() );
        locations.back() = tempPdbAtomLocation; // Replace the latest location with this more precise version.                        
    }
}

/// create upper-case four-character PDB-style atom name, given a more free-form
/// atom name
/* static */ String PdbAtom::canonicalizeAtomName(const String& casualName) 
{
    String canonicalName(casualName);

    // Remove everything before a "/" character in Compound::AtomName
    if (canonicalName.find_first_of("/") != String::npos)
    {
        int slashPos = canonicalName.find_last_of("/");
        canonicalName = canonicalName.substr(slashPos + 1);
    }

    // Convert to upper case
    std::transform(
        canonicalName.begin(), canonicalName.end(), 
        canonicalName.begin(), (int(*)(int)) std::toupper);

    // Perhaps 4 characters with letter in second position is already canonical
    if (   (canonicalName.length() == 4) && 
           (canonicalName[1] >= 'A') &&
           (canonicalName[1] <= 'Z')
       )
        return canonicalName;


    // There are some 4-character hydrogen atom names, so don't prepend
    // a space to those
    if (   ('H' == canonicalName[0]) 
        && (canonicalName.length() >= 4)
        && (' ' != canonicalName[3]) )
    {} // don't prepend a space in the next step

    // If first character is a letter, prepend a space
    // unless it is one of those four character hydrogens above
    else if ( ('A' <= canonicalName[0]) 
           && ('Z' >= canonicalName[0]) )
    {
           canonicalName = " " + canonicalName;
    }

    // Truncate to 4 characters
    if (canonicalName.size() > 4) canonicalName = canonicalName.substr(0, 4);

    // Pad to 4 characters
    while (canonicalName.size() < 4) canonicalName += " ";

    return canonicalName;
}

/// Try to be smart about guessing correct atom name for names that are not 4 characters long
/* static */ std::vector<SimTK::String> PdbAtom::generatePossibleAtomNames(SimTK::String name)
{
    AtomNameList possibleNames;

    // 1) convert to upper case
    std::transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::toupper);

    if (4 == name.length()) 
        possibleNames.push_back(name); // easy case

    else if (4 < name.length())
    {
        // truncate name
        possibleNames.push_back(name.substr(0, 4));
    }

    else
    {
        // 1) append spaces to make four characters
        SimTK::String name1 = name;
        while (name1.length() < 4) name1 += " ";
        possibleNames.push_back(name1);

        // 2) prepend one space, append rest
        SimTK::String name2 = String(" ") + name;
        while (name2.length() < 4) name2 += " ";
        possibleNames.push_back(name2);
    }

    return possibleNames;
}

PdbResidue::PdbResidue(String name, PdbResidueId id)
    : residueId(id)
{
    residueName[0] = name[0];
    residueName[1] = name[1];
    residueName[2] = name[2];
    residueName[3] = 0;
}

PdbResidue::PdbResidue(const Compound& compound, int resNum, const Transform& transform)
    : residueId(resNum)
{
    const String& name = compound.getPdbResidueName();

    residueName[0] = name[0];
    residueName[1] = name[1];
    residueName[2] = name[2];
    residueName[3] = 0;

    for (Compound::AtomIndex aIx(0); aIx < compound.getNumAtoms(); ++aIx)
    {
        Compound::AtomName atomName = compound.getAtomName(aIx);
        addAtom(PdbAtom(compound, atomName, transform));
    }
}

PdbResidue::PdbResidue(
    const State& state, 
    const Compound& compound, 
    int resNum, 
    const Transform& transform)
    : residueId(resNum)
{
    const String& name = compound.getPdbResidueName();

    residueName[0] = name[0];
    residueName[1] = name[1];
    residueName[2] = name[2];
    residueName[3] = 0;

    for (Compound::AtomIndex aIx(0); aIx < compound.getNumAtoms(); ++aIx)
    {
        Compound::AtomName atomName = compound.getAtomName(aIx);
        addAtom(PdbAtom(state, compound, atomName, transform));
    }
}

std::ostream& PdbResidue::write(std::ostream& os, int& nextAtomSerialNumber,String chainId, const Transform& transform) const 
{
    Atoms::const_iterator atomI;
    for (atomI = atoms.begin(); atomI != atoms.end(); ++atomI) 
    {
        atomI->write(os, nextAtomSerialNumber, residueName, residueId, chainId, transform);
    }

    return os;
}

bool PdbResidue::hasAtom(SimTK::String argName) const 
{
    // try variations of atom name spelling
    PdbAtom::AtomNameList possibleNames = PdbAtom::generatePossibleAtomNames(argName);
    PdbAtom::AtomNameList::const_iterator nameI;
    for (nameI = possibleNames.begin(); nameI != possibleNames.end(); ++nameI)
    {
        if (atomIndicesByName.find(*nameI) != atomIndicesByName.end()) return true;
    }

    return false;
}

const PdbAtom& PdbResidue::getAtom(String argName) const 
{
    assert(hasAtom(argName));

    // try variations of atom name spelling
    PdbAtom::AtomNameList possibleNames = PdbAtom::generatePossibleAtomNames(argName);
    PdbAtom::AtomNameList::const_iterator nameI;
    int atomIndex = -1;
    for (nameI = possibleNames.begin(); nameI != possibleNames.end(); ++nameI)
    {
        if (atomIndicesByName.find(*nameI) != atomIndicesByName.end()) {
            atomIndex = atomIndicesByName.find(*nameI)->second;
            break;
        }
    }
    return atoms[atomIndex];
}

PdbAtom& PdbResidue::updAtom(String argName) 
{
    assert(hasAtom(argName));

    // try variations of atom name spelling
    PdbAtom::AtomNameList possibleNames = PdbAtom::generatePossibleAtomNames(argName);
    PdbAtom::AtomNameList::const_iterator nameI;
    int atomIndex = -1;
    for (nameI = possibleNames.begin(); nameI != possibleNames.end(); ++nameI)
    {
        if (atomIndicesByName.find(*nameI) != atomIndicesByName.end()) {
            atomIndex = atomIndicesByName.find(*nameI)->second;
            break;
        }
    }
    return atoms[atomIndex];
}

void PdbResidue::parsePdbLine(const String& line)
{
    if ( (line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM") )  
    {
        String atomName = line.substr(12, 4);

        if (! hasAtom(atomName)) 
        {
            size_t firstLetterPos = atomName.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
            assert(firstLetterPos != std::string::npos); // No letter in atom name (caps only)

            // If the letter appears in the first position of 4, we assume it is a two-letter
            // element symbol, in which case we want the second letter to be lower case.
            String elementSymbol;
            if (firstLetterPos == 0) {
                elementSymbol = atomName.substr(0,2);
                elementSymbol[1] = tolower(elementSymbol[1]);
            } else {
                // Only room for a one-element name now.
                elementSymbol = atomName.substr(firstLetterPos,1);
            }

            // Four character atom names starting with "H"
            // are probably hydrogen, even if they start with "HO", which should
            // properly be Holmium.
            if ( ('H' == atomName[0]) && (std::string::npos == atomName.find_first_of(" ")) )
                elementSymbol = "H";

            const Element& element = Element::getBySymbol(elementSymbol);

            addAtom(PdbAtom(atomName, element));
        }

        atoms[atomIndicesByName[atomName]].parsePdbLine(line);
    } else if (line.substr(0,19) == "REMARK-SIMTK-COORDS"){
        // Assume the last atom added to the "atoms" vector corresponds to the 
        // one from the last ATOM record read. This should be justified since 
        // the addAtom method above contains a push_back call. This could fail 
        // if the addAtom method changes. The only solution to that would be to 
        // add atom names (and probably residue numbers) to the
        // REMARK-SIMTK-COORDS record .. which seems like a pain. (scf)
        atoms.back().parsePdbLine(line);
    }
}

void PdbResidue::addAtom(const PdbAtom& atom)
{
    const String& atomName = atom.getName();
    atomIndicesByName[atomName] = atoms.size();
    atoms.push_back(atom);
}

PdbChain::PdbChain(const Compound& compound,
        const Transform& transform)
    : chainId( compound.getPdbChainId() )
{
    int defaultNextResidueNumber = 1;
    std::cout<<__FILE__<<":"<<__LINE__<<" residues.size() : "<<residues.size()<<std::endl;
    for (const auto& residue : residues) {
        std::cout<<__FILE__<<":"<<__LINE__<<" "<<residue.getResidueId().residueNumber<<std::endl;
    } 
    if (residues.size() > 0)
        defaultNextResidueNumber = residues.back().getResidueId().residueNumber + 1;

    // Delegate chain-level PdbStructure formation to Compound
    // because chain-construction is polymorphic on Compound type
    // e.g. small molecules and residues differ from Proteins and RNAs.
    compound.populateDefaultPdbChain(*this, defaultNextResidueNumber, transform);
}

PdbChain::PdbChain(
        const State& state, 
        const Compound& compound,
        const Transform& transform)
    : chainId( compound.getPdbChainId() )
{
    int defaultNextResidueNumber = 1;

// recall from Pdb.h:
// 370     typedef std::vector<PdbResidue> Residues;
// 371     Residues residues;
// 372     std::map<PdbResidueId, size_t> residueIndicesById;


    if (residues.size() > 0)
        defaultNextResidueNumber = residues.back().getResidueId().residueNumber + 1;

    // Delegate chain-level PdbStructure formation to Compound
    // because chain-construction is polymorphic on Compound type
    // e.g. small molecules and residues differ from Proteins and RNAs.
    compound.populatePdbChain(state, *this, defaultNextResidueNumber, transform);
}

std::ostream& PdbChain::write(std::ostream& os, int& nextAtomSerialNumber, const Transform& transform) const 
{

    if (getChainId().length() > 1) {
	// we have a long chain ID, so need to override the 1-character limit of the PDB ATOM record format.
	os << "REMARK-SimTK-long-chainId "<<getChainId() <<std::endl;
    }

    Residues::const_iterator residueI;
    for (residueI = residues.begin(); residueI != residues.end(); ++residueI) 
    {
        residueI->write(os, nextAtomSerialNumber, getChainId(), transform);
    }
    // S. Flores added TER tags:
    os << "TER" << std::endl;

    if (getChainId().length() > 1) {
	    //  We need to turn off the long chain ID override :
	    os << "REMARK-SimTK-long-chainId"<<std::endl;
    }

    return os;
}

bool PdbChain::hasResidue(PdbResidueId pdbResidueId) const {
    return residueIndicesById.find(pdbResidueId) != residueIndicesById.end();
}

bool PdbChain::hasAtom(String atomName, PdbResidueId residueId) const {
    if (!hasResidue(residueId)) return false;
    return residues[residueIndicesById.find(residueId)->second].hasAtom(atomName);
}

const PdbAtom& PdbChain::getAtom(String atomName, PdbResidueId residueId) const {
    assert(hasAtom(atomName, residueId));
    return residues[residueIndicesById.find(residueId)->second].getAtom(atomName);
}

PdbAtom& PdbChain::updAtom(String atomName, PdbResidueId residueId) {
    assert(hasAtom(atomName, residueId));
    return residues[residueIndicesById.find(residueId)->second].updAtom(atomName);
}

void PdbChain::parsePdbLine(const String& line)
{
    if ( (line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM") ) 
    {
        int residueNumber;
        std::istringstream(line.substr(22, 24)) >> residueNumber;
        char insertionCode = line[26];
        PdbResidueId residueId(residueNumber, insertionCode);

        String residueName = line.substr(17, 3);

        if (! hasResidue(residueId)) {
            residueIndicesById[residueId] = residues.size();
            residues.push_back(PdbResidue(residueName, residueId));
        }

        residues[residueIndicesById[residueId]].parsePdbLine(line);
    } else if ( (line.substr(0,19 ) == "REMARK-SIMTK-COORDS")) { 
        // We assume the residue corresponding to the last ATOM record read is also the 
        // last residue in the "residues" vector. If the residues.push_back call above 
        // changes, this would have to be changed. (scf)
        residues.back().parsePdbLine(line);
    }
}

PdbModel::PdbModel(const Compound& compound, int number,
        const Transform& transform)
    : modelNumber(number)
{
    String chainId = compound.getPdbChainId();

    assert( 0 == chains.size() );
    chainIndicesById[chainId] = chains.size();
    chains.push_back(PdbChain(compound, transform));
}

PdbModel::PdbModel(
        const State& state, 
        const Compound& compound, 
        int number,
        const Transform& transform)
    : modelNumber(number)
{
    String chainId = compound.getPdbChainId();

    assert( 0 == chains.size() );
    chainIndicesById[chainId] = chains.size();
    chains.push_back(PdbChain(state, compound, transform));
}

std::ostream& PdbModel::write(std::ostream& os, const Transform& transform) const 
{
    int nextAtomSerialNumber = 1; // resets to one for each MODEL
    Chains::const_iterator chainI;
    for (chainI = chains.begin(); chainI != chains.end(); ++chainI) 
    {
        chainI->write(os, nextAtomSerialNumber, transform);
    }

    return os;
}

bool PdbModel::hasChain(String id) const {
    return chainIndicesById.find(id) != chainIndicesById.end();
}

bool PdbModel::hasAtom(String atomName, PdbResidueId residueId, String chainId) const {
    if ( !hasChain(chainId) ) return false;
    return chains[chainIndicesById.find(chainId)->second].hasAtom(atomName, residueId);
}

const PdbAtom& PdbModel::getAtom(String atomName, PdbResidueId residueId, String chainId) const {
    assert(hasAtom(atomName, residueId, chainId));
    return chains[chainIndicesById.find(chainId)->second].getAtom(atomName, residueId);
}

PdbAtom& PdbModel::updAtom(String atomName, PdbResidueId residueId, String chainId) {
    assert(hasAtom(atomName, residueId, chainId));
    return chains[chainIndicesById.find(chainId)->second].updAtom(atomName, residueId);
}

//void PdbModel::parsePdbLine(const String& line)
void PdbModel::parsePdbLine(const String& line, String longChainId = String(" "))
{
    if ( (line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM") )
    {
    String chainId = line.substr(21,1);
    //String chainId = line[21];
    if (longChainId.compare(" ") != 0) {
	    if (chainId.compare(" ") != 0) { // if long chain ID's are in use, column 22 is expected to contain a blank.  Anything different is an indicator that the user doesn't understand how to use long chain ID's.
		    std::cout<<__FILE__":"<<__LINE__<<" Bad syntax! We have earlier read a \"REMARK-SimTK-long-chainId\" record which set a long chain ID of "<<longChainId<<". But the other part of this is that all subsequent ATOM or HETATM records should have a blank (\" \") in column 22. You have : "<<line[21]<<".  Please correct this."<<std::endl;
		    exit(1);
	    }
	    chainId = longChainId;
    }; 
      
        updOrCreateChain(chainId).parsePdbLine(line);
    } else if (line.substr(0,19) == "REMARK-SIMTK-COORDS") {
        // We assume the chain corresponding to the last ATOM record read is also the 
        // last chain in the "chains" vector. If the chains.push_back call in 
        // updOrCreateChain changes, this would have to be changed. (scf)
        chains.back().parsePdbLine(line);
    }
}

PdbChain& PdbModel::updOrCreateChain(String chainId)
{
    if (! hasChain(chainId) ) {
        chainIndicesById[chainId] = chains.size();
        chains.push_back(PdbChain(chainId));
    }

    return chains[chainIndicesById[chainId]];
}

const PdbChain& PdbModel::getChain(String chainId) const
{
    if (! hasChain(chainId) )
        SimTK_THROW1(Exception::UndefinedPdbChainId, chainId);

    assert( chainIndicesById.find(chainId) != chainIndicesById.end() );
    int chainIndex = chainIndicesById.find(chainId)->second;
    assert( int(chains.size()) > chainIndex );

    return chains[chainIndex];
}

PdbStructure::PdbStructure(const Compound& compound,
        const Transform& transform) {
    models.push_back(PdbModel(compound, 1, transform));
}

PdbStructure::PdbStructure(
        const State& state, 
        const Compound& compound,
        const Transform& transform) {
    models.push_back(PdbModel(state, compound, 1, transform));
}

PdbStructure::PdbStructure(std::istream& pdbFile) 
{
    char lineBuffer[300];

    //scf added
    SimTK_ASSERT_ALWAYS(pdbFile.good(),
        "PdbStructure::ctor(): can't read pdb file -- was it missing?");

    String longChainId (" ");
    while (pdbFile.good()) 
    {


        pdbFile.getline(lineBuffer, 200);
        String line(lineBuffer);

        // Create a new PdbModel if MODEL card is seen
        // (ignore ENDMDL records)
        if (line.substr(0, 6) == "MODEL")
        {
            // If previous model is still empty, don't create a new one
            if ( (models.size() > 0) && (models.back().chains.size() == 0) ); 

            // Otherwise, create a new model
            else models.push_back( PdbModel(models.size() + 1) );
        }

        // Stop reading if "END" record is seen (but not "ENDMDL").
        else if ( (line.substr(0, 3) == "END") && (line.substr(0, 6) != "ENDMDL") ) 
        {
            break;
        }


       	else if  (line.substr(0,25) == "REMARK-SimTK-long-chainId") {
            std::stringstream u;
            u.str(line);
	    std::vector <String> mystring;
            mystring.clear();
            //char * temp;
            //u.getline(temp,100);
            //std::cout<<__FILE__<<":"<<__LINE__<<" : "<<temp<<std::endl; 
            mystring.push_back(""); // must create three vector elements.
            mystring.push_back("");
            mystring.push_back("");
	    u>>mystring[0]>>mystring[1]>>mystring[2]	;
	    if (mystring[2].length()>0) { // There should only be one argument following "REMARK-SimTK-long-chainId"
		std::cout<<__FILE__<<":"<<__LINE__<<" Bad syntax!  REMARK-SimTK-long-chainId should be followed by either nothing, or a whitespace and then a multi-character chain ID.  If you provide a multi-character chain ID, this overrides the chain ID specified in column 22 of subsequent ATOM records.  If you do not provide a multi-character chain ID, then you turn off this override feature, and go back to reading the chain ID from column 22 of the ATOM record. "<<std::endl; 
		std::cout<<__FILE__<<":"<<__LINE__<<" In this case you have provided an additional argument, "<<mystring[2]<<", which is not permitted."<<std::endl;
		exit(1);
	    }
            
	    if (mystring[1].length() == 1) { // a long chain ID should not be exactly one character long.
		std::cout<<__FILE__<<":"<<__LINE__<<" Bad syntax!  REMARK-SimTK-long-chainId should be followed by either nothing, or a whitespace and then a multi-character chain ID.  If you provide a multi-character chain ID, this overrides the chain ID specified in column 22 of subsequent ATOM records.  If you do not provide a multi-character chain ID, then you turn off this override feature, and go back to reading the chain ID from column 22 of the ATOM record. "<<std::endl; 
		std::cout<<__FILE__<<":"<<__LINE__<<" In this case you have provided a single-character chain ID, "<<mystring[1]<<", which is not permitted. You should just use column 22 as per the PDB ATOM record format. Long chain ID's should be 2 characters or longer."<<std::endl;
		exit(1);
	    }
       
            else if (mystring[1].length() == 0) { // no arguments means set long chain ID back to default value of " ".
                longChainId = " ";
            } else if (mystring[1].length() > 1) { //everything is OK!
	    } else {std::cout<<__FILE__<<":"<<__LINE__<<" Unexplained error!"<<std::endl; exit(1);} // just being paranoid..
            // If we got this far, the long chain ID is kosher.
	    longChainId = mystring[1];
        }


        else if ( (line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM") ||
                  (line.substr(0,19) == "REMARK-SIMTK-COORDS"))                        
        {
            if (models.size() < 1) models.push_back(PdbModel(1));
            models.back().parsePdbLine(line,longChainId); // back() returns a ref to the last element in the vector
        }
    }
}

bool PdbStructure::hasAtom(String atomName, PdbResidueId residueId, String chainId) const {
    if (models.size() < 1) return false;
    return models[0].hasAtom(atomName, residueId, chainId);
}

const PdbAtom& PdbStructure::getAtom(String atomName, PdbResidueId residueId, String chainId) const {
    assert(models.size() > 0);
    return models[0].getAtom(atomName, residueId, chainId);
}

PdbAtom& PdbStructure::updAtom(String atomName, PdbResidueId residueId, String chainId) {
    assert(models.size() > 0);
    return models[0].updAtom(atomName, residueId, chainId);
}

std::ostream& PdbStructure::write(std::ostream& os, Transform transform) const 
{
    Models::const_iterator modelI;
    for (modelI = models.begin(); modelI != models.end(); ++modelI) 
    {
        if (models.size() > 1) 
        {
            os << "MODEL     ";
            os << std::setw(4) << modelI->modelNumber; // right justified (default)
            os << std::endl;
        }

        modelI->write(os, transform);

        if (models.size() > 1) os << "ENDMDL" << std::endl;
    }
    os << "END" << std::endl;

    return os;
}

} // namespace SimTK

std::ostream& operator<<(std::ostream& os, const SimTK::PdbModel& pdbModel) {
	pdbModel.write(os, SimTK::Transform());
	return os;
}

std::ostream& operator<<(std::ostream& os, const SimTK::PdbStructure& pdbStructure) {
	pdbStructure.write(os, SimTK::Transform());
	return os;
}
