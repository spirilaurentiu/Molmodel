#ifndef MOLMODEL_PDB_H_
#define MOLMODEL_PDB_H_

#include "SimTKsimbody.h"
#include <cctype>
#include "molmodel/internal/Compound.h"
#include <map>

namespace SimTK {

namespace Pdb {
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(ResidueIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(ChainIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(ModelIndex);
}

/// Composite key for residue within a chain, 
/// composed of residue number and insertion code
class SimTK_MOLMODEL_EXPORT PdbResidueId {
public:
    explicit PdbResidueId(int num, char iCode = ' ') : residueNumber(num), insertionCode(iCode) {}

    /// < operator is required for use as a key in a hash table
    bool operator<(const PdbResidueId& other) const;

    int residueNumber;
    char insertionCode;
};

namespace Exception 
{
    class UndefinedPdbChainId : public Exception::Base
    {
    public:
        UndefinedPdbChainId(const char* fn, int ln, String chainId)
            : Exception::Base(fn, ln)
        {
            String chainString = " ";
            chainString = chainId;
            // was : chainString[0] = chainId;
            setMessage("Undefined PDB chain Id: '" + chainString + "'");
        }
        
    };

    class DuplicatePdbResidue : public Exception::Base
    {
    public:
        DuplicatePdbResidue(const char* fn, int ln, PdbResidueId pdbResidueId)
            : Exception::Base(fn, ln)
        {
            String message("Duplicate PDB Residue number ");
            message += String(pdbResidueId.residueNumber);

            setMessage(message);
        }
        
    };
    
    class UndefinedAminoAcidResidue : public Exception::Base
    {
    public:
        UndefinedAminoAcidResidue(const char* fn, int ln, String residueName)
            : Exception::Base(fn, ln)
        {
            String message("Unknown amino acid residue name '");
            message += residueName;
            message += "'";

            setMessage(message);
        }
        
    };
}

/// Location information for a PdbAtom, corresponding to one altLoc for a PdbAtom
class SimTK_MOLMODEL_EXPORT PdbAtomLocation {
public:
    explicit PdbAtomLocation(SimTK::Vec3 coords, char altLoc = ' ', SimTK::Real tFac = 0.00, SimTK::Real occ = 1.0) 
        : coordinates(coords), alternateLocationIndicator(altLoc), temperatureFactor(tFac), occupancy(occ)
    {}

    /// write columns 31 to 66 of a PDB ATOM/HETATM record at this location
    std::ostream& writePdb(std::ostream& os, const Transform& transform) const;

    const Vec3& getCoordinates() const {return coordinates;}

    char getAlternateLocationIndicator() const {return alternateLocationIndicator;}

    SimTK::Real getTemperatureFactor() const {return temperatureFactor;}

    SimTK::Real getOccupancy() const {return occupancy;}

private:

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    SimTK::Vec3 coordinates;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

    char alternateLocationIndicator;
    SimTK::Real temperatureFactor;
    SimTK::Real occupancy;
};

/// One atom, which may have more than one location, in the case of static 
/// or dynamic disorder in an X-ray structure.
class SimTK_MOLMODEL_EXPORT PdbAtom {
    friend class PdbResidue;
    typedef std::vector<SimTK::String> AtomNameList;
public:
    explicit PdbAtom(const SimTK::String& name, const Element& e);

    explicit PdbAtom( 
        const class Compound& compound, 
        const String& atomName, 
        const Transform& transform = Transform() );

    explicit PdbAtom( 
        const State& state, 
        const class Compound& compound, 
        const String& atomName, 
        const Transform& transform = Transform() );

    /// If you set this static flag true, then every ATOM line written out
    /// to a PDB file will be follwed by an ignorable REMARK line that
    /// contains the atom location to full precision. The Molmodel Pdb
    /// reader will use that line if it is present to read in more precise
    /// atom locations than is possible in a standard PDB file, which is
    /// limited to 0.001 Angstrom precision. The REMARK line has syntax
    /// REMARK-SIMTK-COORDS X Y Z.
    static void setWriteFullPrecisionLocation(bool saveFullPrecision);

    /// Return the current value of the "write full precision location" 
    /// flag.
    static bool getWriteFullPrecisionLocation();

    std::ostream& write(
        std::ostream& os, 
        int& nextAtomSerialNumber, 
        const char residueName[4], 
        PdbResidueId residueId, 
        String chainId, 
        const Transform& transform) const;

    bool hasLocation(char altLoc) const;

    bool hasLocation() const;

    Vec3 getLocation(char altLoc) const;

    Vec3 getLocation() const;

    PdbAtom& setLocation(const PdbAtomLocation& loc) {
        char alt_loc = loc.getAlternateLocationIndicator();
        if (locationIndicesById.find(alt_loc) == locationIndicesById.end()) {
            locationIndicesById[alt_loc] = locations.size();
            locations.push_back(loc);
        }
        else {
            locations[locationIndicesById.find(alt_loc)->second] = loc;
        }

        return *this;
    }

    const String& getName() const {return atomName;}

    const PdbAtomLocation& getPdbAtomLocation() const {
        assert(hasLocation());
        return locations[0];
    }

    const Vec3& getCoordinates() const {
        return getPdbAtomLocation().getCoordinates();
    }

    char getAlternateLocationIndicator() const {
        return getPdbAtomLocation().getAlternateLocationIndicator();
    }

    SimTK::Real getTemperatureFactor() const {
        return getPdbAtomLocation().getTemperatureFactor();
    }

    SimTK::Real getOccupancy() const {
        return getPdbAtomLocation().getOccupancy();
    }


protected:
    void parsePdbLine(const String& line);

    /// create upper-case four-character PDB-style atom name, given a more free-form
    /// atom name
    static String canonicalizeAtomName(const String& casualName);

    /// Try to be smart about guessing correct atom name for names that are not 4 characters long
    static std::vector<SimTK::String> generatePossibleAtomNames(SimTK::String name);

private:
    Element element;

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    SimTK::String atomName;
    typedef std::vector<PdbAtomLocation> Locations;
    Locations locations;
    std::map<char, int> locationIndicesById;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

    // Pdb writes will be followed by a REMARK that has full precision for
    // atom locations if this is true.
    static bool writeExtraPrecision;
};

/// One residue in a protein or nucleic acid, or a single molecule in the case of non-polymer structures
class SimTK_MOLMODEL_EXPORT PdbResidue {
    friend class PdbChain;
public:
    explicit PdbResidue(String name, PdbResidueId id);

    explicit PdbResidue(
        const class Compound& compound, 
        int residueNumber, 
        const Transform& transform = Transform() );

    explicit PdbResidue(
        const State& state, 
        const class Compound& compound, 
        int residueNumber, 
        const Transform& transform = Transform() );

    std::ostream& write(std::ostream& os, int& nextAtomSerialNumber, String chainId, const Transform& transform) const;

    bool hasAtom(SimTK::String argName) const; 

    const PdbAtom& getAtom(String argName) const;

    PdbAtom& updAtom(String argName);

    const PdbResidueId& getResidueId() const {return residueId;}
    const char* getName() const {return residueName;}
    int getPdbResidueNumber() const {return residueId.residueNumber;}
    char getInsertionCode() const {return residueId.insertionCode;}

    size_t getNumAtoms() const { return atoms.size(); }

    const PdbAtom& getAtom(Pdb::AtomIndex atomIx) const {
        return atoms[atomIx];
    }

    void addAtom(const PdbAtom& atom);

protected:
    void parsePdbLine(const String& line);

private:
    char residueName[4];
    PdbResidueId residueId;

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    typedef std::vector<PdbAtom> Atoms;
    Atoms atoms;
    std::map<SimTK::String, int> atomIndicesByName;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

private:
    // OBSOLETE: use getNumAtoms()
    size_t getNAtoms() const {return getNumAtoms();}
};

class Compound;

/// One molecule in a PDB structure.  One exception is that all of the water molecules in a structure
/// may share a single chainId
class SimTK_MOLMODEL_EXPORT PdbChain {
    friend class PdbModel;
public:
    explicit PdbChain(String id = " ") : chainId(id) {}

    explicit PdbChain(
        const Compound& compound,
        const Transform& transform = Transform());

    explicit PdbChain(
        const State& state, 
        const Compound& compound,
        const Transform& transform = Transform());

    std::ostream& write(std::ostream& os, int& nextAtomSerialNumber, const Transform& transform = Transform()) const;

    bool hasResidue(PdbResidueId pdbResidueId) const;

    bool hasAtom(String atomName, PdbResidueId residueId) const;

    const PdbAtom& getAtom(String atomName, PdbResidueId residueId) const;

    PdbAtom& updAtom(String atomName, PdbResidueId residueId);

    PdbResidue& appendResidue( const PdbResidue& residue )
    {
            //std::cout<<__FILE__<<":"<<__LINE__<<" You have tried to add a residue with ID : "<<residue.getResidueId().residueNumber << residue.getResidueId().insertionCode <<std::endl;//" which is incompatible with : "<<(residueIndicesById.end()).residueNumber<< std::endl;
            /*for (size_t i = residueIndicesById.begin()->second ; i <= residueIndicesById.find(residue.getResidueId())->second; i++) {
		std::cout<<__FILE__<<":"<<__LINE__<< " residue with index "<< i << " is "<< ((residueIndicesById.find(residue.getResidueId()))->first   ).residueNumber << ((residueIndicesById.find(residue.getResidueId()))->first   ).insertionCode  <<std::endl;
            }*/
        if ( residueIndicesById.find(residue.getResidueId()) != residueIndicesById.end() ) {
            std::cout<<__FILE__<<":"<<__LINE__<<" But there is an existing residue : "<<((residueIndicesById.find(residue.getResidueId()))->first   ).residueNumber << ((residueIndicesById.find(residue.getResidueId()))->first   ).insertionCode <<std::endl; //<<residueIndicesById.find(residue.getResidueId()).insertionCode<<std::endl;
            SimTK_THROW1( Exception::DuplicatePdbResidue, residue.getResidueId() );
        }

        residues.push_back( residue );
        residueIndicesById[residue.getResidueId()] = residues.size() - 1;

        return residues.back();
    }

    size_t getNumResidues() const { return residues.size(); }
    const PdbResidue& getResidue(Pdb::ResidueIndex resIx) const {
        return residues[resIx];
    }
    
    size_t getNumAtoms() const { 
        size_t atomCount = 0;
        for (Pdb::ResidueIndex r(0); r < getNumResidues(); ++r) {
            atomCount += getResidue(r).getNumAtoms();
        }

        return atomCount;
    }

    String getChainId() const {return chainId;}


protected:
    void parsePdbLine(const String& line);

private:
    String chainId;

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    typedef std::vector<PdbResidue> Residues;
    Residues residues;
    std::map<PdbResidueId, size_t> residueIndicesById;
    
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

private:
    // OBSOLETE: use getNumResidues()
    size_t getNResidues() const { return getNumResidues(); }
    // OBSOLETE: use getNumAtoms()
    size_t getNAtoms() const { return getNumAtoms(); }
};

/// PDB structure containing one or more molecules (chains), corresponding to one member of an
/// NMR ensemble of alternate structures, or to one frame of a molecular dynamics simulation.
class SimTK_MOLMODEL_EXPORT PdbModel {
    friend class PdbStructure;
public:
    explicit PdbModel(
        const Compound& compound, 
        int number = 1,
        const Transform& transform = Transform());

    explicit PdbModel(
        const State& state, 
        const Compound& compound, int number = 1,
        const Transform& transform = Transform());

    explicit PdbModel(int number) : modelNumber(number) {}

    std::ostream& write(std::ostream& os, const Transform& transform) const;

    bool hasChain(String id) const;

    bool hasAtom(String atomName, PdbResidueId residueId, String chainId) const;

    const PdbAtom& getAtom(String atomName, PdbResidueId residueId, String chainId) const;

    PdbAtom& updAtom(String atomName, PdbResidueId residueId, String chainId);

    size_t getNumChains() const {return chains.size();}
    const PdbChain& getChain(Pdb::ChainIndex chainIx) const {
        return chains[chainIx];
    }
    const PdbChain& getChain(String chainId) const;

    // OBSOLETE; TODO: remove in SimTK 2.0
    size_t getNChains() const {return getNumChains();}

protected:
    void parsePdbLine(const String& line, String longChainId  );

    PdbChain& updOrCreateChain(String chainId);

private:
    int modelNumber;

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    typedef std::vector<PdbChain> Chains;
    Chains chains;
    std::map<String, int> chainIndicesById;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

/// Complete PDB file, possibly including multiple MODELS, in the case of NMR structures or 
/// molecular dynamics trajectories
class SimTK_MOLMODEL_EXPORT PdbStructure {
public:

    /// Construct PdbStructure based on initial default configuration of a Compound
    explicit PdbStructure(
        const Compound& compound,
        const Transform& transform = Transform());

    /// Constructure PdbStructure for a Compound using a particular State
    explicit PdbStructure(
        const State& state, 
        const Compound& compound,
        const Transform& transform = Transform());

    explicit PdbStructure(std::istream& pdbFile);

    bool hasAtom(String atomName, PdbResidueId residueId, String chainId) const;

    const PdbAtom& getAtom(String atomName, PdbResidueId residueId, String chainId) const;

    PdbAtom& updAtom(String atomName, PdbResidueId residueId, String chainId);

    std::ostream& write(std::ostream& os, Transform transform = Transform()) const;

    size_t getNumModels() const {return models.size();}
    const PdbModel& getModel(Pdb::ModelIndex modelIx) const {
        return models[modelIx];
    }

private:

    // avoid dll export warnings for these private types
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif

    typedef std::vector<PdbModel> Models;
    Models models;

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

private:
    // OBSOLETE; use getNumModels() instead
    size_t getNModels() const {return getNumModels();}
};

} // namespace SimTK

SimTK_MOLMODEL_EXPORT std::ostream& operator<<(std::ostream& os, const SimTK::PdbModel& pdbModel);
SimTK_MOLMODEL_EXPORT std::ostream& operator<<(std::ostream& os, const SimTK::PdbStructure& pdbStructure);

#endif // MOLMODEL_PDB_H_
