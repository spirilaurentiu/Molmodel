#ifndef SimTK_MOLMODEL_IONS_H_
#define SimTK_MOLMODEL_IONS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Christopher Bruns                                *
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


#include "SimTKsimbody.h"

#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"

#include <iosfwd> // declare ostream without all the definitions

namespace SimTK {
    
/// Li<sup>+</sup> lithium ion with +1 charge
class LithiumIon : public Molecule {
public:
    LithiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("LI ");

        setBaseAtom( "Li+", Biotype::get("Lithium Ion", "Li+") );

        setCompoundName("Lithium Ion");
    }
	
	static void instantiateBiotypes() {
        if (! Biotype::exists("Lithium Ion", "Li+") )
            Biotype::defineBiotype(Element::Lithium(), 0, "Lithium Ion", "Li+");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Li+ Lithium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Lithium().getAtomicNumber(),
				0, // no covalent bonds
				1.137, // radius
				0.0183 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Lithium Ion", "Li+").getIndex() );
   }
};

/// Na<sup>+</sup> sodium ion with +1 charge
class SodiumIon : public Molecule {
public:
    SodiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("NA ");

        setBaseAtom( "Na+", Biotype::get("Sodium Ion", "Na+") );

        setCompoundName("Sodium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Sodium Ion", "Na+") )
            Biotype::defineBiotype(Element::Sodium(), 0, "Sodium Ion", "Na+");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		if (! dumm.hasAtomClass("Na+ Sodium Ion") ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				"Na+ Sodium Ion",
				Element::Sodium().getAtomicNumber(),
				0, // no covalent bonds
				1.8680, // radius
				0.00277 // well depth
				);
		}

		if (! dumm.hasChargedAtomType("Na+ Sodium Ion") ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				"Na+ Sodium Ion", 
				dumm.getAtomClassIndex("Na+ Sodium Ion"),
				1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex("Na+ Sodium Ion"), Biotype::get("Sodium Ion", "Na+").getIndex() );
   }

};

/// K<sup>+</sup> potassium ion with +1 charge
class PotassiumIon : public Molecule {
public:
    PotassiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("K  ");

        setBaseAtom( "K+", Biotype::get("Potassium Ion", "K+") );

        setCompoundName("Potassium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Potassium Ion", "K+") )
            Biotype::defineBiotype(Element::Potassium(), 0, "Potassium Ion", "K+");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "K+ Potassium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Potassium().getAtomicNumber(),
				0, // no covalent bonds
				2.6580, // radius
				0.000328 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Potassium Ion", "K+").getIndex() );
   }
};

/// Rb<sup>+</sup> rubidium ion with +1 charge
class RubidiumIon : public Molecule {
public:
    RubidiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("RB ");

        setBaseAtom( "Rb+", Biotype::get("Rubidium Ion", "Rb+") );

        setCompoundName("Rubidium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Rubidium Ion", "Rb+") )
            Biotype::defineBiotype(Element::Rubidium(), 0, "Rubidium Ion", "Rb+");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Rb+ Rubidium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Rubidium().getAtomicNumber(),
				0, // no covalent bonds
				2.956, // radius
				0.00017 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Rubidium Ion", "Rb+").getIndex() );
   }
};

/// Cs<sup>+</sup> cesium ion with +1 charge
class CesiumIon : public Molecule {
public:
    CesiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("CS ");

        setBaseAtom( "Cs+", Biotype::get("Cesium Ion", "Cs+") );

        setCompoundName("Cesium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Cesium Ion", "Cs+") )
            Biotype::defineBiotype(Element::Cesium(), 0, "Cesium Ion", "Cs+");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Cs+ Cesium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Cesium().getAtomicNumber(),
				0, // no covalent bonds
				3.390, // radius
				0.0000806 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Cesium Ion", "Cs+").getIndex() );
   }
};

/// Mg<sup>2+</sup> magnesium ion with +2 charge
class MagnesiumIon : public Molecule {
public:
    MagnesiumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("MG ");

        setBaseAtom( "Mg+2", Biotype::get("Magnesium Ion", "Mg+2") );

        setCompoundName("Magnesium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Magnesium Ion", "Mg+2") )
            Biotype::defineBiotype(Element::Magnesium(), 0, "Magnesium Ion", "Mg+2");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Mg+2 Magnesium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Magnesium().getAtomicNumber(),
				0, // no covalent bonds
				0.7926, // radius
				0.8947 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				2.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Magnesium Ion", "Mg+2").getIndex() );
   }
};

/// Ca<sup>2+</sup> calcium ion with +2 charge
class CalciumIon : public Molecule {
public:
    CalciumIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("CA ");

        setBaseAtom( "Ca+2", Biotype::get("Calcium Ion", "Ca+2") );

        setCompoundName("Calcium Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Calcium Ion", "Ca+2") )
            Biotype::defineBiotype(Element::Calcium(), 0, "Calcium Ion", "Ca+2");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Ca+2 Calcium Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Calcium().getAtomicNumber(),
				0, // no covalent bonds
				1.7131, // radius
				0.459789 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				2.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Calcium Ion", "Ca+2").getIndex() );
   }
};

/// Zn<sup>2+</sup> zinc ion with +2 charge
class ZincIon : public Molecule {
public:
    ZincIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("ZN ");

        setBaseAtom( "Zn+2", Biotype::get("Zinc Ion", "Zn+2") );

        setCompoundName("Zinc Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Zinc Ion", "Zn+2") )
            Biotype::defineBiotype(Element::Zinc(), 0, "Zinc Ion", "Zn+2");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Zn+2 Zinc Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Zinc().getAtomicNumber(),
				0, // no covalent bonds
				1.1000, // radius
				0.0125 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				2.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Zinc Ion", "Zn+2").getIndex() );
   }
};

/// Cl<sup>-</sup> chloride ion with -1 charge
class ChlorideIon : public Molecule {
public:
    ChlorideIon() 
    {
        instantiateBiotypes();

        setPdbResidueName("CL ");

        setBaseAtom( "Cl-", Biotype::get("Chloride Ion", "Cl-") );

        setCompoundName("Chloride Ion");
    }

	static void instantiateBiotypes() {
        if (! Biotype::exists("Chloride Ion", "Cl-") )
            Biotype::defineBiotype(Element::Chlorine(), 0, "Chloride Ion", "Cl-");
	}

    // create charged atom types
    // ensure that charges sum to zero, unless molecule has a formal charge
    static void setAmberLikeParameters(DuMMForceFieldSubsystem& dumm)
    {
        instantiateBiotypes();

		const char* ionName = "Cl- Chloride Ion";

		if (! dumm.hasAtomClass(ionName) ) {
			dumm.defineAtomClass_KA(
				dumm.getNextUnusedAtomClassIndex(),
				ionName,
				Element::Chlorine().getAtomicNumber(),
				0, // no covalent bonds
				2.4700, // radius
				0.1000 // well depth
				);
		}

		if (! dumm.hasChargedAtomType(ionName) ) {
			dumm.defineChargedAtomType(
				dumm.getNextUnusedChargedAtomTypeIndex(),
				ionName, 
				dumm.getAtomClassIndex(ionName),
				-1.00
				);
		}

        dumm.setBiotypeChargedAtomType( dumm.getChargedAtomTypeIndex(ionName), Biotype::get("Chloride Ion", "Cl-").getIndex() );
   }
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_IONS_H_
