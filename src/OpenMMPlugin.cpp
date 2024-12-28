/* -------------------------------------------------------------------------- *
 *                     SimTK Molmodel(tm): OpenMM Plugin                      *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

#include "molmodel/internal/OpenMMPlugin.h"
#include "DuMMForceFieldSubsystemRep.h"

#if OPENMM_PLATFORM_CPU
    #include "CpuPlatform.h"
#elif OPENMM_PLATFORM_CUDA
    #include "CudaPlatform.h"
#elif OPENMM_PLATFORM_OPENCL
        #include "OpenCLPlatform.h"
#endif



// -------------------------------------------- 
// Helper functions
// --------------------------------------------

/*!
 * <!-- Minimum number of bits needed to represent a number in binary -->
*/
int requiredBits(int number) {
    if (number == 0) {
        return 1;
    }
    return std::floor(std::log2(number)) + 1;
}

/*!
 * <!-- Integer to binary string -->
*/
std::string toBinary(int number) {
    int numBits = requiredBits(number);
    std::string binary = "";
    for (int i = numBits - 1; i >= 0; --i) {
        int bit = (number >> i) & 1; 
        binary += std::to_string(bit);
    }
    return binary;
}

/*!
 * <!-- Print nonbonded particle information added to OpenMM -->
*/
void stdcout_OpenmmNonbParticle(
    int nax, double sqrtCoulombScale, double charge,
    double sigma, double vdwGlobalScaleFactor, double wellDepth)
{
    std::cout << "OpenMM added particle "
        << nax <<" " << sqrtCoulombScale <<" "
        << charge <<" " << sigma <<" " 
        << vdwGlobalScaleFactor <<" " << wellDepth <<" "
        << std::endl;
}

/*!
 * <!-- Print bond added to OpenMM -->
*/
void stdcout_OpenmmBond(int a1num, int a2num, double bondStretch_d0, double bondStretch_k_t2)
{
    std::cout << "OMMPlugin addBond: " << a1num <<" " << a2num <<" "
        << bondStretch_d0 <<" " << bondStretch_k_t2
    << std::endl;
}

/*!
 * <!-- Print angle added to OpenMM -->
*/
void stdcout_OpenmmAngle(int a1num, int a2num, int a3num, double theta0, double K)
{
    std::cout << "OMMPlugin addAngle: " << a1num <<" " << a2num <<" "<< a3num
        <<" "<< theta0 <<" "<< K
        << std::endl;
}


/*!
 * <!-- Print torsion angle added to OpenMM -->
*/
void stdcout_OpenmmTorsion(int a1num, int a2num, int a3num, int a4num, double periodicity, double theta0, double amplitude)
{
    std::cout << "OMMPlugin addTorsion: " << a1num <<" " << a2num <<" "<< a3num <<" "<< a4num
        <<" "<< periodicity <<" "<< theta0 <<" "<< amplitude
        << std::endl;
}


// -------------------------------------------- 
// Main functions
// --------------------------------------------

/*! <!--
 * Bonded energy calculation only works for fully flexible setup !!
 * Initialize (/ allocate):
 *     - positions cache
 *     - OpenMM forces
 *     - OpenMM system
 * Set:
 *     - temperature and thermostat
 *     - nonbonded method and cutoff distance
 * Check: OpenMM requirements
 * Add to OpenMM:
 *     - atoms
 *     - bonds
 * Set GBSA parameters
 * Add bonds, angles and torsions
 * -->
*/
std::string OpenMMPluginInterface::initializeOpenMM(bool allowReferencePlatform, const SimTK::DuMMForceFieldSubsystemRep* inDumm)
{

    // Pass to internal DuMM pointer
    dumm = inDumm;

    // Allocate positions cache
    // TODO OMM These are memory consuming because repeatedly calling delete becomes slow due to memory fragmentation
    NonbondAtomsPositionsCache = std::vector<OpenMM::Vec3>(dumm->getNumNonbondAtoms());
    PositionsCache = std::vector<OpenMM::Vec3>(dumm->getNumAtoms());

    // Allocate OpenMM forces
    auto ommNonbondedForce = std::make_unique<OpenMM::NonbondedForce>();
    auto ommGBSAOBCForce = std::make_unique<OpenMM::GBSAOBCForce>();
    auto ommHarmonicBondStretch = std::make_unique<OpenMM::HarmonicBondForce>();
    auto ommHarmonicAngleForce = std::make_unique<OpenMM::HarmonicAngleForce>();
    auto ommPeriodicTorsionForce = std::make_unique<OpenMM::PeriodicTorsionForce>();

    // Instantiate the thermostat with adjusted temperature
    Real temperature = 300.0;
    if(dumm->wantOpenMMIntegration){temperature = dumm->temperature;}
    //openMMThermostat = std::make_unique<OpenMM::AndersenThermostat>(temperature, 1);
    openMMThermostat = new OpenMM::AndersenThermostat(temperature, 1);
    openMMThermostat->setRandomNumberSeed(seed);

    // Determine whether OpenMM supports all the features we've asked for.
    // OpenMM does not support 1-2, 1-3, or 1-5 scaling.
    if (   dumm->vdwScale12!=0 || dumm->coulombScale12!=0 
        || dumm->vdwScale13!=0 || dumm->coulombScale13!=0
        || dumm->vdwScale15!=1 || dumm->coulombScale15!=1)
    {
        std::cout << "WARNING: Can't use OpenMM: unsupported vdW or Coulomb scaling required." << std::endl;
        return "";
    }

    // Currently OpenMM supports only the Lorentz-Berthelot L-J combining rule as used by AMBER and CHARMM.
    if (   dumm->vdwGlobalScaleFactor != 0 
        && dumm->vdwMixingRule != DuMMForceFieldSubsystem::LorentzBerthelot) 
    {
        std::cout << "WARNING: Can't use OpenMM: only the Lorentz-Berthelot (AMBER) Lennard Jones mixing rule is supported." << std::endl;
        return "";
    }

    // Allocate OpenMM system and add particles to it
    openMMSystem = std::make_unique<OpenMM::System>();
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax) {

        //openMMSystem->addParticle(masses[nax]);

        const Element& element = Element::getByAtomicNumber(dumm->getAtomElementNum(dumm->getAtomIndexOfNonbondAtom(nax)));
        openMMSystem->addParticle(element.getMass());
        //std::cout << "OMMBUG OpenMM System added particle with mass " <<" "<< element.getMass() << std::endl; std::cout<<std::flush;

    }

    // Nonbonded forces
    if (dumm->coulombGlobalScaleFactor!=0 || dumm->vdwGlobalScaleFactor!=0) {
        ommNonbondedForce->setNonbondedMethod( OpenMM::NonbondedForce::NonbondedMethod( dumm->nonbondedMethod ) );
        ommNonbondedForce->setCutoffDistance( dumm->nonbondedCutoff );
        // nonbondedForce->setUseSwitchingFunction( 0 );

        // Scale charges by sqrt of scale factor so that products of charges 
        // scale linearly.
        const Real sqrtCoulombScale = std::sqrt(dumm->coulombGlobalScaleFactor);

        // Here we'll define all the OpenMM particles, one per DuMM nonbond
        // atom. We'll also build up the list of all 1-2 bonds between nonbond
        // atoms which will be used by OpenMM as an exceptions list, with
        // nonbond interactions excluded for 1-2 and 1-3 connections, and
        // scaled down for 1-4 connections. Since OpenMM doesn't know about 
        // bodies, we can't used the stripped-down cross-body bond lists here.
        // We will look at all the 1-2 bonds for each nonbond atom and keep 
        // those that connect to another nonbond atom.
        //
        // When it is called for the i'th time, it specifies the parameters for the i'th particle.
        std::vector< std::pair<int, int> > ommBonds;
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); 
                                                                        ++nax)
        {   const DuMMAtom&        dummAtom = dumm->getAtom(dumm->getAtomIndexOfNonbondAtom(nax));
            const ChargedAtomType& atype  = dumm->chargedAtomTypes[dummAtom.chargedAtomTypeIndex];
            const AtomClass&       aclass = dumm->atomClasses[atype.atomClassIx];
            const Real             charge = atype.partialCharge;
            const Real             sigma  = 2*aclass.vdwRadius*DuMM::Radius2Sigma;
            const Real             wellDepth = aclass.vdwWellDepth;

            const DuMM::IncludedAtomIndex& iax = dummAtom.getIncludedAtomIndex();
            const DuMM::AtomIndex& dax = dummAtom.atomIndex;

            // Define particle; particle number will be the same as our
            // nonbond index number.
            ommNonbondedForce->addParticle(sqrtCoulombScale*charge, sigma, 
                                        dumm->vdwGlobalScaleFactor*wellDepth);

            // stdcout_OpenmmNonbParticle(nax, sqrtCoulombScale, charge, sigma, dumm->vdwGlobalScaleFactor, wellDepth);
            #ifdef __DRILLING__
                stdcout_OpenmmNonbParticle(nax, sqrtCoulombScale, charge, sigma, dumm->vdwGlobalScaleFactor, wellDepth);
            #endif

            // Collect 1-2 bonds to other nonbond atoms. Note that we 
            // don't care about bodies here -- every atom is considered
            // independent.
            for (unsigned short bix=0; bix < dummAtom.bond12.size(); ++bix) {
                const DuMMAtom& bondedAtom = dumm->getAtom(dummAtom.bond12[bix]);
                if (!bondedAtom.nonbondAtomIndex.isValid()){
                    continue;
                }
                ommBonds.emplace_back(std::make_pair(nax, bondedAtom.nonbondAtomIndex));
            }
        }

        // Register all the 1-2 bonds between nonbond atoms for scaling.
        ommNonbondedForce->createExceptionsFromBonds(ommBonds, dumm->coulombScale14, dumm->vdwScale14);

        // System takes over heap ownership of the force.
        openMMSystem->addForce(ommNonbondedForce.get());
        ommNonbondedForce.release();
    }

    // GBSA
    // When it is called for the i'th time, it specifies the parameters for the i'th particle.
    if (dumm->gbsaGlobalScaleFactor != 0) {
        ommGBSAOBCForce->setSolventDielectric(dumm->gbsaSolventDielectric);
        ommGBSAOBCForce->setSoluteDielectric(dumm->gbsaSoluteDielectric);

        // Watch the units here. OpenMM works exclusively in MD (nm, kJ/mol). 
        // CPU GBSA uses Angstrom, kCal/mol.
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax) 
        {
            ommGBSAOBCForce->addParticle(dumm->gbsaAtomicPartialCharges[nax],
                                      dumm->gbsaRadii[nax]*OpenMM::NmPerAngstrom,
                                      dumm->gbsaObcScaleFactors[nax]); 
        }

        // System takes over heap ownership of the force.
        openMMSystem->addForce(ommGBSAOBCForce.get());
        ommGBSAOBCForce.release();
        
        std::cout << "OpenMMPlugin added GBSA scaled at" << dumm->gbsaGlobalScaleFactor << std::endl;
    }

    // Add bonded forces
    // TODO: As it is now, it should work only with a fully flexible setup (nonbonded index)
    if( ! dumm->wantOpenMMCalcOnlyNonBonded ){

        for (DuMMIncludedBodyIndex incBodyIx(0);
             incBodyIx < dumm->getNumIncludedBodies(); ++incBodyIx) {

            const IncludedBody &inclBody = dumm->includedBodies[incBodyIx];
            assert(inclBody.isValid());

            for (DuMMBondStarterIndex bsx = inclBody.beginBondStarterAtoms;
                 bsx != inclBody.endBondStarterAtoms; ++bsx) {

                const DuMM::IncludedAtomIndex a1num = dumm->bondStarterAtoms[bsx];
                const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                // ADD BONDS (1-2)
                if ((dumm->bondStretchGlobalScaleFactor != 0) ||
                    (dumm->customBondStretchGlobalScaleFactor != 0)) {

                    for (DuMM::IncludedAtomIndex b12(0); b12 < a1.force12.size(); ++b12) {

                        const DuMM::IncludedAtomIndex a2num = a1.force12[b12];
                        const BondStretch &bondStretch = *a1.stretch[b12];

                        if (bondStretch.hasBuiltinTerm()) {

                            // TODO check units: Ang and KcalPerAngstrom2
                            ommHarmonicBondStretch->addBond(a1num, a2num,
                                                 bondStretch.d0, // * OpenMM::NmPerAngstrom,
                                                 bondStretch.k * 2.0); // * OpenMM::KJPerKcal OR * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                            // stdcout_OpenmmBond(a1num, a2num, bondStretch.d0, bondStretch.k * 2.0); std::cout<<std::flush;
                            
                            #ifdef __DRILLING__
                                stdcout_OpenmmBond(a1num, a2num, bondStretch.d0, bondStretch.k * 2.0);
                            #endif

                        }
                    }
                }

                // ADD ANGLES (1-2-3)
                if (dumm->bondBendGlobalScaleFactor != 0
                    || dumm->customBondBendGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                    for (int b13=0; b13 < (int)a1.force13.size(); ++b13) {
                        const DuMM::IncludedAtomIndex a2num = a1.force13[b13][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force13[b13][1];

                        const BondBend& bb = *a1.bend[b13];

                        if (bb.hasBuiltinTerm()) {

                            // TODO: check atom order !!! which should be the central atom?
                            // TODO: check units: degreess and kcal/rad2 ??
                            ommHarmonicAngleForce->addAngle(a1num, a2num, a3num,
                                               bb.theta0,
                                               bb.k * 2 ); // * OpenMM::KJPerKcal);

                            // stdcout_OpenmmAngle(a1num, a2num, a3num, bb.theta0, bb.k); std::cout<<std::flush;

                            #ifdef __DRILLING__
                                PrintOpenMMAngle(a1num, a2num, a3num, bb.theta0, bb.k);
                            #endif
                        }
                    }
                }

                // ADD DIHEDRALS (1-2-3-4)
                if (dumm->bondTorsionGlobalScaleFactor != 0
                    || dumm->customBondTorsionGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                    for (int b14=0; b14 < (int)a1.force14.size(); ++b14) {
                        const DuMM::IncludedAtomIndex a2num = a1.force14[b14][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force14[b14][1];
                        const DuMM::IncludedAtomIndex a4num = a1.force14[b14][2];

                        const BondTorsion& bt = *a1.torsion[b14];

                        if (bt.hasBuiltinTerm()) {
                            for ( int i=0; i < (int) bt.terms.size(); ++i)
                            {
                                ommPeriodicTorsionForce->addTorsion(a1num, a2num, a3num, a4num,
                                                        bt.terms[i].periodicity,
                                                        bt.terms[i].theta0,
                                                        bt.terms[i].amplitude);
                                #ifdef __DRILLING__
                                    stdcout_OpenmmTorsion(a1num, a2num, a3num, a4num, bt.terms[i].periodicity, bt.terms[i].theta0, bt.terms[i].amplitude);
                                #endif                                                        
                            }
                        }
                    }
                }


                // ADD IMPROPERS (1-2-3-4)
                if (dumm->amberImproperTorsionGlobalScaleFactor != 0) {
                    const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                    // TODO: a1num is actually the 3rd one... check openmm order
                    // for (int b14 = 0; b14 < (int)a1.forceImproper14.size(); ++b14)
                    if ((int)a1.forceImproper14.size()>0)
                    {
                        for (int b14 = 0; b14 < 1; ++b14)
                        {
                            //printf("(size: %d)entered a1.forceImproper14: ", (int)a1.forceImproper14.size());
                            const DuMM::IncludedAtomIndex a2num = a1.forceImproper14[b14][0];
                            const DuMM::IncludedAtomIndex a3num = a1.forceImproper14[b14][1];
                            const DuMM::IncludedAtomIndex a4num = a1.forceImproper14[b14][2];

                            const BondTorsion& bt = *a1.aImproperTorsion[b14];

                            if (bt.hasBuiltinTerm()) {
                                for ( int i=0; i < (int) bt.terms.size(); ++i)
                                {
                                    //printf("%d %d %d %d\n", a2num, a3num, a1num, a4num);
                                    ommPeriodicTorsionForce->addTorsion(a2num, a3num, a1num, a4num,
                                                            bt.terms[i].periodicity,
                                                            bt.terms[i].theta0,
                                                            bt.terms[i].amplitude);
                                    #ifdef __DRILLING__
                                        stdcout_OpenmmTorsion(a1num, a2num, a3num, a4num, bt.terms[i].periodicity, bt.terms[i].theta0, bt.terms[i].amplitude);
                                    #endif                                                             
                                }
                            }
                        }
                    }
                }
            }
        }

        openMMSystem->addForce(ommHarmonicBondStretch.get());
        ommHarmonicBondStretch.release();

        openMMSystem->addForce(ommHarmonicAngleForce.get());
        ommHarmonicAngleForce.release();

        openMMSystem->addForce(ommPeriodicTorsionForce.get());
        ommPeriodicTorsionForce.release();
    }

    // Get the thermostat
    //openMMSystem->addForce(openMMThermostat.get());
    // //openMMThermostat.release();
    openMMSystem->addForce(openMMThermostat);

    //std::cout << "OpenMM System number of forces " <<  openMMSystem->getNumForces() << std::endl;

    // Get the integrator
    Real stepsize = 0;
    if (dumm->wantOpenMMIntegration){stepsize = dumm->stepsize;}
    //std::cout << "OMMBUG OpenMM stepsize " << stepsize << std::endl; std::cout<<std::flush;
    openMMIntegrator = std::make_unique<OpenMM::VerletIntegrator>(stepsize); // TODO should release?
    
    // Get the platform
    // By default, OpenMM builds a .so for each platform (CPU, OpenCL and CUDA)
    // When loading that .so, two functions get called
    // 1. registerPlatform() which does what you see below
    // 2. registerKernelFactories() which is used for Drude, Pme, Rpmd and other plugins (which we do not need as of right now)
#if OPENMM_PLATFORM_CPU
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::CpuPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "CPU";

#elif OPENMM_PLATFORM_CUDA
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::CudaPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "CUDA";

#elif OPENMM_PLATFORM_OPENCL
    // if(OpenMM::Platform::getNumPlatforms() == 1)
    {
        platform = std::make_unique<OpenMM::OpenCLPlatform>();
        OpenMM::Platform::registerPlatform(platform.get());
        platform.release();
    }
    constexpr auto PLATFORM_NAME = "OpenCL";
    
#endif

    // Create the context
    try {
        auto& platform = OpenMM::Platform::getPlatformByName(PLATFORM_NAME);
        openMMContext = std::make_unique<OpenMM::Context>(*openMMSystem, *openMMIntegrator, platform);
        const double speed = openMMContext->getPlatform().getSpeed();

        if (speed <= 1 && !allowReferencePlatform) {
            std::cout << "ERROR: OpenMM not used: best available platform was " << PLATFORM_NAME << " with relative speed " << speed << std::endl;
            std::cout << "ERROR: Call setAllowOpenMMReference() if you want to use this anyway." << std::endl;
            return "";
        }

        std::cout << "NOTE: Created OpenMM context with " << PLATFORM_NAME << " platform with relative speed " << speed << std::endl;


    } catch (const std::exception& e) {
        // Could not create this platform so log and try the next one
        std::cout << "ERROR: OpenMM error during initialization: " << e.what() << std::endl;
        return "";
    }

    std::cout << "OpenMMPluginInterface loaded "
        << openMMContext->getPlatform().getName()
        << std::endl;

    // // Prepare positions cache
    // openMMState = openMMContext->getState(OpenMM::State::Positions);
    // const auto numAtoms = openMMState.getPositions().size();

    // atomLocationsCache.reserve(numAtoms);
    // for (const auto& atom : openMMState.getPositions()) {
    //     atomLocationsCache.push_back(SimTK::Vec3(atom[0], atom[1], atom[2]));
    // }


    // ----------------------------------------------
    // PBC - Periodic Boundary Conditions __begin__
    // ----------------------------------------------
    #ifdef __PBC__ // _pbc_

        OpenMM::Vec3 pbcVector_X(1, 0, 0);
        OpenMM::Vec3 pbcVector_Y(0, 1, 0);
        OpenMM::Vec3 pbcVector_Z(0, 0, 1);

        openMMSystem->setDefaultPeriodicBoxVectors(pbcVector_X, pbcVector_Y, pbcVector_Z);
        //ommNonbondedForce->setPMEParameters(double alpha, int nx, int ny, int nz); // _pbc_

        std::cout << "Set OpenMM System PBC " << std::endl;
    
    # endif
    // ----------------------------------------------
    // PBC - Periodic Boundary Conditions __end__
    // ----------------------------------------------


    return openMMContext->getPlatform().getName();
}


/*! <!-- Calculates forces and/or energy and *adds* them into the output
 * parameters -->
*/
void OpenMMPluginInterface::calcOpenMMEnergyAndForces
   (const Vector_<Vec3>&    includedAtomStation_G,
    const Vector_<Vec3>&    includedAtomPos_G,
    bool                    wantForces,
    bool                    wantEnergy,
    Vector_<SpatialVec>&    includedBodyForces_G,
    Real&                   energy) const
{

    if (!(wantForces || wantEnergy))
        return;

    if (!openMMContext) 
        throw std::runtime_error("ERROR: calcOpenMMNonbondedAndGBSAForces(): OpenMM has not been initialized.");

    
    setOpenMMPositions(includedAtomPos_G);

    int openMMStateDataTypes_Drill = 0;

    // #ifdef __DRILLING__
    //     int openMMStateDataTypes = openMMState.getDataTypes();
    //     //std::string openMMStateDataTypes_Str = toBinary(openMMStateDataTypes);
    //     //std::cout << "[OPENMM_DATA_TYPES]: in binary" <<" " << openMMStateDataTypes <<" " << openMMStateDataTypes_Str << std::endl;
    //     openMMStateDataTypes_Drill = ((wantEnergy?OpenMM::State::Forces_drl_bon:0)
    //                                 | (wantEnergy?OpenMM::State::Forces_drl_ang:0)
    //                                 | (wantEnergy?OpenMM::State::Forces_drl_tor:0)
    //                                 | (wantEnergy?OpenMM::State::Forces_drl_n14:0)
    //                                 | (wantEnergy?OpenMM::State::Forces_drl_vdw:0)
    //                                 | (wantEnergy?OpenMM::State::Forces_drl_cou:0));
    // #endif

    // Ask for energy, forces, or both.
    openMMState = openMMContext->getState(
        (wantForces?OpenMM::State::Forces:0) | (wantEnergy?OpenMM::State::Energy:0)
        | openMMStateDataTypes_Drill
    );

    // std::cout << "Energy: " << openMMState.getPotentialEnergy() << std::endl;

    if (wantForces) {
        const std::vector<OpenMM::Vec3>& openMMForces = openMMState.getForces();
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax)
        {
            const DuMM::IncludedAtomIndex iax = 
                dumm->getIncludedAtomIndexOfNonbondAtom(nax);          
            const IncludedAtom& includedAtom = dumm->getIncludedAtom(iax);
            const DuMMIncludedBodyIndex ibx = includedAtom.inclBodyIndex;
            const OpenMM::Vec3& ommForce = openMMForces[nax];
            const Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
            includedBodyForces_G[ibx] += SpatialVec(includedAtomStation_G[iax] % simForce, simForce);

            // Print
            // const DuMM::AtomIndex& dAIx = includedAtom.atomIndex;
            // const DuMMAtom& dummAtom = dumm->getAtom(dAIx);

        }

    }

    if (wantEnergy)
        {energy += openMMState.getPotentialEnergy();}

    TRACE_OPENMM(("OpenMM_Energy\t" +
        std::to_string(openMMState.getPotentialEnergy()) + 
        "\n").c_str());
}




/*!
 * <!-- Integrate trajectory using OpenMM -->
*/
void OpenMMPluginInterface::integrateTrajectory(int steps)
{
    // // Print coordinates after integration
    // std::cout << "After integration:" << std::endl;
    //stdcout_OpenmmPositions("OMMposs");
    // std::cout<<std::flush;

    // for (int i = 0; i < 16; i++) {
    //     openMMSystem->setParticleMass(i, 0);
    // }
    // openMMSystem->setParticleMass(0, 0);
    // openMMContext->reinitialize(true);
    // openMMContext->setPositions(getPositions());

    std::cout << std::setprecision(6) << std::fixed
        << "OMMDEBUG OpenMMPluginInterface::integrateTrajectory for"
        <<" "<<steps <<" x "
        <<" "<<openMMIntegrator->getStepSize()
        <<" at default thermostat T "<< openMMThermostat->getDefaultTemperature()
        << std::endl << std::endl;

    openMMIntegrator->step(steps);

    // // Print coordinates after integration
    // std::cout << "After integration:" << std::endl;
    //stdcout_OpenmmPositions("OMMposs");
    // std::cout<<std::flush;
}



// -------------------------------------------- 
// Interface
// --------------------------------------------

/*!
 * <!-- Set NonbondAtomsPositionsCache and OpenMM positions -->
*/                   
void OpenMMPluginInterface::setOpenMMPositions(
    const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G ) const
{
    assert(NonbondAtomsPositionsCache.size() == dumm->getNumNonbondAtoms());
    assert(includedAtomPos_G.size() == dumm->getNumIncludedAtoms());

    // Positions arrive in an array of all included atoms. Compress that down
    // to just nonbond atoms and convert to OpenMM Vec3 type.
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax)
    {
        const auto& pos_G = includedAtomPos_G[dumm->getIncludedAtomIndexOfNonbondAtom(nax)];
        NonbondAtomsPositionsCache[nax] = OpenMM::Vec3(pos_G[0], pos_G[1], pos_G[2]);
    }

    // Pass the converted positions to OpenMM
    openMMContext->setPositions(NonbondAtomsPositionsCache);
}

/*!
 * <!-- Set PositionsCache and OpenMM positions -->
*/
void OpenMMPluginInterface::setOpenMMPositions(
    const std::vector<SimTK::Vec3>& positions)
{
    // Make sure the memory cache is good
    assert(PositionsCache.size() == positions.size());

    // Copying means converting from one vector type to another
    for (std::size_t i = 0; i < positions.size(); i++) {
        const auto& p = positions[i];
        PositionsCache[i] = OpenMM::Vec3(p[0], p[1], p[2]);
    }

    // Pass the converted positions to OpenMM
    openMMContext->setPositions(PositionsCache);
}



/*!
 * <!--  -->
*/
const std::vector<OpenMM::Vec3>& OpenMMPluginInterface::getPositions() const
{
    openMMState = openMMContext->getState(OpenMM::State::Positions);
    return openMMState.getPositions();
}


/*!
 * <!--  -->
*/
void OpenMMPluginInterface::updateAtomLocationsCache()
{
    // openMMState = openMMContext->getState(OpenMM::State::Positions);
    // return openMMState.getPositions();

    openMMState = openMMContext->getState(OpenMM::State::Positions);
    const auto numAtoms = openMMState.getPositions().size();

    for (size_t i = 0; i < numAtoms; i++) {
        const auto& atom = openMMState.getPositions()[i];
        atomLocationsCache[i] = SimTK::Vec3(atom[0], atom[1], atom[2]);
    }
}



/*!
 * <!--  -->
*/
SimTK::Vec3 OpenMMPluginInterface::getAtomPosition( int dummAtomIndex ) const
{
    SimTK::DuMM::AtomIndex dummAtomIndex_ai(dummAtomIndex);
    SimTK::DuMM::NonbondAtomIndex nonbondedAtomIndex = dumm->getAtom(dummAtomIndex_ai).getNonbondAtomIndex();

    openMMState = openMMContext->getState( OpenMM::State::Positions );
    OpenMM::Vec3 position = openMMState.getPositions() [nonbondedAtomIndex]; 

    return SimTK::Vec3 ( position[0], position[1], position[2] );

    //return atomLocationsCache[nonbondedAtomIndex];

    // openMMState = openMMContext->getState( OpenMM::State::Positions );
    // OpenMM::Vec3 position = openMMState.getPositions() [nonbondedAtomIndex]; 

    // return SimTK::Vec3 ( position[0], position[1], position[2] );

}


/*!
 * <!--  -->
*/
Real OpenMMPluginInterface::calcPotentialEnergy() const
{
    openMMState = openMMContext->getState(OpenMM::State::Energy);
    return openMMState.getPotentialEnergy();
}


/*!
 * <!--  -->
*/
Real OpenMMPluginInterface::calcKineticEnergy() const
{
    openMMState = openMMContext->getState(OpenMM::State::Energy);
    return openMMState.getKineticEnergy();
}

/*!
 * <!--  -->
*/
void OpenMMPluginInterface::setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed) {
    // TODO why check
    // std::cout << "setVelocitiesToTemperature " << temperature << " " << seed << std::endl;
    if (openMMContext){
        openMMThermostat->setDefaultTemperature(temperature);
        openMMContext->setVelocitiesToTemperature(temperature, seed);
    }
}

/*!
 * <!--  -->
*/
void OpenMMPluginInterface::setParticleMass(int index, SimTK::Real mass) {
    openMMSystem->setParticleMass(index, mass);
}

/*!
 * <!--  -->
*/
void OpenMMPluginInterface::setOpenMMMasses(const std::vector<SimTK::Real>& argMasses) {
    this->masses = argMasses;
}


void OpenMMPluginInterface::setSeed(uint32_t seed) {
    this->seed = seed;
}

/*!
 * <!--  -->
*/
void OpenMMPluginInterface::setTimestep(Real stepsize) {
    openMMIntegrator->setStepSize(stepsize);
}


// -------------------------------------------- 
// Debugging information
// --------------------------------------------

/*!
 * <!-- Print OpenMM positions -->
*/
void OpenMMPluginInterface::stdcout_OpenmmPositions(const std::string& header__ )
{
    openMMState = openMMContext->getState(OpenMM::State::Positions);
    const std::vector<OpenMM::Vec3>& omm_positions = openMMState.getPositions();

    for (int oaix = 0; oaix < omm_positions.size(); oaix++) {
        std::cout << header__
            <<" "<< omm_positions[oaix][0]
            <<" "<< omm_positions[oaix][1]
            <<" "<< omm_positions[oaix][2]
            // <<" "<< openMMSystem->getParticleMass(oaix)
            << std::endl;
    }
}




// void OpenMMPluginInterface::setNonbondedCutoff (SimTK::Real cutoff) {
//     assert(!"Not implemented!");
// }
// void OpenMMPluginInterface::setOpenMMPlatform (std::string platform) {
//     assert(!"Not implemented!");
// }
// void OpenMMPluginInterface::setGPUindex(std::string GPUindex) {
//     assert(!"Not implemented!");
// }
// SimTK::Real OpenMMPluginInterface::getNonbondedCutoff() const {
//     assert(!"Not implemented!");
//     return -1;
// }
// std::string OpenMMPluginInterface::getOpenMMPlatform() const {
//     assert(!"Not implemented!");
//     return "";
// }
// std::string OpenMMPluginInterface::getGPUindex() const {
//     assert(!"Not implemented!");
//     return "";
// }


