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


std::string OpenMMPluginInterface::initializeOpenMM(bool allowReferencePlatform, const SimTK::DuMMForceFieldSubsystemRep* inDumm)
{

    // read atoms from dumm
    dumm = inDumm;

    // These are memory because repeatedly calling delete becomes slow due to memory fragmentation
    NonbondAtomsPositionsCache = std::vector<OpenMM::Vec3>(dumm->getNumNonbondAtoms());
    PositionsCache = std::vector<OpenMM::Vec3>(dumm->getNumAtoms());

    // instantiate forces
    auto ommNonbondedForce = std::make_unique<OpenMM::NonbondedForce>();
    auto ommGBSAOBCForce = std::make_unique<OpenMM::GBSAOBCForce>();
    auto ommHarmonicBondStretch = std::make_unique<OpenMM::HarmonicBondForce>();
    auto ommHarmonicAngleForce = std::make_unique<OpenMM::HarmonicAngleForce>();
    auto ommPeriodicTorsionForce = std::make_unique<OpenMM::PeriodicTorsionForce>();

    // instantiate the thermostat with adjusted temperature
    Real temperature = 300;
    if(dumm->wantOpenMMIntegration)
        temperature = dumm->temperature;
    auto openMMThermostat = std::make_unique<OpenMM::AndersenThermostat>(temperature, 1);
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

    // OpenMM system
    openMMSystem = std::make_unique<OpenMM::System>();
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax) {
        const Element& e = Element::getByAtomicNumber(dumm->getAtomElementNum(dumm->getAtomIndexOfNonbondAtom(nax)));
        openMMSystem->addParticle(e.getMass());
    }

    // nonbonded forces
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
        std::vector< std::pair<int, int> > ommBonds;
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); 
                                                                        ++nax) 
        {   const DuMMAtom&        dummAtom = dumm->getAtom(dumm->getAtomIndexOfNonbondAtom(nax));
            const ChargedAtomType& atype  = dumm->chargedAtomTypes[dummAtom.chargedAtomTypeIndex];
            const AtomClass&       aclass = dumm->atomClasses[atype.atomClassIx];
            const Real             charge = atype.partialCharge;
            const Real             sigma  = 2*aclass.vdwRadius*DuMM::Radius2Sigma;
            const Real             wellDepth = aclass.vdwWellDepth;

            // Define particle; particle number will be the same as our
            // nonbond index number.
            ommNonbondedForce->addParticle(sqrtCoulombScale*charge, sigma, 
                                        dumm->vdwGlobalScaleFactor*wellDepth);

            // std::cout << "SP_NEW_LABommp "
            //     << nax <<" " << sqrtCoulombScale <<" "
            //     << charge <<" " << sigma <<" " 
            //     << dumm->vdwGlobalScaleFactor <<" " << wellDepth <<" "
            //     << std::endl;

            // Collect 1-2 bonds to other nonbond atoms. Note that we 
            // don't care about bodies here -- every atom is considered
            // independent.
            for (unsigned short i=0; i < dummAtom.bond12.size(); ++i) {
                const DuMMAtom& b = dumm->getAtom(dummAtom.bond12[i]);
                if (!b.nonbondAtomIndex.isValid())
                    continue;
                ommBonds.emplace_back(std::make_pair(nax, b.nonbondAtomIndex));
            }
        }

        // Register all the 1-2 bonds between nonbond atoms for scaling.
        ommNonbondedForce->createExceptionsFromBonds(ommBonds, dumm->coulombScale14, dumm->vdwScale14);

        // System takes over heap ownership of the force.
        openMMSystem->addForce(ommNonbondedForce.get());
        ommNonbondedForce.release();
    }

    // GBSA
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
    }

    // Bonded
    if( ! dumm->wantOpenMMCalcOnlyNonBonded ){

        //std::cout << "OpenMMPlugin calculate bonded too " << std::endl;

        // TODO !!!!!
        // Be sure that nonbonded index order is equivalent...and all bonded atoms were added as particles
        // As it is now, it should work only with a fully flexible setup

        for (DuMMIncludedBodyIndex incBodyIx(0);
             incBodyIx < dumm->getNumIncludedBodies(); ++incBodyIx) {

            const IncludedBody &inclBody = dumm->includedBodies[incBodyIx];
            assert(inclBody.isValid());

            for (DuMMBondStarterIndex bsx = inclBody.beginBondStarterAtoms;
                 bsx != inclBody.endBondStarterAtoms; ++bsx) {

                const DuMM::IncludedAtomIndex a1num = dumm->bondStarterAtoms[bsx];
                const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);


                // ADD BONDED STRETCHES (1-2)
                if ((dumm->bondStretchGlobalScaleFactor != 0) ||
                    (dumm->customBondStretchGlobalScaleFactor != 0)) {

                    for (DuMM::IncludedAtomIndex b12(0); b12 < a1.force12.size(); ++b12) {

                        const DuMM::IncludedAtomIndex a2num = a1.force12[b12];
                        const BondStretch &bondStretch = *a1.stretch[b12];

                        if (bondStretch.hasBuiltinTerm()) {

                            // TODO check units: Ang and KcalPerAngstrom2 ??
                            ommHarmonicBondStretch->addBond(a1num, a2num,
                                                 bondStretch.d0,
//                                                * OpenMM::NmPerAngstrom,
                                                 bondStretch.k * 2.0);
//                                                * OpenMM::KJPerKcal
//                                                * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);

                            std::cout << "OMMPlug addBond: " << a1num <<" " << a2num <<" "
                               << bondStretch.d0 <<" " << bondStretch.k * 2.0 << std::endl;

                        }
                    }
                }

                // ADD BONDED BEND (1-2-3)
                if (dumm->bondBendGlobalScaleFactor != 0
                    || dumm->customBondBendGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                    for (int b13=0; b13 < (int)a1.force13.size(); ++b13) {
                        const DuMM::IncludedAtomIndex a2num = a1.force13[b13][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force13[b13][1];

                        printf("OMMPlug addAngle: a1num %d, a2num %d, a3num %d\n",
                                a1num, a2num, a3num);

                        const BondBend& bb = *a1.bend[b13];

                        if (bb.hasBuiltinTerm()) {

                            // TODO: check atom order !!! which should be the central atom?
                            // TODO: check units: degreess and kcal/rad2 ??
                            ommHarmonicAngleForce->addAngle(a1num, a2num, a3num,
                                               bb.theta0,
                                               bb.k * 2 );
                            // * OpenMM::KJPerKcal);
                        }
                    }
                }

                // ADD BONDED DIHEDRALS (1-2-3-4)
                if (dumm->bondTorsionGlobalScaleFactor != 0
                    || dumm->customBondTorsionGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm->getIncludedAtom(a1num);

                    for (int b14=0; b14 < (int)a1.force14.size(); ++b14) {
                        const DuMM::IncludedAtomIndex a2num = a1.force14[b14][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force14[b14][1];
                        const DuMM::IncludedAtomIndex a4num = a1.force14[b14][2];

                        printf("OMMPlug addTorsion: a1num %d, a2num %d, a3num %d, a4num %d\n",
                                a1num, a2num, a3num, a4num);

                        const BondTorsion& bt = *a1.torsion[b14];

                        if (bt.hasBuiltinTerm()) {
                            for ( int i=0; i < (int) bt.terms.size(); ++i)
                            {
                                ommPeriodicTorsionForce->addTorsion(a1num, a2num, a3num, a4num,
                                                        bt.terms[i].periodicity,
                                                        bt.terms[i].theta0,
                                                        bt.terms[i].amplitude);
                            }
                        }
                    }
                }


                // ADD BONDED IMPROPERS (1-2-3-4)
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

                            printf("OMMPlug addImproper: a1num %d, a2num %d, a3num %d, a4num %d\n",
                                    a1num, a2num, a3num, a4num);

                            const BondTorsion& bt = *a1.aImproperTorsion[b14];

                            if (bt.hasBuiltinTerm()) {
                                for ( int i=0; i < (int) bt.terms.size(); ++i)
                                {
                                    //printf("%d %d %d %d\n", a2num, a3num, a1num, a4num);
                                    ommPeriodicTorsionForce->addTorsion(a2num, a3num, a1num, a4num,
                                                            bt.terms[i].periodicity,
                                                            bt.terms[i].theta0,
                                                            bt.terms[i].amplitude);
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
    openMMSystem->addForce(openMMThermostat.get());
    openMMThermostat.release();

    //std::cout << "OpenMM System numbe rof forces " <<  openMMSystem->getNumForces() << std::endl;

    // Get the integrator
    Real stepsize = 0.001;
    if (dumm->wantOpenMMIntegration)
        stepsize = dumm->stepsize;
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

    return openMMContext->getPlatform().getName();
}




//-----------------------------------------------------------------------------
//                    updateCoordInOpenMM
//-----------------------------------------------------------------------------
void OpenMMPluginInterface::setOpenMMPositions(
    const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G ) const
{
    assert(NonbondAtomsPositionsCache.size() == dumm->getNumNonbondAtoms());
    assert(includedAtomPos_G.size() == dumm->getNumIncludedAtoms());

    // Positions arrive in an array of all included atoms. Compress that down
    // to just nonbond atoms and convert to OpenMM Vec3 type.
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm->getNumNonbondAtoms(); ++nax)
    {
        const auto& pos_G =
            includedAtomPos_G[dumm->getIncludedAtomIndexOfNonbondAtom(nax)];
        NonbondAtomsPositionsCache[nax] =
            OpenMM::Vec3(pos_G[0], pos_G[1], pos_G[2]);
    }

    // Pass the converted positions to OpenMM
    openMMContext->setPositions(NonbondAtomsPositionsCache);
}

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


//-----------------------------------------------------------------------------
//                    getPositions
//-----------------------------------------------------------------------------
const std::vector<OpenMM::Vec3>& OpenMMPluginInterface::getPositions() const
{
    openMMState = openMMContext->getState(OpenMM::State::Positions);
    return openMMState.getPositions();
}


//-----------------------------------------------------------------------------
//                    getAtomPosition
//-----------------------------------------------------------------------------
SimTK::Vec3 OpenMMPluginInterface::getAtomPosition( int dummAtomIndex ) const
{
    SimTK::DuMM::AtomIndex dummAtomIndex_ai(dummAtomIndex);
    SimTK::DuMM::NonbondAtomIndex nonbondedAtomIndex = dumm->getAtom(dummAtomIndex_ai).getNonbondAtomIndex();

    openMMState = openMMContext->getState( OpenMM::State::Positions );
    OpenMM::Vec3 position = openMMState.getPositions() [nonbondedAtomIndex]; 

    return SimTK::Vec3 ( position[0], position[1], position[2] );

}


//-----------------------------------------------------------------------------
//                    calcPotentialEnergy
//-----------------------------------------------------------------------------
Real OpenMMPluginInterface::calcPotentialEnergy() const
{
    openMMState = openMMContext->getState(OpenMM::State::Energy);
    return openMMState.getPotentialEnergy();
}

//-----------------------------------------------------------------------------
//                    calcKineticEnergy
//-----------------------------------------------------------------------------
Real OpenMMPluginInterface::calcKineticEnergy() const
{
    openMMState = openMMContext->getState(OpenMM::State::Energy);
    return openMMState.getKineticEnergy();
}

//-----------------------------------------------------------------------------
//                    integrateTrajectory
//-----------------------------------------------------------------------------
void OpenMMPluginInterface::integrateTrajectory(int steps)
{
    openMMIntegrator->step(steps);
}



void OpenMMPluginInterface::setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed) {
    // TODO why check
    if (openMMContext)
        openMMContext->setVelocitiesToTemperature(temperature, seed);
}

void OpenMMPluginInterface::setParticleMass(int index, SimTK::Real mass) {
    openMMSystem->setParticleMass(index, mass);
}


//-----------------------------------------------------------------------------
//                    calcOpenMMNonbondedAndGBSAForces
//-----------------------------------------------------------------------------
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

    // Ask for energy, forces, or both.
    openMMState = openMMContext->getState(
        (wantForces?OpenMM::State::Forces:0) | 
        (wantEnergy?OpenMM::State::Energy:0)
        // | (wantEnergy?OpenMM::State::Forces_drl_bon:0)
        // | (wantEnergy?OpenMM::State::Forces_drl_ang:0)
        // | (wantEnergy?OpenMM::State::Forces_drl_tor:0)
    );

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

        //drl BEGIN
        // const std::vector<OpenMM::Vec3>& drl_bon_Forces = openMMState.getForces_drl_bon();
        // printf("drl OpenMMPluginInterface::calcOpenMMEnergyAndForces\n");
        // for (int fIx = 0; fIx < dumm->getNumNonbondAtoms(); ++fIx){
        //     const OpenMM::Vec3& ommForce = drl_bon_Forces[fIx];
        //     const Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        //     printf("drl OMMPlug bon %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
        // }
        // const std::vector<OpenMM::Vec3>& drl_ang_Forces = openMMState.getForces_drl_ang();
        // printf("drl OpenMMPluginInterface::calcOpenMMEnergyAndForces\n");
        // for (int fIx = 0; fIx < dumm->getNumNonbondAtoms(); ++fIx){
        //     const OpenMM::Vec3& ommForce = drl_ang_Forces[fIx];
        //     const Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        //     printf("drl OMMPlug ang %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
        // }
        // const std::vector<OpenMM::Vec3>& drl_tor_Forces = openMMState.getForces_drl_tor();
        // printf("drl OpenMMPluginInterface::calcOpenMMEnergyAndForces\n");
        // for (int fIx = 0; fIx < dumm->getNumNonbondAtoms(); ++fIx){
        //     const OpenMM::Vec3& ommForce = drl_tor_Forces[fIx];
        //     const Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        //     printf("drl OMMPlug tor %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
        // }
        //drl END

    }

    if (wantEnergy)
        energy += openMMState.getPotentialEnergy();

    TRACE_OPENMM(("OpenMM_Energy\t" +
        std::to_string(openMMState.getPotentialEnergy()) + 
        "\n").c_str());
}

void OpenMMPluginInterface::setSeed(uint32_t seed) {
    this->seed = seed;
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


