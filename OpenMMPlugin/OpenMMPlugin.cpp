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

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
    #pragma warning(disable:4251)   // no dll interface for some classes
    #define EXPORT __declspec(dllexport)
#else
    #define EXPORT
#endif

#include "SimTKcommon.h"
#include "../src/OpenMMPlugin.h"
#include "DuMMForceFieldSubsystemRep.h"

#include "OpenMM.h"

#include <string>
#include <vector>
#include <exception>
#include <cassert>
#include <algorithm>

using namespace SimTK;

#define STRINGIZE(var) #var
#define MAKE_VERSION_STRING(maj,min,build)  STRINGIZE(maj.min.build)

#define TRACE_OPENMM(STR) ;
//#define TRACE_OPENMM(STR) printf("%s", STR);

/**
 * This is a concrete implementation of the OpenMMPluginInterface class defined
 * by Molmodel.
 */

class OpenMMInterface : public OpenMMPluginInterface {
public:
    OpenMMInterface(const DuMMForceFieldSubsystemRep& dumm) 
    : dumm(dumm), openMMSystem(0), openMMContext(0), openMMIntegrator(0) {}

    ~OpenMMInterface() {deleteOpenMM();}

    // This reports the version of Molmodel which was current at the time
    // this plugin was compiled.
    std::string getMolmodelVersion() const {
        return MAKE_VERSION_STRING(SimTK_MOLMODEL_MAJOR_VERSION,
                                   SimTK_MOLMODEL_MINOR_VERSION,
                                   SimTK_MOLMODEL_PATCH_VERSION);
    }

    // Call this during Molmodel's realizeTopology() method. Return value
    // is the selected OpenMM Platform name.
    std::string initializeOpenMM(
        bool allowReferencePlatform, 
        std::vector<std::string>& logMessages ) throw();
        //const std::vector<SimTK::Real>& lambda_sterics,
        //const std::vector<SimTK::Real>& lambda_electrostatics) throw();


    // Calculates forces and/or energy and *adds* them into the output
    // parameters.
    void calcOpenMMEnergyAndForces
       (const Vector_<Vec3>&    includedAtomStation_G,
        const Vector_<Vec3>&    includedAtomPos_G,
        bool                    wantForces,
        bool                    wantEnergy,
        Vector_<SpatialVec>&    includedBodyForce_G,
        Real&                   energy) const;


    void setNonbondedCutoff (Real cutoff) ;            /// Set NonbondedCutoff for OpenMM
    void setOpenMMPlatform (std::string platform) ;    /// Set Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    void setGPUindex (std::string GPUindex) ;          /// Set GPU index (if Platform CUDA/OpenCL). Values:"0"/"1"/"0,1"

    Real getNonbondedCutoff () const;                  /// Get NonbondedCutoff for OpenMM
    std::string getOpenMMPlatform () const;            /// Get Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    std::string getGPUindex () const;                  /// Get GPU index. Values: "0"/"1"/"0,1"

    void updLambdaGlobal (Real lambda);
    


private:
    // Put this object back into its just-constructed condition.
    void deleteOpenMM() {
        delete openMMIntegrator; openMMIntegrator=0;
        delete openMMContext;    openMMContext=0;
        delete openMMSystem;     openMMSystem=0;
    }

private:
    const DuMMForceFieldSubsystemRep& dumm;

    OpenMM::System*             openMMSystem;
    OpenMM::Context*            openMMContext;
    OpenMM::Integrator*         openMMIntegrator; // dummy
    int                         cBondForceIx=-1;
    int                         cNonBondForceIx=-1;

};

//-----------------------------------------------------------------------------
//                    SimTK_createOpenMMPluginInterface
//-----------------------------------------------------------------------------
// This is the exported, non-name-mangled symbol that is called from the loading
// process to get access to this plugin's implementation of the OpenMM plugin
// interface defined by DuMM.

extern "C" EXPORT OpenMMPluginInterface*
SimTK_createOpenMMPluginInterface(const DuMMForceFieldSubsystemRep& dumm) {
    return new OpenMMInterface(dumm);
}

//-----------------------------------------------------------------------------
//                              initializeOpenMM
//-----------------------------------------------------------------------------
// sterge c_ dupa ce citesti din fisier
std::string OpenMMInterface::
initializeOpenMM(bool allowReferencePlatform, 
        std::vector<std::string>& logMessages) throw()
        //const std::vector<SimTK::Real>& c_lambda_sterics,
        //const std::vector<SimTK::Real>& c_lambda_electrostatics) throw()
{
    logMessages.clear();

    // Determine whether OpenMM supports all the features we've asked for.

    // OpenMM does not support 1-2, 1-3, or 1-5 scaling.
    if (   dumm.vdwScale12!=0 || dumm.coulombScale12!=0 
        || dumm.vdwScale13!=0 || dumm.coulombScale13!=0
        || dumm.vdwScale15!=1 || dumm.coulombScale15!=1)
    {
        logMessages.push_back(
            "WARNING: Can't use OpenMM: unsupported vdW or Coulomb scaling required.\n");
        return "";
    }

    // Currently OpenMM supports only the Lorentz-Berthelot L-J combining rule
    // as used by AMBER and CHARMM.
    if (   dumm.vdwGlobalScaleFactor != 0 
        && dumm.vdwMixingRule != DuMMForceFieldSubsystem::LorentzBerthelot) 
    {
        logMessages.push_back(
            "WARNING: Can't use OpenMM: only the Lorentz-Berthelot"
            " (AMBER) Lennard Jones mixing rule is supported.\n");
        return "";
    }

try {

        // OpenMM SYSTEM //

    openMMSystem = new OpenMM::System();
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); ++nax) {
        const Element& e = Element::getByAtomicNumber
           (dumm.getAtomElementNum(dumm.getAtomIndexOfNonbondAtom(nax)));
        openMMSystem->addParticle(e.getMass());
    }

        // CLASSIC NONBONDED FORCE //

    if (dumm.coulombGlobalScaleFactor!=0 || dumm.vdwGlobalScaleFactor!=0) {
        OpenMM::NonbondedForce* regularNonbondedForce = new OpenMM::NonbondedForce();

        // Scale charges by sqrt of scale factor so that products of charges 
        // scale linearly.
        const Real sqrtCoulombScale = std::sqrt(dumm.coulombGlobalScaleFactor);


        // Here we'll define all the OpenMM particles, one per DuMM nonbond
        // atom. We'll also build up the list of all 1-2 bonds between nonbond
        // atoms which will be used by OpenMM as an exceptions list, with
        // nonbond interactions excluded for 1-2 and 1-3 connections, and
        // scaled down for 1-4 connections. Since OpenMM doesn't know about 
        // bodies, we can't used the stripped-down cross-body bond lists here.
        // We will look at all the 1-2 bonds for each nonbond atom and keep 
        // those that connect to another nonbond atom.
        std::vector< std::pair<int, int> > ommBonds;
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); 
                                                                        ++nax) 
        {   const DuMMAtom&        a      = dumm.getAtom(dumm.getAtomIndexOfNonbondAtom(nax));
            const ChargedAtomType& atype  = dumm.chargedAtomTypes[a.chargedAtomTypeIndex];
            const AtomClass&       aclass = dumm.atomClasses[atype.atomClassIx];
            const Real             charge = atype.partialCharge;
            const Real             sigma  = 2*aclass.vdwRadius*DuMM::Radius2Sigma;
            const Real             wellDepth = aclass.vdwWellDepth;

            // Define particle; particle number will be the same as our
            // nonbond index number.
            regularNonbondedForce->addParticle(sqrtCoulombScale*charge, sigma, 
                                        dumm.vdwGlobalScaleFactor*wellDepth);

            // Collect 1-2 bonds to other nonbond atoms. Note that we 
            // don't care about bodies here -- every atom is considered
            // independent.
            for (unsigned short i=0; i < a.bond12.size(); ++i) {
                const DuMMAtom& b = dumm.getAtom(a.bond12[i]);
                if (!b.nonbondAtomIndex.isValid()) continue;
                ommBonds.push_back(std::pair<int,int>(nax,b.nonbondAtomIndex));
            }
        }

        // Register all the 1-2 bonds between nonbond atoms for scaling.
        regularNonbondedForce->createExceptionsFromBonds
                            (ommBonds, dumm.coulombScale14, dumm.vdwScale14);
        //logMessages.push_back("Regular NB force initialized!");

    
        // CUSTOM NONBONDED FORCES //

        std::string forceExpression = "U_sterics+U_electrostatics;";
        forceExpression += "U_electrostatics=((lambda_electrostatics)^1)*ONE_4PI_EPS0*chargeprod/reff_electrostatics;";
        forceExpression += "reff_electrostatics=sigma*((0*(1.0-(lambda_electrostatics))^1+(r/sigma)^2))^(1/2);";
        forceExpression += "U_sterics=((lambda_sterics)^1)*4*epsilon*x*(x-1.0);";
        forceExpression += "x=(sigma/reff_sterics)^6;reff_sterics=sigma*((0.5*(1.0-(lambda_sterics))^1+(r/sigma)^6))^(1/6);";
        forceExpression += "chargeprod=charge1*charge2;epsilon=sqrt(epsilon1*epsilon2);sigma=0.5*(sigma1+sigma2);";
        forceExpression += "lambda_sterics=max(lambda_sterics1,lambda_sterics2);";
        forceExpression += "lambda_electrostatics=max(lambda_electrostatics1,lambda_electrostatics2)";

        
        OpenMM::CustomNonbondedForce* customNonbondedForce = new OpenMM::CustomNonbondedForce(forceExpression);
         

        // Attempt to reproduce standard OpenMM NB energy, with code borrowed from
        // https://github.com/openmm/openmm/issues/3281
            
        // Define per particle parameters
        customNonbondedForce->addPerParticleParameter("charge");
        customNonbondedForce->addPerParticleParameter("epsilon");
        customNonbondedForce->addPerParticleParameter("sigma");
        customNonbondedForce->addPerParticleParameter("lambda_sterics");
        customNonbondedForce->addPerParticleParameter("lambda_electrostatics");
        customNonbondedForce->addGlobalParameter("ONE_4PI_EPS0", 138.935456);

        std::vector<double> NBparams(5);      

        // Scale charges by sqrt of scale factor so that products of charges 
        // scale linearly.

        std::vector<SimTK::Real> lambda_sterics, lambda_electrostatics; // sterge asta dupa ce citesti din fisier

        for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); ++nax) 

        {   const DuMMAtom&        a      = dumm.getAtom(dumm.getAtomIndexOfNonbondAtom(nax));
            const ChargedAtomType& atype  = dumm.chargedAtomTypes[a.chargedAtomTypeIndex];
            const AtomClass&       aclass = dumm.atomClasses[atype.atomClassIx];
            const Real             charge = atype.partialCharge;
            const Real             sigma  = 2*aclass.vdwRadius*DuMM::Radius2Sigma;
            Real                   lambda = dumm.lambdaGlobal;
            // const Real             lambda_sterics = 1;
            // const Real             lambda_electrostatics = 1;
            const Real             wellDepth = aclass.vdwWellDepth;

            // Define particle; particle number will be the same as our
            // nonbond index number.

            lambda_sterics.push_back(lambda);
            lambda_electrostatics.push_back(lambda);
            
            NBparams[0] = sqrtCoulombScale*charge;
            NBparams[1] = dumm.vdwGlobalScaleFactor*wellDepth;
            NBparams[2] = sigma;
            //NBparams[3] = lambda_sterics[nax];
            NBparams[3] = lambda;
            NBparams[4] = lambda;
            logMessages.push_back("NOTE: Lambda_sterics has a value of " + String(NBparams[3]) + " for particle number " + String(nax));
            //NBparams[4] = lambda_electrostatics[nax];
            logMessages.push_back("NOTE: Lambda_electrostatics has a value of " + String(NBparams[4]) + " for particle number " + String(nax));

            customNonbondedForce->addParticle(NBparams);

            // Collect 1-2 bonds to other nonbond atoms. Note that we 
            // don't care about bodies here -- every atom is considered
            // independent.
            for (unsigned short i=0; i < a.bond12.size(); ++i) {
                const DuMMAtom& b = dumm.getAtom(a.bond12[i]);
                if (!b.nonbondAtomIndex.isValid()) continue;
                ommBonds.push_back(std::pair<int,int>(nax,b.nonbondAtomIndex));
            }
        }            
        
        // Register all the 1-2 bonds between nonbond atoms for scaling.
        //customNonbondedForce->createExclusionsFromBonds
        //                    (ommBonds, 2);      

        // Using the "real" nonbonded force, we get a list of exceptions, which we turn
        // into exclusions. Then we create a customBondForce object which will include
        // all exceptions, scaled accordingly

        // We add a customBondForce to handle all the exceptions generated by NonbondedForce
        // It has the same expression as the CustomNonbondedForce, as we want it to simulate
        // that type of interaction

        std::string bondForceExpression = "U_sterics+U_electrostatics;";
        bondForceExpression += "U_electrostatics=((lambda_electrostatics)^1)*ONE_4PI_EPS0*chargeprod/reff_electrostatics;";
        bondForceExpression += "reff_electrostatics=sigma*((0*(1.0-(lambda_electrostatics))^1+(r/sigma)^2))^(1/2);";
        bondForceExpression += "U_sterics=((lambda_sterics)^1)*4*epsilon*x*(x-1.0);";
        bondForceExpression += "x=(sigma/reff_sterics)^6;reff_sterics=sigma*((0.5*(1.0-(lambda_sterics))^1+(r/sigma)^6))^(1/6);";
        // bondforceExpression += "lambda_sterics=max(lambda_sterics1,lambda_sterics2);";
        // bondforceExpression += "lambda_electrostatics=max(lambda_electrostatics1,lambda_electrostatics2)";

        OpenMM::CustomBondForce* customBondForce = new OpenMM::CustomBondForce("4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r;"
        "ONE_4PI_EPS0 = 138.935456");

        // Define per bond parameters
        customBondForce->addPerBondParameter("chargeprod");
        customBondForce->addPerBondParameter("epsilon");
        customBondForce->addPerBondParameter("sigma");
        customBondForce->addPerBondParameter("lambda_sterics");
        customBondForce->addPerBondParameter("lambda_electrostatics");

        std::vector<double> bondParams(5);

        //logMessages.push_back("Original NB Force has a number of " + String(regularNonbondedForce -> getNumExceptions()) + " nonbonded exceptions.");
        for (int i=0; i < regularNonbondedForce -> getNumExceptions(); ++i){
            int part1Exc=0;
            int part2Exc=0;
            double chargeProdExc=0.0;
            double sigmaExc=0.0;
            double epsilonExc=0.0;

            regularNonbondedForce -> getExceptionParameters(i, part1Exc, part2Exc,
                                            chargeProdExc, sigmaExc, epsilonExc);

            //logMessages.push_back("Exclusion added for particles " + String(part1Exc) + " and " + String(part2Exc) + ".");
            customNonbondedForce -> addExclusion(part1Exc, part2Exc);
            
            bondParams[0] = chargeProdExc;
            bondParams[1] = epsilonExc;
            bondParams[2] = sigmaExc;
            bondParams[3] = std::max(lambda_sterics[part1Exc], lambda_sterics[part2Exc]);
            bondParams[4] = std::max(lambda_electrostatics[part1Exc], lambda_electrostatics[part2Exc]);
            //bondParams[3] = std::max(lambda[part1Exc], lambda[part2Exc]);
            //bondParams[4] = std::max(lambda[part1Exc], lambda[part2Exc]);
            customBondForce -> addBond(part1Exc, part2Exc, bondParams);

            std::string out = "bond sterics/electrostatics: ";
            out += std::to_string(bondParams[3]);
            out += " ";
            out += std::to_string(bondParams[4]);

            //logMessages.push_back(out.c_str());
        }

        // System takes over heap ownership of the force.
        openMMSystem->addForce(customNonbondedForce);
        openMMSystem->addForce(customBondForce);

        cBondForceIx    = openMMSystem->getNumForces()-1;
        cNonBondForceIx = openMMSystem->getNumForces()-2;
        logMessages.push_back ("################## " +  String(cBondForceIx) + "##############");
        //std::cout << cNonBondForceIx << std::endl;


    }
    

        // GBSA //

    if (dumm.gbsaGlobalScaleFactor != 0) {
        OpenMM::GBSAOBCForce* GBSAOBCForce   = new OpenMM::GBSAOBCForce();
        GBSAOBCForce->setSolventDielectric(dumm.gbsaSolventDielectric);
        GBSAOBCForce->setSoluteDielectric(dumm.gbsaSoluteDielectric);

        // Watch the units here. OpenMM works exclusively in MD (nm, kJ/mol). 
        // CPU GBSA uses Angstrom, kCal/mol.
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); 
                                                                        ++nax) 
        {   GBSAOBCForce->addParticle(dumm.gbsaAtomicPartialCharges[nax],
                                      dumm.gbsaRadii[nax]*OpenMM::NmPerAngstrom,
                                      dumm.gbsaObcScaleFactors[nax]); 
        }

        // System takes over heap ownership of the force.
        openMMSystem->addForce(GBSAOBCForce);
    }




    // BONDED

    if( ! dumm.wantOpenMMCalcOnlyNonBonded ){

        // TODO !!!!!
        // Be sure that nonbonded index order is equivalent...and all bonded atoms were added as particles
        // As it is now, it should work only with a fully flexible setup

        OpenMM::HarmonicBondForce *bondStretch = new OpenMM::HarmonicBondForce();
        OpenMM::HarmonicAngleForce     *bondBend    = new OpenMM::HarmonicAngleForce();
        OpenMM::PeriodicTorsionForce   *bondTorsion = new OpenMM::PeriodicTorsionForce();


        for (DuMMIncludedBodyIndex incBodyIx(0);
             incBodyIx < dumm.getNumIncludedBodies(); ++incBodyIx) {

            const IncludedBody &inclBody = dumm.includedBodies[incBodyIx];
            assert(inclBody.isValid());

            for (DuMMBondStarterIndex bsx = inclBody.beginBondStarterAtoms;
                 bsx != inclBody.endBondStarterAtoms; ++bsx) {

                const DuMM::IncludedAtomIndex a1num = dumm.bondStarterAtoms[bsx];
                const IncludedAtom &a1 = dumm.getIncludedAtom(a1num);


                // ADD BONDED STRETCHES (1-2)
                if (dumm.bondStretchGlobalScaleFactor != 0
                    || dumm.customBondStretchGlobalScaleFactor != 0) {

                    for (DuMM::IncludedAtomIndex b12(0); b12 < a1.force12.size(); ++b12) {

                        const DuMM::IncludedAtomIndex a2num = a1.force12[b12];
                        const BondStretch &bs = *a1.stretch[b12];

                        if (bs.hasBuiltinTerm()) {

                            // TODO check units: Ang and KcalPerAngstrom2 ??
                            bondStretch->addBond(a1num, a2num,
                                                 bs.d0,
//                                                * OpenMM::NmPerAngstrom,
                                                 bs.k * 2.0);
//                                                * OpenMM::KJPerKcal
//                                                * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm);
                        }
                    }
                }

                // ADD BONDED BEND (1-2-3)
                if (dumm.bondBendGlobalScaleFactor != 0
                    || dumm.customBondBendGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm.getIncludedAtom(a1num);

                    for (int b13=0; b13 < (int)a1.force13.size(); ++b13) {
                        const DuMM::IncludedAtomIndex a2num = a1.force13[b13][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force13[b13][1];

                        const BondBend& bb = *a1.bend[b13];

                        if (bb.hasBuiltinTerm()) {

                            // TODO: check atom order !!! which should be the central atom?
                            // TODO: check units: degreess and kcal/rad2 ??
                            bondBend->addAngle(a1num, a2num, a3num,
                                               bb.theta0,
                                               bb.k * 2 );
                            // * OpenMM::KJPerKcal);
                        }
                    }
                }

                // ADD BONDED DIHEDRALS (1-2-3-4)
                if (dumm.bondTorsionGlobalScaleFactor != 0
                    || dumm.customBondTorsionGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm.getIncludedAtom(a1num);

                    for (int b14=0; b14 < (int)a1.force14.size(); ++b14) {
                        const DuMM::IncludedAtomIndex a2num = a1.force14[b14][0];
                        const DuMM::IncludedAtomIndex a3num = a1.force14[b14][1];
                        const DuMM::IncludedAtomIndex a4num = a1.force14[b14][2];

                        const BondTorsion& bt = *a1.torsion[b14];

                        if (bt.hasBuiltinTerm()) {
                            for ( int i=0; i < (int) bt.terms.size(); ++i)
                            {
                                bondTorsion->addTorsion(a1num, a2num, a3num, a4num,
                                                        bt.terms[i].periodicity,
                                                        bt.terms[i].theta0,
                                                        bt.terms[i].amplitude);
                            }
                        }
                    }
                }


                // ADD BONDED IMPROPERS (1-2-3-4)
                if (dumm.amberImproperTorsionGlobalScaleFactor != 0) {

                    const IncludedAtom &a1 = dumm.getIncludedAtom(a1num);

                    // TODO: a1num is actually the 3rd one... check openmm order
                    for (int b14=0; b14 < (int)a1.forceImproper14.size(); ++b14) {
                        const DuMM::IncludedAtomIndex a2num = a1.forceImproper14[b14][0];
                        const DuMM::IncludedAtomIndex a3num = a1.forceImproper14[b14][1];
                        const DuMM::IncludedAtomIndex a4num = a1.forceImproper14[b14][2];

                        const BondTorsion& bt = *a1.aImproperTorsion[b14];

                        if (bt.hasBuiltinTerm()) {
                            for ( int i=0; i < (int) bt.terms.size(); ++i)
                            {
                                bondTorsion->addTorsion(a1num, a2num, a3num, a4num,
                                                        bt.terms[i].periodicity,
                                                        bt.terms[i].theta0,
                                                        bt.terms[i].amplitude);
                            }
                        }
                    }
                }


            }
        }

        openMMSystem->addForce(bondStretch);
        openMMSystem->addForce(bondBend);
        openMMSystem->addForce(bondTorsion);

    }



        // OpenMM CONTEXT //

    const std::vector<std::string> pluginsLoaded = 
        OpenMM::Platform::loadPluginsFromDirectory(OpenMM::Platform::getDefaultPluginsDirectory());
    logMessages.push_back("NOTE: Loaded " + String(pluginsLoaded.size()) + " OpenMM plugins:");
    for (unsigned i=0; i < pluginsLoaded.size(); ++i)
        logMessages.back() += " " + pluginsLoaded[i];

    const int nPlatforms = OpenMM::Platform::getNumPlatforms();
    logMessages.push_back("NOTE: OpenMM has " + String(nPlatforms) + " Platforms registered: ");
    for (int i = 0; i < nPlatforms; ++i) {
        const OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
        logMessages.back() += " " + platform.getName();
    }

    // This is just a dummy to keep OpenMM happy; we're not using it for anything
    // so it doesn't matter what kind of integrator we pick.
    openMMIntegrator = new OpenMM::VerletIntegrator(0.1);
    openMMContext = new OpenMM::Context(*openMMSystem, *openMMIntegrator);
    const std::string pname = openMMContext->getPlatform().getName();
    const double speed = openMMContext->getPlatform().getSpeed();

    if (speed <= 1 && !allowReferencePlatform) {
        logMessages.push_back(
            "WARNING: DuMM: OpenMM not used: best available platform was "
                  + pname + " with relative speed=" + String(speed)
                  + ".\nCall setAllowOpenMMReference() if you want to use this anyway.\n");
        deleteOpenMM();
        return "";
    }

    return openMMContext->getPlatform().getName();
} 
catch (const std::exception& e) {
    logMessages.push_back(std::string("ERROR: OpenMM error during initialization: ") + e.what());
    deleteOpenMM();
    return "";
} 
catch (...) {
    logMessages.push_back("ERROR: Unknown exception during OpenMM initialization.");
    deleteOpenMM();
    return "";
}
}

//-----------------------------------------------------------------------------
//                    calcOpenMMNonbondedAndGBSAForces
//-----------------------------------------------------------------------------
void OpenMMInterface::calcOpenMMEnergyAndForces
   (const Vector_<Vec3>&    includedAtomStation_G,
    const Vector_<Vec3>&    includedAtomPos_G,
    bool                    wantForces,
    bool                    wantEnergy,
    Vector_<SpatialVec>&    includedBodyForces_G,
    Real&                   energy) const
{
    assert(includedAtomStation_G.size() == dumm.getNumIncludedAtoms());
    assert(includedAtomPos_G.size()     == dumm.getNumIncludedAtoms());
    assert(includedBodyForces_G.size()  == dumm.getNumIncludedBodies());

    if (!(wantForces || wantEnergy))
        return;

    if (!openMMContext) 
        throw std::runtime_error("ERROR: calcOpenMMNonbondedAndGBSAForces(): OpenMM has not been initialized.");

    // Positions arrive in an array of all included atoms. Compress that down
    // to just nonbond atoms and convert to OpenMM Vec3 type.
    std::vector<OpenMM::Vec3> positions(dumm.getNumNonbondAtoms());
    for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); ++nax)
    {   const Vec3& pos = 
            includedAtomPos_G[dumm.getIncludedAtomIndexOfNonbondAtom(nax)];
        positions[nax] = OpenMM::Vec3(pos[0], pos[1], pos[2]); }

    openMMContext->setPositions(positions);

    // Ask for energy, forces, or both.
    const OpenMM::State openMMState = 
        openMMContext->getState(  (wantForces?OpenMM::State::Forces:0)
                                | (wantEnergy?OpenMM::State::Energy:0));

    if (wantForces) {
        const std::vector<OpenMM::Vec3>& openMMForces = openMMState.getForces();
        for (DuMM::NonbondAtomIndex nax(0); nax < dumm.getNumNonbondAtoms(); 
                                                                        ++nax)
        {   const DuMM::IncludedAtomIndex iax = 
                dumm.getIncludedAtomIndexOfNonbondAtom(nax);          
            const IncludedAtom& a = dumm.getIncludedAtom(iax);
            const DuMMIncludedBodyIndex ibx = a.inclBodyIndex;
            const OpenMM::Vec3& ommForce = openMMForces[nax];
            const Vec3 f(ommForce[0], ommForce[1], ommForce[2]);
            includedBodyForces_G[ibx] += 
                SpatialVec(includedAtomStation_G[iax] % f, f);
        }
    }

    if (wantEnergy)
        energy += openMMState.getPotentialEnergy();
    TRACE_OPENMM(("OpenMM_Energy\t" + std::to_string(openMMState.getPotentialEnergy()) +  "\n").c_str());
}

void OpenMMInterface::updLambdaGlobal(Real lambda){
    ;
}


