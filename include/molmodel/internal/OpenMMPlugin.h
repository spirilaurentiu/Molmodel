#ifndef SimTK_MOLMODEL_OPENMM_PLUGIN_H_
#define SimTK_MOLMODEL_OPENMM_PLUGIN_H_

#include "SimTKcommon.h"
#include "OpenMM.h"

#include <string>
#include <vector>
#include <string>
#include <fstream>
#include <streambuf>
#include <exception>
#include <cassert>

#ifndef __DRILLING__
#define __DRILLING__
#endif

// Minimum number of bits needed to represent a number in binary
int requiredBits(int number) ;

// Integer to binary string
std::string toBinary(int number) ;

namespace SimTK {
    class DuMMForceFieldSubsystemRep;
}

/**
 * This is the interface which must be implemented by a DLL which wants
 * to serve as an OpenMM Plugin to Molmodel. The DLL should define a concrete
 * class deriving from this one which implements all the virtual methods.
 */
class OpenMMPluginInterface {
public:

    // Call this during Molmodel's realizeTopology() method. Return value
    // is the selected OpenMM Platform name.
    std::string initializeOpenMM(bool allowReferencePlatform, const SimTK::DuMMForceFieldSubsystemRep* inDumm);

    // Calculates forces and/or energy and *adds* them into the output
    // parameters.
    void calcOpenMMEnergyAndForces
       (const SimTK::Vector_<SimTK::Vec3>&    includedAtomStation_G,
        const SimTK::Vector_<SimTK::Vec3>&    includedAtomPos_G,
        bool                    wantForces,
        bool                    wantEnergy,
        SimTK::Vector_<SimTK::SpatialVec>&    includedBodyForce_G,
        SimTK::Real&                   energy) const;

    void setOpenMMPositions(const std::vector<SimTK::Vec3>& positions);
    const std::vector<OpenMM::Vec3>& getPositions() const;

    SimTK::Vec3 getAtomPosition( int dummAtomIndex ) const;
    SimTK::Real calcPotentialEnergy() const;
    SimTK::Real calcKineticEnergy() const;
    void integrateTrajectory(int steps);
    void setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed);
    void setParticleMass(int index, SimTK::Real mass);
    void setSeed(uint32_t seed);
    void setTimestep(SimTK::Real timestep);

    // void setNonbondedCutoff (SimTK::Real cutoff) ;            /// Set NonbondedCutoff for OpenMM
    // void setOpenMMPlatform (std::string platform) ;    /// Set Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    // void setGPUindex (std::string GPUindex) ;          /// Set GPU index (if Platform CUDA/OpenCL). Values:"0"/"1"/"0,1"

    // SimTK::Real getNonbondedCutoff () const;                  /// Get NonbondedCutoff for OpenMM
    // std::string getOpenMMPlatform () const;            /// Get Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    // std::string getGPUindex () const;                  /// Get GPU index. Values: "0"/"1"/"0,1"

        //drl BEGIN

        const std::vector<std::vector<double>>& getEnergies_drl_bon(){return openMMState.getEnergies_drl_bon();}
        const std::vector<std::vector<double>>& getEnergies_drl_ang(){return openMMState.getEnergies_drl_ang();}
        const std::vector<std::vector<double>>& getEnergies_drl_tor(){return openMMState.getEnergies_drl_tor();}
        const std::vector<std::vector<double>>& getEnergies_drl_n14(){return openMMState.getEnergies_drl_n14();}           
        const std::vector<std::vector<double>>& getEnergies_drl_vdw(){return openMMState.getEnergies_drl_vdw();}
        const std::vector<std::vector<double>>& getEnergies_drl_cou(){return openMMState.getEnergies_drl_cou();}           
        const std::vector<OpenMM::Vec3>& getForces_drl_bon(){return openMMState.getForces_drl_bon();}
        const std::vector<OpenMM::Vec3>& getForces_drl_ang(){return openMMState.getForces_drl_ang();}
        const std::vector<OpenMM::Vec3>& getForces_drl_tor(){return openMMState.getForces_drl_tor();}
        const std::vector<OpenMM::Vec3>& getForces_drl_n14(){return openMMState.getForces_drl_n14();}

        // const std::vector<std::vector<double>>& getEnergies_drl_bon(){}
        // const std::vector<std::vector<double>>& getEnergies_drl_ang(){}
        // const std::vector<std::vector<double>>& getEnergies_drl_tor(){}
        // const std::vector<std::vector<double>>& getEnergies_drl_n14(){}           
        // const std::vector<std::vector<double>>& getEnergies_drl_vdw(){}
        // const std::vector<std::vector<double>>& getEnergies_drl_cou(){}           
        // const std::vector<OpenMM::Vec3>& getForces_drl_bon(){}
        // const std::vector<OpenMM::Vec3>& getForces_drl_ang(){}
        // const std::vector<OpenMM::Vec3>& getForces_drl_tor(){}
        // const std::vector<OpenMM::Vec3>& getForces_drl_n14(){}  

        //drl END    

private:
    void setOpenMMPositions(const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G) const;

    // std::string OpenMMPlatform;
    // std::string OpenMMGPUindex;

    // mutable is ugly, but we are forced to use it as some functions must be declared const
    mutable std::vector<OpenMM::Vec3> NonbondAtomsPositionsCache;
    mutable OpenMM::State openMMState;
    std::vector<OpenMM::Vec3> PositionsCache;

    // These must be destroyed in reverse order (from context to system)
    std::unique_ptr<OpenMM::Platform> platform;
    std::unique_ptr<OpenMM::System> openMMSystem;
    std::unique_ptr<OpenMM::Integrator> openMMIntegrator;
    std::unique_ptr<OpenMM::Context> openMMContext;

    const SimTK::DuMMForceFieldSubsystemRep* dumm = nullptr;

    uint32_t seed = 0;
};

#endif //SimTK_MOLMODEL_OPENMM_PLUGIN_H_



