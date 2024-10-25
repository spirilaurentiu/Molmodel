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

// #ifndef __DRILLING__
// #define __DRILLING__
// #endif

#ifndef __PBC__ // _pbc_
#define __PBC__
#endif


// Minimum number of bits needed to represent a number in binary
int requiredBits(int number) ;

// Integer to binary string
std::string toBinary(int number) ;

namespace SimTK {
    class DuMMForceFieldSubsystemRep;
}

//==============================================================================
//                           CLASS OpenMMPluginInterface
//==============================================================================
/** 
 * This is the interface which must be implemented by a DLL which wants
 * to serve as an OpenMM Plugin to Molmodel. The DLL should define a concrete
 * class deriving from this one which implements all the virtual methods.
**/
class OpenMMPluginInterface {
public:

    /** @name Main functions.**/
    /**@{**/

	/**	
	* @brief Call this during Molmodel's realizeTopology() method. Return value
    * is the selected OpenMM Platform name.
	* @param allowReferencePlatform Allow OpenMM reference platform
	* @param inDumm DuMM force field
	* @return OpenMM context's platform's name
	*/    
    std::string initializeOpenMM(bool allowReferencePlatform,
        const SimTK::DuMMForceFieldSubsystemRep* inDumm);

	/**	
	* @brief Calculates forces and/or energy and *adds* them into the output
    * parameters.
	* @param includedAtomStation_G,includedAtomPos_G input
	* @param includedBodyForce_G, energy output
	*/    
    void calcOpenMMEnergyAndForces
       (const SimTK::Vector_<SimTK::Vec3>&    includedAtomStation_G,
        const SimTK::Vector_<SimTK::Vec3>&    includedAtomPos_G,
        bool                                  wantForces,
        bool                                  wantEnergy,
        SimTK::Vector_<SimTK::SpatialVec>&    includedBodyForce_G,
        SimTK::Real&                          energy) const;

	/**	
	* @brief Integrate trajectory using OpenMM
	* @param steps nof steps
	*/ 
    void integrateTrajectory(int steps);
    
	/**@}**/

    /** @name Interface.**/
    /**@{**/

	/**	
	* @brief Set PositionsCache and OpenMM positions
	* @param positions std::vector of SimTK::Vec3s
	*/ 
    void setOpenMMPositions(const std::vector<SimTK::Vec3>& positions);

	/**	
	* @brief Set NonbondAtomsPositionsCache and OpenMM positions
	* @param positions std::vector of SimTK::Vec3s
	*/ 
    void setOpenMMPositions(const SimTK::Vector_<SimTK::Vec3>& includedAtomPos_G) const;

    const std::vector<OpenMM::Vec3>& getPositions() const;
    
    void updateAtomLocationsCache();

    SimTK::Vec3 getAtomPosition( int dummAtomIndex ) const;

    SimTK::Real calcPotentialEnergy() const;
    SimTK::Real calcKineticEnergy() const;

    void setParticleMass(int index, SimTK::Real mass);
    void setOpenMMMasses(const std::vector<SimTK::Real>& masses);

    void setSeed(uint32_t seed);
    void setTimestep(SimTK::Real timestep);

    void setVelocitiesToTemperature(SimTK::Real temperature, uint32_t seed);

    // void setNonbondedCutoff (SimTK::Real cutoff) ;     /// Set NonbondedCutoff for OpenMM
    // void setOpenMMPlatform (std::string platform) ;    /// Set Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    // void setGPUindex (std::string GPUindex) ;          /// Set GPU index (if Platform CUDA/OpenCL). Values:"0"/"1"/"0,1"
    // SimTK::Real getNonbondedCutoff () const;           /// Get NonbondedCutoff for OpenMM
    // std::string getOpenMMPlatform () const;            /// Get Platform to use for OpenMM ('CPU', 'CUDA', 'OpenCL')
    // std::string getGPUindex () const;                  /// Get GPU index. Values: "0"/"1"/"0,1"

	/**@}**/

    /** @name Debugging information.**/
    /**@{**/

    /**	
	* @brief Print OpenMM positions
	* @param header__ word to be in front of every line
	*/ 
    void stdcout_OpenmmPositions(const std::string& header__ );

    /**@{**/


    /** @name Drilling drl.**/
    /**@{**/

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

	/**@}**/

private:


    std::vector<SimTK::Vec3> atomLocationsCache;

    // std::string OpenMMPlatform;
    // std::string OpenMMGPUindex;

    // mutable is ugly, but we are forced to use it as some functions must be declared const
    mutable std::vector<OpenMM::Vec3> NonbondAtomsPositionsCache;
    mutable OpenMM::State openMMState;
    std::vector<OpenMM::Vec3> PositionsCache;
    std::vector<SimTK::Real> masses;

    // These must be destroyed in reverse order (from context to system)
    std::unique_ptr<OpenMM::Platform> platform;
    std::unique_ptr<OpenMM::System> openMMSystem;
    std::unique_ptr<OpenMM::Integrator> openMMIntegrator;
    std::unique_ptr<OpenMM::Context> openMMContext;

    //std::unique_ptr<OpenMM::AndersenThermostat> openMMThermostat;
    OpenMM::AndersenThermostat* openMMThermostat;

    const SimTK::DuMMForceFieldSubsystemRep* dumm = nullptr;

    uint32_t seed = 0;
};

#endif //SimTK_MOLMODEL_OPENMM_PLUGIN_H_



