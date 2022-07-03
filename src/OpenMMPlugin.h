#ifndef SimTK_MOLMODEL_OPENMM_PLUGIN_H_
#define SimTK_MOLMODEL_OPENMM_PLUGIN_H_

#include "SimTKcommon.h"
#include <string>
#include <vector>
#include <string>
#include <fstream>
#include <streambuf>

/**
 * This is the interface which must be implemented by a DLL which wants
 * to serve as an OpenMM Plugin to Molmodel. The DLL should define a concrete
 * class deriving from this one which implements all the virtual methods.
 */
class OpenMMPluginInterface {
public:
    virtual ~OpenMMPluginInterface() {}

    // This is the Molmodel version number at the time the
    // OpenMMPlugin was built.
    virtual std::string getMolmodelVersion() const = 0;

    // Call this during Molmodel's realizeTopology() method. Return value
    // is the selected OpenMM Platform name, or the empty string if OpenMM
    // was not initialized. Log messages may be returned whether we succeed
    // or fail -- they are intended to provide a human some hint as to what
    // happened and why.
    // This method will never throw an exception because failure just means
    // you shouldn't use OpenMM.
    virtual std::string initializeOpenMM
       (bool allowReferencePlatform, 
        std::vector<std::string>& logMessages) throw() = 0;
        //const std::vector<SimTK::Real>& lambda_sterics,
        //const std::vector<SimTK::Real>& lambda_electrostatics) throw() = 0;

    // Calculates forces and/or energy and *adds* them into the output
    // parameters. This will throw an exception if something goes wrong.
    virtual void calcOpenMMEnergyAndForces
       (const SimTK::Vector_<SimTK::Vec3>&  atomStation_G,
        const SimTK::Vector_<SimTK::Vec3>&  atomPos_G,
        bool                                wantForces,
        bool                                wantEnergy,
        SimTK::Vector_<SimTK::SpatialVec>&  forces,
        SimTK::Real&                        energy) const = 0;

    //virtual void updLambdaGlobalIFC (SimTK::Real& lambda);
    virtual void updLambdaGlobalIFC
       (std::vector<SimTK::Real> lambdaPair) const = 0;
    //virtual void updLambdaGlobalIFC ();


private:
    std::string OpenMMPlatform;
    std::string OpenMMGPUindex;
};

namespace SimTK {
    class DuMMForceFieldSubsystemRep;
}

/**
 * This class defines the OpenMM Plugin library (dll, so, dylib) and 
 * encapsulates the policy we use to locate the library if we aren't
 * given an absolute pathname.
 *
 * To use this, declare a variable of type OpenMMPlugin, then when
 * you are ready to load the library and check whether it succeeds,
 * call the load() method:
 * <pre>
 *      OpenMMPlugin openMM;
 *      if (!openMM.load())
 *          failed to load -- see Plugin for how to get error message
 *      otherwise it loaded and you can use it
 * </pre>
 * It is harmless to call load() multiple times; if the library is
 * already loaded it will just return true. The library will unload
 * when the OpenMMPlugin object goes out of scope, or when unload()
 * is explicitly called.
 */

class OpenMMPlugin : public SimTK::Plugin {
public:
    // Constructor supplies default base name for this Plugin and sets
    // up search rule, which in this case is the lib/plugins directory
    // of the SimTK installation directory.
    OpenMMPlugin() : SimTK::Plugin("OpenMMPlugin") {
        // // Install directory is either the contents of this environment
        // // variable if it exists, or if not then it is the system default
        // // installation directory with /SimTK appended. Then /lib/plugins
        // // is added to that.


        char * envvar = std::getenv("OpenMMPlugin_PATH");
        if ( envvar != NULL ) { addSearchDirectory( envvar ); }

        // Quick & dirty fix for Laurentiu
        //addSearchDirectory("/home/pcuser/git4/Robosample/build/release/Molmodel/");
        //addSearchDirectory("/home/teo/2022/SNEED/Robosample/build/release/Molmodel/");

        // addSearchDirectory(SimTK::Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK")
        //                     + "lib/plugins");
        //addSearchDirectory("/usr/local/lib/plugins/");
        addSearchDirectory("/home/teo/2022/SNEED/Robosample/install/release/Simbody01/lib/plugins/");

        std::cout << "BLABLA1" << getSearchPath();
        /*std::cout << A;
        for (auto a : A) {
            std::cout << "SEARCH PATH: " << a << std::endl;
        }*/

        // The above works for a standard installation. We don't do that here.
        // When installing locally, we must look for the library locally.
        // The build script writes a file containing the path to the plugin.
        // We load that file and search the lib at the specified path.
        //std::ifstream fin(SimTK::Pathname::getThisExecutableDirectory() + "/openmmplugin");
        //std::string dir((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
        //addSearchDirectory(dir);
    }

    // This is the only defined exported method for this kind of plugin.
    // Call it to get an interface object you can use to communicate with
    // OpenMM. (Don't forget to delete that object when you're done.)
    SimTK_PLUGIN_DEFINE_FUNCTION1(OpenMMPluginInterface*,
                                  SimTK_createOpenMMPluginInterface,
                                  const SimTK::DuMMForceFieldSubsystemRep&);

    std::string OpenMMPluginPATH;

private:
};

#endif //SimTK_MOLMODEL_OPENMM_PLUGIN_H_



