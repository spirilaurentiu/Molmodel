#include "molmodel/internal/PeriodicVmdReporter.h"

#include "molmodel/internal/Compound.h"
#include "molmodel/internal/VmdConnection.h"
#include "molmodel/internal/CompoundSystem.h"
#include <iostream>
#include <iomanip>
#include <vector>

// Microsoft build environment
#ifdef _MSC_VER
#include <Windows.h>
#define sleep(x) Sleep(1000*x)
#else
//HOREA
#include <unistd.h>
#endif

namespace SimTK {

void PeriodicVmdReporter::handleEvent(const State& state) const 
{
    // checkForConnection() can be slow when there is none
    // so only check every so often
    static int timeStepNumber = 0;
    int checkFrequency = 20;
    
    // Try our best to connect to vmd
    if (! vmdConnection.clientIsConnected() )
    {
        if (blockWaitingForVmdConnection)
        {
            while ( ! vmdConnection.clientIsConnected() ) 
            {
                std::cerr << "Waiting for vmd IMD connection on port " << vmdConnection.getSocketNumber() << "..." << std::endl;

                vmdConnection.checkForConnection();

                sleep(1);

            }  
        }
        else
        {
            if ( 0 == (timeStepNumber % checkFrequency) )
            vmdConnection.checkForConnection();
        } 
    }

    // Write the coordinates, if we are connected to vmd
    if ( vmdConnection.clientIsConnected() )
    {
        system.realize(state, Stage::Position);
        
        // To ensure that the atom coordinates are in the same order as PDB output,
        // create PDB output, and parse that to generate the coordinates
        std::stringstream pdbString;
        
        // 1) Write pdb coordinates into the string
        int nextAtomSerialNumber = 1; // atom serial number for each compound picks up where previous compound left off
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
            system.getCompound(c).writePdb(state, pdbString, nextAtomSerialNumber);
        
        // 2) Read pdb coordinates from the string
        std::vector<VmdFloat3> vmdCoordinates;
        char lineBuffer[102];
        while ( ! pdbString.eof() )
        {
            pdbString.getline(lineBuffer, 100);
            std::string line(lineBuffer);
            
            if ( (line.length() > 55) && ((line.substr(0, 6) == "ATOM  ") || (line.substr(0, 6) == "HETATM")) )
            {
                float x, y, z;
                std::istringstream(line.substr(30, 8)) >> x;
                std::istringstream(line.substr(38, 8)) >> y;
                std::istringstream(line.substr(46, 8)) >> z;
                vmdCoordinates.push_back(VmdFloat3(x, y, z));
            }
        }
        
        vmdConnection.sendCoordinates(vmdCoordinates);
    }

    ++timeStepNumber;
    
}

} // namespace SimTK

