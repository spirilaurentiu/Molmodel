#include "molmodel/internal/VmdConnection.h"

extern "C" {
#include "vmdsock.h"
}
#include "imd.h"
#include <iostream>
#include <vector>

// Microsoft build environment
#ifdef _MSC_VER
#include <Windows.h>
#define sleep(x) Sleep(1000*x)
#else
//HOREA
#include <unistd.h>
#endif

using namespace std;

static bool debugVmdConnection = true;

///////////////////////////////
//// VmdConnection methods ////
///////////////////////////////

namespace SimTK {

VmdConnection::VmdConnection(int portNumber)
:   socket(NULL),
    clientSocket(NULL),
    socketNumber(portNumber)
{
    if (debugVmdConnection)
        cerr << "Opening socket connection to VMD on port " << portNumber << "..." << endl;
    
    vmdsock_init();
    socket = vmdsock_create(); // uses malloc...
    vmdsock_bind(socket, portNumber);
    vmdsock_listen(socket);
}

VmdConnection::~VmdConnection()
{
    if (debugVmdConnection)
        cerr << "Closing socket connection to VMD..." << endl;
    
    if ( clientIsConnected() )
    {
        closeClientConnection();
    }
    
    if (socket != NULL)
    {
        vmdsock_shutdown(socket);
        vmdsock_destroy(socket);
    }
}

bool VmdConnection::checkForConnection() 
{
    if ( ! clientIsConnected() )
    {
        if (socket == NULL) 
        {
            cerr << "ERROR: socket is NULL" << endl;
            return false;
        }
        if (vmdsock_selread(socket, 0) > 0)
        {
            clientSocket = vmdsock_accept(socket);
            if ( imd_handshake(clientSocket) )
            {
                clientSocket = NULL;
            }
        }
    }

    if ( clientIsConnected() )
    {
        sleep(1); // wait a second for VMD to respond

        int length;
        if ( (vmdsock_selread(clientSocket, 0) != 1) ||
             (imd_recv_header(clientSocket, &length) != IMD_GO) ) 
        {
            closeClientConnection();
        }
    }
    
    if ( clientIsConnected() ) 
        return true;
    else
        return false;
}

bool VmdConnection::clientIsConnected() const 
{
    return (NULL != clientSocket);
}

void VmdConnection::sendCoordinates(std::vector<VmdFloat3> coords)
{
    if ( imd_send_fcoords(clientSocket, coords.size(), &coords[0][0]) )
    {
        closeClientConnection();
    }
}

void VmdConnection::closeClientConnection() 
{
    imd_disconnect(clientSocket);
    vmdsock_shutdown(clientSocket);
    vmdsock_destroy(clientSocket);
    clientSocket = NULL;
}

int VmdConnection::getSocketNumber() const {return socketNumber;}

} // namespace SimTK


