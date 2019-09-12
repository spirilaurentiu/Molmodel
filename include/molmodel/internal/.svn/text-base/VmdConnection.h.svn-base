#ifndef SimTK_MOLMODEL_VMDCONNECTION_H_
#define SimTK_MOLMODEL_VMDCONNECTION_H_

#include "molmodel/internal/common.h"
#include <vector>

namespace SimTK {
    
class SimTK_MOLMODEL_EXPORT VmdFloat3 
{
public:
    VmdFloat3(float x, float y, float z)
    {
        d[0] = x;
        d[1] = y;
        d[2] = z;
    }
        
    float& operator[](int i) {return d[i];}
    const float& operator[](int i) const {return d[i];}

private:
    float d[3];
};

class SimTK_MOLMODEL_EXPORT VmdConnection
{
public:
    
    explicit VmdConnection(int portNumber);
    
    ~VmdConnection();
    
    bool checkForConnection() ;
    
    bool clientIsConnected() const;
    
    void sendCoordinates(std::vector<VmdFloat3> coords);

    void closeClientConnection();

    int getSocketNumber() const;

private:
    int socketNumber;
    void *socket;
    void *clientSocket;
};

} // namespace SimTK

#endif // SimTK_MOLMODEL_VMDCONNECTION_H_
