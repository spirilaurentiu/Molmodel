#ifndef VOXEL_HASH_HPP_
#define VOXEL_HASH_HPP_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include <vector>
#include <map>
// #include "md_units.hpp"

#ifdef __GNUC__
#include <unordered_map>
#define HASH_MAP_NAMESPACE __gnu_cxx
#else
#include <unordered_map>
#define HASH_MAP_NAMESPACE stdext
#endif

// Subystem for managing atom locations.
// Depends on matter subsystem

namespace SimTK {

// Elide use of boost::units for now
namespace units { namespace md 
{
    typedef Real dimensionless_t;
    typedef Real length_t;
    typedef Real mass_t;
    typedef Real area_t;
    typedef Real inverse_volume_t;
    typedef Vec3 location_t;

    static const Real nanometers = 1.0;
    static const Real daltons = 1.0;
    static const Real cubic_nanometer = 1.0;
    static const Real square_nanometers = 1.0;
}}

template< class T >
class SimTK_MOLMODEL_EXPORT VoxelHash 
{
public:
    typedef std::pair<SimTK::Vec3, T> VoxelItem;
    typedef std::vector< VoxelItem > Voxel;

    class VoxelIndex {
    public:
        VoxelIndex(int x, int y, int z, size_t numBuckets) 
            : x(x), y(y), z(z), numBuckets(numBuckets) 
        {}

        // operator<() needed for map
        bool operator<(const VoxelIndex& other) const {
            if      (x < other.x) return true;
            else if (x > other.x) return false;
            else if (y < other.y) return true;
            else if (y > other.y) return false;
            else if (z < other.z) return true;
            else return false;
        }
        bool operator==(const VoxelIndex& other) const {
            if      (x != other.x) return false;
            else if (y != other.y) return false;
            else if (z != other.z) return false;
            else return true;
        }

        // size_t conversion for use by hash_map
        operator size_t() const {
            return hash_value();
        }
        
        size_t hash_value() const {
            // multiply each index by a different large prime number
            size_t n = (size_t)(x * 0x8da6b343 + y * 0xd8163841 + z * 0xcb1ab31f);
            n = n % numBuckets;
            if (n < 0) n += numBuckets;
            return n;            
        }

        int x;
        int y;
        int z;
        size_t numBuckets;
    };

private:
    VoxelIndex getVoxelIndex(const Vec3& location) const 
    {
        SimTK::Real voxSize = voxelSize;
        int x = int(floor(location[0] / voxSize));
        int y = int(floor(location[1] / voxSize));
        int z = int(floor(location[2] / voxSize));
        return VoxelIndex(x, y, z, numBuckets);
    }

    SimTK::units::md::length_t voxelSize;

    // TODO - replace std::map with hash_map or unordered map
    // std::map<VoxelIndex, Voxel> voxelMap;
    std::map<VoxelIndex, Voxel> voxelMap;
    // HASH_MAP_NAMESPACE::hash_map<VoxelIndex, Voxel> voxelMap;
    size_t numBuckets;
    
public:
    VoxelHash(SimTK::units::md::length_t voxelSize, int numBuckets) 
        : voxelSize(voxelSize), numBuckets(numBuckets)
    {}

    void insert(const T& item, const SimTK::Vec3& location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if ( voxelMap.find(voxelIndex) == voxelMap.end() ) voxelMap[voxelIndex] = Voxel();
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back( VoxelItem(location, item) );
    }

    void findNeighbors(std::vector<T>& neighbors, const SimTK::Vec3& locationI, SimTK::units::md::length_t maxDistance) const
    {
        SimTK::units::md::area_t maxDistanceSquared(maxDistance * maxDistance);

        int dIndex = int(maxDistance / voxelSize) + 1; // How may voxels away do we have to look?
        VoxelIndex centerVoxelIndex(getVoxelIndex(locationI));

        for (int x = centerVoxelIndex.x - dIndex; x <= centerVoxelIndex.x + dIndex; ++x)
            for (int y = centerVoxelIndex.y - dIndex; y <= centerVoxelIndex.y + dIndex; ++y)
                for (int z = centerVoxelIndex.z - dIndex; z <= centerVoxelIndex.z + dIndex; ++z)
                {
                    const VoxelIndex voxelIndex(x, y, z, numBuckets);

                    // TODO - store a list of neighbors for each voxel
                    if ( voxelMap.find(voxelIndex) == voxelMap.end() ) continue; // no such voxel; skip

                    const Voxel& voxel = voxelMap.find(voxelIndex)->second;
                    typename std::vector<VoxelItem>::const_iterator itemIter;
                    for (itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter)
                    {
                        const SimTK::Vec3 r = locationI - itemIter->first;
                        SimTK::units::md::area_t dSquared(dot(r, r) * units::md::square_nanometers);
                        if (dSquared > maxDistanceSquared) continue; // beyond cutoff

                        neighbors.push_back(itemIter->second); // store neighbor
                    }
                }
    }

};

} // namespace SimTK


#endif /* VOXEL_HASH_HPP_ */
