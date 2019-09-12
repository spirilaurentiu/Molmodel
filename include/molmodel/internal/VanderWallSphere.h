#ifndef SimTK_MOLMODEL_VANDERWALLSPHERE_H_
#define SimTK_MOLMODEL_VANDERWALLSPHERE_H_

#include "SimTKsimbody.h"
#include "molmodel/internal/common.h"
#include "molmodel/internal/DuMMForceFieldSubsystem.h"

namespace SimTK {

class SimTK_MOLMODEL_EXPORT VanderWallSphere : public Force::Custom {
public:
	VanderWallSphere(GeneralForceSubsystem& forces, DuMMForceFieldSubsystem& dumm, Vec3 center, Real radius, Real vdwRadius, Real wellDepth);
};

} // namespace SimTK

#endif

