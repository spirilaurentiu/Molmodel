#ifndef SimTK_TINKER_AMBER99_PARAMS_H_
#define SimTK_TINKER_AMBER99_PARAMS_H_

#include "molmodel/internal/common.h"

namespace SimTK {

class DuMMForceFieldSubsystem;

void SimTK_MOLMODEL_EXPORT populateAmber99Params(DuMMForceFieldSubsystem& dumm);

} // namespace SimTK

#endif // SimTK_TINKER_AMBER99_PARAMS_H_
