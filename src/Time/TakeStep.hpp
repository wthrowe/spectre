// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Time/Actions/ChangeStepSize.hpp"
#include "Time/Actions/RecordTimeStepperData.hpp"
#include "Time/Actions/UpdateU.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace Parallel::Tags {
struct Metavariables;
}  // namespace Parallel::Tags
/// \endcond

/// Bundled method for recording the current system state in the history, and
/// updating the evolved variables and step size.
///
/// This function is used to encapsulate any needed logic for updating the
/// system.  In LTS mode, this includes adjusting the size of the next step.
/// If \p allow_step_rejection is true, current step size may also be rejected
/// and adjusted.  The \p allow_step_rejection argument is ignored in GTS mode.
template <typename System, bool LocalTimeStepping,
          typename StepChoosersToUse = AllStepChoosers, typename DbTags>
void take_step(const gsl::not_null<db::DataBox<DbTags>*> box,
               const bool allow_step_rejection) {
  record_time_stepper_data<System>(box);
  if constexpr (LocalTimeStepping) {
    do {
      update_u<System>(box);
    } while (
        not change_step_size<StepChoosersToUse>(box, allow_step_rejection));
  } else {
    update_u<System>(box);
  }
}
