// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeSteppers/HeunImex.hpp"

#include "Time/EvolutionOrdering.hpp"
#include "Time/History.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace TimeSteppers {

size_t HeunImex::order() const { return 2; }

size_t HeunImex::error_estimate_order() const { return 1; }

uint64_t HeunImex::number_of_substeps() const { return 2; }

uint64_t HeunImex::number_of_substeps_for_error() const {
  return number_of_substeps();
}

size_t HeunImex::number_of_past_steps() const { return 0; }

double HeunImex::stable_step() const { return 1.0; }

TimeStepId HeunImex::next_time_id(const TimeStepId& current_id,
                                  const TimeDelta& time_step) const {
  switch (current_id.substep()) {
    case 0:
      ASSERT(current_id.substep_time() == current_id.step_time(),
             "Wrong substep time");
      return {current_id.time_runs_forward(), current_id.slab_number(),
              current_id.step_time(), 1, current_id.step_time() + time_step};
    case 1:
      ASSERT(current_id.substep_time() == current_id.step_time() + time_step,
             "Wrong substep time");
      return {current_id.time_runs_forward(), current_id.slab_number(),
              current_id.step_time() + time_step};
    default:
      ERROR("Bad substep value: " << current_id.substep());
  }
}

TimeStepId HeunImex::next_time_id_for_error(const TimeStepId& current_id,
                                            const TimeDelta& time_step) const {
  return next_time_id(current_id, time_step);
}

template <typename T>
void HeunImex::update_u_impl(const gsl::not_null<T*> u,
                             const gsl::not_null<UntypedHistory<T>*> history,
                             const TimeDelta& time_step) const {
  ASSERT(history->integration_order() == 2,
         "Fixed-order stepper cannot run at order "
         << history->integration_order());
  const size_t substep = (history->end() - 1).time_step_id().substep();

  // Clean up old history
  if (substep == 0) {
    history->mark_unneeded(history->end() - 1);
  }

  switch (substep) {
    case 0: {
      *u += time_step.value() * *history->begin().derivative();
      break;
    }
    case 1: {
      *u += 0.5 * time_step.value() *
            (*(history->begin() + 1).derivative() -
             *history->begin().derivative());
      break;
    }
    default:
      ERROR("Bad substep value: " << substep);
  }
}

template <typename T>
bool HeunImex::update_u_impl(const gsl::not_null<T*> u,
                             const gsl::not_null<T*> u_error,
                             const gsl::not_null<UntypedHistory<T>*> history,
                             const TimeDelta& time_step) const {
  *u_error = *u;
  update_u_impl(u, history, time_step);
  // error estimate is only available when completing a full step
  if ((history->end() - 1).time_step_id().substep() == 1) {
    *u_error -= *u;
    return true;
  }
  return false;
}

template <typename T>
bool HeunImex::dense_update_u_impl(const gsl::not_null<T*> u,
                                   const UntypedHistory<T>& history,
                                   const double time) const {
  if ((history.end() - 1).time_step_id().substep() != 0) {
    return false;
  }
  const double step_start = history.front().value();
  const double step_end = history.back().value();
  if (time == step_end) {
    // Special case necessary for dense output at the initial time,
    // before taking a step.
    return true;
  }
  const evolution_less<double> before{step_end > step_start};
  if (history.size() == 1 or before(step_end, time)) {
    return false;
  }

  const double time_step = step_end - step_start;
  const double output_fraction = (time - step_start) / time_step;
  ASSERT(output_fraction >= 0, "Attempting dense output at time " << time
         << ", but already progressed past " << step_start);
  ASSERT(output_fraction <= 1,
         "Requested time (" << time << " not within step [" << step_start
         << ", " << step_end << "]");

  // We only need linear interpolation to get the correct convergence
  // order.  We have enough data to produce a smooth extension,
  // though, and that's nice, even if it doesn't formally improve the
  // accuracy.
  *u -= time_step * (1.0 - output_fraction) *
        ((1.0 - output_fraction) *
             (0.5 * *history.begin().derivative() +
              (0.5 + output_fraction) * *(history.begin() + 1).derivative()) +
         square(output_fraction) * *(history.begin() + 2).derivative());
  return true;
}

template <typename T>
bool HeunImex::can_change_step_size_impl(
    const TimeStepId& time_id, const UntypedHistory<T>& /*history*/) const {
  return time_id.substep() == 0;
}

TIME_STEPPER_DEFINE_OVERLOADS(HeunImex)

bool operator==(const HeunImex& /*lhs*/, const HeunImex& /*rhs*/) {
  return true;
}

bool operator!=(const HeunImex& lhs, const HeunImex& rhs) {
  return not(rhs == lhs);
}
}  // namespace TimeSteppers

PUP::able::PUP_ID TimeSteppers::HeunImex::my_PUP_ID = 0;  // NOLINT
