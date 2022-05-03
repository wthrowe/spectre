// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <cstdint>
#include <pup.h>

#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"  // IWYU pragma: keep
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class TimeDelta;
struct TimeStepId;
namespace TimeSteppers {
template <typename T>
class UntypedHistory;
}  // namespace TimeSteppers
/// \endcond

namespace TimeSteppers {

/*!
 * \ingroup TimeSteppersGroup
 *
 * Heun's method, a second order Runge-Kutta method.
 *
 * The CFL factor/stable step size is 1.0, the same as Euler's method.
 */
class HeunImex : public TimeStepper {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help = {"Heun's method."};

  HeunImex() = default;
  HeunImex(const HeunImex&) = default;
  HeunImex& operator=(const HeunImex&) = default;
  HeunImex(HeunImex&&) = default;
  HeunImex& operator=(HeunImex&&) = default;
  ~HeunImex() override = default;

  size_t order() const override;

  size_t error_estimate_order() const override;

  uint64_t number_of_substeps() const override;

  uint64_t number_of_substeps_for_error() const override;

  size_t number_of_past_steps() const override;

  double stable_step() const override;

  TimeStepId next_time_id(const TimeStepId& current_id,
                          const TimeDelta& time_step) const override;

  TimeStepId next_time_id_for_error(const TimeStepId& current_id,
                                    const TimeDelta& time_step) const override;

  WRAPPED_PUPable_decl_template(HeunImex);  // NOLINT

  explicit HeunImex(CkMigrateMessage* /*unused*/) {}

 private:
  template <typename T>
  void update_u_impl(gsl::not_null<T*> u,
                     gsl::not_null<UntypedHistory<T>*> history,
                     const TimeDelta& time_step) const;

  template <typename T>
  bool update_u_impl(gsl::not_null<T*> u, gsl::not_null<T*> u_error,
                     gsl::not_null<UntypedHistory<T>*> history,
                     const TimeDelta& time_step) const;

  template <typename T>
  bool dense_update_u_impl(gsl::not_null<T*> u,
                           const UntypedHistory<T>& history, double time) const;

  template <typename T>
  bool can_change_step_size_impl(const TimeStepId& time_id,
                                 const UntypedHistory<T>& history) const;

  TIME_STEPPER_DECLARE_OVERLOADS
};

bool operator==(const HeunImex& lhs, const HeunImex& rhs);

bool operator!=(const HeunImex& lhs, const HeunImex& rhs);
}  // namespace TimeSteppers
