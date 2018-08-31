// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class TimeId.

#pragma once

#include <cstddef>
#include <functional>
#include <iosfwd>

#include "ErrorHandling/Assert.hpp"
#include "Time/Time.hpp"

namespace PUP {
class er;
}  // namespace PUP

/// \ingroup TimeGroup
///
/// A unique identifier for the temporal state of an integrated
/// system.
class TimeId {
 public:
  TimeId() = default;
  /// Create a TimeId at the start of a step
  TimeId(const bool time_runs_forward, const size_t slab_number,
         const Time& time) noexcept
      : time_runs_forward_(time_runs_forward),
        slab_number_(slab_number),
        step_time_(time),
        substep_(0),
        time_(time) {}
  /// Create a TimeId at a substep at time `time` in a step starting
  /// at time `step_time`.
  TimeId(const bool time_runs_forward, const size_t slab_number,
         const Time& step_time, const size_t substep, const Time& time) noexcept
      : time_runs_forward_(time_runs_forward),
        slab_number_(slab_number),
        step_time_(step_time),
        substep_(substep),
        time_(time) {
    ASSERT(substep_ != 0 or step_time_ == time_,
           "Initial substep must align with the step.");
  }

  bool time_runs_forward() const noexcept { return time_runs_forward_; }
  size_t slab_number() const noexcept { return slab_number_; }
  /// Time at the start of the current step
  const Time& step_time() const noexcept { return step_time_; }
  size_t substep() const noexcept { return substep_; }
  /// Time of the current substep
  const Time& time() const noexcept { return time_; }

  bool is_at_slab_boundary() const noexcept {
    return substep_ == 0 and time_.is_at_slab_boundary();
  }

  /// If this object represents the (evolution-direction) end of a
  /// slab advance the slab, increase the slab number, and transfer
  /// the time to the new slab.
  void canonicalize() noexcept;

  // clang-tidy: google-runtime-references
  void pup(PUP::er& p) noexcept;  // NOLINT

 private:
  bool time_runs_forward_{false};
  size_t slab_number_{0};
  Time step_time_{};
  size_t substep_{0};
  Time time_{};
};

bool operator==(const TimeId& a, const TimeId& b) noexcept;
bool operator!=(const TimeId& a, const TimeId& b) noexcept;
bool operator<(const TimeId& a, const TimeId& b) noexcept;
bool operator<=(const TimeId& a, const TimeId& b) noexcept;
bool operator>(const TimeId& a, const TimeId& b) noexcept;
bool operator>=(const TimeId& a, const TimeId& b) noexcept;

std::ostream& operator<<(std::ostream& s, const TimeId& id) noexcept;

size_t hash_value(const TimeId& id) noexcept;

namespace std {
template <>
struct hash<TimeId> {
  size_t operator()(const TimeId& id) const noexcept;
};
}  // namespace std
