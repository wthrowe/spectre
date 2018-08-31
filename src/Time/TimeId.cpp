// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/TimeId.hpp"

#include <boost/functional/hash.hpp>
#include <ostream>
#include <pup.h>

void TimeId::canonicalize() noexcept {
  if (is_at_slab_boundary() and
      (time_runs_forward_ ? time_.is_at_slab_end()
                          : time_.is_at_slab_start())) {
    const Slab new_slab =
        time_runs_forward_ ? time_.slab().advance() : time_.slab().retreat();
    ++slab_number_;
    time_ = time_.with_slab(new_slab);
    step_time_ = time_;
  }
}

void TimeId::pup(PUP::er& p) noexcept {
  p | time_runs_forward_;
  p | slab_number_;
  p | step_time_;
  p | substep_;
  p | time_;
}

bool operator==(const TimeId& a, const TimeId& b) noexcept {
  ASSERT(a.time_runs_forward() == b.time_runs_forward(),
         "Time is not running in a consistent direction");
  const bool equal = a.slab_number() == b.slab_number() and
                     a.step_time() == b.step_time() and
                     a.substep() == b.substep();
  // This could happen if we have a local-time-stepping substep
  // method.  If we implement any of those the comparison operators
  // for TimeId will have to be revisited.
  ASSERT(not equal or a.time() == b.time(),
         "IDs at same step and substep but different times");
  return equal;
}

bool operator!=(const TimeId& a, const TimeId& b) noexcept {
  return not(a == b);
}

bool operator<(const TimeId& a, const TimeId& b) noexcept {
  ASSERT(a.time_runs_forward() == b.time_runs_forward(),
         "Time is not running in a consistent direction");
  return a.slab_number() < b.slab_number() or
         (a.slab_number() == b.slab_number() and
          ((a.time_runs_forward() ? a.step_time() < b.step_time()
                                  : a.step_time() > b.step_time()) or
           (a.step_time() == b.step_time() and a.substep() < b.substep())));
}
bool operator<=(const TimeId& a, const TimeId& b) noexcept {
  return not(b < a);
}
bool operator>(const TimeId& a, const TimeId& b) noexcept {
  return b < a;
}
bool operator>=(const TimeId& a, const TimeId& b) noexcept {
  return not(a < b);
}

std::ostream& operator<<(std::ostream& s, const TimeId& id) noexcept {
  return s << id.slab_number() << ':' << id.step_time() << ':' << id.substep()
           << ':' << id.time();
}

size_t hash_value(const TimeId& id) noexcept {
  size_t h = 0;
  boost::hash_combine(h, id.slab_number());
  boost::hash_combine(h, id.step_time());
  boost::hash_combine(h, id.substep());
  boost::hash_combine(h, id.time());
  return h;
}

// clang-tidy: do not modify std namespace (okay for hash)
namespace std {  // NOLINT
size_t hash<TimeId>::operator()(const TimeId& id) const noexcept {
  return boost::hash<TimeId>{}(id);
}
}  // namespace std
