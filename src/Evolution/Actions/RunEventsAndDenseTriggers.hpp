// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "Evolution/EventsAndDenseTriggers/Tags.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace Tags {
struct TimeStep;
}  // namespace Tags
/// \endcond

namespace evolution::Actions {
/// \ingroup ActionsGroup
/// \ingroup EventsAndTriggersGroup
/// \brief Run the events and dense triggers
///
/// Uses:
/// - DataBox: EventsAndDenseTriggers, as required by events and triggers
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies: nothing
template <typename PrimFromCon = void>
struct RunEventsAndDenseTriggers {
 private:
  // RAII object to restore the time and variables changed by dense
  // output.
  template <typename DbTags, typename Tag>
  class StateRestorer {
   public:
    StateRestorer(const gsl::not_null<db::DataBox<DbTags>*> box) noexcept
        : box_(box) {}

    void save() noexcept {
      // Only store the value the first time, because after that we
      // are seeing the value after the previous change instead of the
      // original.
      if (not value_.has_value()) {
        value_ = db::get<Tag>(*box_);
      }
    }

    ~StateRestorer() {
      if (value_.has_value()) {
        db::mutate<Tag>(
            box_,
            [this](const gsl::not_null<typename Tag::type*> value) noexcept {
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 11
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 11
              *value = *value_;
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 11
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 11
            });
      }
    }

   private:
    gsl::not_null<db::DataBox<DbTags>*> box_ = nullptr;
    std::optional<typename Tag::type> value_{};
  };

 public:
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTags>&&, Parallel::AlgorithmExecution> apply(
      db::DataBox<DbTags>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const component) noexcept {
    using system = typename Metavariables::system;
    using variables_tag = typename system::variables_tag;

    const auto& time_step_id = db::get<::Tags::TimeStepId>(box);
    if (time_step_id.slab_number() < 0) {
      // Skip dense output during self-start
      return {std::move(box), Parallel::AlgorithmExecution::Continue};
    }

    auto& events_and_dense_triggers =
        db::get_mutable_reference<::evolution::Tags::EventsAndDenseTriggers>(
            make_not_null(&box));

    const auto step_end =
        time_step_id.step_time() + db::get<::Tags::TimeStep>(box);
    const evolution_less<double> before{time_step_id.time_runs_forward()};

    StateRestorer<DbTags, ::Tags::Time> time_restorer(make_not_null(&box));
    StateRestorer<DbTags, variables_tag> variables_restorer(
        make_not_null(&box));
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
    auto primitives_restorer = [&box]() noexcept {
      if constexpr (system::has_primitive_and_conservative_vars) {
        return StateRestorer<DbTags, typename system::primitive_variables_tag>(
            make_not_null(&box));
      } else {
        (void)box;
        return 0;
      }
    }();
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
    for (;;) {
      const double next_trigger = events_and_dense_triggers.next_trigger(box);
      if (before(step_end.value(), next_trigger)) {
        return {std::move(box), Parallel::AlgorithmExecution::Continue};
      }

      if (db::get<::Tags::Time>(box) != next_trigger) {
        time_restorer.save();
        db::mutate<::Tags::Time>(
            make_not_null(&box),
            [&next_trigger](const gsl::not_null<double*> time) noexcept {
              *time = next_trigger;
            });
      }

      const auto triggered = events_and_dense_triggers.is_ready(
          box, cache, array_index, component);
      using TriggeringState = std::decay_t<decltype(triggered)>;
      switch (triggered) {
        case TriggeringState::NotReady:
          return {std::move(box), Parallel::AlgorithmExecution::Retry};
        case TriggeringState::NeedsEvolvedVariables:
          {
            using history_tag = ::Tags::HistoryEvolvedVariables<variables_tag>;
            bool dense_output_succeeded = false;
            variables_restorer.save();
            db::mutate<variables_tag>(
                make_not_null(&box),
                [&dense_output_succeeded, &next_trigger](
                    gsl::not_null<typename variables_tag::type*> vars,
                    const TimeStepper& stepper,
                    const typename history_tag::type& history) noexcept {
                  dense_output_succeeded =
                      stepper.dense_update_u(vars, history, next_trigger);
                },
                db::get<::Tags::TimeStepper<>>(box), db::get<history_tag>(box));
            if (not dense_output_succeeded) {
              // With LTS, we may need to wait for neighbor data
              // rather than take another step.
              static_assert(not Metavariables::local_time_stepping,
                            "LTS support for RunDenseEventsAndTriggers is not "
                            "implemented.");

              // Need to take another time step
              return {std::move(box), Parallel::AlgorithmExecution::Continue};
            }

            static_assert(system::has_primitive_and_conservative_vars !=
                              std::is_same_v<PrimFromCon, void>,
                          "Primitive update scheme not provided.");
            if constexpr (system::has_primitive_and_conservative_vars) {
              primitives_restorer.save();
              db::mutate_apply<PrimFromCon>(make_not_null(&box));
            }
          }
          [[fallthrough]];
        default:
          break;
      }

      events_and_dense_triggers.run_events(box, cache, array_index, component);
    }
  }
};

struct InitializeRunEventsAndDenseTriggers {
  using initialization_tags =
      tmpl::list<evolution::Tags::EventsAndDenseTriggers>;
  using initialization_tags_to_keep = initialization_tags;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTags>& box,
                    tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*component*/) noexcept {
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace evolution::Actions
