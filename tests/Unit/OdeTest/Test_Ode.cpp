// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Factory.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Times.hpp"
#include "Evolution/EventsAndDenseTriggers/Tags.hpp"
#include "Evolution/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "Framework/ActionTesting.hpp"
#include "Options/Comparator.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Actions/RecordTimeStepperData.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/Actions/UpdateU.hpp"
#include "Time/History.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Time/Slab.hpp"
#include "Time/StepChoosers/ErrorControl.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/RungeKutta3.hpp"
#include "Utilities/MakeVector.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename Metavariables>
struct ComputeTimeDerivativeFIXME {
  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    using system = typename Metavariables::system;
    using variables_tag = typename system::variables_tag;
    using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
    db::mutate<dt_variables_tag>(make_not_null(&box), [](
        const gsl::not_null<typename dt_variables_tag::type*> dt_vars) {
      dt_vars->initialize(dt_vars->number_of_grid_points(), 1.0);
    });
    return {std::move(box)};
  }
};

struct Var : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct System {
  using variables_tag = Tags::Variables<tmpl::list<Var>>;
  static constexpr bool has_primitive_and_conservative_vars = false;
};

class RecordTime : public Event {
 public:
  RecordTime() = default;

  explicit RecordTime(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
  WRAPPED_PUPable_decl_template(RecordTime);  // NOLINT
#pragma GCC diagnostic pop

  using compute_tags_for_observation_box = tmpl::list<>;

  using argument_tags = tmpl::list<Tags::Time>;

  template <typename Metavariables, typename ArrayIndex,
            typename ParallelComponent>
  void operator()(const double time,
                  const Parallel::GlobalCache<Metavariables>& /*cache*/,
                  const ArrayIndex& /*array_index*/,
                  const ParallelComponent* const /*meta*/) const {
    recorded_time = time;
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex,
            typename ParallelComponent>
  bool is_ready(const Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const ParallelComponent* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return false; }

  static double recorded_time;
};

double RecordTime::recorded_time = 0.0;
PUP::able::PUP_ID RecordTime::my_PUP_ID = 0;

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using Phase = typename metavariables::Phase;
  // Avoid `variables_tag` to make sure it's not accidentally used
  // from somewhere else.
  using system_variables_tag = typename metavariables::system::variables_tag;

  using const_global_cache_tags = tmpl::list<Tags::TimeStepper<TimeStepper>>;
  using simple_tags = tmpl::flatten<db::AddSimpleTags<
      system_variables_tag, db::add_tag_prefix<Tags::dt, system_variables_tag>,
      Tags::HistoryEvolvedVariables<system_variables_tag>, Tags::TimeStepId,
      Tags::Next<Tags::TimeStepId>, Tags::TimeStep, Tags::Next<Tags::TimeStep>,
      Tags::Time, Tags::NeverUsingTimeSteppingErrorControl,
      evolution::Tags::EventsAndDenseTriggers,
      evolution::Tags::PreviousTriggerTime>>;
  using compute_tags = db::AddComputeTags<>;

  using step_actions = tmpl::list<
      ComputeTimeDerivativeFIXME<metavariables>,
      Actions::RecordTimeStepperData<>,
      evolution::Actions::RunEventsAndDenseTriggers<>,
      Actions::UpdateU<>>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Phase, Phase::Initialization,
          tmpl::list<
              ActionTesting::InitializeDataBox<simple_tags, compute_tags>,
              Actions::SetupDataBox>>,

      Parallel::PhaseActions<
          Phase, Phase::InitializeTimeStepperHistory,
          SelfStart::self_start_procedure<step_actions,
                                          typename metavariables::system>>,

      Parallel::PhaseActions<
          Phase, Phase::Evolve,
          tmpl::list<Actions::RunEventsAndTriggers, step_actions,
                     Actions::AdvanceTime>>>;
};

struct Metavariables {
  using system = System;
  static constexpr bool local_time_stepping = false;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DenseTrigger, DenseTriggers::standard_dense_triggers>,
        tmpl::pair<Event, tmpl::list<Events::Completion, RecordTime>>,
        tmpl::pair<TimeSequence<double>,
                   TimeSequences::all_time_sequences<double>>,
        tmpl::pair<TimeSequence<std::uint64_t>,
                   TimeSequences::all_time_sequences<std::uint64_t>>,
        tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                         Triggers::time_triggers>>>;
  };

  using component_list = tmpl::list<Component<Metavariables>>;
  enum class Phase {
    Initialization,
    InitializeTimeStepperHistory,
    Evolve,
    Exit
  };
};

void test(const bool forward_in_time, const bool use_dense_trigger) {
  using component = Component<Metavariables>;
  using Vars = System::variables_tag::type;
  using DtVars = db::add_tag_prefix<Tags::dt, System::variables_tag>::type;

  const double final_trigger_time = forward_in_time ? 5.3 : -4.5;

  const Slab initial_slab(1.0, 2.0);
  const auto initial_time =
      forward_in_time ? initial_slab.start() : initial_slab.end();
  const auto initial_time_step =
      (forward_in_time ? 1 : -1) * initial_slab.duration();
  Vars initial_value(1);
  get(get<Var>(initial_value)) = 1.0;

  TimeSteppers::RungeKutta3 time_stepper{};
  EventsAndTriggers::Storage events_and_triggers;
  evolution::EventsAndDenseTriggers::ConstructionType events_and_dense_triggers;
  auto events = make_vector<std::unique_ptr<Event>>(
      std::make_unique<Events::Completion>(), std::make_unique<RecordTime>());
  if (use_dense_trigger) {
    events_and_dense_triggers.emplace(
        std::make_unique<DenseTriggers::Times>(
            std::make_unique<TimeSequences::Specified<double>>(
                std::vector{final_trigger_time})),
        std::move(events));
  } else {
    events_and_triggers.emplace(
        std::make_unique<Triggers::TimeCompares>(
            forward_in_time
                ? Options::Comparator::Comparison::GreaterThanOrEqualTo
                : Options::Comparator::Comparison::LessThanOrEqualTo,
            final_trigger_time),
        std::move(events));
  }
  ActionTesting::MockRuntimeSystem<Metavariables> runner{
      {std::make_unique<decltype(time_stepper)>(time_stepper),
       EventsAndTriggers(std::move(events_and_triggers))}};
  ActionTesting::emplace_component_and_initialize<component>(
      make_not_null(&runner), 0,
      {initial_value, DtVars(1),
       TimeSteppers::History<Vars>(time_stepper.order()),
       TimeStepId{},
       TimeStepId(forward_in_time,
                  -static_cast<int64_t>(time_stepper.number_of_past_steps()),
                  initial_time),
       initial_time_step, initial_time_step, initial_time.value(), false,
       evolution::EventsAndDenseTriggers(std::move(events_and_dense_triggers)),
       std::nullopt});

  ActionTesting::next_action<component>(make_not_null(&runner), 0);

  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::InitializeTimeStepperHistory);

  while (not ActionTesting::get_terminate<component>(runner, 0)) {
    ActionTesting::next_action<component>(make_not_null(&runner), 0);
  }

  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::Evolve);

  while (not ActionTesting::get_terminate<component>(runner, 0)) {
    ActionTesting::next_action<component>(make_not_null(&runner), 0);
  }

  const double final_time =
      ActionTesting::get_databox_tag<component, Tags::Time>(runner, 0);
  if (use_dense_trigger) {
    CHECK(RecordTime::recorded_time == final_trigger_time);
  } else {
    const evolution_less<double> less{forward_in_time};
    CHECK(not less(final_time, final_trigger_time));
    CHECK(less(final_time, final_trigger_time + initial_time_step.value()));
  }
  CHECK_ITERABLE_APPROX(
      get(ActionTesting::get_databox_tag<component, Var>(runner, 0)),
      get(get<Var>(initial_value)) + (final_time - initial_time.value()));
}
}  // namespace

SPECTRE_TEST_CASE("OdeTest", "[Unit]") {
  Parallel::register_classes_with_charm<TimeSteppers::RungeKutta3>();
  Parallel::register_factory_classes_with_charm<Metavariables>();
  test(true, true);
  test(true, false);
  test(false, true);
  test(false, false);
}
