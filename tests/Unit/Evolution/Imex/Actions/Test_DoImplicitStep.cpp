// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Imex/Actions/DoImplicitStep.hpp"
#include "Evolution/Imex/GuessResult.hpp"
#include "Evolution/Imex/Mode.hpp"
#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Evolution/Imex/Tags/Mode.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/Heun2.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename Var>
struct Sector : tt::ConformsTo<imex::protocols::ImplicitSector> {
  using tensors = tmpl::list<Var>;

  struct source {
    using return_tags = tmpl::list<Tags::Source<Var>>;
    using argument_tags = tmpl::list<Var>;
    static void apply(const gsl::not_null<Scalar<DataVector>*> source,
                      const Scalar<DataVector>& var) {
      get(*source) = -get(var);
    }
  };

  using jacobian = imex::NoJacobianBecauseSolutionIsAnalytic;

  struct initial_guess {
    using return_tags = tmpl::list<Var>;
    using argument_tags = tmpl::list<>;
    static imex::GuessResult apply(
        const gsl::not_null<Scalar<DataVector>*> var,
        const Variables<tmpl::list<Var>>& inhomogeneous_terms,
        const double implicit_weight) {
      get(*var) = get(get<Var>(inhomogeneous_terms)) / (1.0 + implicit_weight);
      return imex::GuessResult::ExactSolution;
    }
  };

  using tags_from_evolution = tmpl::list<>;
  using simple_tags = tmpl::list<>;
  using compute_tags = tmpl::list<>;
  using source_prep = tmpl::list<>;
  using jacobian_prep = tmpl::list<>;
  using initial_guess_prep = tmpl::list<>;
};

struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct Var2 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct System : tt::ConformsTo<imex::protocols::ImexSystem> {
  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2>>;
  using implicit_sectors = tmpl::list<Sector<Var1>, Sector<Var2>>;
};

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using simple_tags =
      db::AddSimpleTags<Tags::TimeStepId, Tags::TimeStep, System::variables_tag,
                        imex::Tags::Mode,
                        Tags::TimeStepper<TimeSteppers::Heun2>,
                        imex::Tags::ImplicitHistory<Sector<Var1>>,
                        imex::Tags::ImplicitHistory<Sector<Var2>>>;

  using compute_tags = db::AddComputeTags<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 simple_tags, compute_tags>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             tmpl::list<imex::Actions::DoImplicitStep>>>;
};

struct Metavariables {
  using system = System;
  using component_list = tmpl::list<Component<Metavariables>>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Imex.Actions.DoImplicitStep",
                  "[Unit][Evolution][Actions]") {
  Parallel::register_classes_with_charm<TimeSteppers::Heun2>();
  using component = Component<Metavariables>;

  const size_t number_of_grid_points = 5;

  const Slab slab(1.0, 3.0);
  const TimeStepId initialize_time_step_id(true, 0, slab.start());
  const TimeStepId time_step_id(true, 0, slab.start(), 1, slab.end().value());
  const auto time_step = slab.duration();

  System::variables_tag::type initial_vars(number_of_grid_points);
  get(get<Var1>(initial_vars)) = 2.0;
  get(get<Var2>(initial_vars)) = 3.0;

  imex::Tags::ImplicitHistory<Sector<Var1>>::type history1(2);
  history1.insert(initialize_time_step_id, decltype(history1)::no_value,
                  -get(get<Var1>(initial_vars)));
  imex::Tags::ImplicitHistory<Sector<Var2>>::type history2(2);
  history2.insert(initialize_time_step_id, decltype(history2)::no_value,
                  -get(get<Var2>(initial_vars)));

  ActionTesting::MockRuntimeSystem<Metavariables> runner{{}};

  ActionTesting::emplace_component_and_initialize<component>(
      &runner, 0,
      {time_step_id, time_step, initial_vars, imex::Mode::Implicit,
       std::make_unique<TimeSteppers::Heun2>(), std::move(history1),
       std::move(history2)});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  runner.next_action<component>(0);

  const auto& box = ActionTesting::get_databox<component>(runner, 0);
  const auto& final_vars = db::get<System::variables_tag>(box);

  const double dt = time_step.value();
  const double step_factor = (1.0 - 0.5 * dt) / (1.0 + 0.5 * dt);

  CHECK(get(get<Var1>(final_vars)) ==
        step_factor * get(get<Var1>(initial_vars)));
  CHECK(get(get<Var2>(final_vars)) ==
        step_factor * get(get<Var2>(initial_vars)));
}
