// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <random>
#include <type_traits>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Imex/GuessResult.hpp"
#include "Evolution/Imex/Mode.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/SolveImplicitSector.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Evolution/Imex/Tags/Jacobian.hpp"
#include "Evolution/Imex/Tags/Mode.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Evolution/Imex/TestSector.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/Heun2.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {
// Set temporarily to verify that the solver correctly skips most of
// the work when the step is explicit.
bool performing_step_with_no_implicit_term = false;

struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct Var2 : db::SimpleTag {
  using type = tnsr::II<DataVector, 2>;
};

struct Var3 : db::SimpleTag {
  using type = tnsr::I<DataVector, 2>;
};

struct NonTensor : db::SimpleTag {
  using type = double;
};

// These next several tags aren't used in the calculation, just for
// testing DataBox handling.
struct TensorFromEvolution : db::SimpleTag {
  using type = Scalar<DataVector>;
};

using VariablesFromEvolution = Tags::Variables<tmpl::list<TensorFromEvolution>>;

struct TensorTemporary : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct VariablesTemporary : db::SimpleTag {
  using type = Variables<tmpl::list<TensorTemporary>>;
};

struct SomeComputeTagBase : db::SimpleTag {
  using type = double;
};

struct SomeComputeTag : SomeComputeTagBase, db::ComputeTag {
  using base = SomeComputeTagBase;
  using argument_tags =
      tmpl::list<Var1, TensorFromEvolution, VariablesTemporary>;
  static void function(
      const gsl::not_null<double*> result, const Scalar<DataVector>& var1,
      const Scalar<DataVector>& from_evolution,
      const Variables<tmpl::list<TensorTemporary>>& temporary) {
    // Check the initialization of the temporary Variables in the
    // solver DataBox.  None of the mutators modify the object, so it
    // should always have that state.
    CHECK(temporary.number_of_grid_points() == 1);
    CHECK(get(get<TensorTemporary>(temporary))[0] == 0.0);

    // Check slicing
    CHECK(get(from_evolution).size() == 1);
    CHECK(get(from_evolution)[0] == 2.0 * get(var1)[0]);

    *result = get(var1)[0] + 1.0;
  }
};

enum class PrepId { InitialGuess, Source, Jacobian, Shared };

struct RecordPreparersForTest : db::SimpleTag {
  using type = std::array<std::pair<Var2::type, Var3::type>, 4>;
};

template <PrepId Prep>
struct Preparer {
  using return_tags = tmpl::list<RecordPreparersForTest>;
  using argument_tags = tmpl::list<Var2, Var3>;

  static void apply(
      const gsl::not_null<RecordPreparersForTest::type*> prep_run_values,
      const Var2::type& var2, const Var3::type& var3) {
    CHECK(not performing_step_with_no_implicit_term);

    std::pair current_values{var2, var3};
    CHECK((*prep_run_values)[static_cast<size_t>(Prep)] != current_values);
    (*prep_run_values)[static_cast<size_t>(Prep)] = std::move(current_values);
  }
};
// End stuff only used for DataBox handling

// [initial_guess]
struct AnalyticSolution {
  using return_tags = tmpl::list<Var2, Var3>;
  using argument_tags = tmpl::list<Var1, NonTensor>;
  static imex::GuessResult apply(
      const gsl::not_null<tnsr::II<DataVector, 2>*> var2,
      const gsl::not_null<tnsr::I<DataVector, 2>*> var3,
      const Scalar<DataVector>& var1, const double non_tensor,
      const Variables<tmpl::list<Var2, Var3>>& inhomogeneous_terms,
      const double implicit_weight) {
    // Solution for source terms
    // S[v2^ij] = v3^i v3^j - nt v2^ij
    // S[v3^i] = -v1 v3^i

    // Solving  v3^i = X - w v1 v3^i  gives  v3^i = X / (1 + w v1)
    tenex::evaluate<ti::I>(var3, get<Var3>(inhomogeneous_terms)(ti::I) /
                                     (1.0 + implicit_weight * var1()));
    tenex::evaluate<ti::I, ti::J>(
        var2, (get<Var2>(inhomogeneous_terms)(ti::I, ti::J) +
               implicit_weight * (*var3)(ti::I) * (*var3)(ti::J)) /
                  (1.0 + implicit_weight * non_tensor));
    return imex::GuessResult::ExactSolution;
  }
};
// [initial_guess]

struct InitialGuess {
  using return_tags = tmpl::list<Var2, Var3>;
  using argument_tags = tmpl::list<>;
  static imex::GuessResult apply(
      const gsl::not_null<tnsr::II<DataVector, 2>*> var2,
      const gsl::not_null<tnsr::I<DataVector, 2>*> var3,
      const Variables<tmpl::list<Var2, Var3>>& /*inhomogeneous_terms*/,
      const double /*implicit_weight*/) {
    CHECK(not performing_step_with_no_implicit_term);

    for (auto& component : *var2) {
      component *= 2.0;
    }
    for (auto& component : *var3) {
      component *= 3.0;
    }
    return imex::GuessResult::InitialGuess;
  }
};

struct Source {
  using return_tags = tmpl::list<::Tags::Source<Var2>, ::Tags::Source<Var3>>;
  using argument_tags = tmpl::list<Var1, Var2, Var3, NonTensor,
                                   RecordPreparersForTest, SomeComputeTagBase>;

  static void apply(const gsl::not_null<tnsr::II<DataVector, 2>*> source_var2,
                    const gsl::not_null<tnsr::I<DataVector, 2>*> source_var3,
                    const Scalar<DataVector>& var1,
                    const tnsr::II<DataVector, 2>& var2,
                    const tnsr::I<DataVector, 2>& var3, const double non_tensor,
                    const RecordPreparersForTest::type& prep_run_values,
                    const double compute_tag_value) {
    CHECK(not performing_step_with_no_implicit_term);

    const std::pair current_values{var2, var3};
    CHECK(prep_run_values[static_cast<size_t>(PrepId::Shared)] ==
          current_values);
    CHECK(prep_run_values[static_cast<size_t>(PrepId::Source)] ==
          current_values);

    CHECK(compute_tag_value == get(var1)[0] + 1.0);

    work(source_var2, source_var3, var1, var2, var3, non_tensor);
  }

  // Used in the test below.  Not part of the IMEX interface.
  static void work(const gsl::not_null<tnsr::II<DataVector, 2>*> source_var2,
                   const gsl::not_null<tnsr::I<DataVector, 2>*> source_var3,
                   const Scalar<DataVector>& var1,
                   const tnsr::II<DataVector, 2>& var2,
                   const tnsr::I<DataVector, 2>& var3,
                   const double non_tensor) {
    tenex::evaluate<ti::I, ti::J>(
        source_var2,
        var3(ti::I) * var3(ti::J) - non_tensor * var2(ti::I, ti::J));
    tenex::evaluate<ti::I>(source_var3, -var1() * var3(ti::I));
  }
};

struct Jacobian {
  using return_tags =
      tmpl::list<imex::Tags::Jacobian<Var2, ::Tags::Source<Var2>>,
                 imex::Tags::Jacobian<Var3, ::Tags::Source<Var2>>,
                 imex::Tags::Jacobian<Var3, ::Tags::Source<Var3>>>;
  using argument_tags =
      tmpl::list<Var1, Var3, NonTensor, RecordPreparersForTest>;

  static void apply(const gsl::not_null<tnsr::iiJJ<DataVector, 2>*> dvar2_dvar2,
                    const gsl::not_null<tnsr::iJJ<DataVector, 2>*> dvar2_dvar3,
                    const gsl::not_null<tnsr::iJ<DataVector, 2>*> dvar3_dvar3,
                    const Scalar<DataVector>& var1,
                    const tnsr::I<DataVector, 2>& var3, const double non_tensor,
                    const RecordPreparersForTest::type& prep_run_values) {
    CHECK(not performing_step_with_no_implicit_term);

    // We don't need var2 for anything else in this function.  Hard to
    // imagine a way not taking one of the variables as an argument
    // could break anything, but easy to test so we do it.
    CHECK(prep_run_values[static_cast<size_t>(PrepId::Shared)].second == var3);
    CHECK(prep_run_values[static_cast<size_t>(PrepId::Jacobian)].second ==
          var3);

    std::fill(dvar2_dvar3->begin(), dvar2_dvar3->end(), 0.0);
    for (size_t i = 0; i < 2; ++i) {
      dvar2_dvar2->get(i, i, i, i) = -non_tensor;
      dvar2_dvar3->get(i, i, i) = 2.0 * var3.get(i);
      dvar3_dvar3->get(i, i) = -get(var1);
      for (size_t j = 0; j < i; ++j) {
        dvar2_dvar2->get(i, j, i, j) = -non_tensor;
        dvar2_dvar3->get(i, i, j) += var3.get(j);
        dvar2_dvar3->get(j, i, j) += var3.get(i);
      }
    }
  }
};

// [ImplicitSector]
template <bool TestWithAnalyticSolution>
struct ImplicitSector : tt::ConformsTo<imex::protocols::ImplicitSector> {
  using tensors = tmpl::list<Var2, Var3>;

  using tags_from_evolution =
      tmpl::list<Var1, NonTensor, VariablesFromEvolution>;
  using simple_tags = tmpl::list<RecordPreparersForTest, VariablesTemporary>;
  using compute_tags = tmpl::list<SomeComputeTag>;

  using initial_guess_prep = tmpl::list<Preparer<PrepId::InitialGuess>>;
  using source_prep =
      tmpl::list<Preparer<PrepId::Shared>, Preparer<PrepId::Source>>;
  using jacobian_prep =
      tmpl::list<Preparer<PrepId::Shared>, Preparer<PrepId::Jacobian>>;

  using initial_guess = tmpl::conditional_t<TestWithAnalyticSolution,
                                            AnalyticSolution, InitialGuess>;

  using source = Source;
  using jacobian =
      tmpl::conditional_t<TestWithAnalyticSolution,
                          imex::NoJacobianBecauseSolutionIsAnalytic, Jacobian>;
};
// [ImplicitSector]

// ::tensors doesn't depend on the template parameter
using sector_variables_tag = Tags::Variables<ImplicitSector<false>::tensors>;
using SectorVariables = sector_variables_tag::type;

tuples::TaggedTuple<sector_variables_tag, Var1, NonTensor,
                    VariablesFromEvolution>
arbitrary_test_values() {
  SectorVariables explicit_values(1);
  tnsr::II<DataVector, 2>& var2 = get<Var2>(explicit_values);
  get<0, 0>(var2) = 3.0;
  get<0, 1>(var2) = 4.0;
  get<1, 1>(var2) = 5.0;
  tnsr::I<DataVector, 2>& var3 = get<Var3>(explicit_values);
  get<0>(var3) = 6.0;
  get<1>(var3) = 7.0;

  Scalar<DataVector> var1{};
  get(var1) = DataVector{8.0};
  const double non_tensor = 9.0;
  Variables<tmpl::list<TensorFromEvolution>> test_variables(1);
  get<TensorFromEvolution>(test_variables)[0] = 2.0 * get(var1)[0];
  return {std::move(explicit_values), std::move(var1), non_tensor,
          std::move(test_variables)};
}

template <bool TestWithAnalyticSolution>
void test_test_sector() {
  using sector = ImplicitSector<TestWithAnalyticSolution>;
  auto values = arbitrary_test_values();
  TestHelpers::imex::test_sector<sector>(
      1.0e-1, 1.0e-12, std::move(get<sector_variables_tag>(values)),
      {std::move(get<Var1>(values)), get<NonTensor>(values),
       std::move(get<VariablesFromEvolution>(values))});
}

void test_internal_jacobian_ordering() {
  // This test doesn't make sense on the analytic solution version.
  using sector = ImplicitSector<false>;

  const Slab slab(0.0, 1.0);
  TimeSteppers::History<SectorVariables> history(2);
  history.insert(TimeStepId(true, 0, slab.start()), decltype(history)::no_value,
                 db::prefix_variables<Tags::dt, SectorVariables>(1, 3.0));

  auto values = arbitrary_test_values();
  auto evolution_box = db::create<db::AddSimpleTags<
      Tags::TimeStepper<TimeSteppers::Heun2>, Tags::TimeStep,
      imex::Tags::ImplicitHistory<sector>, sector_variables_tag, Var1,
      NonTensor, VariablesFromEvolution>>(
      std::make_unique<TimeSteppers::Heun2>(), slab.duration(),
      std::move(history), std::move(get<sector_variables_tag>(values)),
      std::move(get<Var1>(values)), get<NonTensor>(values),
      std::move(get<VariablesFromEvolution>(values)));

  imex::solve_implicit_sector_detail::ImplicitSolver<
      sector, std::decay_t<decltype(evolution_box)>>
      solver(evolution_box);
  solver.set_index(
      make_not_null(
          &db::get_mutable_reference<imex::Tags::ImplicitHistory<sector>>(
              make_not_null(&evolution_box))),
      0);
  solver.compute_initial_guess();
  const auto initial_guess = solver.initial_guess();
  const auto jacobian = solver.jacobian(initial_guess);

  auto deriv_approx = Approx::custom().epsilon(1.0e-12);
  // gsl_multiroot wants jacobian[i][j] = dfi/dxj
  for (size_t j = 0; j < initial_guess.size(); ++j) {
    const auto derivative =
        numerical_derivative(solver, initial_guess, j, 1.0e-1);
    for (size_t i = 0; i < initial_guess.size(); ++i) {
      CHECK(jacobian[i][j] == deriv_approx(derivative[i]));
    }
  }
}

template <bool TestWithAnalyticSolution>
void test_solve_implicit_sector(const imex::Mode solve_mode) {
  using sector = ImplicitSector<TestWithAnalyticSolution>;
  // No solve is done with an analytic solution.
  const bool doing_semi_implicit_solve =
      solve_mode == imex::Mode::SemiImplicit and not TestWithAnalyticSolution;
  // We handle v1 entirely explicitly and v2, v3 entirely implicitly.
  // The evolution equations for the latter two (coded in `Source`
  // above) are
  // d/dt[v2^ij] = v3^i v3^j - nt v2^ij
  // d/dt[v3^i] = -v1 v3^i

  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2, Var3>>;
  using implicit_variables_source_tag =
      Tags::Variables<tmpl::list<::Tags::Source<Var2>, ::Tags::Source<Var3>>>;
  using DtImplicitVariables =
      Variables<tmpl::list<::Tags::dt<Var2>, ::Tags::dt<Var3>>>;
  using history_tag = imex::Tags::ImplicitHistory<sector>;

  const size_t number_of_grid_points = 5;
  const Slab slab(3.0, 5.0);
  const TimeStepId initial_time_step_id(true, 0, slab.start());
  const auto time_step = Slab(3.0, 5.0).duration() / 3;

  MAKE_GENERATOR(gen);
  // Keep values positive to prevent the denominators in the analytic
  // solution from becoming small.
  std::uniform_real_distribution<double> dist(0.0, 5.0);
  const auto non_tensor = make_with_random_values<double>(make_not_null(&gen),
                                                          make_not_null(&dist));
  const auto initial_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), make_not_null(&dist), number_of_grid_points);
  auto box = db::create<
      db::AddSimpleTags<variables_tag, NonTensor, VariablesFromEvolution,
                        Tags::TimeStepper<TimeSteppers::Heun2>, Tags::TimeStep,
                        history_tag, imex::Tags::Mode>>(
      initial_vars, non_tensor, VariablesFromEvolution::type{},
      std::make_unique<TimeSteppers::Heun2>(), time_step,
      typename history_tag::type{2}, solve_mode);

  // Perform updates as if taking an explicit step.
  const auto simulate_explicit_step = [&dist, &gen, &initial_vars](
                                          const auto box,
                                          const TimeStepId& time_step_id) {
    db::mutate<history_tag, Var1, Var2, Var3, VariablesFromEvolution>(
        box,
        [&dist, &gen, &initial_vars, &time_step_id](
            const gsl::not_null<typename history_tag::type*> history,
            const gsl::not_null<Var1::type*> var1,
            const gsl::not_null<Var2::type*> var2,
            const gsl::not_null<Var3::type*> var3,
            const gsl::not_null<VariablesFromEvolution::type*> test_variables,
            const NonTensor::type& non_tensor) {
          implicit_variables_source_tag::type source_vars(number_of_grid_points,
                                                          0.0);
          Source::work(&get<Tags::Source<Var2>>(source_vars),
                       &get<Tags::Source<Var3>>(source_vars), *var1, *var2,
                       *var3, non_tensor);

          history->insert(
              time_step_id, history_tag::type::no_value,
              source_vars
                  .reference_with_different_prefixes<DtImplicitVariables>());
          // Update the explicitly evolved variable.
          fill_with_random_values(var1, make_not_null(&gen),
                                  make_not_null(&dist));
          // The explicit time derivative for var2 and var3 is
          // zero, so the explicit integration will consider them
          // constant and reset them to the initial value.
          *var2 = get<Var2>(initial_vars);
          *var3 = get<Var3>(initial_vars);
          // This isn't evolved but we test obtaining it from the
          // evolution box by checking for this value.
          test_variables->initialize(get(*var1).size());
          get(get<TensorFromEvolution>(*test_variables)) = 2.0 * get(*var1);
        },
        db::get<NonTensor>(*box));
  };

  simulate_explicit_step(make_not_null(&box), initial_time_step_id);

  auto guess = initial_vars;
  InitialGuess::apply(&get<Var2>(guess), &get<Var3>(guess), {}, {});

  imex::solve_implicit_sector<sector>(make_not_null(&box));

  const double dt = time_step.value();
  const auto final_vars = db::get<variables_tag>(box);
  Var2::type expected_var2{};
  Var3::type expected_var3{};
  Var2::type expected_var2_final{};
  if (not doing_semi_implicit_solve) {
    // The first implicit substep for the Heun stepper is
    // y(dt) = y(0) + dt/2 (d/d[y(0)] + d/dt[y(dt)])

    // The analytic solution for the result of the first substep is
    // v3^i(dt) = v3^i(0) (1 - dt/2 v1(0)) / (1 + dt/2 v1(dt))
    // v2^ij(dt) = (v2^ij(0) (1 - dt/2 nt) +
    //              + dt/2 (v3^i(0) v3^j(0) + v3^i(dt) v3^j(dt)))
    //             / (1 + dt/2 nt)

    tenex::evaluate<ti::I>(make_not_null(&expected_var3),
                           (1.0 - 0.5 * dt * get<Var1>(initial_vars)()) /
                               (1.0 + 0.5 * dt * get<Var1>(final_vars)()) *
                               get<Var3>(initial_vars)(ti::I));
    tenex::evaluate<ti::I, ti::J>(
        make_not_null(&expected_var2),
        ((1.0 - 0.5 * dt * non_tensor) * get<Var2>(initial_vars)(ti::I, ti::J) +
         0.5 * dt *
             (get<Var3>(initial_vars)(ti::I) * get<Var3>(initial_vars)(ti::J) +
              expected_var3(ti::I) * expected_var3(ti::J))) /
            (1.0 + 0.5 * dt * non_tensor));

    // The second implicit substep is simpler, since it isn't actually
    // implicit, and is in fact the same as the first substep.
    expected_var2_final = expected_var2;
  } else {
    // For a semi-implicit method, we expand the source term around
    // the initial guess:
    //
    // S(y(dt)) = S(guess + (y(dt) - guess))
    //   ~= S(guess) + J(guess) (y(dt) - guess)
    //
    // For Heun's method
    //
    // y(dt) = y(0) + dt/2 (d/d[y(0)] + d/dt[y(dt)])
    //
    // this gives the equation
    //
    // [1 - dt/2 J(guess)] (y(dt) - guess)
    //     = y(0) - guess + dt/2 (d/d[y(0)] + S(guess))
    //
    // with
    //
    // S[v2^ij] = v3^i v3^j - nt v2^ij
    // S[v3^i] = -v1 v3^i
    // J[v2^ij] : dS[v2^ij]/dv2^kl = -nt delta^i_k delta^j_l
    //            dS[v2^ij]/dv3^k = v3^i delta^j_k + v3^j delta^i_k
    // J[v3^i] : dS[v3^i]/dv2^jk = 0
    //           dS[v3^i]/dv3^j = -v1 delta^i_j
    //
    // This is fairly easy to solve as (denoting the guess as g2 and g3):
    //
    // v3^i(dt) = (1 - dt/2 v1(0)) / (1 + dt/2 v1(dt)) v3^i(0)
    //     (solved exactly as the source is linear)
    //
    // v2^ij(dt) =
    //    (1 - dt/2 nt) / (1 + dt/2 nt) v2^ij(0)
    //    + dt/2 / (1 + dt/2 nt) [
    //         v3^i(0) v3^j(0) + g3^i v3^j(dt) + v3^i(dt) g3^j - g3^i g3^j]

    tenex::evaluate<ti::I>(make_not_null(&expected_var3),
                           (1.0 - 0.5 * dt * get<Var1>(initial_vars)()) /
                               (1.0 + 0.5 * dt * get<Var1>(final_vars)()) *
                               get<Var3>(initial_vars)(ti::I));
    tenex::evaluate<ti::I, ti::J>(
        make_not_null(&expected_var2),
        ((1.0 - 0.5 * dt * non_tensor) * get<Var2>(initial_vars)(ti::I, ti::J) +
         0.5 * dt *
             (get<Var3>(initial_vars)(ti::I) * get<Var3>(initial_vars)(ti::J) +
              get<Var3>(guess)(ti::I) * expected_var3(ti::J) +
              expected_var3(ti::I) * get<Var3>(guess)(ti::J) -
              get<Var3>(guess)(ti::I) * get<Var3>(guess)(ti::J))) /
            (1.0 + 0.5 * dt * non_tensor));

    // The second substep is explicit, so doesn't do a semi-implicit
    // solve, but still uses the results from the first substep.  Var3
    // is again exact.
    tenex::evaluate<ti::I, ti::J>(
        make_not_null(&expected_var2_final),
        (1.0 - 0.5 * dt * non_tensor) * get<Var2>(initial_vars)(ti::I, ti::J) +
            0.5 * dt *
                (get<Var3>(initial_vars)(ti::I) *
                     get<Var3>(initial_vars)(ti::J) +
                 expected_var3(ti::I) * expected_var3(ti::J) -
                 non_tensor * expected_var2(ti::I, ti::J)));
  }
  CHECK_ITERABLE_APPROX(get<Var2>(final_vars), expected_var2);
  CHECK_ITERABLE_APPROX(get<Var3>(final_vars), expected_var3);

  CHECK(db::get<history_tag>(box).size() == 1);
  CHECK(db::get<history_tag>(box).substeps().empty());

  simulate_explicit_step(make_not_null(&box),
                         initial_time_step_id.next_substep(time_step, 1.0));
  performing_step_with_no_implicit_term = true;
  imex::solve_implicit_sector<sector>(make_not_null(&box));
  performing_step_with_no_implicit_term = false;
  CHECK_ITERABLE_APPROX(get<Var2>(db::get<variables_tag>(box)),
                        expected_var2_final);
  CHECK_ITERABLE_APPROX(get<Var3>(db::get<variables_tag>(box)), expected_var3);

  CHECK(db::get<history_tag>(box).size() == 1);
  CHECK(db::get<history_tag>(box).substeps().size() == 1);

  // Take another substep just to test the history cleanup.
  simulate_explicit_step(make_not_null(&box),
                         initial_time_step_id.next_step(time_step));
  imex::solve_implicit_sector<sector>(make_not_null(&box));

  CHECK(db::get<history_tag>(box).size() == 1);
  CHECK(db::get<history_tag>(box).substeps().empty());
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Imex.solve_implicit_sector",
                  "[Unit][Evolution]") {
  test_test_sector<false>();
  test_test_sector<true>();
  test_internal_jacobian_ordering();
  test_solve_implicit_sector<false>(imex::Mode::Implicit);
  test_solve_implicit_sector<true>(imex::Mode::Implicit);
  test_solve_implicit_sector<false>(imex::Mode::SemiImplicit);
  test_solve_implicit_sector<true>(imex::Mode::SemiImplicit);
}
