// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <random>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/Protocols/ImplicitSource.hpp"
#include "Evolution/Imex/Protocols/ImplicitSourceJacobian.hpp"
#include "Evolution/Imex/SolveImplicitSector.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Evolution/Imex/Tags/Jacobian.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Evolution/Imex/TestSector.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/HeunImex.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Protocols/StaticReturnApplyable.hpp"
#include "Utilities/TMPL.hpp"

namespace {
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

// [ImplicitSector]
struct Sector : tt::ConformsTo<imex::protocols::ImplicitSector> {
  using tensors = tmpl::list<Var2, Var3>;

  struct source : tt::ConformsTo<imex::protocols::ImplicitSource>,
                  tt::ConformsTo<protocols::StaticReturnApplyable> {
    using return_tags = tmpl::list<::Tags::Source<Var2>, ::Tags::Source<Var3>>;
    using argument_tags = tmpl::list<Var1, Var2, Var3, NonTensor>;

    static void apply(const gsl::not_null<tnsr::II<DataVector, 2>*> source_var2,
                      const gsl::not_null<tnsr::I<DataVector, 2>*> source_var3,
                      const Scalar<DataVector>& var1,
                      const tnsr::II<DataVector, 2>& var2,
                      const tnsr::I<DataVector, 2>& var3,
                      const double non_tensor);
  };

  struct source_jacobian
      : tt::ConformsTo<imex::protocols::ImplicitSourceJacobian>,
        tt::ConformsTo<protocols::StaticReturnApplyable> {
    using return_tags =
        tmpl::list<imex::Tags::Jacobian<Var2, ::Tags::Source<Var2>>,
                   imex::Tags::Jacobian<Var3, ::Tags::Source<Var2>>,
                   imex::Tags::Jacobian<Var3, ::Tags::Source<Var3>>>;
    using argument_tags = tmpl::list<Var1, Var3, NonTensor>;

    static void apply(
        const gsl::not_null<tnsr::iiJJ<DataVector, 2>*> dvar2_dvar2,
        const gsl::not_null<tnsr::iJJ<DataVector, 2>*> dvar2_dvar3,
        const gsl::not_null<tnsr::iJ<DataVector, 2>*> dvar3_dvar3,
        const Scalar<DataVector>& var1, const tnsr::I<DataVector, 2>& var3,
        const double non_tensor);
  };
};
// [ImplicitSector]

void Sector::source::apply(
    const gsl::not_null<tnsr::II<DataVector, 2>*> source_var2,
    const gsl::not_null<tnsr::I<DataVector, 2>*> source_var3,
    const Scalar<DataVector>& var1, const tnsr::II<DataVector, 2>& var2,
    const tnsr::I<DataVector, 2>& var3, const double non_tensor) {
  for (size_t i = 0; i < 2; ++i) {
    source_var3->get(i) = -get(var1) * var3.get(i);
    for (size_t j = 0; j <= i; ++j) {
      source_var2->get(i, j) =
          var3.get(i) * var3.get(j) - non_tensor * var2.get(i, j);
    }
  }
}

void Sector::source_jacobian::apply(
    const gsl::not_null<tnsr::iiJJ<DataVector, 2>*> dvar2_dvar2,
    const gsl::not_null<tnsr::iJJ<DataVector, 2>*> dvar2_dvar3,
    const gsl::not_null<tnsr::iJ<DataVector, 2>*> dvar3_dvar3,
    const Scalar<DataVector>& var1, const tnsr::I<DataVector, 2>& var3,
    const double non_tensor) {
  std::fill(dvar2_dvar2->begin(), dvar2_dvar2->end(), 0.0);
  std::fill(dvar2_dvar3->begin(), dvar2_dvar3->end(), 0.0);
  std::fill(dvar3_dvar3->begin(), dvar3_dvar3->end(), 0.0);
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

}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Imex.solve_implicit_sector",
                  "[Unit][Evolution]") {
  {
    Scalar<DataVector> var1{};
    get(var1) = DataVector{3.0};
    tnsr::II<DataVector, 2> var2{};
    get<0, 0>(var2) = DataVector{4.0};
    get<0, 1>(var2) = DataVector{5.0};
    get<1, 1>(var2) = DataVector{6.0};
    tnsr::I<DataVector, 2> var3{};
    get<0>(var3) = DataVector{7.0};
    get<1>(var3) = DataVector{8.0};
    const double non_tensor = 9.0;
    TestHelpers::imex::test_sector<Sector, Var1, Var2, Var3, NonTensor>(
        {std::move(var1), std::move(var2), std::move(var3), non_tensor});
  }

  // Heun first substep:
  // y(dt) = y(0) + dt/2 (d/d[y(0)] + d/dt[y(dt)])

  // For simplicity, we handle v1 entirely explicitly and v2, v3
  // entirely implicitly.
  // d/dt[v2_ij] = v3_i v3_j - nt v2_ij
  // d/dt[v3_i] = -v1 v3_i

  // Analytic solution:
  // v3_i(dt) = v3_i(0) (1 - dt/2 v1(0)) / (1 + dt/2 v1(dt))
  // v2_ij(dt) = (v2_ij(0) (1 - dt/2 nt) +
  //              + dt/2 (v3_i(0) v3_j(0) + v3_i(dt) v3_j(dt))) / (1 + dt/2 nt)

  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2, Var3>>;
  using implicit_variables_source_tag =
      Tags::Variables<tmpl::list<::Tags::Source<Var2>, ::Tags::Source<Var3>>>;
  using DtImplicitVariables =
      Variables<tmpl::list<::Tags::dt<Var2>, ::Tags::dt<Var3>>>;

  const size_t number_of_grid_points = 5;
  const auto time_step = Slab(3.0, 5.0).duration() / 3;

  MAKE_GENERATOR(gen);
  // Keep values positive to prevent the denominators in the analytic
  // solution from becoming small.
  std::uniform_real_distribution<double> dist(0.0, 5.0);
  const auto non_tensor = make_with_random_values<double>(make_not_null(&gen),
                                                          make_not_null(&dist));
  const auto initial_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), make_not_null(&dist), number_of_grid_points);
  // We overwrite these values below.  This doesn't need to be in the
  // DataBox for the test, but putting it there lets us use
  // db::mutate_apply to calculate it.
  auto source_vars =
      make_with_random_values<implicit_variables_source_tag::type>(
          make_not_null(&gen), make_not_null(&dist), number_of_grid_points);
  auto box = db::create<
      db::AddSimpleTags<variables_tag, implicit_variables_source_tag, NonTensor,
                        Tags::TimeStepper<TimeSteppers::HeunImex>,
                        Tags::TimeStep, imex::Tags::ImplicitHistory<Sector>>>(
      initial_vars, std::move(source_vars), non_tensor,
      std::make_unique<TimeSteppers::HeunImex>(), time_step,
      imex::Tags::ImplicitHistory<Sector>::type{2});
  db::mutate_apply<Sector::source>(make_not_null(&box));
  db::mutate<imex::Tags::ImplicitHistory<Sector>, Var1>(
      make_not_null(&box),
      [&dist, &gen, &time_step](
          const gsl::not_null<imex::Tags::ImplicitHistory<Sector>::type*>
              history,
          const gsl::not_null<Var1::type*> var1,
          const implicit_variables_source_tag::type& implicit_vars_source) {
        history->insert(
            TimeStepId(true, 0, time_step.slab().start()),
            implicit_vars_source
                .reference_with_different_prefixes<DtImplicitVariables>());
        fill_with_random_values(var1, make_not_null(&gen),
                                make_not_null(&dist));
      },
      db::get<implicit_variables_source_tag>(box));

  imex::solve_implicit_sector<Sector>(make_not_null(&box));

  const double dt = time_step.value();
  const auto final_vars = db::get<variables_tag>(box);
  Var3::type expected_var3(number_of_grid_points);
  tenex::evaluate<ti::I>(make_not_null(&expected_var3),
                         (1.0 - 0.5 * dt * get<Var1>(initial_vars)()) /
                             (1.0 + 0.5 * dt * get<Var1>(final_vars)()) *
                             get<Var3>(initial_vars)(ti::I));
  CHECK_ITERABLE_APPROX(get<Var3>(final_vars), expected_var3);
  Var2::type expected_var2(number_of_grid_points);
  tenex::evaluate<ti::I, ti::J>(
      make_not_null(&expected_var2),
      ((1.0 - 0.5 * dt * non_tensor) * get<Var2>(initial_vars)(ti::I, ti::J) +
       0.5 * dt *
           (get<Var3>(initial_vars)(ti::I) * get<Var3>(initial_vars)(ti::J) +
            expected_var3(ti::I) * expected_var3(ti::J))) /
          (1.0 + 0.5 * dt * non_tensor));
  CHECK_ITERABLE_APPROX(get<Var2>(final_vars), expected_var2);
}
