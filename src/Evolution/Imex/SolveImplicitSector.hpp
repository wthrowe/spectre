// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/ExtractPoint.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tags.hpp"
#include "DataStructures/Tensor/Metafunctions.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/Mode.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "NumericalAlgorithms/LinearSolver/Lapack.hpp"
#include "NumericalAlgorithms/RootFinding/GslMultiRoot.hpp"
#include "Time/History.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/DoNotMove.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
namespace imex::Tags {
struct Mode;
}  // namespace imex::Tags
/// \endcond

namespace imex {
namespace solve_implicit_sector_detail {
template <typename ImplicitSector, typename DbTags>
class ImplicitSolve {
  static_assert(tt::conforms_to_v<ImplicitSector, protocols::ImplicitSector>);

  //FIXME add helper to protocol, rename, make required or void, test
  // verify doesn't modify input.
  template <typename LocalSector, typename = std::void_t<>>
  struct get_helper {
    struct type {
      using return_tags = tmpl::list<>;
      using argument_tags = tmpl::list<>;
      static void apply() {}
    };
  };

  template <typename LocalSector>
  struct get_helper<LocalSector, std::void_t<typename LocalSector::helper>> {
    using type = tmpl::conditional_t<std::is_same_v<typename LocalSector::helper, void>, typename get_helper<void>::type, typename LocalSector::helper>;
  };

  using helper = typename get_helper<ImplicitSector>::type;

  using non_sector_argument_list = tmpl::list_difference<
      tmpl::remove_duplicates<tmpl::append<
          typename ImplicitSector::source::argument_tags,
          typename ImplicitSector::source_jacobian::argument_tags,
          typename helper::argument_tags>>,
      typename ImplicitSector::tensors>;

  using non_sector_tensor_argument_list =
      tmpl::filter<non_sector_argument_list,
                   tt::is_a<Tensor, tmpl::bind<tmpl::type_from, tmpl::_1>>>;

  // The tricky part of supporting this would be handling requests for
  // the evolved variables themselves, because the implicit solve is
  // updating the values of some of the tensors.
  static_assert(
      tmpl::none<
          non_sector_argument_list,
          tt::is_a<Variables, tmpl::bind<tmpl::type_from, tmpl::_1>>>::value,
      "Slicing Variables for implicit sector sources is not implemented.  "
      "Take the individual Tensors as arguments instead.");

 public:
  using Vars = Variables<typename ImplicitSector::tensors>;
  using DtVars = db::prefix_variables<::Tags::dt, Vars>;
  using SourceVars = db::prefix_variables<::Tags::Source, Vars>;
  using History = TimeSteppers::History<DtVars>;

  static constexpr size_t dimension = Vars::number_of_independent_components;

  ImplicitSolve(const gsl::not_null<History*> history,
                const db::DataBox<DbTags>& box, const Vars& explicit_value,
                const size_t index)
      : history_(history),
        box_(&box),
        explicit_value_(&explicit_value),
        index_(index),
        implicit_weight_(db::get<::Tags::TimeStepper<>>(*box_).implicit_weight(
            *history_, db::get<::Tags::TimeStep>(*box_))),
        non_sector_tensor_arguments_(1) {
    tmpl::for_each<non_sector_tensor_argument_list>([&](auto tag) {
      using Tag = tmpl::type_from<decltype(tag)>;
      overwrite_point(make_not_null(&get<Tag>(non_sector_tensor_arguments_)),
                      extract_point(get<Tag>(*box_), index_), 0);
    });
  }

  std::array<double, dimension> operator()(
      std::array<double, dimension> value) const {
    Vars value_variables(value.data(), value.size());
    call_helper(value_variables, typename helper::return_tags{},
                typename helper::argument_tags{});
    const auto source = call_source(
        value_variables, typename ImplicitSector::source::return_tags{},
        typename ImplicitSector::source::argument_tags{});
    value_variables = *explicit_value_ - value_variables;
    db::get<::Tags::TimeStepper<>>(*box_).update_u_implicit(
        make_not_null(&value_variables), history_,
        source.template reference_with_different_prefixes<DtVars>(),
        db::get<::Tags::TimeStep>(*box_));
    return value;
  }

  std::array<std::array<double, dimension>, dimension> jacobian(
      std::array<double, dimension> value) const {
    const Vars value_variables(value.data(), value.size());
    call_helper(value_variables, typename helper::return_tags{},
                typename helper::argument_tags{});
    auto result =
        implicit_weight_ *
        call_jacobian(
            value_variables,
            typename ImplicitSector::source_jacobian::return_tags{},
            typename ImplicitSector::source_jacobian::argument_tags{});
    for (size_t i = 0; i < dimension; ++i) {
      result[i][i] -= 1.0;
    }
    return result;
  }

 private:
  template <typename Tag>
  const auto& from_vars_or_box(const Vars& vars) const {
    if constexpr (tmpl::list_contains_v<typename Vars::tags_list, Tag>) {
      return get<Tag>(vars);
    } else if constexpr (tmpl::list_contains_v<decltype(helper_results_), Tag>) {
      //FIXME rename func or something?
      return get<Tag>(helper_results_);
    } else if constexpr (tt::is_a_v<Tensor, typename Tag::type>) {
      return get<Tag>(non_sector_tensor_arguments_);
    } else {
      return get<Tag>(*box_);
    }
  }

  template <typename... ReturnTags, typename... ArgumentTags>
  void call_helper(const Vars& implicit_vars,
                   tmpl::list<ReturnTags...> /*meta*/,
                   tmpl::list<ArgumentTags...> /*meta*/) const {
    //FIXME avoid recomputation
    //FIXME copy from box?  Reuse old value?
    helper_results_ =
        make_with_value<decltype(helper_results_)>(size_t{1}, 0.0);
    helper::apply(make_not_null(&get<ReturnTags>(helper_results_))...,
                  from_vars_or_box<ArgumentTags>(implicit_vars)...);
  }

  template <typename... ReturnTags, typename... ArgumentTags>
  SourceVars call_source(const Vars& implicit_vars,
                         tmpl::list<ReturnTags...> /*meta*/,
                         tmpl::list<ArgumentTags...> /*meta*/) const {
    SourceVars result(1);
    ImplicitSector::source::apply(
        make_not_null(&get<ReturnTags>(result))...,
        from_vars_or_box<ArgumentTags>(implicit_vars)...);
    return result;
  }

  template <typename ImplicitTags, typename... ReturnTags,
            typename... ArgumentTags>
  std::array<std::array<double, dimension>, dimension> call_jacobian(
      const Variables<ImplicitTags>& implicit_vars,
      tmpl::list<ReturnTags...> /*meta*/,
      tmpl::list<ArgumentTags...> /*meta*/) const {
    Variables<tmpl::list<ReturnTags...>> result_variables(1);
    ImplicitSector::source_jacobian::apply(
        make_not_null(&get<ReturnTags>(result_variables))...,
        from_vars_or_box<ArgumentTags>(implicit_vars)...);

    std::array<std::array<double, dimension>, dimension> result{};
    // The storage order for the tensors does not match the order for
    // the returned array, so we have to copy components individually.
    const auto copy_result = [&](auto return_tag) {
      using ReturnTag = decltype(return_tag);
      using dependent_tag = typename ReturnTag::dependent;
      using independent_tag = typename ReturnTag::independent;
      // Despite repeated references to it, the result of this is
      // independent of implicit_vars.  The variable is just used for
      // calculating offsets.
      const auto& dependent_for_offsets =
          get<typename dependent_tag::tag>(implicit_vars);
      const auto& independent_for_offsets = get<independent_tag>(implicit_vars);
      for (size_t dependent_component = 0;
           dependent_component < dependent_for_offsets.size();
           ++dependent_component) {
        const auto dependent_index =
            dependent_for_offsets.get_tensor_index(dependent_component);
        auto& result_row = result[static_cast<size_t>(
            dependent_for_offsets[dependent_component].data() -
            implicit_vars.data())];
        for (size_t independent_component = 0;
             independent_component < independent_for_offsets.size();
             ++independent_component) {
          const auto independent_index =
              independent_for_offsets.get_tensor_index(independent_component);
          result_row[static_cast<size_t>(
              independent_for_offsets[independent_component].data() -
              implicit_vars.data())] =
              get<ReturnTag>(result_variables)
                  .get(concatenate(dependent_index, independent_index))[0];
        }
      }
      return 0;
    };
    expand_pack(copy_result(ReturnTags{})...);
    return result;
  }

  gsl::not_null<History*> history_;
  gsl::not_null<const db::DataBox<DbTags>*> box_;
  gsl::not_null<const Vars*> explicit_value_;
  size_t index_;

  double implicit_weight_;
  Variables<non_sector_tensor_argument_list> non_sector_tensor_arguments_;
  //FIXME more efficient storage, mutable yuck
  mutable tuples::tagged_tuple_from_typelist<typename helper::return_tags>
      helper_results_;
};
}  // namespace solve_implicit_sector_detail

//FIXME implicit_weight = 0 should be handled specially somewhere.

//FIXME where should the new history value be calculated?  Probably
//not here, so primitives are guaranteed to be consistent.

/// Perform the implicit solve for one implicit sector.
///
/// This will update the tensors in the implicit sector and clean up
/// the corresponding time stepper history.  A new history entry is
/// not added, because that should be done with the same values of the
/// variables used for the explicit portion of the time derivative,
/// which may still undergo variable-fixing-like corrections. FIXME or
/// has that already happened?
template <typename ImplicitSector, typename DbTags>
void solve_implicit_sector(const gsl::not_null<db::DataBox<DbTags>*> box) {
  using ImplicitSolve =
      solve_implicit_sector_detail::ImplicitSolve<ImplicitSector, DbTags>;
  using ImplicitVars = Variables<typename ImplicitSector::tensors>;

  // Copy the variables instead of updating them in place to avoid
  // reevaluating compute items after each point.  Compute items
  // depending on the implicit variables won't work anyway, because
  // they will not be updated after each iteration of the root find.
  auto vars = db::apply<typename ImplicitVars::tags_list>(
      [](const auto&... tensors) {
        return tmpl::as_pack<typename ImplicitVars::tags_list>(
            [&](auto... tensor_tags) {
              Variables<typename ImplicitVars::tags_list> result(
                  get_first_argument(tensors...)[0].size());
              expand_pack((get<tmpl::type_from<decltype(tensor_tags)>>(result) =
                               tensors)...);
              return result;
            });
      },
      *box);
  auto& implicit_history =
      db::get_mutable_reference<Tags::ImplicitHistory<ImplicitSector>>(box);

  Matrix semi_implicit_jacobian;
  typename ImplicitSolve::History pointwise_history(
      implicit_history.integration_order());
  for (size_t point = 0; point < vars.number_of_grid_points(); ++point) {
    pointwise_history.mark_unneeded(pointwise_history.end());
    for (auto volume_derivative = implicit_history.derivatives_begin();
         volume_derivative != implicit_history.derivatives_end();
         ++volume_derivative) {
      pointwise_history.insert(volume_derivative.time_step_id(),
                               extract_point(*volume_derivative, point));
    }
    std::array<double, ImplicitVars::number_of_independent_components>
        pointwise_array;
    ImplicitVars pointwise_vars(pointwise_array.data(), pointwise_array.size());
    pointwise_vars = do_not_move(extract_point(vars, point));

    switch (db::get<Tags::Mode>(*box)) {
      case Mode::Implicit:
        {
          const double tolerance = 1.0e-10;  // FIXME
          const size_t max_iterations = 100;
          pointwise_array = RootFinder::gsl_multiroot(
              ImplicitSolve(&pointwise_history, *box, pointwise_vars, point),
              pointwise_array,
              RootFinder::StoppingConditions::Residual(tolerance),
              max_iterations);
          break;
        }
      case Mode::SemiImplicit:
        {
          const ImplicitSolve solve(&pointwise_history, *box, pointwise_vars,
                                    point);
          std::array<double, ImplicitVars::number_of_independent_components>
              implicit_correction_array = solve(pointwise_array);
          DataVector implicit_correction_dv(implicit_correction_array.data(),
                                            implicit_correction_array.size());
          implicit_correction_dv *= -1.0;
          semi_implicit_jacobian = solve.jacobian(pointwise_array);
          //FIXME check for errors
          lapack::general_matrix_linear_solve(&implicit_correction_dv,
                                              &semi_implicit_jacobian);
          pointwise_array += implicit_correction_array;
          break;
        }
      case Mode::Explicit:
        {
          // FIXME don't do pointwise
          db::get<::Tags::TimeStepper<>>(*box).update_u(
              make_not_null(&pointwise_vars), &pointwise_history,
              db::get<::Tags::TimeStep>(*box));
          break;
        }
      default:
        ERROR("Invalid mode");
    }

    overwrite_point(make_not_null(&vars), pointwise_vars, point);
  }

  db::mutate_apply<typename ImplicitVars::tags_list, tmpl::list<>>(
      [&](const auto... tensors) {
        tmpl::as_pack<typename ImplicitVars::tags_list>(
            [&](auto... tensor_tags) {
              expand_pack(
                  (*tensors =
                       get<tmpl::type_from<decltype(tensor_tags)>>(vars))...);
            });
      },
      box);

  // Copy the cleanup to the history for the full element.
  implicit_history.mark_unneeded(implicit_history.end() -
                                 static_cast<int>(pointwise_history.size()));
}
}  // namespace imex
