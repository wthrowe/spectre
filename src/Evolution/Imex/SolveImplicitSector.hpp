// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/ExtractPoint.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Imex/GuessResult.hpp"
#include "Evolution/Imex/Mode.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Evolution/Imex/Tags/Jacobian.hpp"
#include "Evolution/Imex/Tags/Mode.hpp"
#include "NumericalAlgorithms/LinearSolver/Lapack.hpp"
#include "NumericalAlgorithms/RootFinding/GslMultiRoot.hpp"
#include "Time/History.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/ImexTimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace imex {
namespace solve_implicit_sector_detail {
template <typename ImplicitSector, typename EvolutionBox>
class ImplicitSolver {
  static_assert(
      tt::assert_conforms_to_v<ImplicitSector, protocols::ImplicitSector>);

  using system_variables_tag =
      db::creation_tag<tmpl::front<typename ImplicitSector::tensors>,
                       EvolutionBox>;

  // Retrieving compute tags based on the system variables will cause
  // full recomputation in the volume after the solve of each point.
  // Retrieving the system variables themselves, is fine, of course.
  static_assert(
      tmpl::all<
          typename ImplicitSector::tags_from_evolution,
          tmpl::or_<
              tmpl::bind<tmpl::list_contains,
                         tmpl::pin<typename system_variables_tag::tags_list>,
                         tmpl::_1>,
              tmpl::not_<tmpl::bind<db::tag_depends_on,
                                    tmpl::pin<system_variables_tag>, tmpl::_1,
                                    tmpl::pin<EvolutionBox>>>>>::value,
      "Do not include tags computed from the system variables in "
      "`tags_from_evolution`.  This would result in unnecessary "
      "computations.  Instead, include the tag in the sector's "
      "`compute_tags`.");

  using sector_variables_tag =
      ::Tags::Variables<typename ImplicitSector::tensors>;
  using SectorVariables = typename sector_variables_tag::type;
  static constexpr size_t solve_dimension =
      SectorVariables::number_of_independent_components;

  struct EvolutionBoxTag : db::SimpleTag {
    using type = const EvolutionBox*;
  };

  struct SolverPointIndex : db::SimpleTag {
    using type = size_t;
  };

  template <typename Tag, typename = typename Tag::type>
  struct FromEvolution : Tag, db::ReferenceTag {
    using base = Tag;
    using parent_tag = EvolutionBoxTag;
    using argument_tags = tmpl::list<parent_tag>;
    static const typename base::type& get(
        const EvolutionBox* const evolution_box) {
      return db::get<Tag>(*evolution_box);
    }
  };

  template <typename Tag, typename VariablesTags>
  struct FromEvolution<Tag, Variables<VariablesTags>> : Tag, db::ComputeTag {
    using base = Tag;
    using argument_tags = tmpl::list<EvolutionBoxTag, SolverPointIndex>;
    static constexpr auto function(
        const gsl::not_null<typename base::type*> result,
        const EvolutionBox* const evolution_box, const size_t index) {
      result->initialize(1);
      extract_point(result, db::get<Tag>(*evolution_box), index);
    }
  };

  // Tensor<DataVector> always allocates, so instead of creating one
  // we create a one-tensor Variables, and the DataBox will allow
  // access to the tensor transparently.  We could instead manually
  // set all the DataVectors as non-owning, pointing at individual
  // doubles, but that's more work and it's not clear it gains us
  // anything.
  template <typename Tag, typename Symm, typename IndexList>
  struct FromEvolution<Tag, Tensor<DataVector, Symm, IndexList>>
      : ::Tags::Variables<tmpl::list<Tag>>, db::ComputeTag {
    using base = ::Tags::Variables<tmpl::list<Tag>>;
    using argument_tags = tmpl::list<EvolutionBoxTag, SolverPointIndex>;
    static constexpr auto function(
        const gsl::not_null<typename base::type*> result,
        const EvolutionBox* const evolution_box, const size_t index) {
      result->initialize(1);
      extract_point(make_not_null(&get<Tag>(*result)),
                    db::get<Tag>(*evolution_box), index);
    }
  };

  struct ExplicitValue : db::SimpleTag {
    using type = SectorVariables;
  };

  struct ExplicitValueCompute : ExplicitValue, db::ComputeTag {
    using base = ExplicitValue;
    using argument_tags = tmpl::list<EvolutionBoxTag, SolverPointIndex>;
    static constexpr auto function(
        const gsl::not_null<typename base::type*> result,
        const EvolutionBox* const evolution_box, const size_t index) {
      result->initialize(1);
      tmpl::for_each<typename ImplicitSector::tensors>([&](auto tag) {
        using Tag = tmpl::type_from<decltype(tag)>;
        extract_point(make_not_null(&get<Tag>(*result)),
                      db::get<Tag>(*evolution_box), index);
      });
    }
  };

  using all_mutators = tmpl::remove_duplicates<
      tmpl::append<tmpl::list<typename ImplicitSector::source,
                              typename ImplicitSector::jacobian>,
                   typename ImplicitSector::source_prep,
                   typename ImplicitSector::jacobian_prep,
                   typename ImplicitSector::initial_guess_prep>>;

  using source_tag = db::add_tag_prefix<::Tags::Source, sector_variables_tag>;
  using jacobian_tag =
      ::Tags::Variables<jacobian_tags<typename ImplicitSector::tensors,
                                      typename source_tag::type::tags_list>>;

  using internal_simple_tags =
      tmpl::list<EvolutionBoxTag, SolverPointIndex, sector_variables_tag,
                 source_tag, jacobian_tag>;
  using internal_compute_tags = tmpl::list<ExplicitValueCompute>;
  using tags_from_evolution =
      tmpl::transform<typename ImplicitSector::tags_from_evolution,
                      tmpl::bind<FromEvolution, tmpl::_1>>;

  using simple_tags =
      tmpl::append<internal_simple_tags, typename ImplicitSector::simple_tags>;
  using compute_tags = tmpl::append<internal_compute_tags, tags_from_evolution,
                                    typename ImplicitSector::compute_tags>;

  using SolveBox =
      db::compute_databox_type<tmpl::append<simple_tags, compute_tags>>;

 public:
  ImplicitSolver(const EvolutionBox& evolution_box)
      : solve_box_(db::create<simple_tags, compute_tags>()),
        // The implicit weight depends on on the step pattern, not on
        // any of the values in the history.
        implicit_weight_(
            db::get<::Tags::TimeStepper<>>(evolution_box)
                .implicit_weight(db::get<Tags::ImplicitHistory<ImplicitSector>>(
                                     evolution_box),
                                 db::get<::Tags::TimeStep>(evolution_box))) {
    db::mutate_apply<
        tmpl::push_front<
            tmpl::filter<
                simple_tags,
                tt::is_a<Variables, tmpl::bind<tmpl::type_from, tmpl::_1>>>,
            EvolutionBoxTag>,
        tmpl::list<>>(
        [&evolution_box](
            const gsl::not_null<const EvolutionBox**> evolution_box_pointer,
            const auto... vars) {
          *evolution_box_pointer = &evolution_box;
          expand_pack((vars->initialize(1, 0.0), 0)...);
        },
        make_not_null(&solve_box_));
  }

  void set_index(
      const gsl::not_null<TimeSteppers::History<SectorVariables>*> history,
      const size_t index) {
    db::mutate<SolverPointIndex>(
        make_not_null(&solve_box_),
        [&index](const gsl::not_null<size_t*> box_index) {
          *box_index = index;
        });

    const auto& evolution_box = *db::get<EvolutionBoxTag>(solve_box_);
    inhomogeneous_terms_ = db::get<ExplicitValue>(solve_box_);
    const ImexTimeStepper& time_stepper =
        db::get<::Tags::TimeStepper<>>(evolution_box);
    const TimeDelta& time_step = db::get<::Tags::TimeStep>(evolution_box);
    time_stepper.add_inhomogeneous_implicit_terms(
        make_not_null(&inhomogeneous_terms_), history, time_step);
  }

  std::array<double, solve_dimension> operator()(
      const std::array<double, solve_dimension>& sector_variables_array) const {
    ASSERT(implicit_weight_ != 0.0,
           "Should not be performing solves on explicit substeps");
    set_sector_variables(sector_variables_array);
    run_mutators<tmpl::push_back<typename ImplicitSector::source_prep,
                                 typename ImplicitSector::source>>();
    std::array<double, solve_dimension> residual_array{};
    SectorVariables residual(residual_array.data(), residual_array.size());
    residual =
        inhomogeneous_terms_ - db::get<sector_variables_tag>(solve_box_) +
        implicit_weight_ *
            db::get<db::add_tag_prefix<::Tags::Source, sector_variables_tag>>(
                solve_box_);
    return residual_array;
  }

  std::array<std::array<double, solve_dimension>, solve_dimension> jacobian(
      const std::array<double, solve_dimension>& sector_variables_array) const {
    ASSERT(implicit_weight_ != 0.0,
           "Should not be performing solves on explicit substeps");
    set_sector_variables(sector_variables_array);
    run_mutators<tmpl::push_back<typename ImplicitSector::jacobian_prep,
                                 typename ImplicitSector::jacobian>>();

    std::array<std::array<double, solve_dimension>, solve_dimension>
        jacobian_array{};
    // The storage order for the tensors does not match the required
    // order for the returned array, so we have to copy components
    // individually.
    //
    // Despite repeated references to then, the result of this is
    // independent of the *_for_offsets variables.  They are only used
    // for calculating offsets into the returned array.
    const auto& variables_for_offsets =
        db::get<sector_variables_tag>(solve_box_);
    tmpl::for_each<typename SectorVariables::tags_list>(
        [&](auto dependent_tag_v) {
          using dependent_tag = tmpl::type_from<decltype(dependent_tag_v)>;
          const auto& dependent_for_offsets =
              get<dependent_tag>(variables_for_offsets);
          for (size_t dependent_component = 0;
               dependent_component < dependent_for_offsets.size();
               ++dependent_component) {
            const auto dependent_index =
                dependent_for_offsets.get_tensor_index(dependent_component);
            auto& result_row = jacobian_array[static_cast<size_t>(
                dependent_for_offsets[dependent_component].data() -
                variables_for_offsets.data())];
            tmpl::for_each<typename SectorVariables::tags_list>(
                [&](auto independent_tag_v) {
                  using independent_tag =
                      tmpl::type_from<decltype(independent_tag_v)>;
                  using jacobian_component_tag =
                      imex::Tags::Jacobian<independent_tag,
                                           ::Tags::Source<dependent_tag>>;
                  const auto& independent_for_offsets =
                      get<independent_tag>(variables_for_offsets);

                  for (size_t independent_component = 0;
                       independent_component < independent_for_offsets.size();
                       ++independent_component) {
                    const auto independent_index =
                        independent_for_offsets.get_tensor_index(
                            independent_component);
                    result_row[static_cast<size_t>(
                        independent_for_offsets[independent_component].data() -
                        variables_for_offsets.data())] =
                        get<jacobian_component_tag>(solve_box_)
                            .get(concatenate(independent_index,
                                             dependent_index))[0];
                  }
                });
          }
        });

    jacobian_array *= implicit_weight_;

    for (size_t i = 0; i < solve_dimension; ++i) {
      jacobian_array[i][i] -= 1.0;
    }
    return jacobian_array;
  }

  GuessResult compute_initial_guess() {
    if (implicit_weight_ == 0.0) {
      db::mutate<sector_variables_tag>(
          make_not_null(&solve_box_),
          [this](const gsl::not_null<SectorVariables*> sector_variables) {
            *sector_variables = inhomogeneous_terms_;
          });
      return GuessResult::ExactSolution;
    }

    set_sector_variables(db::get<ExplicitValue>(solve_box_));
    run_mutators<typename ImplicitSector::initial_guess_prep>();
    return db::mutate_apply<typename ImplicitSector::initial_guess>(
        make_not_null(&solve_box_), inhomogeneous_terms_, implicit_weight_);
  }

  std::array<double, solve_dimension> initial_guess() const {
    std::array<double, solve_dimension> guess_array{};
    SectorVariables guess(guess_array.data(), guess_array.size());
    // The variables were modified to the initial guess in
    // compute_initial_guess().
    guess = db::get<sector_variables_tag>(solve_box_);
    return guess_array;
  }

 private:
  void set_sector_variables(
      std::array<double, solve_dimension> sector_variables_array) const {
    const SectorVariables sector_variables(sector_variables_array.data(),
                                           sector_variables_array.size());
    set_sector_variables(sector_variables);
  }

  void set_sector_variables(const SectorVariables& sector_variables) const {
    if (sector_variables == most_recent_sector_variables_) {
      return;
    }
    most_recent_sector_variables_ = sector_variables;
    db::mutate<sector_variables_tag>(
        make_not_null(&solve_box_),
        [&sector_variables](const gsl::not_null<SectorVariables*> vars) {
          *vars = sector_variables;
        });
    completed_mutators_ = decltype(completed_mutators_){};
  }

  template <typename Mutators>
  void run_mutators() const {
    tmpl::for_each<Mutators>([this](auto mutator_v) {
      using mutator = tmpl::type_from<decltype(mutator_v)>;
      if (not get<RanMutator<mutator>>(completed_mutators_)) {
        db::mutate_apply<mutator>(make_not_null(&solve_box_));
        get<RanMutator<mutator>>(completed_mutators_) = true;
      }
    });
  }

  template <typename Mutator>
  struct RanMutator {
    using type = bool;
  };

  // Re mutables: This struct is only used locally in serial
  // single-threaded implicit solves.  The gsl_multiroot interface
  // takes a const solver object, but we want to be able to share
  // calculations between the source and jacobian calculations.
  // NOLINTNEXTLINE(spectre-mutable)
  mutable SolveBox solve_box_;
  double implicit_weight_;
  SectorVariables inhomogeneous_terms_{};
  // NOLINTNEXTLINE(spectre-mutable)
  mutable SectorVariables most_recent_sector_variables_{};
  // NOLINTNEXTLINE(spectre-mutable)
  mutable tuples::tagged_tuple_from_typelist<
      tmpl::transform<all_mutators, tmpl::bind<RanMutator, tmpl::_1>>>
      completed_mutators_{};
};
}  // namespace solve_implicit_sector_detail

/// Perform the implicit solve for one implicit sector.
///
/// This will update the tensors in the implicit sector and clean up
/// the corresponding time stepper history.  A new history entry is
/// not added, because that should be done with the same values of the
/// variables used for the explicit portion of the time derivative,
/// which may still undergo variable-fixing-like corrections.
template <typename ImplicitSector, typename DbTags>
void solve_implicit_sector(const gsl::not_null<db::DataBox<DbTags>*> box) {
  using ImplicitVars = Variables<typename ImplicitSector::tensors>;
  solve_implicit_sector_detail::ImplicitSolver<ImplicitSector,
                                               db::DataBox<DbTags>>
      solver(*box);

  Matrix semi_implicit_jacobian{};

  const size_t number_of_grid_points =
      db::get<tmpl::front<typename ImplicitSector::tensors>>(*box)
          .begin()
          ->size();
  auto& implicit_history =
      db::get_mutable_reference<imex::Tags::ImplicitHistory<ImplicitSector>>(
          box);
  using History = std::decay_t<decltype(implicit_history)>;
  History pointwise_history{};
  for (size_t point = 0; point < number_of_grid_points; ++point) {
    transform(make_not_null(&pointwise_history), implicit_history,
              [&](const auto& v) { return extract_point(v, point); });

    std::array<double, ImplicitVars::number_of_independent_components>
        pointwise_vars_array;
    ImplicitVars pointwise_vars(pointwise_vars_array.data(),
                                pointwise_vars_array.size());
    solver.set_index(make_not_null(&pointwise_history), point);
    if (solver.compute_initial_guess() == GuessResult::ExactSolution) {
      pointwise_vars_array = solver.initial_guess();
    } else {
      switch (db::get<Tags::Mode>(*box)) {
        case Mode::Implicit: {
          // FIXME where should these be specified?
          const double tolerance = 1.0e-10;
          const size_t max_iterations = 100;
          pointwise_vars_array = RootFinder::gsl_multiroot(
              solver, solver.initial_guess(),
              RootFinder::StoppingConditions::Residual(tolerance),
              max_iterations);
          break;
        }
        case Mode::SemiImplicit: {
          const auto initial_guess = solver.initial_guess();
          auto correction_array = solver(initial_guess);
          DataVector correction(correction_array.data(),
                                correction_array.size());
          correction *= -1.0;
          semi_implicit_jacobian = solver.jacobian(initial_guess);
          const int lapack_info =
              lapack::general_matrix_linear_solve(&correction,
                                                  &semi_implicit_jacobian);
          if (lapack_info != 0) {
            if (lapack_info < 0) {
              ERROR("LAPACK invalid argument: " << -lapack_info);
            } else {
              ERROR("Semi-implicit inversion was singular at\n"
                    << pointwise_vars);
            }
          }
          pointwise_vars_array = initial_guess + correction_array;
          break;
        }
        default:
          ERROR("Invalid implicit mode");
      }
    }

    // Write the result into the evolution variables.
    db::mutate_apply<typename ImplicitVars::tags_list, tmpl::list<>>(
        [&](const auto... tensors) {
          tmpl::as_pack<typename ImplicitVars::tags_list>(
              [&](auto... tensor_tags) {
                expand_pack((
                    overwrite_point(tensors,
                                    get<tmpl::type_from<decltype(tensor_tags)>>(
                                        pointwise_vars),
                                    point),
                    0)...);
              });
        },
        box);
  }

  // Copy the cleanup to the history for the full element.
  if (pointwise_history.substeps().empty()) {
    implicit_history.clear_substeps();
  }
  while (implicit_history.size() > pointwise_history.size()) {
    implicit_history.pop_front();
  }
}
}  // namespace imex
