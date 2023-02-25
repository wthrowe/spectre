// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
class DataVector;
template <typename X, typename Symm, typename IndexList>
class Tensor;
namespace Tags {
template <typename Tag>
struct Source;
}  // namespace Tags
/// \endcond

namespace imex::protocols {
/// Protocol for an implicit sector of an IMEX system.
///
/// An implicit sector describes the sources for one implicit solve
/// performed during IMEX evolution.  A system may have multiple
/// implicit sectors, but they must be independent, i.e., their
/// sources must not depend on any tensors in other sectors.
///
/// Classes implementing this protocol must define:
/// * a `tensors` type alias of tags for the variables to be solved for
/// * lists used to construct a DataBox during the pointwise implicit solve:
///   * `tags_from_evolution` for tags to be made available from the evolution
///      DataBox.  Volume quantities will be reduced to have one grid point,
///      with the appropriate value for the point being solved for.
///   * `simple_tags` for temporaries (e.g., primitives)
///   * `compute_tags`
/// * a `source` type to be passed to `db::mutate_apply` to compute the
///   sources.
/// * a `jacobian` type to be passed to `db::mutate_apply` to compute the
///   source jacobian.  If the implicit equation can always be solved
///   analytically for the sector, the jacobian is not required and this may
///   be the type `imex::NoJacobianBecauseSolutionIsAnalytic`.
/// * an `initial_guess` type to be passed to `db::mutate_apply`, taking
///   additional arguments for the inhomogeneous terms \f$X\f$ and implicit
///   weight \f$w\f$ in the equation to be solved: \f$u = X + w S(u)\f$.  (See
///   example below.)  It must return a `GuessResult` indicating whether the
///   implicit equation has been solved analytically or whether the numerical
///   solve should continue.  This mutator will not be called if the implicit
///   weight is zero, as a system-independent analytic solution is available
///   in that case.  If using the value of the explicit step as an initial
///   guess is acceptable, this can be the type `imex::GuessExplicitResult`.
/// * lists `source_prep`, `jacobian_prep`, and `initial_guess_prep` that will
///   be called before the corresponding main mutator, e.g., for computing
///   primitives.  Mutators appearing in multiple lists, as well as the
///   `source` and `jacobian` mutators, will be skipped if they have already
///   been applied for the current point.  Note that the `source_prep` mutators
///   are only used during the implicit solve, and any preparation needed
///   before the `source` call in the main action loop to record the history
///   is the responsibility of the user.
///
/// All `Variables` in the DataBox, including the sources and source
/// jacobian, will be initialized to zero with a single grid point.
struct ImplicitSector {
  template <typename ConformingType>
  struct test {
    using tensors = typename ConformingType::tensors;
    using source = typename ConformingType::source;
    using jacobian = typename ConformingType::jacobian;
    using initial_guess = typename ConformingType::initial_guess;

    using tags_from_evolution = typename ConformingType::tags_from_evolution;
    using simple_tags = typename ConformingType::simple_tags;
    using compute_tags = typename ConformingType::compute_tags;

    using source_prep = typename ConformingType::source_prep;
    using jacobian_prep = typename ConformingType::jacobian_prep;
    using initial_guess_prep = typename ConformingType::initial_guess_prep;

    static_assert(tt::is_a_v<tmpl::list, tensors>);
    static_assert(tt::is_a_v<tmpl::list, tags_from_evolution>);
    static_assert(tt::is_a_v<tmpl::list, simple_tags>);
    static_assert(tt::is_a_v<tmpl::list, compute_tags>);
    static_assert(tt::is_a_v<tmpl::list, source_prep>);
    static_assert(tt::is_a_v<tmpl::list, jacobian_prep>);
    static_assert(tt::is_a_v<tmpl::list, initial_guess_prep>);

    static_assert(
        tmpl::all<
            tensors,
            tt::is_a<Tensor, tmpl::bind<tmpl::type_from, tmpl::_1>>>::value);

    using source_tensors =
        tmpl::transform<tensors, tmpl::bind<::Tags::Source, tmpl::_1>>;

    static_assert(
        std::is_same_v<
            tmpl::list_difference<source_tensors, typename source::return_tags>,
            tmpl::list<>> and
        std::is_same_v<
            tmpl::list_difference<typename source::return_tags, source_tensors>,
            tmpl::list<>>,
        "Implicit source must provide sources for the entire sector.");

    template <typename T>
    struct is_a_tensor_of_data_vector : std::false_type {};

    template <typename Symm, typename IndexList>
    struct is_a_tensor_of_data_vector<Tensor<DataVector, Symm, IndexList>>
        : std::true_type {};

    static_assert(tmpl::none<simple_tags,
                             is_a_tensor_of_data_vector<
                                 tmpl::bind<tmpl::type_from, tmpl::_1>>>::value,
                  "Do not include tags for Tensor<DataVector> in simple_tags, "
                  "because they trigger many memory allocations.  Add the "
                  "tensors as part of a Variables instead.");
  };
};
}  // namespace imex::protocols
