// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "Evolution/Imex/Protocols/ImplicitSource.hpp"
#include "Evolution/Imex/Protocols/ImplicitSourceJacobian.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <typename X, typename Symm, typename IndexList>
class Tensor;
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
/// * a `source` type conforming to protocols::ImplicitSource
/// * a `source_jacobian` type conforming to protocols::ImplicitSourceJacobian.
struct ImplicitSector {
  template <typename ConformingType>
  struct test {
    using tensors = typename ConformingType::tensors;
    using source = typename ConformingType::source;
    using source_jacobian = typename ConformingType::source_jacobian;

    static_assert(tt::is_a_v<tmpl::list, tensors>);
    static_assert(
        tmpl::all<
            tensors,
            tt::is_a<Tensor, tmpl::bind<tmpl::type_from, tmpl::_1>>>::value);

    static_assert(tt::assert_conforms_to_v<source, ImplicitSource>);
    static_assert(
        tt::assert_conforms_to_v<source_jacobian, ImplicitSourceJacobian>);

    template <typename T>
    struct get_tag {
      using type = typename T::tag;
    };

    using sourced_tensors =
        tmpl::transform<typename source::return_tags, get_tag<tmpl::_1>>;
    static_assert(
        std::is_same_v<tmpl::list_difference<sourced_tensors, tensors>,
                       tmpl::list<>> and
            std::is_same_v<tmpl::list_difference<tensors, sourced_tensors>,
                           tmpl::list<>>,
        "Implicit source must provide sources for the entire sector.");

    template <typename T>
    struct get_dependent {
      using type = typename T::dependent;
    };
    template <typename T>
    struct get_independent {
      using type = typename T::independent;
    };

    static_assert(
        std::is_same_v<
            tmpl::list_difference<
                tmpl::transform<typename source_jacobian::return_tags,
                                get_independent<tmpl::_1>>,
                tensors>,
            tmpl::list<>>,
        "Source Jacobian independent variables must be tensors in the sector.");
    static_assert(std::is_same_v<
                      tmpl::list_difference<
                          tmpl::transform<typename source_jacobian::return_tags,
                                          get_tag<get_dependent<tmpl::_1>>>,
                          tensors>,
                      tmpl::list<>>,
                  "Source Jacobian dependent variables must be sources of "
                  "tensors in the sector.");
  };
};
}  // namespace imex::protocols
