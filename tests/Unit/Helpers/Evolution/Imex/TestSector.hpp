// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <limits>
#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Imex/Protocols/ImplicitSector.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace TestHelpers::imex {
namespace test_sector_detail {
template <typename Function, typename Results, typename Args>
void apply(const gsl::not_null<Results*> results, const Args& arguments) {
  tmpl::as_pack<typename Function::return_tags>([&](auto... return_tags) {
    tmpl::as_pack<typename Function::argument_tags>([&](auto... argument_tags) {
      Function::apply(
          make_not_null(
              &get<tmpl::type_from<decltype(return_tags)>>(*results))...,
          get<tmpl::type_from<decltype(argument_tags)>>(arguments)...);
    });
  });
}
}  // namespace test_sector_detail

/// Check that an implicit sector conforms to
/// ::imex::protocols::ImplicitSector and that its `source` and
/// `source_jacobian` are consistent with each other.
template <typename Sector, typename... Arguments>
void test_sector(const tuples::TaggedTuple<Arguments...>& arguments) {
  const size_t differentiation_points = 5;  // Must be odd
  const double dx = 1.0e-2;
  const double deriv_tolerance = 1.0e-5;
  static_assert(
      tt::assert_conforms_to_v<Sector, ::imex::protocols::ImplicitSector>);

  auto jacobian = make_with_value<tuples::tagged_tuple_from_typelist<
      typename Sector::source_jacobian::return_tags>>(
      get<tmpl::front<typename Sector::tensors>>(arguments),
      std::numeric_limits<double>::signaling_NaN());
  test_sector_detail::apply<typename Sector::source_jacobian>(
      make_not_null(&jacobian), arguments);

  const auto spectral_points =
      Spectral::collocation_points<Spectral::Basis::Legendre,
                                   Spectral::Quadrature::Gauss>(
          differentiation_points);
  const Matrix& differentiation_matrix =
      Spectral::differentiation_matrix<Spectral::Basis::Legendre,
                                       Spectral::Quadrature::Gauss>(
          differentiation_points);
  DataVector differentiation_weights(differentiation_points);
  for (size_t i = 0; i < differentiation_points; ++i) {
    differentiation_weights[i] =
        differentiation_matrix((differentiation_points - 1) / 2, i);
  }

  tmpl::for_each<typename Sector::tensors>([&](auto independent_v) {
    using independent = tmpl::type_from<decltype(independent_v)>;
    for (size_t independent_component = 0;
         independent_component < independent::type::size();
         ++independent_component) {
      tuples::TaggedTuple<Arguments...> differentiation_arguments{};
      const auto expand_argument = [&](auto arg_tag) {
        using ArgTag = decltype(arg_tag);
        auto& dest = get<ArgTag>(differentiation_arguments);
        const auto& source = get<ArgTag>(arguments);
        if constexpr (tt::is_a_v<Tensor, typename ArgTag::type>) {
          if constexpr (std::is_same_v<typename ArgTag::type::type,
                                       DataVector>) {
            for (size_t i = 0; i < dest.size(); ++i) {
              ASSERT(source[i].size() == 1,
                     "Must supply one-element DataVectors");
              dest[i] = DataVector(differentiation_points, source[i][0]);
            }
          } else {
            dest = source;
          }
        } else {
          dest = source;
        }
        return 0;
      };
      expand_pack(expand_argument(Arguments{})...);
      get<independent>(differentiation_arguments)[independent_component] +=
          dx * spectral_points;

      auto values = make_with_value<tuples::tagged_tuple_from_typelist<
          typename Sector::source::return_tags>>(
          differentiation_points, std::numeric_limits<double>::signaling_NaN());
      test_sector_detail::apply<typename Sector::source>(
          make_not_null(&values), differentiation_arguments);

      tmpl::for_each<typename Sector::tensors>([&](auto unprefixed_dependent) {
        using dependent =
            ::Tags::Source<tmpl::type_from<decltype(unprefixed_dependent)>>;
        using jacobian_tag = ::imex::Tags::Jacobian<independent, dependent>;
        for (size_t dependent_component = 0;
             dependent_component < dependent::type::size();
             ++dependent_component) {
          const auto jacobian_index = concatenate(
              independent::type::get_tensor_index(independent_component),
              dependent::type::get_tensor_index(dependent_component));
          CAPTURE(pretty_type::get_name<jacobian_tag>());
          CAPTURE(jacobian_index);
          if constexpr (tmpl::list_contains_v<
                            typename Sector::source_jacobian::return_tags,
                            jacobian_tag>) {
            const auto& component_values =
                get<dependent>(values)[dependent_component];
            const double scale = max(abs(component_values));
            auto deriv_approx =
                Approx::custom().scale(scale).epsilon(deriv_tolerance);
            const double derivative =
                dot(differentiation_weights, component_values) / dx;
            CHECK(derivative == deriv_approx(get<jacobian_tag>(jacobian).get(
                                    jacobian_index)[0]));
          } else {
            const auto& component_values =
                get<dependent>(values)[dependent_component];
            CHECK(component_values ==
                  DataVector(differentiation_points, component_values[0]));
          }
        }
      });
    }
  });
}
}  // namespace TestHelpers::imex
