// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Metafunctions.hpp"

namespace imex {
namespace Tags {
/// Tag for the derivative of an implicit source (`Dependent`) with
/// respect to an implicit variable (`Independent`).
template <typename Independent, typename Dependent>
struct Jacobian : db::SimpleTag {
  using independent = Independent;
  using dependent = Dependent;
  using type = TensorMetafunctions::concatenate_indices<
      TensorMetafunctions::change_all_valences<typename Independent::type>,
      typename Dependent::type>;
};
}  // namespace Tags

/// Create a list of all jacobian tags for the dependence of \p
/// DependentList on \p IndependentList.
///
/// This can be used to construct the argument for a `Variables` from
/// the `tags_list`s of to other `Variables`.
template <typename IndependentList, typename DependentList>
using jacobian_tags = tmpl::join<tmpl::transform<
    IndependentList,
    tmpl::lazy::transform<
        tmpl::pin<DependentList>,
        tmpl::defer<tmpl::bind<::imex::Tags::Jacobian, tmpl::parent<tmpl::_1>,
                               tmpl::_1>>>>>;
}  // namespace imex
