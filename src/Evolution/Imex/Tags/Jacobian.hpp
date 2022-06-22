// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Metafunctions.hpp"

namespace imex::Tags {
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
}  // namespace imex::Tags
