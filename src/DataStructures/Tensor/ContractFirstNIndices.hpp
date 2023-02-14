// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Expressions/TensorExpression.hpp"
#include "DataStructures/Tensor/Expressions/TensorIndex.hpp"
#include "Utilities/TMPL.hpp"

template <typename X, typename Symm, typename IndexList>
class Tensor;

namespace detail {
// Get the values that encode the generic tensor indices for the first operand,
// second operand, and result tensor so that the first NumIndicesToContract
// indices will contract and the TensorExpression written with them will be
// valid
//
// Note: the result tensor's index ordering is simply the order of the free
// indices of the first operand concatenated with the index ordering of the
// free indices of the second operand. For example, if we want to contract
// the first 2 indices with
// `R(ti::a, ti::b, ti::c, ti::i) * S(ti::A, ti::B, ti::d)`, the result tensor's
// index order will be `ti::c, ti::i, ti::d`.
template <size_t NumIndicesToContract, size_t NumIndices1, size_t NumIndices2>
constexpr auto
get_tensor_index_values_for_tensors_to_contract_and_result_tensor(
    const std::array<bool, NumIndices1>& index_type_is_spacetime1,
    const std::array<bool, NumIndices1>& valence_is_lower1,
    const std::array<bool, NumIndices2>& index_type_is_spacetime2,
    const std::array<bool, NumIndices2>& valence_is_lower2) {
  constexpr size_t num_result_indices =
      NumIndices1 + NumIndices2 - 2 * NumIndicesToContract;

  // (first op index values, second op index values, result tensor index values)
  std::tuple<std::array<size_t, NumIndices1>, std::array<size_t, NumIndices2>,
             std::array<size_t, num_result_indices>>
      tensor_index_values{};

  // the next lower spacetime index value that we have not yet used
  size_t next_lower_spacetime_value = 0;
  // the next lower spatial index value that we have not yet used
  size_t next_lower_spatial_value = tenex::TensorIndex_detail::spatial_sentinel;

  // assign first operand's tensor index values
  for (size_t i = 0; i < NumIndices1; i++) {
    const bool index1_is_spacetime = index_type_is_spacetime1[i];
    const bool index1_is_lower = valence_is_lower1[i];
    if (index1_is_spacetime) {
      std::get<0>(tensor_index_values)[i] =
          index1_is_lower ? next_lower_spacetime_value
                          : tenex::get_tensorindex_value_with_opposite_valence(
                                next_lower_spacetime_value);
      next_lower_spacetime_value++;
    } else {
      std::get<0>(tensor_index_values)[i] =
          index1_is_lower ? next_lower_spatial_value
                          : tenex::get_tensorindex_value_with_opposite_valence(
                                next_lower_spatial_value);
      next_lower_spatial_value++;
    }
  }

  // assign the first NumIndicesToContract index values of the second operand so
  // that they will contract with the first NumIndicesToContract indices of the
  // first operand
  for (size_t i = 0; i < NumIndicesToContract; i++) {
    std::get<1>(tensor_index_values)[i] =
        tenex::get_tensorindex_value_with_opposite_valence(
            std::get<0>(tensor_index_values)[i]);
  }

  // assign the remaining index values of the second operand so that they are
  // not duplicates of any previously-used indices
  for (size_t i = NumIndicesToContract; i < NumIndices2; i++) {
    const bool index2_is_spacetime = index_type_is_spacetime2[i];
    const bool index2_is_lower = valence_is_lower2[i];
    if (index2_is_spacetime) {
      std::get<1>(tensor_index_values)[i] =
          index2_is_lower ? next_lower_spacetime_value
                          : tenex::get_tensorindex_value_with_opposite_valence(
                                next_lower_spacetime_value);
      next_lower_spacetime_value++;
    } else {
      std::get<1>(tensor_index_values)[i] =
          index2_is_lower ? next_lower_spatial_value
                          : tenex::get_tensorindex_value_with_opposite_valence(
                                next_lower_spatial_value);
      next_lower_spatial_value++;
    }
  }

  // assign the free indices of the first operand to the result tensor indices
  for (size_t i = 0; i < NumIndices1 - NumIndicesToContract; i++) {
    std::get<2>(tensor_index_values)[i] =
        std::get<0>(tensor_index_values)[i + NumIndicesToContract];
  }

  // assign the free indices of the second operand to the result tensor indices
  for (size_t i = 0; i < NumIndices2 - NumIndicesToContract; i++) {
    std::get<2>(tensor_index_values)[i + NumIndices1 - NumIndicesToContract] =
        std::get<1>(tensor_index_values)[i + NumIndicesToContract];
  }

  return tensor_index_values;
}

// Note: assumes the index types of the indices to be contracted are the same,
// so you can contract a spacetime index with a spacetime index and a spatial
// with a spatial, but if you ask to contract a spacetime with a spatial,
// it doesn't guess that you want to sum over spatial indices, there will just
// be a compile error. Support for this, of course, can be added.
template <size_t NumIndicesToContract, typename T1, typename T2,
          size_t NumIndices1 = tmpl::size<typename T1::symmetry>::value,
          size_t NumIndices2 = tmpl::size<typename T2::symmetry>::value,
          size_t NumResultIndices =
              NumIndices1 + NumIndices2 - 2 * NumIndicesToContract>
struct contract_first_n_indices_impl;

template <size_t NumIndicesToContract, typename X1, typename Symm1,
          typename... Indices1, typename X2, typename Symm2,
          typename... Indices2, size_t NumIndices1, size_t NumIndices2,
          size_t NumResultIndices>
struct contract_first_n_indices_impl<
    NumIndicesToContract, Tensor<X1, Symm1, tmpl::list<Indices1...>>,
    Tensor<X2, Symm2, tmpl::list<Indices2...>>, NumIndices1, NumIndices2,
    NumResultIndices> {
  // whether the indices of the operands are upper, lower, spatial, and/or
  // spacetime indices
  static constexpr std::array<bool, NumIndices1> valence_is_lower1 = {
      {(Indices1::ul == UpLo::Lo)...}};
  static constexpr std::array<bool, NumIndices2> valence_is_lower2 = {
      {(Indices2::ul == UpLo::Lo)...}};
  static constexpr std::array<bool, NumIndices1> index_type_is_spacetime1 = {
      {(Indices1::index_type == IndexType::Spacetime)...}};
  static constexpr std::array<bool, NumIndices2> index_type_is_spacetime2 = {
      {(Indices2::index_type == IndexType::Spacetime)...}};

  // the values of the generic indices (i.e. `TensorIndex::value`s) that
  // uniquely identify different generic indices

  // tuple of first operand's index values, second operand's index values, and
  // the result tensor's index values
  static constexpr auto tensor_index_values =
      get_tensor_index_values_for_tensors_to_contract_and_result_tensor<
          NumIndicesToContract>(index_type_is_spacetime1, valence_is_lower1,
                                index_type_is_spacetime2, valence_is_lower2);
  static constexpr std::array<size_t, NumIndices1> tensor_index_values1 =
      std::get<0>(tensor_index_values);
  static constexpr std::array<size_t, NumIndices2> tensor_index_values2 =
      std::get<1>(tensor_index_values);
  static constexpr std::array<size_t, NumResultIndices>
      result_tensor_index_values = std::get<2>(tensor_index_values);

  // contract first N indices by evaluating `TensorExpression` with the computed
  // `TensorIndex` values
  template <typename T1, typename T2, size_t... Ints1, size_t... Ints2,
            size_t... ResultInts>
  static auto apply(const T1& tensor1, const T2& tensor2,
                    const std::index_sequence<Ints1...>& /*seq1*/,
                    const std::index_sequence<Ints2...>& /*seq2*/,
                    const std::index_sequence<ResultInts...>& /*result_seq*/) {
    return tenex::evaluate_with_lhs_tensorindex_types<
        TensorIndex<result_tensor_index_values[ResultInts]>...>(
        tensor1(TensorIndex<tensor_index_values1[Ints1]>{}...) *
        tensor2(TensorIndex<tensor_index_values2[Ints2]>{}...));
  }
};
}  // namespace detail

template <size_t NumIndicesToContract, typename T1, typename T2>
auto contract_first_n_indices(const T1& tensor1, const T2& tensor2) {
  static constexpr size_t num_indices1 =
      tmpl::size<typename T1::symmetry>::value;
  static constexpr size_t num_indices2 =
      tmpl::size<typename T2::symmetry>::value;
  static constexpr size_t result_num_indices =
      num_indices1 + num_indices2 - 2 * NumIndicesToContract;

  static_assert(
      NumIndicesToContract <= num_indices1 and
          NumIndicesToContract <= num_indices2,
      "Cannot request to contract more indices than one of the tensors has.");

  return detail::contract_first_n_indices_impl<NumIndicesToContract, T1, T2>::
      apply(tensor1, tensor2, std::make_index_sequence<num_indices1>{},
            std::make_index_sequence<num_indices2>{},
            std::make_index_sequence<result_num_indices>{});
}
