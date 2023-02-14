// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <climits>
#include <cstddef>
#include <random>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/ContractFirstNIndices.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace {

template <typename Generator, typename DataType>
void test(const gsl::not_null<Generator*> generator,
          const DataType& used_for_size) {
  std::uniform_real_distribution<> distribution(-1.0, 1.0);

  const auto R =
      make_with_random_values<tnsr::abc<DataType, 3, Frame::Inertial>>(
          generator, distribution, used_for_size);
  const auto S =
      make_with_random_values<tnsr::ABc<DataType, 3, Frame::Inertial>>(
          generator, distribution, used_for_size);

  // contract first index (spacetime)
  // tnsr::abCd
  const Tensor<DataType, Symmetry<4, 3, 2, 1>,
               index_list<SpacetimeIndex<3, UpLo::Lo, Frame::Inertial>,
                          SpacetimeIndex<3, UpLo::Lo, Frame::Inertial>,
                          SpacetimeIndex<3, UpLo::Up, Frame::Inertial>,
                          SpacetimeIndex<3, UpLo::Lo, Frame::Inertial>>>
      RS_contract_first_index = contract_first_n_indices<1>(R, S);

  for (size_t b = 0; b < 4; b++) {
    for (size_t c = 0; c < 4; c++) {
      for (size_t d = 0; d < 4; d++) {
        for (size_t e = 0; e < 4; e++) {
          DataType expected_sum = make_with_value<DataType>(used_for_size, 0.0);
          for (size_t a = 0; a < 4; a++) {
            expected_sum += R.get(a, b, c) * S.get(a, d, e);
          }
          CHECK_ITERABLE_APPROX(RS_contract_first_index.get(b, c, d, e),
                                expected_sum);
        }
      }
    }
  }

  // contract first two indices (both spacetime)
  const tnsr::ab<DataType, 3, Frame::Inertial> RS_contract_first_2_indices =
      contract_first_n_indices<2>(R, S);

  for (size_t c = 0; c < 4; c++) {
    for (size_t d = 0; d < 4; d++) {
      DataType expected_sum = make_with_value<DataType>(used_for_size, 0.0);
      for (size_t a = 0; a < 4; a++) {
        for (size_t b = 0; b < 4; b++) {
          expected_sum += R.get(a, b, c) * S.get(a, b, d);
        }
      }
      CHECK_ITERABLE_APPROX(RS_contract_first_2_indices.get(c, d),
                            expected_sum);
    }
  }

  const auto G =
      make_with_random_values<tnsr::Ijaa<DataType, 3, Frame::Inertial>>(
          generator, distribution, used_for_size);
  // tnsr::iJA
  const auto H = make_with_random_values<
      Tensor<DataType, Symmetry<3, 2, 1>,
             index_list<SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                        SpatialIndex<3, UpLo::Up, Frame::Inertial>,
                        SpacetimeIndex<3, UpLo::Up, Frame::Inertial>>>>(
      generator, distribution, used_for_size);

  // contract first two indices (both spatial) of two tensors of different rank
  const tnsr::aaB<DataType, 3, Frame::Inertial> GH =
      contract_first_n_indices<2>(G, H);
  // for checking that having the smaller rank as the first arg gives us the
  // "same" result mathematically (though the LHS index order will be different)
  const tnsr::Abb<DataType, 3, Frame::Inertial> HG =
      contract_first_n_indices<2>(H, G);

  for (size_t a = 0; a < 4; a++) {
    for (size_t b = 0; b < 4; b++) {
      for (size_t c = 0; c < 4; c++) {
        DataType expected_sum = make_with_value<DataType>(used_for_size, 0.0);
        for (size_t i = 0; i < 3; i++) {
          for (size_t j = 0; j < 3; j++) {
            expected_sum += G.get(i, j, a, b) * H.get(i, j, c);
          }
        }
        CHECK_ITERABLE_APPROX(GH.get(a, b, c), expected_sum);
        CHECK_ITERABLE_APPROX(HG.get(c, a, b), expected_sum);
      }
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.DataStructures.Tensor.ContractFirstNIndices",
                  "[DataStructures][Unit]") {
  MAKE_GENERATOR(generator);

  test(make_not_null(&generator), std::numeric_limits<double>::signaling_NaN());
  test(make_not_null(&generator),
       DataVector(5, std::numeric_limits<double>::signaling_NaN()));
}
