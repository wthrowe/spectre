// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/ExtractPoint.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct ScalarTag : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct TensorTag : db::SimpleTag {
  using type = tnsr::ii<DataVector, 2>;
};

template <typename... Structure>
Tensor<DataVector, Structure...> promote_to_data_vectors(
    const Tensor<double, Structure...>& tensor) {
  Tensor<DataVector, Structure...> result(1_st);
  for (size_t i = 0; i < result.size(); ++i) {
    result[i] = tensor[i];
  }
  return result;
}
}  // namespace

SPECTRE_TEST_CASE("Unit.DataStructures.ExtractPoint",
                  "[DataStructures][Unit]") {
  // Test with owning Tensors as well as Variables to make sure there
  // are no assumptions about memory layout in the tensor versions.
  const Scalar<DataVector> scalar{DataVector{3.0, 4.0}};
  tnsr::ii<DataVector, 2> tensor;
  get<0, 0>(tensor) = DataVector{1.0, 2.0};
  get<0, 1>(tensor) = DataVector{3.0, 4.0};
  get<1, 1>(tensor) = DataVector{5.0, 6.0};
  Variables<tmpl::list<ScalarTag, TensorTag>> variables(2);
  get<ScalarTag>(variables) = scalar;
  get<TensorTag>(variables) = tensor;

  Scalar<DataVector> reconstructed_scalar{DataVector(2)};
  tnsr::ii<DataVector, 2> reconstructed_tensor(DataVector(2));
  Variables<tmpl::list<ScalarTag, TensorTag>> reconstructed_variables(2);

  const auto check_point = [&](const size_t index,
                               const double scalar_component,
                               const std::array<double, 3>& tensor_components) {
    const Scalar<double> expected_scalar{scalar_component};
    tnsr::ii<double, 2> expected_tensor;
    get<0, 0>(expected_tensor) = tensor_components[0];
    get<0, 1>(expected_tensor) = tensor_components[1];
    get<1, 1>(expected_tensor) = tensor_components[2];
    Variables<tmpl::list<ScalarTag, TensorTag>> expected_variables(1);
    get<ScalarTag>(expected_variables) =
        promote_to_data_vectors(expected_scalar);
    get<TensorTag>(expected_variables) =
        promote_to_data_vectors(expected_tensor);
    CHECK(extract_point(scalar, index) == expected_scalar);
    CHECK(extract_point(tensor, index) == expected_tensor);
    CHECK(extract_point(variables, index) == expected_variables);
    overwrite_point(make_not_null(&reconstructed_scalar), expected_scalar,
                    index);
    overwrite_point(make_not_null(&reconstructed_tensor), expected_tensor,
                    index);
    overwrite_point(make_not_null(&reconstructed_variables), expected_variables,
                    index);
  };
  check_point(0, 3.0, {{1.0, 3.0, 5.0}});
  check_point(1, 4.0, {{2.0, 4.0, 6.0}});

  CHECK(reconstructed_scalar == scalar);
  CHECK(reconstructed_tensor == tensor);
  CHECK(reconstructed_variables == variables);
}
