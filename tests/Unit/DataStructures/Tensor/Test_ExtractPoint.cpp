// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/ExtractPoint.hpp"
#include "DataStructures/Tensor/Tensor.hpp"


SPECTRE_TEST_CASE("Unit.DataStructures.Tensor.ExtaractPoint",
                  "[DataStructures][Unit]") {
  const Scalar<DataVector> scalar{DataVector{3.0, 4.0, 5.0}};
  CHECK(extract_point(scalar, 0) == Scalar<double>{3.0});
  CHECK(extract_point(scalar, 1) == Scalar<double>{4.0});
  CHECK(extract_point(scalar, 2) == Scalar<double>{5.0});

  tnsr::ii<DataVector, 2> tensor;
  get<0, 0>(tensor) = DataVector{1.0, 2.0};
  get<0, 1>(tensor) = DataVector{3.0, 4.0};
  get<1, 1>(tensor) = DataVector{5.0, 6.0};
  {
    tnsr::ii<double, 2> expected;
    get<0, 0>(tensor) = 1.0;
    get<0, 1>(tensor) = 3.0;
    get<1, 1>(tensor) = 5.0;
    CHECK(extract_point(tensor, 0) == expected);
  }
  {
    tnsr::ii<double, 2> expected;
    get<0, 0>(tensor) = 2.0;
    get<0, 1>(tensor) = 4.0;
    get<1, 1>(tensor) = 6.0;
    CHECK(extract_point(tensor, 1) == expected);
  }
}
}  // namespace
