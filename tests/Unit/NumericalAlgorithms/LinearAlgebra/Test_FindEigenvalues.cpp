// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <blaze/math/DynamicMatrix.h>
#include <complex>
#include <cstddef>
#include <random>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearAlgebra/FindEigenvalues.hpp"
#include "Utilities/Gsl.hpp"

namespace {
void test_result() {
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> matrix(2, 2);
  matrix(0, 0).real(dist(gen));
  matrix(0, 0).imag(dist(gen));
  matrix(0, 1).real(dist(gen));
  matrix(0, 1).imag(dist(gen));
  matrix(1, 0).real(dist(gen));
  matrix(1, 0).imag(dist(gen));
  matrix(1, 1).real(dist(gen));
  matrix(1, 1).imag(dist(gen));

  const auto expected_sum = matrix(0, 0) + matrix(1, 1);
  const auto expected_product =
      matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0);

  const auto allocating_eigenvalues = find_eigenvalues(matrix);
  REQUIRE(allocating_eigenvalues.size() == 2);
  const auto result_sum = sum(allocating_eigenvalues);
  const auto result_product = prod(allocating_eigenvalues);
  CHECK(result_sum.real() == approx(expected_sum.real()));
  CHECK(result_sum.imag() == approx(expected_sum.imag()));
  CHECK(result_product.real() == approx(expected_product.real()));
  CHECK(result_product.imag() == approx(expected_product.imag()));

  ComplexDataVector nonallocating_eigenvalues(2);
  DataVector scratch{};
  {
    auto matrix_copy = matrix;
    find_eigenvalues(&nonallocating_eigenvalues, &matrix_copy, &scratch);
    CHECK(nonallocating_eigenvalues == allocating_eigenvalues);
    CHECK(scratch.size() != 0);
  }

  const auto needed_scratch = scratch.size();
  scratch.destructive_resize(needed_scratch + 1);

  find_eigenvalues(&nonallocating_eigenvalues, &matrix, &scratch);
  CHECK(nonallocating_eigenvalues == allocating_eigenvalues);
  CHECK(scratch.size() == needed_scratch + 1);
}

void test_errors() {
#ifdef SPECTRE_DEBUG
  CHECK_THROWS_WITH(
      find_eigenvalues({2, 3}),
      Catch::Matchers::ContainsSubstring("Matrix should be square"));
#endif
}

SPECTRE_TEST_CASE("Unit.Numerical.LinearAlgebra.Eigenvalue",
                  "[NumericalAlgorithms][LinearAlgebra][Unit]") {
  test_result();
  test_errors();
}
}  // namespace
