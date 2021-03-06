// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <catch.hpp>

#include "tests/Unit/TestHelpers.hpp"
#include "tests/Unit/Utilities/Test_Blas.hpp"

// [[OutputRegex, TRANSA must be upper or lower case N, T, or C. See the BLAS
// documentation for help.]]
[[noreturn]] TEST_CASE("Unit.Utilities.Blas.dgemm_error_transa_false",
                       "[Unit][Utilities]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  test_blas_asserts_for_bad_char::dgemm_error_transa_false();
  ERROR("Assertion test failed incorrectly");
#endif
}

// [[OutputRegex, TRANSB must be upper or lower case N, T, or C. See the BLAS
// documentation for help.]]
[[noreturn]] TEST_CASE("Unit.Utilities.Blas.dgemm_error_transb_false",
                       "[Unit][Utilities]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  test_blas_asserts_for_bad_char::dgemm_error_transb_false();
  ERROR("Assertion test failed incorrectly");
#endif
}

// [[OutputRegex, TRANSA must be upper or lower case N, T, or C. See the BLAS
// documentation for help.]]
[[noreturn]] TEST_CASE("Unit.Utilities.Blas.dgemm_error_transa_true",
                       "[Unit][Utilities]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  test_blas_asserts_for_bad_char::dgemm_error_transa_true();
  ERROR("Assertion test failed incorrectly");
#endif
}

// [[OutputRegex, TRANSB must be upper or lower case N, T, or C. See the BLAS
// documentation for help.]]
[[noreturn]] TEST_CASE("Unit.Utilities.Blas.dgemm_error_transb_true",
                       "[Unit][Utilities]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  test_blas_asserts_for_bad_char::dgemm_error_transb_true();
  ERROR("Assertion test failed incorrectly");
#endif
}

// [[OutputRegex, TRANS must be upper or lower case N, T, or C. See the BLAS
// documentation for help.]]
[[noreturn]] TEST_CASE("Unit.Utilities.Blas.dgemv_error_trans",
                       "[Unit][Utilities]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  test_blas_asserts_for_bad_char::dgemv_error_trans();
  ERROR("Assertion test failed incorrectly");
#endif
}
