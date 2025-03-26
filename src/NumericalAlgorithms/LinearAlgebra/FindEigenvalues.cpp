// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearAlgebra/FindEigenvalues.hpp"

#include <blaze/math/DynamicMatrix.h>
#include <complex>
#include <cstddef>
#include <limits>
#include <ostream>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

// LAPACK routine to do the eigenvalue problem
extern "C" {
// The final two arguments are the "hidden" lengths of the first two.
// https://gcc.gnu.org/onlinedocs/gfortran/Argument-passing-conventions.html
extern void zgeev_(char*, char*, int*, double*, int*, double*, double*, int*,
                   double*, int*, double*, int*, double*, int*, size_t, size_t);
}

void find_eigenvalues(
    const gsl::not_null<ComplexDataVector*> eigenvalues,
    const gsl::not_null<
        blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>*>
        matrix,
    const gsl::not_null<DataVector*> scratch) {
  // Sanity checks on the sizes of the vectors and matrices
  const size_t number_of_rows = matrix->rows();
  ASSERT(number_of_rows == matrix->columns(),
         "Matrix should be square, but has " << matrix->rows() << " rows and "
                                             << matrix->columns()
                                             << " columns.");
  eigenvalues->destructive_resize(number_of_rows);

  // Set up parameters for the lapack call
  // Lapack uses chars to decide whether to compute the left eigenvectors,
  // the right eigenvectors, both, or neither. 'N' means do not compute,
  // 'V' means do compute. Note: not const because lapack does not want this
  // option const.
  char compute_left_eigenvectors = 'N';
  char compute_right_eigenvectors = 'N';

  // Lapack expects the sizes to be ints, not size_t.
  // NOTE: not const because lapack function arguments are not const.
  auto matrix_and_vector_size = static_cast<int>(number_of_rows);

  // Size of the unused eigenvector output.  Must be positive
  // according to the docs, even though unused.
  int unused_eigenvector_size = 1;

  // Lapack uses a work vector.  -1 indicates a query for the optimal
  // size.
  int work_size = -1;
  std::complex<double> optimal_work_size(
      std::numeric_limits<double>::signaling_NaN(),
      std::numeric_limits<double>::signaling_NaN());

  //  Lapack uses an integer called info to return its status
  //  info = 0 : success
  //  info = -i: ith argument had bad value
  //  info > 0: some other failure
  int info = 0;

  int matrix_spacing = matrix->spacing();

  // Query for the work size
  zgeev_(&compute_left_eigenvectors, &compute_right_eigenvectors,
         &matrix_and_vector_size, reinterpret_cast<double*>(matrix->data()),
         &matrix_spacing, reinterpret_cast<double*>(eigenvalues->data()),
         nullptr, &unused_eigenvector_size, nullptr, &unused_eigenvector_size,
         reinterpret_cast<double*>(&optimal_work_size), &work_size, nullptr,
         &info, 1, 1);

  if (UNLIKELY(info != 0)) {
    ERROR("Lapack failed to compute workspace requirements. Lapack's INFO = "
          << info);
  }

  ASSERT(optimal_work_size.imag() == 0.0, "Lapack requests complex work");
  work_size = static_cast<int>(optimal_work_size.real());
  const size_t needed_scratch =
      2 * static_cast<size_t>(work_size) + 2 * number_of_rows;
  if (scratch->size() < needed_scratch) {
    scratch->destructive_resize(needed_scratch);
  }

  zgeev_(&compute_left_eigenvectors, &compute_right_eigenvectors,
         &matrix_and_vector_size, reinterpret_cast<double*>(matrix->data()),
         &matrix_spacing, reinterpret_cast<double*>(eigenvalues->data()),
         nullptr, &unused_eigenvector_size, nullptr, &unused_eigenvector_size,
         reinterpret_cast<double*>(scratch->data()), &work_size,
         scratch->data() + 2 * work_size, &info, 1, 1);

  if (UNLIKELY(info != 0)) {
    ERROR("Lapack failed to compute eigenvectors. Lapack's zgeev INFO = "
          << info);
  }
}

ComplexDataVector find_eigenvalues(
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> matrix) {
  ComplexDataVector eigenvalues{};
  DataVector scratch{};
  find_eigenvalues(&eigenvalues, &matrix, &scratch);
  return eigenvalues;
}
