// Distributed under the MIT License.
// See LICENSE.txt for details.

// FIXME
#include <blaze/math/DynamicMatrix.h>
//#include <boost/program_options.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/LiftFlux.hpp"
#include "NumericalAlgorithms/LinearAlgebra/FindEigenvalues.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/DifferentiationMatrix.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

// Charm looks for this function but since we build without a main function or
// main module we just have it be empty
extern "C" void CkRegisterMainModule(void) {}

namespace {
struct DummyField : db::SimpleTag {
  using type = Scalar<DataVector>;
};

// FIXME name, factor
ComplexDataVector eigenvalues_for_phase(const size_t num_points,
                                        const double element_phase) {
  blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> derivative =
      Spectral::differentiation_matrix<Spectral::Basis::Legendre,
                                       Spectral::Quadrature::GaussLobatto>(
                                           num_points);
  // Scale element to unit length.
  derivative *= -2.0;

  // LGL
  Variables<tmpl::list<DummyField>> boundary_correction(1, 1.0);
  // Scale element to unit length
  const Scalar<DataVector> magnitude_of_face_normal(1_st, 2.0);
  dg::lift_flux(make_not_null(&boundary_correction), num_points,
                magnitude_of_face_normal);
  // For upwind LGL, only one endpoint is affected
  // FIXME signs?
  derivative(0, 0) += get(get<DummyField>(boundary_correction))[0];
  derivative(0, num_points - 1) -=
      std::polar(1.0, element_phase) *
      get(get<DummyField>(boundary_correction))[0];

  return find_eigenvalues(std::move(derivative));

  // LG
  // dg::lift_boundary_terms_gauss_points(...)
}

// Matrix time_stepper_matrix(const TimeStepper& stepper) {
//   Matrix stepper_matrix{};
// }
}  // namespace

int main(const int argc, char** const argv) {
  //try {
    const size_t phase_steps = 100;
    for (size_t i = 0; i < phase_steps; ++i) {
      const double phase = 2 * M_PI / phase_steps * i;
      const auto eigenvalues = eigenvalues_for_phase(7, phase);
      for (const auto& eigenvalue : eigenvalues) {
        Parallel::printf("%.18g\t%.18g\n", eigenvalue.real(),
                         eigenvalue.imag());
      }
    }
  // } catch (const bpo::error& e) {
  //   ERROR_NO_TRACE(e.what());
  // }
  return 0;
}
