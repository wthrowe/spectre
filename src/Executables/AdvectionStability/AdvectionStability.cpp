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
#include "Time/History.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
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

using ComplexMatrix =
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;

//using options = tmpl::list<OptionTags::TimeStepper>;

// FIXME name, factor
ComplexDataVector eigenvalues_for_phase(const size_t num_points,
                                        const double element_phase) {
  ComplexMatrix derivative =
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

double time_stepper_amplification(
    const TimeStepper& stepper, const double time_step,
    const std::complex<double>& deriv_eigenvalue) {
  const size_t hist_size = stepper.number_of_past_steps() + 1;
  ComplexMatrix stepper_operator(hist_size, hist_size, 0.0);
  // History aging entries
  for (size_t i = 1; i < hist_size; ++i) {
    stepper_operator(i, i - 1) = 1.0;
  }

  for (size_t history_entry = 0; history_entry < hist_size; ++history_entry) {
    const Slab slab(0.0, 20.0 * time_step);
    const auto step = slab.duration() / 20;
    TimeSteppers::History<std::complex<double>> history(
        variants::get<TimeSteppers::Tags::FixedOrder>(stepper.order()));
    for (size_t i = 0; i < hist_size; ++i) {
      const std::complex<double> init_value = i == history_entry ? 1.0 : 0.0;
      history.insert_initial(TimeStepId(true, 0, slab.end() - step * (i + 1)),
                             init_value, init_value * deriv_eigenvalue);
    }

    std::complex<double> value = 1.0;
    TimeStepId time_step_id(true, 0, slab.end() - step);
    do {
      stepper.update_u(make_not_null(&value), history, step);
      stepper.clean_history(make_not_null(&history));
      time_step_id = stepper.next_time_id(time_step_id, step);
      history.insert(time_step_id, value, value * deriv_eigenvalue);
    } while (time_step_id.substep() != 0);
    stepper_operator(0, history_entry) = value;
  }

  return max(abs(find_eigenvalues(std::move(stepper_operator))));
}

double largest_amplification(const TimeStepper& stepper, const double time_step,
                             const size_t num_points) {
  const size_t phase_steps = 1000;
  double max_amplification = 0.0;
  for (size_t i = 0; i < phase_steps; ++i) {
    const double phase = 2 * M_PI / phase_steps * i;
    const auto eigenvalues = eigenvalues_for_phase(num_points, phase);
    for (const auto& eigenvalue : eigenvalues) {
      const double amplification = time_stepper_amplification(
          stepper, time_step, eigenvalue);
      max_amplification = std::max(amplification, max_amplification);
    }
  }
  return max_amplification;
}
}  // namespace

int main(const int argc, char** const argv) {
  //try {
  // const double time_step = 1.0e-1;
  // const TimeSteppers::Rk3HesthavenSsp stepper{};
    // const size_t phase_steps = 100;
    // for (size_t i = 0; i < phase_steps; ++i) {
    //   const double phase = 2 * M_PI / phase_steps * i;
    //   const auto eigenvalues = eigenvalues_for_phase(7, phase);
    //   for (const auto& eigenvalue : eigenvalues) {
    //     const double amplification = time_stepper_amplification(
    //         stepper, time_step, eigenvalue);
    //     Parallel::printf("%.18g\t%.18g\t%.18g\n", eigenvalue.real(),
    //                      eigenvalue.imag(), amplification);
    //   }
    // }

  //const TimeSteppers::Rk3HesthavenSsp stepper{};
  //const TimeSteppers::DormandPrince5 stepper{};
  const TimeSteppers::AdamsBashforth stepper{2};
  const size_t num_points = 7;
  const double min_step = 1.0e-3;
  const double max_step = 1.0;
  const int samples = 301;
  for (int sample = 0; sample <= samples; ++sample) {
    const double sample_fraction = static_cast<double>(sample) / samples;
    const double time_step = std::pow(min_step, 1.0 - sample_fraction) *
                             std::pow(max_step, sample_fraction);
    Parallel::printf("%.18g\t%.18g\n", time_step,
                     largest_amplification(stepper, time_step, num_points));
  }

  // } catch (const bpo::error& e) {
  //   ERROR_NO_TRACE(e.what());
  // }
  return 0;
}
