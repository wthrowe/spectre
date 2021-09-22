// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/VariableFixing/FixToAtmosphere.hpp"

#include <pup.h>  // IWYU pragma: keep

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"

// IWYU pragma: no_include <array>

// IWYU pragma: no_forward_declare EquationsOfState::EquationOfState
// IWYU pragma: no_forward_declare Tensor

namespace VariableFixing {

template <size_t Dim>
FixToAtmosphere<Dim>::FixToAtmosphere(const double density_of_atmosphere,
                                      const double density_cutoff,
                                      const double transition_density_cutoff,
                                      const double max_velocity_magnitude,
                                      const Options::Context& context)
    : density_of_atmosphere_(density_of_atmosphere),
      density_cutoff_(density_cutoff),
      transition_density_cutoff_(transition_density_cutoff),
      max_velocity_magnitude_(max_velocity_magnitude) {
  if (density_of_atmosphere_ > density_cutoff_) {
    PARSE_ERROR(context, "The cutoff density ("
                             << density_cutoff_
                             << ") must be greater than or equal to the "
                                "density value in the atmosphere ("
                             << density_of_atmosphere_ << ')');
  }
  if (transition_density_cutoff_ < density_of_atmosphere_ or
      transition_density_cutoff_ > 10.0 * density_of_atmosphere_) {
    PARSE_ERROR(context, "The transition density must be in ["
                             << density_of_atmosphere_ << ", "
                             << 10 * density_of_atmosphere_ << "], but is "
                             << transition_density_cutoff_);
  }
  if (transition_density_cutoff_ <= density_cutoff_) {
    PARSE_ERROR(context, "The transition density cutoff ("
                             << transition_density_cutoff_
                             << ") must be bigger than the density cutoff ("
                             << density_cutoff_ << ")");
  }
}

template <size_t Dim>
// NOLINTNEXTLINE(google-runtime-references)
void FixToAtmosphere<Dim>::pup(PUP::er& p) {
  p | density_of_atmosphere_;
  p | density_cutoff_;
  p | transition_density_cutoff_;
  p | max_velocity_magnitude_;
}

template <size_t Dim>
template <size_t ThermodynamicDim>
void FixToAtmosphere<Dim>::operator()(
    const gsl::not_null<Scalar<DataVector>*> rest_mass_density,
    const gsl::not_null<Scalar<DataVector>*> specific_internal_energy,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        spatial_velocity,
    const gsl::not_null<Scalar<DataVector>*> lorentz_factor,
    const gsl::not_null<Scalar<DataVector>*> pressure,
    const gsl::not_null<Scalar<DataVector>*> specific_enthalpy,
    const tnsr::ii<DataVector, Dim, Frame::Inertial>& spatial_metric,
    const EquationsOfState::EquationOfState<true, ThermodynamicDim>&
        equation_of_state) const {
  for (size_t i = 0; i < rest_mass_density->get().size(); i++) {
    if (UNLIKELY(rest_mass_density->get()[i] < density_cutoff_)) {
      set_density_to_atmosphere(rest_mass_density, specific_internal_energy,
                                pressure, specific_enthalpy, equation_of_state,
                                i);
      for (size_t d = 0; d < Dim; ++d) {
        spatial_velocity->get(d)[i] = 0.0;
      }
      lorentz_factor->get()[i] = 1.0;
    } else if (UNLIKELY(rest_mass_density->get()[i] <
                        transition_density_cutoff_)) {
      set_to_magnetic_free_transition(rest_mass_density, spatial_velocity,
                                      lorentz_factor, spatial_metric, i);
    }
  }
}

namespace {
template <typename T0, typename... Ts>
bool all_same(const T0& t0, const Ts&... ts) {
  return (... and (t0 == ts));
}

template <typename T0, typename... Ts>
std::string all_same_error(const T0& t0, const Ts&... ts) {
  std::ostringstream ss;
  ss << " Values: " << t0;
  (..., (ss << " " << ts));
  return ss.str();
}
}  // namespace

template <size_t Dim>
template <size_t ThermodynamicDim>
void FixToAtmosphere<Dim>::fix_ghost_data(
    gsl::not_null<Scalar<DataVector>*> rest_mass_density,
    gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        lorentz_factor_times_spatial_velocity,
    gsl::not_null<Scalar<DataVector>*> pressure,
    const tnsr::ii<DataVector, Dim, Frame::Inertial>& spatial_metric,
    const EquationsOfState::EquationOfState<true, ThermodynamicDim>&
        equation_of_state) const {
  ASSERT(
      all_same(get(*rest_mass_density).size(),
               get<0>(*lorentz_factor_times_spatial_velocity).size(),
               get(*pressure).size(), get<0, 0>(spatial_metric).size()),
      all_same_error(get(*rest_mass_density).size(),
                     get<0>(*lorentz_factor_times_spatial_velocity).size(),
                     get(*pressure).size(), get<0, 0>(spatial_metric).size()));
  for (size_t i = 0; i < rest_mass_density->get().size(); i++) {
    if (UNLIKELY(rest_mass_density->get()[i] < density_cutoff_)) {
      set_density_to_atmosphere(rest_mass_density, {}, pressure, {},
                                equation_of_state, i);
      for (size_t d = 0; d < Dim; ++d) {
        lorentz_factor_times_spatial_velocity->get(d)[i] = 0.0;
      }
    } else if (UNLIKELY(rest_mass_density->get()[i] <
                        transition_density_cutoff_)) {
      set_to_magnetic_free_transition(rest_mass_density,
                                      lorentz_factor_times_spatial_velocity,
                                      spatial_metric, i);
    }
  }
}

template <size_t Dim>
template <size_t ThermodynamicDim>
void FixToAtmosphere<Dim>::set_density_to_atmosphere(
    const gsl::not_null<Scalar<DataVector>*> rest_mass_density,
    const std::optional<gsl::not_null<Scalar<DataVector>*>>
        specific_internal_energy,
    const gsl::not_null<Scalar<DataVector>*> pressure,
    const std::optional<gsl::not_null<Scalar<DataVector>*>> specific_enthalpy,
    const EquationsOfState::EquationOfState<true, ThermodynamicDim>&
        equation_of_state,
    const size_t grid_index) const {
  rest_mass_density->get()[grid_index] = density_of_atmosphere_;
  Scalar<double> atmosphere_density{density_of_atmosphere_};
  if constexpr (ThermodynamicDim == 1) {
    pressure->get()[grid_index] =
        get(equation_of_state.pressure_from_density(atmosphere_density));
    if (specific_internal_energy.has_value()) {
      (*specific_internal_energy)->get()[grid_index] =
          get(equation_of_state.specific_internal_energy_from_density(
              atmosphere_density));
      (*specific_enthalpy)->get()[grid_index] = get(
          equation_of_state.specific_enthalpy_from_density(atmosphere_density));
    }
  } else if constexpr (ThermodynamicDim == 2) {
    Scalar<double> atmosphere_energy{0.0};
    pressure->get()[grid_index] =
        get(equation_of_state.pressure_from_density_and_energy(
            atmosphere_density, atmosphere_energy));
    if (specific_internal_energy.has_value()) {
      (*specific_internal_energy)->get()[grid_index] = get(atmosphere_energy);
      (*specific_enthalpy)->get()[grid_index] =
          get(equation_of_state.specific_enthalpy_from_density_and_energy(
              atmosphere_density, atmosphere_energy));
    }
  }
}

namespace {
template <size_t Dim>
bool cap_magnitude(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> vector,
    const tnsr::ii<DataVector, Dim, Frame::Inertial>& spatial_metric,
    const double maximum_magnitude, const size_t grid_index) {
  double magnitude = 0.0;
  for (size_t j = 0; j < Dim; ++j) {
    magnitude += vector->get(j)[grid_index] * vector->get(j)[grid_index] *
                 spatial_metric.get(j, j)[grid_index];
    for (size_t k = j + 1; k < Dim; ++k) {
      magnitude += 2.0 * vector->get(j)[grid_index] *
                   vector->get(k)[grid_index] *
                   spatial_metric.get(j, k)[grid_index];
    }
  }
  magnitude = sqrt(magnitude);
  const bool needs_fixing = magnitude > maximum_magnitude;
  if (needs_fixing) {
    for (size_t j = 0; j < Dim; ++j) {
      vector->get(j)[grid_index] *= maximum_magnitude / magnitude;
    }
  }
  return needs_fixing;
}
}  // namespace

template <size_t Dim>
void FixToAtmosphere<Dim>::set_to_magnetic_free_transition(
    const gsl::not_null<Scalar<DataVector>*> rest_mass_density,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        spatial_velocity,
    const gsl::not_null<Scalar<DataVector>*> lorentz_factor,
    const tnsr::ii<DataVector, Dim, Frame::Inertial>& spatial_metric,
    const size_t grid_index) const {
  const double scale_factor =
      (get(*rest_mass_density)[grid_index] - density_cutoff_) /
      (transition_density_cutoff_ - density_cutoff_);
  const double max_mag_of_velocity = scale_factor * max_velocity_magnitude_;
  if (cap_magnitude(spatial_velocity, spatial_metric, max_mag_of_velocity,
                    grid_index)) {
    get(*lorentz_factor)[grid_index] =
        1.0 / sqrt(1.0 - max_mag_of_velocity * max_mag_of_velocity);
  }
}

template <size_t Dim>
void FixToAtmosphere<Dim>::set_to_magnetic_free_transition(
    const gsl::not_null<Scalar<DataVector>*> rest_mass_density,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        lorentz_factor_times_spatial_velocity,
    const tnsr::ii<DataVector, Dim, Frame::Inertial>& spatial_metric,
    const size_t grid_index) const {
  const double scale_factor =
      (get(*rest_mass_density)[grid_index] - density_cutoff_) /
      (transition_density_cutoff_ - density_cutoff_);
  const double max_mag_of_velocity = scale_factor * max_velocity_magnitude_;
  const double max_lorentz_factor_times_velocity =
      max_mag_of_velocity /
      sqrt(1.0 - max_mag_of_velocity * max_mag_of_velocity);
  cap_magnitude(lorentz_factor_times_spatial_velocity, spatial_metric,
                max_lorentz_factor_times_velocity, grid_index);
}

template <size_t Dim>
bool operator==(const FixToAtmosphere<Dim>& lhs,
                const FixToAtmosphere<Dim>& rhs) {
  return lhs.density_of_atmosphere_ == rhs.density_of_atmosphere_ and
         lhs.density_cutoff_ == rhs.density_cutoff_ and
         lhs.transition_density_cutoff_ == rhs.transition_density_cutoff_ and
         lhs.max_velocity_magnitude_ == rhs.max_velocity_magnitude_;
}

template <size_t Dim>
bool operator!=(const FixToAtmosphere<Dim>& lhs,
                const FixToAtmosphere<Dim>& rhs) {
  return not(lhs == rhs);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define THERMO_DIM(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATION(r, data)                                     \
  template class FixToAtmosphere<DIM(data)>;                       \
  template bool operator==(const FixToAtmosphere<DIM(data)>& lhs,  \
                           const FixToAtmosphere<DIM(data)>& rhs); \
  template bool operator!=(const FixToAtmosphere<DIM(data)>& lhs,  \
                           const FixToAtmosphere<DIM(data)>& rhs);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION

#define INSTANTIATION(r, data)                                                \
  template void FixToAtmosphere<DIM(data)>::operator()(                       \
      const gsl::not_null<Scalar<DataVector>*> rest_mass_density,             \
      const gsl::not_null<Scalar<DataVector>*> specific_internal_energy,      \
      const gsl::not_null<tnsr::I<DataVector, DIM(data), Frame::Inertial>*>   \
          spatial_velocity,                                                   \
      const gsl::not_null<Scalar<DataVector>*> lorentz_factor,                \
      const gsl::not_null<Scalar<DataVector>*> pressure,                      \
      const gsl::not_null<Scalar<DataVector>*> specific_enthalpy,             \
      const tnsr::ii<DataVector, DIM(data), Frame::Inertial>& spatial_metric, \
      const EquationsOfState::EquationOfState<true, THERMO_DIM(data)>&        \
          equation_of_state) const;                                           \
  template void FixToAtmosphere<DIM(data)>::fix_ghost_data(                   \
      gsl::not_null<Scalar<DataVector>*> rest_mass_density,                   \
      gsl::not_null<tnsr::I<DataVector, DIM(data), Frame::Inertial>*>         \
          lorentz_factor_times_spatial_velocity,                              \
      gsl::not_null<Scalar<DataVector>*> pressure,                            \
      const tnsr::ii<DataVector, DIM(data), Frame::Inertial>& spatial_metric, \
      const EquationsOfState::EquationOfState<true, THERMO_DIM(data)>&        \
          equation_of_state) const;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3), (1, 2))

#undef DIM
#undef THERMO_DIM
#undef INSTANTIATION

}  // namespace VariableFixing
