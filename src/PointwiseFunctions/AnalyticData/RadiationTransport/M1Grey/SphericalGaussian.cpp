// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/RotationallySymmetric.tpp"

namespace RadiationTransport::M1Grey::AnalyticData {
template <size_t Dim>
SphericalGaussian<Dim>::SphericalGaussian(
    const double radius, const double emissivity_and_opacity_amplitude,
    const double outer_radius, const double outer_opacity)
    : radius_(radius),
      emissivity_and_opacity_amplitude_(emissivity_and_opacity_amplitude),
      outer_radius_(outer_radius),
      outer_opacity_(outer_opacity) {}

template <size_t Dim>
void SphericalGaussian<Dim>::pup(PUP::er& p) {
  p | radius_;
  p | emissivity_and_opacity_amplitude_;
  p | outer_radius_;
  p | outer_opacity_;
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::energy_profile(const DataVector& r) const {
  // Just a smooth function.  Not particularly close to the actual
  // solution.
  return exp(-square(1.0 / radius_) * square(r));
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::emission_profile(const DataVector& r) const {
  return emissivity_and_opacity_amplitude_ *
         exp(-square(1.0 / radius_) * square(r));
}

template <size_t Dim>
DataVector SphericalGaussian<Dim>::absorption_profile(const DataVector& r) const {
  return emissivity_and_opacity_amplitude_ *
             exp(-square(1.0 / radius_) * square(r)) +
         outer_opacity_ * step_function(r - outer_radius_);
}

template class SphericalGaussian<2>;
template class SphericalGaussian<3>;
INSTANTIATE_ROTATIONALLY_SYMMETRIC(SphericalGaussian<2>)
INSTANTIATE_ROTATIONALLY_SYMMETRIC(SphericalGaussian<3>)
}  // namespace RadiationTransport::M1Grey::AnalyticData
