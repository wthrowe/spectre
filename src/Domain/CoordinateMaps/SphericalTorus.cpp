// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/SphericalTorus.hpp"

#include <pup.h>

#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace domain::CoordinateMaps {
SphericalTorus::SphericalTorus(double r_min, double r_max, double phi_max,
                               double fraction_of_torus,
                               const Options::Context& context)
    : r_min_(r_min),
      r_max_(r_max),
      phi_max_(phi_max),
      fraction_of_torus_(fraction_of_torus) {
  if (r_min_ <= 0.0) {
    PARSE_ERROR(context, "Minimum radius must be positive.");
  }
  if (r_max_ <= r_min_) {
    PARSE_ERROR(context, "Maximum radius must be greater than minimum radius.");
  }
  if (phi_max_ <= 0.0) {
    PARSE_ERROR(context, "Polar extent must be positive.");
  }
  if (phi_max_ >= 0.5 * M_PI) {
    PARSE_ERROR(context, "Cannot cover z-axis with a torus.");
  }
  if (fraction_of_torus_ <= 0.0) {
    PARSE_ERROR(context, "Fraction of torus included must be positive.");
  }
  ASSERT(fraction_of_torus_ <= 1.0,
         "Fraction of torus included must be at most 1.");
}

template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, 3> SphericalTorus::operator()(
    const std::array<T, 3>& source_coords) const noexcept {
  const auto r = radius(source_coords[0]);
  const auto theta = M_PI * fraction_of_torus_ * source_coords[2];
  const auto phi = phi_max_ * source_coords[1];

  return {{r * cos(theta) * cos(phi), r * sin(theta) * cos(phi), r * sin(phi)}};
}

std::optional<std::array<double, 3>> SphericalTorus::inverse(
    const std::array<double, 3>& target_coords) const noexcept {
  const double r =
      std::hypot(target_coords[0], target_coords[1], target_coords[2]);
  const double theta = std::atan2(target_coords[1], target_coords[0]);
  const double phi = std::atan2(target_coords[2],
                                std::hypot(target_coords[0], target_coords[1]));

  return {
      {radius_inverse(r), phi / phi_max_, theta / (M_PI * fraction_of_torus_)}};
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
SphericalTorus::jacobian(const std::array<T, 3>& source_coords) const noexcept {
  using UnwrappedT = tt::remove_cvref_wrap_t<T>;
  auto jacobian = make_with_value<tnsr::Ij<UnwrappedT, 3, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);

  const auto r = radius(source_coords[0]);
  const auto theta = M_PI * fraction_of_torus_ * source_coords[2];
  const auto phi = phi_max_ * source_coords[1];
  const UnwrappedT cos_theta = cos(theta);
  const UnwrappedT sin_theta = sin(theta);
  const UnwrappedT cos_phi = cos(phi);
  const UnwrappedT sin_phi = sin(phi);

  get<0, 0>(jacobian) = 0.5 * (r_max_ - r_min_) * cos_theta * cos_phi;
  get<0, 1>(jacobian) = -phi_max_ * r * cos_theta * sin_phi;
  get<0, 2>(jacobian) = -M_PI * fraction_of_torus_ * r * sin_theta * cos_phi;
  get<1, 0>(jacobian) = 0.5 * (r_max_ - r_min_) * sin_theta * cos_phi;
  get<1, 1>(jacobian) = -phi_max_ * r * sin_theta * sin_phi;
  get<1, 2>(jacobian) = M_PI * fraction_of_torus_ * r * cos_theta * cos_phi;
  get<2, 0>(jacobian) = 0.5 * (r_max_ - r_min_) * sin_phi;
  get<2, 1>(jacobian) = phi_max_ * r * cos_phi;
  return jacobian;
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
SphericalTorus::inv_jacobian(
    const std::array<T, 3>& source_coords) const noexcept {
  using UnwrappedT = tt::remove_cvref_wrap_t<T>;
  auto inv_jacobian = make_with_value<tnsr::Ij<UnwrappedT, 3, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);

  const auto r = radius(source_coords[0]);
  const auto theta = M_PI * fraction_of_torus_ * source_coords[2];
  const auto phi = phi_max_ * source_coords[1];
  const UnwrappedT cos_theta = cos(theta);
  const UnwrappedT sin_theta = sin(theta);
  const UnwrappedT cos_phi = cos(phi);
  const UnwrappedT sin_phi = sin(phi);

  get<0, 0>(inv_jacobian) = 2.0 / (r_max_ - r_min_) * cos_theta * cos_phi;
  get<0, 1>(inv_jacobian) = 2.0 / (r_max_ - r_min_) * sin_theta * cos_phi;
  get<0, 2>(inv_jacobian) = 2.0 / (r_max_ - r_min_) * sin_phi;
  get<1, 0>(inv_jacobian) = -(1.0 / phi_max_) * cos_theta * sin_phi / r;
  get<1, 1>(inv_jacobian) = -(1.0 / phi_max_) * sin_theta * sin_phi / r;
  get<1, 2>(inv_jacobian) = (1.0 / phi_max_) * cos_phi / r;
  get<2, 0>(inv_jacobian) =
      -(1.0 / (M_PI * fraction_of_torus_)) * sin_theta / (r * cos_phi);
  get<2, 1>(inv_jacobian) =
      (1.0 / (M_PI * fraction_of_torus_)) * cos_theta / (r * cos_phi);
  return inv_jacobian;
}

template <typename T>
tnsr::Ijj<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
SphericalTorus::hessian(const std::array<T, 3>& source_coords) const noexcept {
  using UnwrappedT = tt::remove_cvref_wrap_t<T>;
  auto hessian = make_with_value<tnsr::Ijj<UnwrappedT, 3, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);

  const double r_factor = 0.5 * (r_max_ - r_min_);
  const double theta_factor = M_PI * fraction_of_torus_;
  const double phi_factor = phi_max_;
  const auto r = radius(source_coords[0]);
  const auto theta = theta_factor * source_coords[2];
  const auto phi = phi_factor * source_coords[1];
  const UnwrappedT cos_theta = cos(theta);
  const UnwrappedT sin_theta = sin(theta);
  const UnwrappedT cos_phi = cos(phi);
  const UnwrappedT sin_phi = sin(phi);

  get<0, 0, 1>(hessian) = -r_factor * phi_factor * cos_theta * sin_phi;
  get<0, 0, 2>(hessian) = -r_factor * theta_factor * sin_theta * cos_phi;
  get<0, 1, 1>(hessian) = -square(phi_factor) * r * cos_theta * cos_phi;
  get<0, 1, 2>(hessian) = phi_factor * theta_factor * r * sin_theta * sin_phi;
  get<0, 2, 2>(hessian) = -square(theta_factor) * r * cos_theta * cos_phi;
  get<1, 0, 1>(hessian) = -r_factor * phi_factor * sin_theta * sin_phi;
  get<1, 0, 2>(hessian) = r_factor * theta_factor * cos_theta * cos_phi;
  get<1, 1, 1>(hessian) = -square(phi_factor) * r * sin_theta * cos_phi;
  get<1, 1, 2>(hessian) = -phi_factor * theta_factor * r * cos_theta * sin_phi;
  get<1, 2, 2>(hessian) = -square(theta_factor) * r * sin_theta * cos_phi;
  get<2, 0, 1>(hessian) = r_factor * phi_factor * cos_phi;
  get<2, 1, 1>(hessian) = -square(phi_factor) * r * sin_phi;
  return hessian;
}

template <typename T>
tnsr::Ijk<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
SphericalTorus::derivative_of_inv_jacobian(
    const std::array<T, 3>& source_coords) const noexcept {
  using UnwrappedT = tt::remove_cvref_wrap_t<T>;
  auto result = make_with_value<tnsr::Ijk<UnwrappedT, 3, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);

  const double r_factor = 0.5 * (r_max_ - r_min_);
  const double theta_factor = M_PI * fraction_of_torus_;
  const double phi_factor = phi_max_;
  const auto r = radius(source_coords[0]);
  const auto theta = M_PI * fraction_of_torus_ * source_coords[2];
  const auto phi = phi_max_ * source_coords[1];
  const UnwrappedT cos_theta = cos(theta);
  const UnwrappedT sin_theta = sin(theta);
  const UnwrappedT cos_phi = cos(phi);
  const UnwrappedT sin_phi = sin(phi);

  get<0, 0, 1>(result) = -phi_factor / r_factor * cos_theta * sin_phi;
  get<0, 0, 2>(result) = -theta_factor / r_factor * sin_theta * cos_phi;
  get<0, 1, 1>(result) = -phi_factor / r_factor * sin_theta * sin_phi;
  get<0, 1, 2>(result) = theta_factor / r_factor * cos_theta * cos_phi;
  get<0, 2, 1>(result) = phi_factor / r_factor * cos_phi;
  get<1, 0, 0>(result) =
      r_factor / phi_factor * cos_theta * sin_phi / square(r);
  get<1, 0, 1>(result) = -cos_theta * cos_phi / r;
  get<1, 0, 2>(result) = theta_factor / phi_factor * sin_theta * sin_phi / r;
  get<1, 1, 0>(result) =
      r_factor / phi_factor * sin_theta * sin_phi / square(r);
  get<1, 1, 1>(result) = -sin_theta * cos_phi / r;
  get<1, 1, 2>(result) = -theta_factor / phi_factor * cos_theta * sin_phi / r;
  get<1, 2, 0>(result) = -r_factor / phi_factor * cos_phi / square(r);
  get<1, 2, 1>(result) = -sin_phi / r;
  get<2, 0, 0>(result) =
      r_factor / theta_factor * sin_theta / (square(r) * cos_phi);
  get<2, 0, 1>(result) =
      -phi_factor / theta_factor * sin_theta * sin_phi / (r * square(cos_phi));
  get<2, 0, 2>(result) = -cos_theta / (r * cos_phi);

  get<2, 1, 0>(result) =
      -r_factor / theta_factor * cos_theta / (square(r) * cos_phi);
  get<2, 1, 1>(result) =
      phi_factor / theta_factor * cos_theta * sin_phi / (r * square(cos_phi));
  get<2, 1, 2>(result) = -sin_theta / (r * cos_phi);
  return result;
}

void SphericalTorus::pup(PUP::er& p) noexcept {
  p | r_min_;
  p | r_max_;
  p | phi_max_;
  p | fraction_of_torus_;
}

template <typename T>
tt::remove_cvref_wrap_t<T> SphericalTorus::radius(const T& x) const noexcept {
  return 0.5 * r_min_ * (1.0 - x) + 0.5 * r_max_ * (1.0 + x);
}

template <typename T>
tt::remove_cvref_wrap_t<T> SphericalTorus::radius_inverse(
    const T& r) const noexcept {
  return ((r - r_min_) - (r_max_ - r)) / (r_max_ - r_min_);
}

bool operator==(const SphericalTorus& lhs, const SphericalTorus& rhs) noexcept {
  return lhs.r_min_ == rhs.r_min_ and lhs.r_max_ == rhs.r_max_ and
         lhs.phi_max_ == rhs.phi_max_ and
         lhs.fraction_of_torus_ == rhs.fraction_of_torus_;
}

bool operator!=(const SphericalTorus& lhs, const SphericalTorus& rhs) noexcept {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                  \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, 3>                \
  SphericalTorus::operator()(const std::array<DTYPE(data), 3>& source_coords) \
      const noexcept;                                                         \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame>  \
  SphericalTorus::jacobian(const std::array<DTYPE(data), 3>& source_coords)   \
      const noexcept;                                                         \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame>  \
  SphericalTorus::inv_jacobian(                                               \
      const std::array<DTYPE(data), 3>& source_coords) const noexcept;        \
  template tnsr::Ijj<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame> \
  SphericalTorus::hessian(const std::array<DTYPE(data), 3>& source_coords)    \
      const noexcept;                                                         \
  template tnsr::Ijk<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame> \
  SphericalTorus::derivative_of_inv_jacobian(                                 \
      const std::array<DTYPE(data), 3>& source_coords) const noexcept;

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector,
                                      std::reference_wrapper<const double>,
                                      std::reference_wrapper<const DataVector>))
#undef INSTANTIATE
#undef DTYPE

}  // namespace domain::CoordinateMaps
