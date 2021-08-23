// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/Options.hpp"
#include "Utilities/TypeTraits/RemoveReferenceWrapper.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace domain::CoordinateMaps {
/*!
 * \ingroup CoordinateMapsGroup
 * Torus made by removing polar cones from a spherical shell, mapping the
 * source coordinates \f$(x, y, z)\f$ as \f$(r, \phi, \theta)\f$, scaled to
 * \f$[-1,1]\f$.
 */
class SphericalTorus {
 public:
  static constexpr size_t dim = 3;

  struct RadialRange {
    using type = std::array<double, 2>;
    static constexpr Options::String help = "Radial extent of the torus.";
  };

  struct PhiMax {
    using type = double;
    static constexpr Options::String help = "Angular height of the torus.";
    static type lower_bound() noexcept { return 0.0; }
    static type upper_bound() noexcept { return 0.5 * M_PI; }
  };

  struct FractionOfTorus {
    using type = double;
    static constexpr Options::String help = "Fraction of orbit covered.";
    static type lower_bound() noexcept { return 0.0; }
    static type upper_bound() noexcept { return 1.0; }
  };

  static constexpr Options::String help =
      "Torus made by removing polar cones from a spherical shell, mapping the\n"
      "source coordinates (x, y, z) as (r, phi, theta), scaled to [-1,1].";

  using options = tmpl::list<RadialRange, PhiMax, FractionOfTorus>;

  SphericalTorus(const std::array<double, 2>& radial_range,
                 const double phi_max, const double fraction_of_torus,
                 const Options::Context& context = {})
      : SphericalTorus(radial_range[0], radial_range[1], phi_max,
                       fraction_of_torus, context) {}

  SphericalTorus(double r_min, double r_max, double phi_max,
                 double fraction_of_torus = 1.0,
                 const Options::Context& context = {});

  SphericalTorus() = default;

  template <typename T>
  std::array<tt::remove_cvref_wrap_t<T>, 3> operator()(
      const std::array<T, 3>& source_coords) const noexcept;

  std::optional<std::array<double, 3>> inverse(
      const std::array<double, 3>& target_coords) const noexcept;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> jacobian(
      const std::array<T, 3>& source_coords) const noexcept;

  template <typename T>
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> inv_jacobian(
      const std::array<T, 3>& source_coords) const noexcept;

  template <typename T>
  tnsr::Ijj<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> hessian(
      const std::array<T, 3>& source_coords) const noexcept;

  template <typename T>
  tnsr::Ijk<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
  derivative_of_inv_jacobian(
      const std::array<T, 3>& source_coords) const noexcept;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) noexcept;

  bool is_identity() const noexcept { return false; }

 private:
  template <typename T>
  tt::remove_cvref_wrap_t<T> radius(const T& x) const noexcept;

  template <typename T>
  tt::remove_cvref_wrap_t<T> radius_inverse(const T& x) const noexcept;

  friend bool operator==(const SphericalTorus& lhs,
                         const SphericalTorus& rhs) noexcept;

  double r_min_ = std::numeric_limits<double>::signaling_NaN();
  double r_max_ = std::numeric_limits<double>::signaling_NaN();
  double phi_max_ = std::numeric_limits<double>::signaling_NaN();
  double fraction_of_torus_ = std::numeric_limits<double>::signaling_NaN();
};

bool operator!=(const SphericalTorus& lhs, const SphericalTorus& rhs) noexcept;
}  // namespace domain::CoordinateMaps
