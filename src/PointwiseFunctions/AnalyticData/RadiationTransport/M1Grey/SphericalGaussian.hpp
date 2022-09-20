// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

//FIXME
#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/RotationallySymmetric.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
class DataVector;
/// \endcond

namespace RadiationTransport::M1Grey::AnalyticData {
// FIXME doc
/// \kappa_0 exp(-(r/r_0)^2)
template <size_t Dim>
class SphericalGaussian : public RotationallySymmetric<SphericalGaussian<Dim>> {
  static_assert(Dim == 2 or Dim == 3);
 public:
  static std::string name() {
    return Dim == 3 ? "SphericalGaussian" : "CylindricalGaussian";
  }

  /// The Gaussian width \f$r_0\f$.
  struct Radius {
    using type = double;
    static constexpr Options::String help = "Gaussian radius";
  };

  /// The emissivity and absorption opacity amplitude \f$\kappa_0\f$.
  struct EmissivityAndOpacityAmplitude {
    using type = double;
    static constexpr Options::String help =
        "Amplitude of emissivity and absorption opacity Gaussian";
  };

  /// The radius of the absorbing boundary region.
  struct OuterRadius {
    using type = double;
    static constexpr Options::String help = "Radius of outer absorption region";
  };

  /// The absorption opacity of the outer region.
  struct OuterOpacity {
    using type = double;
    static constexpr Options::String help =
        "Opacity of outer absorption region";
  };

  static constexpr Options::String help = {
      "A Gaussian profile emitting and absorbing neutrinos."};
  using options = tmpl::list<Radius, EmissivityAndOpacityAmplitude, OuterRadius,
                             OuterOpacity>;

  SphericalGaussian() = default;
  explicit SphericalGaussian(CkMigrateMessage* /*message*/) {}

  SphericalGaussian(double radius, double emissivity_and_opacity_amplitude,
                    double outer_radius, double outer_opacity);

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  friend class RotationallySymmetric<SphericalGaussian<Dim>>;

  static constexpr bool spherical = Dim == 3;

  DataVector energy_profile(const DataVector& r) const;
  DataVector emission_profile(const DataVector& r) const;
  DataVector absorption_profile(const DataVector& r) const;

  double radius_;
  double emissivity_and_opacity_amplitude_;
  double outer_radius_;
  double outer_opacity_;
};
}  // namespace RadiationTransport::M1Grey::AnalyticData
