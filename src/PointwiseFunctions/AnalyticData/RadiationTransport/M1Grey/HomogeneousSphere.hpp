// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
class DataVector;
/// \endcond

namespace RadiationTransport::M1Grey::AnalyticData {
/*!
 * \brief FIXME
 */
template <typename Derived>
class HomogeneousSphereImpl : public MarkAsAnalyticData {
 public:
  /// The sphere radius.
  struct Radius {
    using type = double;
    static constexpr Options::String help = "Sphere radius";
  };

  /// The emissivity and absorption opacity.
  struct EmissivityAndOpacity {
    using type = double;
    static constexpr Options::String help = "Emissivity and absorption opacity";
  };

  /// The radius of the empty sphere
  struct OuterRadius {
    using type = double;
    static constexpr Options::String help = "Radius or outer absorption region";
  };

  /// The absorption opacity of the exterior
  struct OuterOpacity {
    using type = double;
    static constexpr Options::String help =
        "Opacity of outer absorption region";
  };

  using options =
      tmpl::list<Radius, EmissivityAndOpacity, OuterRadius, OuterOpacity>;

  HomogeneousSphereImpl() = default;
  explicit HomogeneousSphereImpl(CkMigrateMessage* /*message*/) {}

  HomogeneousSphereImpl(double radius, double emissivity_and_opacity,
                        double outer_radius, double outer_opacity);

  /// @{
  /// Retrieve fluid and neutrino variables
  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x,
                 tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<
                     Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeE<
          Frame::Inertial, NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x,
                 tmpl::list<RadiationTransport::M1Grey::Tags::TildeS<
                     Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeS<
          Frame::Inertial, NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x,
                 tmpl::list<RadiationTransport::M1Grey::Tags::GreyEmissivity<
                     NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<
          RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(
      const tnsr::I<DataVector, 3>& x,
      tmpl::list<RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<
          NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                                 GreyAbsorptionOpacity<NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(
      const tnsr::I<DataVector, 3>& x,
      tmpl::list<RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<
          NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                                 GreyScatteringOpacity<NeutrinoSpecies>>;

  auto variables(const tnsr::I<DataVector, 3>& x,
                 tmpl::list<hydro::Tags::LorentzFactor<DataVector>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataVector>>;

  auto variables(
      const tnsr::I<DataVector, 3>& x,
      tmpl::list<hydro::Tags::SpatialVelocity<DataVector, 3>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataVector, 3>>;
  /// @}

  /// Retrieve the metric variables
  template <typename Tag>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataVector, 3>& x,
                                     tmpl::list<Tag> /*meta*/) const {
    return gr::Solutions::Minkowski<3>{}.variables(x, 0.0, tmpl::list<Tag>{});
  }

  /// Retrieve a collection of variables
  template <typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataVector, 3>& x,
                                         tmpl::list<Tags...> /*meta*/) const {
    return {get<Tags>(variables(x, tmpl::list<Tags>{}))...};
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  friend bool operator==(const HomogeneousSphereImpl& lhs,
                         const HomogeneousSphereImpl& rhs) {
    return lhs.radius_ == rhs.radius_ and
           lhs.emissivity_and_opacity_ == rhs.emissivity_and_opacity_ and
           lhs.outer_radius_ == rhs.outer_radius_ and
           lhs.outer_opacity_ == rhs.outer_opacity_;
  }

  double radius_ = std::numeric_limits<double>::signaling_NaN();
  double emissivity_and_opacity_ = std::numeric_limits<double>::signaling_NaN();
  double outer_radius_ = std::numeric_limits<double>::signaling_NaN();
  double outer_opacity_ = std::numeric_limits<double>::signaling_NaN();
};

template <typename Derived>
bool operator!=(const HomogeneousSphereImpl<Derived>& lhs,
                const HomogeneousSphereImpl<Derived>& rhs);

// FIXME doc
class HomogeneousCylinder : public HomogeneousSphereImpl<HomogeneousCylinder> {
 public:
  static constexpr Options::String help = {
      "A homogeneous cylinder along the z axis emitting and absorbing "
      "neutrinos."};

  using HomogeneousSphereImpl<HomogeneousCylinder>::HomogeneousSphereImpl;
  // FIXME friend private?
  static auto radius_squared(const tnsr::I<DataVector, 3>& x);
};

// FIXME doc
class HomogeneousSphere : public HomogeneousSphereImpl<HomogeneousSphere> {
 public:
  static constexpr Options::String help = {
      "A homogeneous sphere emitting and absorbing neutrinos."};

  using HomogeneousSphereImpl<HomogeneousSphere>::HomogeneousSphereImpl;
  static auto radius_squared(const tnsr::I<DataVector, 3>& x);
};
}  // namespace RadiationTransport::M1Grey::AnalyticData
