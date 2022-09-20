// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

//FIXME
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
template <typename Derived>
class RotationallySymmetric : public MarkAsAnalyticData {
 public:
  RotationallySymmetric() = default;

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
};
}  // namespace RadiationTransport::M1Grey::AnalyticData
