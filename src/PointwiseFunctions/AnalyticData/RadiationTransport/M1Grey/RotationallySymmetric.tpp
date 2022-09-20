// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/RotationallySymmetric.hpp"

//FIXME
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <cmath>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/RadiationTransport/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/Math.hpp"

namespace RadiationTransport::M1Grey::AnalyticData {
namespace {
template <bool Spherical>
auto radius(const tnsr::I<DataVector, 3>& x) {
  if constexpr (Spherical) {
    return sqrt(square(get<0>(x)) + square(get<1>(x)) + square(get<2>(x)));
  } else {
    return sqrt(square(get<0>(x)) + square(get<1>(x)));
  }
}
}  // namespace

template <typename Derived>
template <typename NeutrinoSpecies>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<
        Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeE<
        Frame::Inertial, NeutrinoSpecies>> {
  return {Scalar<DataVector>{static_cast<const Derived&>(*this).energy_profile(radius<Derived::spherical>(x))}};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::TildeS<
        Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeS<
        Frame::Inertial, NeutrinoSpecies>> {
  return {make_with_value<tnsr::i<DataVector, 3, Frame::Inertial>>(x, 0.0)};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyEmissivity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<
        RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>> {
  return {Scalar<DataVector>{static_cast<const Derived&>(*this).emission_profile(radius<Derived::spherical>(x))}};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                               GreyAbsorptionOpacity<NeutrinoSpecies>> {
  return {Scalar<DataVector>{
      static_cast<const Derived&>(*this).absorption_profile(radius<Derived::spherical>(x))}};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                               GreyScatteringOpacity<NeutrinoSpecies>> {
  return {make_with_value<Scalar<DataVector>>(x, 0.0)};
}

template <typename Derived>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<hydro::Tags::LorentzFactor<DataVector>> /*meta*/) const
    -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataVector>> {
  return {make_with_value<Scalar<DataVector>>(x, 1.0)};
}

template <typename Derived>
auto RotationallySymmetric<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<hydro::Tags::SpatialVelocity<DataVector, 3>> /*meta*/) const
    -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataVector, 3>> {
  return {make_with_value<tnsr::I<DataVector, 3, Frame::Inertial>>(x, 0.0)};
}
}  // namespace RadiationTransport::M1Grey::AnalyticData

#define ROTATIONALLY_SYMMETRIC_DERIVED(data) BOOST_PP_TUPLE_ELEM(0, data)
#define ROTATIONALLY_SYMMETRIC_TAG(data) BOOST_PP_TUPLE_ELEM(1, data)
#define ROTATIONALLY_SYMMETRIC_NTYPE(data) BOOST_PP_TUPLE_ELEM(2, data)
#define ROTATIONALLY_SYMMETRIC_EBIN(data) BOOST_PP_TUPLE_ELEM(3, data)
#define ROTATIONALLY_SYMMETRIC_GENERATE_LIST(z, n, _) BOOST_PP_COMMA_IF(n) n

#define ROTATIONALLY_SYMMETRIC_NEUTRINO_BINS \
  (BOOST_PP_REPEAT(MAX_NUMBER_OF_NEUTRINO_ENERGY_BINS,  \
                   ROTATIONALLY_SYMMETRIC_GENERATE_LIST, _))
#define ROTATIONALLY_SYMMETRIC_INSTANTIATE_M1_FUNCTION_WITH_FRAME(_, data)                            \
  template tuples::TaggedTuple<ROTATIONALLY_SYMMETRIC_TAG(data) < Frame::Inertial,                    \
                               ROTATIONALLY_SYMMETRIC_NTYPE(data) < ROTATIONALLY_SYMMETRIC_EBIN(data)> >>                    \
  RadiationTransport::M1Grey::AnalyticData::RotationallySymmetric<ROTATIONALLY_SYMMETRIC_DERIVED(data)>::variables( \
          const tnsr::I<DataVector, 3>& x,                                     \
          tmpl::list<ROTATIONALLY_SYMMETRIC_TAG(data) < Frame::Inertial, ROTATIONALLY_SYMMETRIC_NTYPE(data) < ROTATIONALLY_SYMMETRIC_EBIN(data)> >> \
          /*meta*/) const;

#define ROTATIONALLY_SYMMETRIC_INSTANTIATE_M1_FUNCTION(_, data)                                \
  template tuples::TaggedTuple<ROTATIONALLY_SYMMETRIC_TAG(data) < ROTATIONALLY_SYMMETRIC_NTYPE(data) < ROTATIONALLY_SYMMETRIC_EBIN(data)> >> \
      RadiationTransport::M1Grey::AnalyticData::RotationallySymmetric<ROTATIONALLY_SYMMETRIC_DERIVED(data)>::variables(                  \
          const tnsr::I<DataVector, 3>& x,                              \
          tmpl::list<ROTATIONALLY_SYMMETRIC_TAG(data) < ROTATIONALLY_SYMMETRIC_NTYPE(data) < ROTATIONALLY_SYMMETRIC_EBIN(data)> >>           \
          /*meta*/) const;

#define INSTANTIATE_ROTATIONALLY_SYMMETRIC(derived) \
GENERATE_INSTANTIATIONS(ROTATIONALLY_SYMMETRIC_INSTANTIATE_M1_FUNCTION_WITH_FRAME, (derived),\
                        (RadiationTransport::M1Grey::Tags::TildeE,\
                         RadiationTransport::M1Grey::Tags::TildeS),\
                        (neutrinos::ElectronNeutrinos,\
                         neutrinos::ElectronAntiNeutrinos,\
                         neutrinos::HeavyLeptonNeutrinos),\
                        ROTATIONALLY_SYMMETRIC_NEUTRINO_BINS)\
GENERATE_INSTANTIATIONS(\
    ROTATIONALLY_SYMMETRIC_INSTANTIATE_M1_FUNCTION, (derived),\
    (RadiationTransport::M1Grey::Tags::GreyEmissivity,\
     RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity,\
     RadiationTransport::M1Grey::Tags::GreyScatteringOpacity),\
    (neutrinos::ElectronNeutrinos, neutrinos::ElectronAntiNeutrinos,\
     neutrinos::HeavyLeptonNeutrinos),\
    ROTATIONALLY_SYMMETRIC_NEUTRINO_BINS) \
template class RadiationTransport::M1Grey::AnalyticData::RotationallySymmetric<derived>;

