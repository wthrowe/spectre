// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/HomogeneousSphere.hpp"

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
template <typename Derived>
HomogeneousSphereImpl<Derived>::HomogeneousSphereImpl(
    const double radius, const double emissivity_and_opacity,
    const double outer_radius, const double outer_opacity)
    : radius_(radius),
      emissivity_and_opacity_(emissivity_and_opacity),
      outer_radius_(outer_radius),
      outer_opacity_(outer_opacity) {}

template <typename Derived>
template <typename NeutrinoSpecies>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<
        Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeE<
        Frame::Inertial, NeutrinoSpecies>> {
  const DataVector r = sqrt(Derived::radius_squared(x));
  return Scalar<DataVector>{
      ((outer_radius_ - r) * step_function(outer_radius_ - r) -
       (radius_ - r) * step_function(radius_ - r)) /
          (outer_radius_ - radius_) +
      1.0e-12};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::TildeS<
        Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeS<
        Frame::Inertial, NeutrinoSpecies>> {
  return {make_with_value<tnsr::i<DataVector, 3, Frame::Inertial>>(x, 0.0)};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyEmissivity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<
        RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>> {
  return {Scalar<DataVector>{
      emissivity_and_opacity_ *
      step_function(square(radius_) - Derived::radius_squared(x))}};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                               GreyAbsorptionOpacity<NeutrinoSpecies>> {
  return {Scalar<DataVector>{
      emissivity_and_opacity_ *
          step_function(square(radius_) - Derived::radius_squared(x)) +
      outer_opacity_ *
          step_function(Derived::radius_squared(x) - square(outer_radius_))}};
}

template <typename Derived>
template <typename NeutrinoSpecies>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<
        NeutrinoSpecies>> /*meta*/) const
    -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                               GreyScatteringOpacity<NeutrinoSpecies>> {
  return {make_with_value<Scalar<DataVector>>(x, 0.0)};
}

template <typename Derived>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<hydro::Tags::LorentzFactor<DataVector>> /*meta*/) const
    -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataVector>> {
  return {make_with_value<Scalar<DataVector>>(x, 1.0)};
}

template <typename Derived>
auto HomogeneousSphereImpl<Derived>::variables(
    const tnsr::I<DataVector, 3>& x,
    tmpl::list<hydro::Tags::SpatialVelocity<DataVector, 3>> /*meta*/) const
    -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataVector, 3>> {
  return {make_with_value<tnsr::I<DataVector, 3, Frame::Inertial>>(x, 0.0)};
}

template <typename Derived>
void HomogeneousSphereImpl<Derived>::pup(PUP::er& p) {
  p | radius_;
  p | emissivity_and_opacity_;
  p | outer_radius_;
  p | outer_opacity_;
}

template <typename Derived>
bool operator!=(const HomogeneousSphereImpl<Derived>& lhs,
                const HomogeneousSphereImpl<Derived>& rhs) {
  return not(lhs == rhs);
}

auto HomogeneousCylinder::radius_squared(const tnsr::I<DataVector, 3>& x) {
  return square(get<0>(x)) + square(get<1>(x));
}

auto HomogeneousSphere::radius_squared(const tnsr::I<DataVector, 3>& x) {
  return square(get<0>(x)) + square(get<1>(x)) + square(get<2>(x));
}

#define DERIVED_CLASSES (HomogeneousCylinder, HomogeneousSphere)

#define DERIVED(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)
#define NTYPE(data) BOOST_PP_TUPLE_ELEM(2, data)
#define EBIN(data) BOOST_PP_TUPLE_ELEM(3, data)
#define GENERATE_LIST(z, n, _) BOOST_PP_COMMA_IF(n) n

#define INSTANTIATE_M1_FUNCTION_WITH_FRAME(_, data)                            \
  template tuples::TaggedTuple<TAG(data) < Frame::Inertial,                    \
                               NTYPE(data) < EBIN(data)> >>                    \
      HomogeneousSphereImpl<DERIVED(data)>::variables(                         \
          const tnsr::I<DataVector, 3>& x,                                     \
          tmpl::list<TAG(data) < Frame::Inertial, NTYPE(data) < EBIN(data)> >> \
          /*meta*/) const;

#define temp_list \
  (BOOST_PP_REPEAT(MAX_NUMBER_OF_NEUTRINO_ENERGY_BINS, GENERATE_LIST, _))

GENERATE_INSTANTIATIONS(INSTANTIATE_M1_FUNCTION_WITH_FRAME, DERIVED_CLASSES,
                        (RadiationTransport::M1Grey::Tags::TildeE,
                         RadiationTransport::M1Grey::Tags::TildeS),
                        (neutrinos::ElectronNeutrinos,
                         neutrinos::ElectronAntiNeutrinos,
                         neutrinos::HeavyLeptonNeutrinos),
                        temp_list)

#undef temp_list
#undef INSTANTIATE_M1_FUNCTION_WITH_FRAME
#undef DERIVED
#undef TAG
#undef NTYPE
#undef EBIN
#undef GENERATE_LIST

#define DERIVED(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)
#define NTYPE(data) BOOST_PP_TUPLE_ELEM(2, data)
#define EBIN(data) BOOST_PP_TUPLE_ELEM(3, data)
#define GENERATE_LIST(z, n, _) BOOST_PP_COMMA_IF(n) n

#define INSTANTIATE_M1_FUNCTION(_, data)                                \
  template tuples::TaggedTuple<TAG(data) < NTYPE(data) < EBIN(data)> >> \
      HomogeneousSphereImpl<DERIVED(data)>::variables(                  \
          const tnsr::I<DataVector, 3>& x,                              \
          tmpl::list<TAG(data) < NTYPE(data) < EBIN(data)> >>           \
          /*meta*/) const;

#define temp_list \
  (BOOST_PP_REPEAT(MAX_NUMBER_OF_NEUTRINO_ENERGY_BINS, GENERATE_LIST, _))

GENERATE_INSTANTIATIONS(
    INSTANTIATE_M1_FUNCTION, DERIVED_CLASSES,
    (RadiationTransport::M1Grey::Tags::GreyEmissivity,
     RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity,
     RadiationTransport::M1Grey::Tags::GreyScatteringOpacity),
    (neutrinos::ElectronNeutrinos, neutrinos::ElectronAntiNeutrinos,
     neutrinos::HeavyLeptonNeutrinos),
    temp_list)

#undef INSTANTIATE_M1_FUNCTION
#undef temp_list
#undef DERIVED
#undef TAG
#undef NTYPE
#undef EBIN
#undef GENERATE_LIST

#define DERIVED(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) \
  template class HomogeneousSphereImpl<DERIVED(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, DERIVED_CLASSES)

#undef INSTANTIATE
#undef DERIVED

#undef DERIVED_CLASSES
}  // namespace RadiationTransport::M1Grey::AnalyticData
