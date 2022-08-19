#include <iostream>// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/RadiationTransport/M1Grey/M1HydroCoupling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

// IWYU pragma: no_forward_declare Tensor

namespace {
struct densitized_eta_minus_kappaJ : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct kappaT_lapse : db::SimpleTag {
  using type = Scalar<DataVector>;
};
}  // namespace

namespace RadiationTransport::M1Grey::detail {

void compute_m1_hydro_coupling_impl(
    const gsl::not_null<Scalar<DataVector>*> source_n,
    const gsl::not_null<tnsr::i<DataVector, 3>*> source_i,
    const Scalar<DataVector>& emissivity,
    const Scalar<DataVector>& absorption_opacity,
    const Scalar<DataVector>& scattering_opacity,
    const Scalar<DataVector>& comoving_energy_density,
    const Scalar<DataVector>& comoving_momentum_density_normal,
    const tnsr::i<DataVector, 3>& comoving_momentum_density_spatial,
    const tnsr::I<DataVector, 3>& fluid_velocity,
    const Scalar<DataVector>& fluid_lorentz_factor,
    const Scalar<DataVector>& lapse,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric) {
  // std::cerr.precision(17);
  // std::cerr << "Coupling: " << comoving_energy_density << "\n";
  Variables<tmpl::list<hydro::Tags::SpatialVelocityOneForm<DataVector, 3>,
                       densitized_eta_minus_kappaJ, kappaT_lapse>>
      temp_tensors(get(lapse).size());
  // Dimension of spatial tensors
  constexpr size_t spatial_dim = 3;

  auto& dens_e_minus_kJ = get<densitized_eta_minus_kappaJ>(temp_tensors);
  get(dens_e_minus_kJ) =
      get(lapse) * get(fluid_lorentz_factor) *
      (get(sqrt_det_spatial_metric) * get(emissivity) -
       get(absorption_opacity) * get(comoving_energy_density));
  auto& kT_lapse = get<kappaT_lapse>(temp_tensors);
  get(kT_lapse) =
      get(lapse) * (get(absorption_opacity) + get(scattering_opacity));
  auto& fluid_velocity_i =
      get<hydro::Tags::SpatialVelocityOneForm<DataVector, 3>>(temp_tensors);
  raise_or_lower_index(make_not_null(&fluid_velocity_i), fluid_velocity,
                       spatial_metric);

  get(*source_n) = get(dens_e_minus_kJ) +
                   get(kT_lapse) * get(comoving_momentum_density_normal);
  for (size_t i = 0; i < spatial_dim; i++) {
    source_i->get(i) = fluid_velocity_i.get(i) * get(dens_e_minus_kJ) -
                       get(kT_lapse) * comoving_momentum_density_spatial.get(i);
  }
}

namespace {
namespace LocalTags {
//FIXME clean up unused
struct DummySpecies;
using TildeE = Tags::TildeE<Frame::Inertial, DummySpecies>;
using TildeHSpatial = Tags::TildeHSpatial<Frame::Inertial, DummySpecies>;
using TildeJ = Tags::TildeJ<DummySpecies>;
using TildeS = Tags::TildeS<Frame::Inertial, DummySpecies>;
using TildeSVector = Tags::TildeSVector<Frame::Inertial>;
}  // namespace LocalTags
}  // namespace

void compute_m1_hydro_coupling_jacobian_impl(
    const gsl::not_null<Scalar<DataVector>*> deriv_e_source_e,
    const gsl::not_null<tnsr::i<DataVector, 3>*> deriv_e_source_s,
    const gsl::not_null<tnsr::I<DataVector, 3>*> deriv_s_source_e,
    const gsl::not_null<tnsr::Ij<DataVector, 3>*> deriv_s_source_s,
    const tnsr::i<DataVector, 3>& tilde_s,
    const Scalar<DataVector>& tilde_e,
    const Scalar<DataVector>& emissivity,
    const Scalar<DataVector>& absorption_opacity,
    const Scalar<DataVector>& scattering_opacity,
    const tnsr::I<DataVector, 3>& fluid_velocity,
    const Scalar<DataVector>& fluid_lorentz_factor,
    const Scalar<DataVector>& closure_factor,
    const Scalar<DataVector>& comoving_energy_density,
    const tnsr::i<DataVector, 3>& comoving_momentum_density_spatial,
    const Scalar<DataVector>& comoving_momentum_density_normal,
    const Scalar<DataVector>& lapse,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric) {
  const double s_squared_floor = 1.0e-150;

  //FIXME reorder
  Variables<tmpl::list<
      imex::Tags::Jacobian<LocalTags::TildeE, LocalTags::TildeJ>,
      imex::Tags::Jacobian<LocalTags::TildeS, LocalTags::TildeJ>,
      imex::Tags::Jacobian<LocalTags::TildeE, LocalTags::TildeHSpatial>,
      imex::Tags::Jacobian<LocalTags::TildeS, LocalTags::TildeHSpatial>,
      LocalTags::TildeSVector,
      ::Tags::TempScalar<0>,
      ::Tags::TempScalar<1>,
      ::Tags::TempScalar<2>,
      ::Tags::TempScalar<3>,
      ::Tags::TempScalar<4>,
      ::Tags::TempScalar<5>,
      ::Tags::TempScalar<6>,
      ::Tags::TempScalar<7>,
      ::Tags::TempScalar<8>,
      ::Tags::TempScalar<9>,
      ::Tags::TempScalar<10>,
      ::Tags::TempScalar<11>,
      ::Tags::TempScalar<12>,
      ::Tags::TempScalar<13>,
      ::Tags::TempScalar<14>,
      ::Tags::TempScalar<15>,
      ::Tags::TempScalar<16>,
      ::Tags::Tempi<0, 3>
>>
      temporaries(get(emissivity).size());

  auto& eddington_factor = get<::Tags::TempScalar<16>>(temporaries);
  tenex::evaluate(make_not_null(&eddington_factor),
                  1.0 / 3.0 + 2.0 / 15.0 * square(closure_factor()) * (3.0 + closure_factor() * (-1.0 + 3.0 * closure_factor())));

  auto& d_thick = get<::Tags::TempScalar<15>>(temporaries);
  tenex::evaluate(make_not_null(&d_thick), 1.5 * (1.0 - eddington_factor()));

  auto& d_thin = get<::Tags::TempScalar<14>>(temporaries);
  tenex::evaluate(make_not_null(&d_thin), 1.0 - d_thick());

  auto& total_opacity = get<::Tags::TempScalar<0>>(temporaries);
  tenex::evaluate(make_not_null(&total_opacity),
                  absorption_opacity() + scattering_opacity());

  auto& fluid_velocity_lower = get<::Tags::Tempi<0, 3>>(temporaries);
  tenex::evaluate<ti::i>(make_not_null(&fluid_velocity_lower),
                         spatial_metric(ti::i, ti::j) * fluid_velocity(ti::J));

  auto& tilde_s_upper = get<LocalTags::TildeSVector>(temporaries);
  tenex::evaluate<ti::I>(make_not_null(&tilde_s_upper),
                         inverse_spatial_metric(ti::I, ti::J) * tilde_s(ti::j));

  auto& inverse_s_norm = get<::Tags::TempScalar<1>>(temporaries);
  tenex::evaluate(
      make_not_null(&inverse_s_norm),
      1.0 / (tilde_s_upper(ti::I) * tilde_s(ti::i) + s_squared_floor));

  auto& fluid_velocity_norm = get<::Tags::TempScalar<2>>(temporaries);
  tenex::evaluate(make_not_null(&fluid_velocity_norm),
                  fluid_velocity(ti::I) * fluid_velocity_lower(ti::i));

  auto& s_dot_fluid_velocity = get<::Tags::TempScalar<3>>(temporaries);
  tenex::evaluate(make_not_null(&s_dot_fluid_velocity),
                  tilde_s(ti::i) * fluid_velocity(ti::I));

  auto& denom = get<::Tags::TempScalar<4>>(temporaries);
  tenex::evaluate(make_not_null(&denom),
                  1.0 / (1.0 + 2.0 * square(fluid_lorentz_factor())));

  auto& deriv_e_h_velocity_coef = get<::Tags::TempScalar<5>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_e_h_velocity_coef),
      -cube(fluid_lorentz_factor()) *
      (d_thin() * (1.0 + square(s_dot_fluid_velocity()) * inverse_s_norm()) +
       4.0 * d_thick() * denom()));

  auto& deriv_e_h_s_coef = get<::Tags::TempScalar<6>>(temporaries);
  tenex::evaluate(make_not_null(&deriv_e_h_s_coef),
                  -d_thin() * fluid_lorentz_factor() * s_dot_fluid_velocity() *
                      inverse_s_norm());

  auto& deriv_s_h_trace_coef = get<::Tags::TempScalar<7>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_s_h_trace_coef),
      fluid_lorentz_factor() *
          (1.0 -
           d_thin() * tilde_e() * s_dot_fluid_velocity() * inverse_s_norm() -
           d_thick() * fluid_velocity_norm()));

  auto& deriv_s_j_velocity_coef = get<::Tags::TempScalar<8>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_s_j_velocity_coef),
      - 2.0 * fluid_lorentz_factor() *
      (deriv_s_h_trace_coef() +
       d_thick() * fluid_lorentz_factor() * fluid_velocity_norm() * denom()));

  auto& deriv_s_h_vv_coef = get<::Tags::TempScalar<9>>(temporaries);
  tenex::evaluate(make_not_null(&deriv_s_h_vv_coef),
                  fluid_lorentz_factor() *
                      (2.0 * fluid_lorentz_factor() * deriv_s_h_trace_coef() -
                       d_thick() * denom()));

  auto& deriv_s_h_sv_coef = get<::Tags::TempScalar<10>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_s_h_sv_coef),
      -d_thin() * fluid_lorentz_factor() * tilde_e() * inverse_s_norm());

  auto& deriv_s_h_ss_coef = get<::Tags::TempScalar<11>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_s_h_ss_coef),
      -2.0 * s_dot_fluid_velocity() * inverse_s_norm() * deriv_s_h_sv_coef());

  auto& deriv_s_j_s_coef = get<::Tags::TempScalar<12>>(temporaries);
  tenex::evaluate(
      make_not_null(&deriv_s_j_s_coef),
      -fluid_lorentz_factor() * s_dot_fluid_velocity() * deriv_s_h_ss_coef());

  auto& deriv_s_h_vs_coef = get<::Tags::TempScalar<13>>(temporaries);
  tenex::evaluate(make_not_null(&deriv_s_h_vs_coef),
                  -fluid_lorentz_factor() * deriv_s_j_s_coef());

  auto& deriv_e_j =
      get<imex::Tags::Jacobian<LocalTags::TildeE, LocalTags::TildeJ>>(
          temporaries);
  tenex::evaluate(make_not_null(&deriv_e_j),
                  square(fluid_lorentz_factor()) *
                      (1.0 + d_thin() * square(s_dot_fluid_velocity()) *
                                 inverse_s_norm() +
                       d_thick() * (3.0 - 2.0 * square(fluid_lorentz_factor())) *
                           fluid_velocity_norm() * denom()));

  auto& deriv_s_j =
      get<imex::Tags::Jacobian<LocalTags::TildeS, LocalTags::TildeJ>>(
          temporaries);
  tenex::evaluate<ti::J>(make_not_null(&deriv_s_j),
                         deriv_s_j_velocity_coef() * fluid_velocity(ti::J) +
                         deriv_s_j_s_coef() * tilde_s_upper(ti::J));

  auto& deriv_e_h =
      get<imex::Tags::Jacobian<LocalTags::TildeE, LocalTags::TildeHSpatial>>(
          temporaries);
  tenex::evaluate<ti::i>(
      make_not_null(&deriv_e_h),
      deriv_e_h_velocity_coef() * fluid_velocity_lower(ti::i) +
      deriv_e_h_s_coef() * tilde_s(ti::i));

  auto& deriv_s_h =
      get<imex::Tags::Jacobian<LocalTags::TildeS, LocalTags::TildeHSpatial>>(
          temporaries);
  tenex::evaluate<ti::J, ti::i>(
      make_not_null(&deriv_s_h),
      deriv_s_h_vv_coef() * fluid_velocity_lower(ti::i) *
          fluid_velocity(ti::J) +
      deriv_s_h_ss_coef() * tilde_s(ti::i) * tilde_s_upper(ti::J) +
      deriv_s_h_vs_coef() * fluid_velocity_lower(ti::i) * tilde_s_upper(ti::J) +
      deriv_s_h_sv_coef() * tilde_s(ti::i) * fluid_velocity(ti::J));
  for (size_t i = 0; i < 3; ++i) {
    deriv_s_h.get(i, i) += get(deriv_s_h_trace_coef);
  }

  Scalar<DataVector> j_difference{};
  tenex::evaluate(
      make_not_null(&j_difference),
      square(fluid_lorentz_factor()) * fluid_velocity_norm() * denom() *
      ((2.0 * square(fluid_lorentz_factor()) - 3.0) * tilde_e() -
       4.0 * square(fluid_lorentz_factor()) * s_dot_fluid_velocity())
      + square(fluid_lorentz_factor()) * tilde_e() *
      square(s_dot_fluid_velocity()) * inverse_s_norm()
                  );
  //FIXME consider using 4-tensor or avoiding Hn
  tnsr::i<DataVector, 3> h_spatial_difference{};//FIXME names?
  Scalar<DataVector> h_normal_difference{};
  {
    Scalar<DataVector> s_coef{};
    tenex::evaluate(
        make_not_null(&s_coef),
        fluid_lorentz_factor() *
        (1.0 - tilde_e() * s_dot_fluid_velocity() * inverse_s_norm()) -
        1.0 / fluid_lorentz_factor());
    Scalar<DataVector> n_coef{}; //FIXME why did I call this n_coef?
    tenex::evaluate(
        make_not_null(&n_coef),
        cube(fluid_lorentz_factor()) * denom() *
        ((2.0 * square(fluid_lorentz_factor()) - 3.0) * tilde_e() -
         (2.0 * square(fluid_lorentz_factor()) - 1.0) * s_dot_fluid_velocity()));
    tenex::evaluate<ti::i>(
        make_not_null(&h_spatial_difference),
        s_coef() * tilde_s(ti::i) +
        (s_coef() * square(fluid_lorentz_factor()) * s_dot_fluid_velocity()
        - n_coef()) * fluid_velocity_lower(ti::i));
    tenex::evaluate(
        make_not_null(&h_normal_difference),
        -s_coef() * square(fluid_lorentz_factor()) * s_dot_fluid_velocity()
        + n_coef() * fluid_velocity_norm());
  }

  Scalar<DataVector> deriv_dthin_prefactor{size_t{1}};//FIXME find bug test case
  tenex::evaluate(
      make_not_null(&deriv_dthin_prefactor),
      1.0 /
      (-comoving_momentum_density_spatial(ti::i) *
       inverse_spatial_metric(ti::I, ti::J) * h_spatial_difference(ti::j) +
       comoving_momentum_density_normal() * h_normal_difference() +
       square(closure_factor()) * comoving_energy_density() * j_difference() +
       5.0 / 3.0 * square(comoving_energy_density()) /
       (2.0 + closure_factor() * (-1.0 + closure_factor() * 4.0))));

  // (H1_i gamma^ij H2_j - H1n H2n)
  // deriv_e_h // OK
  Scalar<DataVector> deriv_e_hn{}; // FIXME existing coeffs?
  // 2 W^2 + 1 = W^2 (3 - v^2)
  tenex::evaluate(
      make_not_null(&deriv_e_hn),
      cube(fluid_lorentz_factor()) *
      (d_thick() * 4.0 * fluid_velocity_norm() * denom()
      + d_thin() *
      (fluid_velocity_norm()
       + square(s_dot_fluid_velocity()) * inverse_s_norm())));
  Scalar<DataVector> deriv_e_dthin{};
  tenex::evaluate(
      make_not_null(&deriv_e_dthin),
      deriv_dthin_prefactor() *
      (comoving_momentum_density_spatial(ti::i) *
       inverse_spatial_metric(ti::I, ti::J) * deriv_e_h(ti::j) -
       comoving_momentum_density_normal() * deriv_e_hn() -
       square(closure_factor()) * comoving_energy_density() * deriv_e_j()));
  Scalar<DataVector> deriv_e_j_fixed{};
  tenex::evaluate(make_not_null(&deriv_e_j_fixed),
                  deriv_e_j() + j_difference() * deriv_e_dthin());
  tnsr::i<DataVector, 3> deriv_e_h_fixed{};
  tenex::evaluate<ti::i>(
      make_not_null(&deriv_e_h_fixed),
      deriv_e_h(ti::i) + h_spatial_difference(ti::i) * deriv_e_dthin());

  tnsr::I<DataVector, 3> deriv_s_hn{}; // FIXME existing coeffs?
  tenex::evaluate<ti::I>(
      make_not_null(&deriv_s_hn),
      fluid_lorentz_factor() * (
          1.0
          - 6.0 * d_thick() * square(fluid_lorentz_factor()) * denom()
          + 2.0 * d_thin() * square(fluid_lorentz_factor()) * (tilde_e() * s_dot_fluid_velocity() * inverse_s_norm() - 1.0)
       ) * fluid_velocity(ti::I) -
        2.0 * d_thin() * cube(fluid_lorentz_factor()) * tilde_e() * square(s_dot_fluid_velocity()) * square(inverse_s_norm()) *
        tilde_s_upper(ti::I));
  tnsr::I<DataVector, 3> deriv_s_dthin{};
  tenex::evaluate<ti::K>(
      make_not_null(&deriv_s_dthin),
      deriv_dthin_prefactor() *
      (comoving_momentum_density_spatial(ti::i) *
       inverse_spatial_metric(ti::I, ti::J) * deriv_s_h(ti::K, ti::j) -
       comoving_momentum_density_normal() * deriv_s_hn(ti::K) -
       square(closure_factor()) * comoving_energy_density() *
       deriv_s_j(ti::K)));
  tnsr::I<DataVector, 3> deriv_s_j_fixed{};
  tenex::evaluate<ti::I>(
      make_not_null(&deriv_s_j_fixed),
      deriv_s_j(ti::I) + j_difference() * deriv_s_dthin(ti::I));
  tnsr::Ij<DataVector, 3> deriv_s_h_fixed{};
  tenex::evaluate<ti::J, ti::i>(
      make_not_null(&deriv_s_h_fixed),
      deriv_s_h(ti::J, ti::i) + h_spatial_difference(ti::i) * deriv_s_dthin(ti::J));

  deriv_e_j = deriv_e_j_fixed;//FIXME
  deriv_e_h = deriv_e_h_fixed;
  deriv_s_j = deriv_s_j_fixed;
  deriv_s_h = deriv_s_h_fixed;

  tenex::evaluate(deriv_e_source_e,
                  -lapse() * fluid_lorentz_factor() *
                      (total_opacity() - scattering_opacity() * deriv_e_j()));
  tenex::evaluate<ti::J>(deriv_s_source_e,
                         lapse() * fluid_lorentz_factor() *
                             (scattering_opacity() * deriv_s_j(ti::J) +
                              total_opacity() * fluid_velocity(ti::J)));
  tenex::evaluate<ti::i>(
      deriv_e_source_s,
      -lapse() * (total_opacity() * deriv_e_h(ti::i) +
                  fluid_lorentz_factor() * absorption_opacity() * deriv_e_j() *
                      fluid_velocity_lower(ti::i)));
  tenex::evaluate<ti::J, ti::i>(
      deriv_s_source_s,
      -lapse() * (total_opacity() * deriv_s_h(ti::J, ti::i) +
                  fluid_lorentz_factor() * absorption_opacity() *
                      fluid_velocity_lower(ti::i) * deriv_s_j(ti::J)));
}
}  // namespace RadiationTransport::M1Grey::detail
