// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <utility>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/System.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "Helpers/Evolution/Imex/TestSector.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct DummySpecies {
  static constexpr size_t energy_bin = 0;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.RadiationTransport.M1Grey.System.Imex",
                  "[Unit][M1Grey]") {
  // FIXME test with two species
  using system =
      RadiationTransport::M1Grey::System<tmpl::list<DummySpecies>>;
  static_assert(tt::assert_conforms_to_v<system, imex::protocols::ImexSystem>);
  using sector = tmpl::front<system::implicit_sectors>;

  using tilde_e_tag =
      RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial, DummySpecies>;
  using tilde_s_tag =
      RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial, DummySpecies>;

  //FIXME better values?
  Variables<tmpl::list<tilde_e_tag, tilde_s_tag>> sector_variables(1);
  auto& tilde_e = get<tilde_e_tag>(sector_variables);
  get(tilde_e) = 10.2;
  auto& tilde_s = get<tilde_s_tag>(sector_variables);
  get<0>(tilde_s) = 0.1;
  get<1>(tilde_s) = 0.2;
  get<2>(tilde_s) = 0.3;

  Scalar<DataVector> lapse{};
  get(lapse) = DataVector{1.0};
  tnsr::ii<DataVector, 3> spatial_metric{};
  get<0, 0>(spatial_metric) = DataVector{0.1};
  get<0, 1>(spatial_metric) = DataVector{0.0};
  get<0, 2>(spatial_metric) = DataVector{0.0};
  get<1, 1>(spatial_metric) = DataVector{0.1};
  get<1, 2>(spatial_metric) = DataVector{0.0};
  get<2, 2>(spatial_metric) = DataVector{0.1};
  tnsr::II<DataVector, 3> inverse_spatial_metric{};
  Scalar<DataVector> sqrt_det_spatial_metric{};
  determinant_and_inverse(make_not_null(&sqrt_det_spatial_metric),
                          make_not_null(&inverse_spatial_metric),
                          spatial_metric);
  get(sqrt_det_spatial_metric) = sqrt(get(sqrt_det_spatial_metric));

  Scalar<DataVector> emissivity{};
  get(emissivity) = DataVector{2.0};
  Scalar<DataVector> absorption_opacity{};
  get(absorption_opacity) = DataVector{4.0};
  Scalar<DataVector> scattering_opacity{};
  get(scattering_opacity) = DataVector{3.0};
  tnsr::I<DataVector, 3> fluid_velocity{};
  get<0>(fluid_velocity) = DataVector{0.1};
  get<1>(fluid_velocity) = DataVector{0.1};
  get<2>(fluid_velocity) = DataVector{0.1};
  Scalar<DataVector> lorentz_factor{1_st};
  tenex::evaluate(
      make_not_null(&lorentz_factor),
      1.0 / sqrt(1.0 - fluid_velocity(ti::I) * fluid_velocity(ti::J) *
                 spatial_metric(ti::i, ti::j)));
  // Internal root-find only solves to 1.0e-6.  Numerical derivatives
  // in the test could amplify that, but it doesn't seem to be an
  // issue in practice.
  const double tolerance = 1.0e-5;
  TestHelpers::imex::test_sector<sector>(
      1.0e-3, tolerance, sector_variables,
      {std::move(lorentz_factor), std::move(fluid_velocity),
       std::move(emissivity),
       std::move(absorption_opacity), std::move(scattering_opacity),
       std::move(lapse),
       std::move(spatial_metric),
       std::move(inverse_spatial_metric),
       std::move(sqrt_det_spatial_metric)
});
}
