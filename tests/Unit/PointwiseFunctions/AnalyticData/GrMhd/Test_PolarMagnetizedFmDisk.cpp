// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <limits>
#include <string>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/SphericalTorus.hpp"
#include "Evolution/TypeTraits.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/PointwiseFunctions/AnalyticData/TestHelpers.hpp"
#include "Helpers/PointwiseFunctions/AnalyticSolutions/GeneralRelativity/VerifyGrSolution.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/MagnetizedFmDisk.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/PolarMagnetizedFmDisk.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace {

static_assert(not evolution::is_analytic_solution_v<
              grmhd::AnalyticData::PolarMagnetizedFmDisk>);
static_assert(
    evolution::is_analytic_data_v<grmhd::AnalyticData::PolarMagnetizedFmDisk>);

void test_create_from_options() noexcept {
  const auto disk =
      TestHelpers::test_creation<grmhd::AnalyticData::PolarMagnetizedFmDisk>(
          "DiskParameters:\n"
          "  BhMass: 1.3\n"
          "  BhDimlessSpin: 0.345\n"
          "  InnerEdgeRadius: 6.123\n"
          "  MaxPressureRadius: 14.2\n"
          "  PolytropicConstant: 0.065\n"
          "  PolytropicExponent: 1.654\n"
          "  ThresholdDensity: 0.42\n"
          "  InversePlasmaBeta: 85.0\n"
          "  BFieldNormGridRes: 6\n"
          "TorusParameters:\n"
          "  RadialRange: [2.5, 3.6]\n"
          "  PhiMax: 0.9\n"
          "  FractionOfTorus: 0.7\n");
  CHECK(disk ==
        grmhd::AnalyticData::PolarMagnetizedFmDisk(
            grmhd::AnalyticData::MagnetizedFmDisk(1.3, 0.345, 6.123, 14.2,
                                                  0.065, 1.654, 0.42, 85.0, 6),
            domain::CoordinateMaps::SphericalTorus(2.5, 3.6, 0.9, 0.7)));
}

struct Wrapper {
  template <typename Tags>
  auto variables(const tnsr::I<double, 3>& x, const double /*t*/,
                 Tags tags) const noexcept {
    return disk_->variables(x, tags);
  }
  const grmhd::AnalyticData::PolarMagnetizedFmDisk* disk_;
};

template <typename DataType>
void test_variables(const DataType& used_for_size) {
  const double bh_mass = 1.12;
  const double bh_dimless_spin = 0.97;
  const double inner_edge_radius = 6.2;
  const double max_pressure_radius = 11.6;
  const double polytropic_constant = 0.034;
  const double polytropic_exponent = 1.65;
  const double threshold_density = 0.14;
  const double inverse_plasma_beta = 0.023;

  const grmhd::AnalyticData::PolarMagnetizedFmDisk disk(
      grmhd::AnalyticData::MagnetizedFmDisk(
          bh_mass, bh_dimless_spin, inner_edge_radius, max_pressure_radius,
          polytropic_constant, polytropic_exponent, threshold_density,
          inverse_plasma_beta),
      domain::CoordinateMaps::SphericalTorus(3.0, 20.0, 1.0, 0.3));

  const auto coords = make_with_value<tnsr::I<DataType, 3>>(used_for_size, 0.5);

  TestHelpers::AnalyticData::test_tag_retrieval(
      disk, coords,
      typename grmhd::AnalyticData::PolarMagnetizedFmDisk::tags<DataType>{});

  if constexpr (std::is_same_v<DataType, double>) {
    TestHelpers::VerifyGrSolution::verify_consistency(Wrapper{&disk}, 0.0,
                                                      coords, 1.0e-2, 1.0e-10);
  }

  //FIXME Is there something easy to test for correctness?
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticData.GrMhd.PolarMagFmDisk",
                  "[Unit][PointwiseFunctions]") {
  test_create_from_options();

  test_variables(std::numeric_limits<double>::signaling_NaN());
  test_variables(DataVector(5));
}
