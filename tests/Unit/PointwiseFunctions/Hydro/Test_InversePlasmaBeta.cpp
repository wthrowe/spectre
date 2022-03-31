// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "PointwiseFunctions/Hydro/InversePlasmaBeta.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace hydro {
namespace {
template <typename DataType>
void test_inverse_plasma_beta(const DataType& used_for_size) {
  const auto comoving_magnetic_field_squared =
      make_with_value<Scalar<DataType>>(used_for_size, 1.5);
  const auto fluid_pressure =
      make_with_value<Scalar<DataType>>(used_for_size, 2.25);
  const auto expected =
      make_with_value<Scalar<DataType>>(used_for_size, 1.5 / (2.0 * 2.25));
  const auto computed =
      inverse_plasma_beta(comoving_magnetic_field_squared, fluid_pressure);
  CHECK_EQUAL_APPROX(computed, expected);

  const auto box = db::create<
      db::AddSimpleTags<Tags::ComovingMagneticFieldSquared, Tags::Pressure>,
      db::AddComputeTags<Tags::InversePlasmaBetaCompute>>(
      comoving_magnetic_field_squared, fluid_pressure);
  CHECK(db::get<Tags::InversePlasmaBeta<DataType>>(box) == computed);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Hydro.InvserePlasmaBeta",
                  "[Unit][Hydro]") {
  GENERATE_UNINITIALIZED_DOUBLE_AND_DATAVECTOR;
  CHECK_FOR_DOUBLES_AND_DATAVECTORS(test_inverse_plasma_beta);
}
}  // namespace hydro
