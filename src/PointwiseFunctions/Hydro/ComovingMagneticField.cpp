// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/Hydro/ComovingMagneticField.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace hydro {

template <typename DataType>
void comoving_magnetic_field_squared(
    gsl::not_null<Scalar<DataType>*> result,
    const tnsr::I<DataType, 3>& magnetic_field,
    const tnsr::I<DataType, 3>& spatial_velocity,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::ii<DataType, 3>& spatial_metric) {
  const auto b_dot_v =
      dot_product(magnetic_field, spatial_velocity, spatial_metric);
  dot_product(result, magnetic_field, magnetic_field, spatial_metric);
  get(*result) /= square(get(lorentz_factor));
  get(*result) += square(get(b_dot_v));
}

template <typename DataType>
Scalar<DataType> comoving_magnetic_field_squared(
    const tnsr::I<DataType, 3>& magnetic_field,
    const tnsr::I<DataType, 3>& spatial_velocity,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::ii<DataType, 3>& spatial_metric) {
  Scalar<DataType> result{};
  comoving_magnetic_field_squared(make_not_null(&result), magnetic_field,
                                  spatial_velocity, lorentz_factor,
                                  spatial_metric);
  return result;
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                    \
  template void comoving_magnetic_field_squared(                \
      gsl::not_null<Scalar<DTYPE(data)>*> result,               \
      const tnsr::I<DTYPE(data), 3>& magnetic_field,            \
      const tnsr::I<DTYPE(data), 3>& spatial_velocity,          \
      const Scalar<DTYPE(data)>& lorentz_factor,                \
      const tnsr::ii<DTYPE(data), 3>& spatial_metric);          \
  template Scalar<DTYPE(data)> comoving_magnetic_field_squared( \
      const tnsr::I<DTYPE(data), 3>& magnetic_field,            \
      const tnsr::I<DTYPE(data), 3>& spatial_velocity,          \
      const Scalar<DTYPE(data)>& lorentz_factor,                \
      const tnsr::ii<DTYPE(data), 3>& spatial_metric);

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector))

#undef INSTANTIATE
#undef DTYPE
}  // namespace hydro
