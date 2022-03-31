// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace hydro {
/// @{
/*!
 * \ingroup EquationsOfStateGroup
 * \brief Computes the inverse plasma beta
 *
 * The inverse plasma beta \f$\beta^{-1} = b^2 / (2 p)\f$, where
 * \f$b^2\f$ is the square of the comoving magnetic field amplitude
 * and \f$p\f$ is the fluid pressure.
 */
template <typename DataType>
void inverse_plasma_beta(
    gsl::not_null<Scalar<DataType>*> result,
    const Scalar<DataType>& comoving_magnetic_field_squared,
    const Scalar<DataType>& fluid_pressure);

template <typename DataType>
Scalar<DataType> inverse_plasma_beta(
    const Scalar<DataType>& comoving_magnetic_field_squared,
    const Scalar<DataType>& fluid_pressure);
/// @}

namespace Tags {
/// Can be retrieved using `hydro::Tags::InversePlasmaBeta`
template <typename DataType>
struct InversePlasmaBetaCompute : InversePlasmaBeta<DataType>, db::ComputeTag {
  using base = InversePlasmaBeta<DataType>;
  using return_type = Scalar<DataType>;

  using argument_tags =
      tmpl::list<ComovingMagneticFieldSquared<DataType>, Pressure<DataType>>;

  static constexpr void (*function)(
      gsl::not_null<Scalar<DataType>*>, const Scalar<DataType>&,
      const Scalar<DataType>&) = inverse_plasma_beta;
};
}  // namespace Tags
}  // namespace hydro
