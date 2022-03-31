// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Inertial;
}  // namespace Frame
/// \endcond

namespace hydro {
/// @{
/*!
 * \ingroup EquationsOfStateGroup
 * \brief Computes the square of the comoving magnetic field
 *
 * The square of the comoving magnetic field is \f$B^2 / W^2 + (B^i
 * v_i)^2)\f$ where \f$B^i\f$ is the magnetic field, \f$v^i\f$ is the
 * spatial velocity of the fluid, and \$W\f$ is the Lorentz factor.
 */
template <typename DataType>
void comoving_magnetic_field_squared(
    gsl::not_null<Scalar<DataType>*> result,
    const tnsr::I<DataType, 3>& magnetic_field,
    const tnsr::I<DataType, 3>& spatial_velocity,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::ii<DataType, 3>& spatial_metric);

template <typename DataType>
Scalar<DataType> comoving_magnetic_field_squared(
    const tnsr::I<DataType, 3>& magnetic_field,
    const tnsr::I<DataType, 3>& spatial_velocity,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::ii<DataType, 3>& spatial_metric);
/// @}

namespace Tags {
/// Can be retrieved using `hydro::Tags::ComovingMagneticFieldSquared`
template <typename DataType>
struct ComovingMagneticFieldSquaredCompute
    : ComovingMagneticFieldSquared<DataType>, db::ComputeTag {
  using base = ComovingMagneticFieldSquared<DataType>;
  using return_type = Scalar<DataType>;

  using argument_tags =
      tmpl::list<MagneticField<DataType, 3>, SpatialVelocity<DataType, 3>,
                 LorentzFactor<DataType>,
                 gr::Tags::SpatialMetric<3, Frame::Inertial, DataType>>;

  static constexpr void (*function)(
      gsl::not_null<Scalar<DataType>*>, const tnsr::I<DataType, 3>&,
      const tnsr::I<DataType, 3>&, const Scalar<DataType>&,
      const tnsr::ii<DataType, 3>&) = comoving_magnetic_field_squared;
};
}  // namespace Tags
}  // namespace hydro
