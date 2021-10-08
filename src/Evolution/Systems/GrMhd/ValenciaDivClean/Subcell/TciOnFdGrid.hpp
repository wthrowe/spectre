// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/Tags/Inactive.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Subcell/TciOptions.hpp"
#include "Evolution/Systems/GrMhd/ValenciaDivClean/Tags.hpp"
#include "Evolution/VariableFixing/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Mesh;
/// \endcond

namespace grmhd::ValenciaDivClean::subcell {
/*!
 * \brief The troubled-cell indicator run on the FD grid to check if the
 * corresponding DG solution is admissible.
 *
 * The following checks are done in the order they are listed:
 *
 * - if `grmhd::ValenciaDivClean::Tags::VariablesNeededFixing` is
 *   `true` and the maximum of \f$\tilde{D}/\sqrt{\gamma}\f$ is
 *   greater than `tci_options.atmosphere_density`, then we remain on
 *   FD.
 * - if `min(tilde_d)` is less than
 *   `tci_options.minimum_rest_mass_density_times_lorentz_factor` or if
 *   `min(tilde_tau)` is less than `tci_options.minimum_tilde_tau` then the we
 *   remain on FD.
 * - apply the Persson TCI to \f$\tilde{D}\f$ and \f$\tilde{\tau}\f$
 */
struct TciOnFdGrid {
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<
      evolution::dg::subcell::Tags::Inactive<
          grmhd::ValenciaDivClean::Tags::TildeD>,
      evolution::dg::subcell::Tags::Inactive<
          grmhd::ValenciaDivClean::Tags::TildeTau>,
      evolution::dg::subcell::Tags::Inactive<
          grmhd::ValenciaDivClean::Tags::TildeB<>>,
      evolution::dg::subcell::Tags::Inactive<gr::Tags::SqrtDetSpatialMetric<>>,
      grmhd::ValenciaDivClean::Tags::VariablesNeededFixing,
      domain::Tags::Mesh<3>, Tags::TciOptions>;
  static bool apply(const Scalar<DataVector>& tilde_d,
                    const Scalar<DataVector>& tilde_tau,
                    const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_b,
                    const Scalar<DataVector>& sqrt_det_spatial_metric,
                    bool vars_needed_fixing, const Mesh<3>& dg_mesh,
                    const TciOptions& tci_options, double persson_exponent);
};
}  // namespace grmhd::ValenciaDivClean::subcell
