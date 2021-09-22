// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GrMhd/ValenciaDivClean/Subcell/PrimitiveGhostData.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/DgSubcell/Projection.hpp"
#include "Evolution/VariableFixing/FixToAtmosphere.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/TMPL.hpp"

namespace grmhd::ValenciaDivClean::subcell {
auto PrimitiveGhostDataOnSubcells::apply(
    const Variables<hydro::grmhd_tags<DataVector>>& prims)
    -> Variables<prims_to_reconstruct_tags> {
  Variables<prims_to_reconstruct_tags> vars_to_reconstruct(
      prims.number_of_grid_points());
  get<hydro::Tags::RestMassDensity<DataVector>>(vars_to_reconstruct) =
      get<hydro::Tags::RestMassDensity<DataVector>>(prims);
  get<hydro::Tags::Pressure<DataVector>>(vars_to_reconstruct) =
      get<hydro::Tags::Pressure<DataVector>>(prims);
  get<hydro::Tags::MagneticField<DataVector, 3>>(vars_to_reconstruct) =
      get<hydro::Tags::MagneticField<DataVector, 3>>(prims);
  get<hydro::Tags::DivergenceCleaningField<DataVector>>(vars_to_reconstruct) =
      get<hydro::Tags::DivergenceCleaningField<DataVector>>(prims);

  auto& lorentz_factor_time_spatial_velocity =
      get<hydro::Tags::LorentzFactorTimesSpatialVelocity<DataVector, 3>>(
          vars_to_reconstruct) =
          get<hydro::Tags::SpatialVelocity<DataVector, 3>>(prims);
  for (size_t i = 0; i < 3; ++i) {
    lorentz_factor_time_spatial_velocity.get(i) *=
        get(get<hydro::Tags::LorentzFactor<DataVector>>(prims));
  }
  return vars_to_reconstruct;
}

template <size_t ThermodynamicDim>
auto PrimitiveGhostDataToSlice::apply(
    const Variables<hydro::grmhd_tags<DataVector>>& prims,
    const Mesh<3>& dg_mesh, const Mesh<3>& subcell_mesh,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const EquationsOfState::EquationOfState<true, ThermodynamicDim>&
        equation_of_state,
    const VariableFixing::FixToAtmosphere<3>& fix_to_atmosphere)
    -> Variables<prims_to_reconstruct_tags> {
  auto projected_prims = evolution::dg::subcell::fd::project(
      PrimitiveGhostDataOnSubcells::apply(prims), dg_mesh,
      subcell_mesh.extents());
  fix_to_atmosphere.fix_ghost_data(
      &get<hydro::Tags::RestMassDensity<DataVector>>(projected_prims),
      &get<hydro::Tags::LorentzFactorTimesSpatialVelocity<DataVector, 3>>(
          projected_prims),
      &get<hydro::Tags::Pressure<DataVector>>(projected_prims), spatial_metric,
      equation_of_state);
  return projected_prims;
}

#define THERMO_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(_, data)                                             \
  template Variables<PrimitiveGhostDataToSlice::prims_to_reconstruct_tags> \
  PrimitiveGhostDataToSlice::apply(                                        \
      const Variables<hydro::grmhd_tags<DataVector>>& prims,               \
      const Mesh<3>& dg_mesh, const Mesh<3>& subcell_mesh,                 \
      const tnsr::ii<DataVector, 3>& spatial_metric,                       \
      const EquationsOfState::EquationOfState<true, THERMO_DIM(data)>&     \
          equation_of_state,                                               \
      const VariableFixing::FixToAtmosphere<3>& fix_to_atmosphere);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2))

#undef INSTANTIATION
#undef THERMO_DIM
}  // namespace grmhd::ValenciaDivClean::subcell
