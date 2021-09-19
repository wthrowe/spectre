// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <type_traits>

#include "DataStructures/Tensor/EagerMath/Determinant.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/SphericalTorus.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/MagnetizedFmDisk.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Solutions.hpp"
#include "PointwiseFunctions/GeneralRelativity/Transform.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace grmhd::AnalyticData {
/*!
 * \brief Magnetized fluid disk orbiting a Kerr black hole. FIXME
 *
 * In the context of simulating accretion disks, this class implements a widely
 * used (e.g. \cite Gammie2003, \cite Porth2016rfi, \cite White2015omx)
 * initial setup for the GRMHD variables, consisting of a Fishbone-Moncrief disk
 * \cite Fishbone1976apj (see also
 * RelativisticEuler::Solutions::FishboneMoncriefDisk),
 * threaded by a weak poloidal magnetic field. The magnetic field is constructed
 * from an axially symmetric toroidal magnetic potential which, in Kerr
 * ("spherical Kerr-Schild") coordinates, has the form
 *
 * \f{align*}
 * A_\phi(r,\theta) \propto \text{max}(\rho(r,\theta) - \rho_\text{thresh}, 0),
 * \f}
 *
 * where \f$\rho_\text{thresh}\f$ is a user-specified threshold density that
 * confines the magnetic flux to exist inside of the fluid disk only. A commonly
 * used value for this parameter is
 * \f$\rho_\text{thresh} = 0.2\rho_\text{max}\f$, where \f$\rho_\text{max}\f$
 * is the maximum value of
 * the rest mass density in the disk. Using this potential, the Eulerian
 * magnetic field takes the form
 *
 * \f{align*}
 * B^r = \frac{F_{\theta\phi}}{\sqrt{\gamma}},\quad
 * B^\theta = \frac{F_{\phi r}}{\sqrt{\gamma}},\quad B^\phi = 0,
 * \f}
 *
 * where \f$F_{ij} = \partial_i A_j - \partial_j A_i\f$ are the spatial
 * components of the Faraday tensor, and \f$\gamma\f$ is the determinant of the
 * spatial metric. The magnetic field is then normalized so that the
 * plasma-\f$\beta\f$ parameter, \f$\beta = 2p/b^2\f$, equals some value
 * specified by the user. Here, \f$p\f$ is the fluid pressure, and
 *
 * \f{align*}
 * b^2 = b^\mu b_\mu = \frac{B_iB^i}{W^2} + (B^iv_i)^2
 * \f}
 *
 * is the norm of the magnetic field in the fluid frame, with \f$v_i\f$ being
 * the spatial velocity, and \f$W\f$ the Lorentz factor.
 */
class PolarMagnetizedFmDisk : public MarkAsAnalyticData,
                              public grmhd::AnalyticDataBase {
 public:
  struct DiskParameters {
    using type = MagnetizedFmDisk;
    static constexpr Options::String help = "Parameters for the disk.";
  };

  struct TorusParameters {
    using type = domain::CoordinateMaps::SphericalTorus;
    static constexpr Options::String help =
        "Parameters for the evolution region.";
  };

  using options = tmpl::list<DiskParameters, TorusParameters>;

  static constexpr Options::String help =
      "Magnetized Fishbone-Moncrief disk in polar coordinates.";

  PolarMagnetizedFmDisk() = default;

  PolarMagnetizedFmDisk(MagnetizedFmDisk fm_disk,
                        domain::CoordinateMaps::SphericalTorus torus_map);

  /// The grmhd variables.
  ///
  /// \note The functions are optimized for retrieving the hydro variables
  /// before the metric variables.
  template <typename DataType, typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataType, 3>& x,
                                         tmpl::list<Tags...> /*meta*/) const {
    // In this function, we label the coordinates this solution works
    // in with Frame::BlockLogical, and the coordinates the wrapped
    // solution uses Frame::Inertial.  This means the input and output
    // have to be converted to the correct label.

    // We have to copy here to match the CoordinateMap interface.
    tnsr::I<DataType, 3, Frame::BlockLogical> x_logical;
    std::array<DataType, 3> x_array;
    for (size_t i = 0; i < 3; ++i) {
      x_logical.get(i) = x.get(i);
      gsl::at(x_array, i) = x.get(i);
    }

    const tnsr::I<DataType, 3> observation_coordinates(torus_map_(x_logical));

    using dependencies = tmpl::map<
        tmpl::pair<gr::AnalyticSolution<3>::DerivShift<DataType>,
                   gr::Tags::Shift<3, Frame::Inertial, DataType>>,
        tmpl::pair<gr::AnalyticSolution<3>::DerivSpatialMetric<DataType>,
                   gr::Tags::SpatialMetric<3, Frame::Inertial, DataType>>>;
    using required_tags = tmpl::remove_duplicates<
        tmpl::remove<tmpl::list<Tags..., tmpl::at<dependencies, Tags>...>,
                     tmpl::no_such_type_>>;

    auto observation_data =
        fm_disk_.variables(observation_coordinates, required_tags{});

    const auto jacobian = torus_map_.jacobian(x_logical);
    const auto inv_jacobian = torus_map_.inv_jacobian(x_logical);

    const auto change_frame = [this, &inv_jacobian, &jacobian, &x_array](
                                  const auto& data, auto tag) {
      using Tag = decltype(tag);
      auto result =
          transform::to_different_frame(get<Tag>(data), jacobian, inv_jacobian);

      if constexpr (std::is_same_v<
                        Tag, gr::AnalyticSolution<3>::DerivShift<DataType>>) {
        const auto deriv_inv_jacobian =
            bare_torus_map_.derivative_of_inv_jacobian(x_array);
        const auto& shift =
            get<gr::Tags::Shift<3, Frame::Inertial, DataType>>(data);
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
              result.get(i, j) +=
                  deriv_inv_jacobian.get(j, k, i) * shift.get(k);
            }
          }
        }
      } else if constexpr (std::is_same_v<
                               Tag, gr::AnalyticSolution<3>::DerivSpatialMetric<
                                        DataType>>) {
        const auto hessian = bare_torus_map_.hessian(x_array);
        const auto& spatial_metric =
            get<gr::Tags::SpatialMetric<3, Frame::Inertial, DataType>>(data);
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            for (size_t k = j; k < 3; ++k) {
              for (size_t l = 0; l < 3; ++l) {
                for (size_t m = 0; m < 3; ++m) {
                  result.get(i, j, k) +=
                      (hessian.get(l, j, i) * jacobian.get(m, k) +
                       hessian.get(l, k, i) * jacobian.get(m, j)) *
                      spatial_metric.get(l, m);
                }
              }
            }
          }
        }
      } else if constexpr (std::is_same_v<
                               Tag, gr::Tags::SqrtDetSpatialMetric<DataType>>) {
        get(result) *= abs(get(determinant(jacobian)));
      }

      typename Tag::type result_with_replaced_frame{};
      std::copy(std::move_iterator(result.begin()),
                std::move_iterator(result.end()),
                result_with_replaced_frame.begin());
      return result_with_replaced_frame;
    };

    return {change_frame(observation_data, Tags{})...};
  }

  using equation_of_state_type = MagnetizedFmDisk::equation_of_state_type;
  const equation_of_state_type& equation_of_state() const {
    return fm_disk_.equation_of_state();
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  friend bool operator==(const PolarMagnetizedFmDisk& lhs,
                         const PolarMagnetizedFmDisk& rhs);

  MagnetizedFmDisk fm_disk_;
  domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial,
                        domain::CoordinateMaps::SphericalTorus>
      torus_map_;
  domain::CoordinateMaps::SphericalTorus bare_torus_map_;
};

bool operator!=(const PolarMagnetizedFmDisk& lhs,
                const PolarMagnetizedFmDisk& rhs);
}  // namespace grmhd::AnalyticData
