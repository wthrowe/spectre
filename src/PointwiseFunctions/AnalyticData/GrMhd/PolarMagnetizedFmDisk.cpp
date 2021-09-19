// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/GrMhd/PolarMagnetizedFmDisk.hpp"

#include <pup.h>
#include <utility>

#include "Domain/CoordinateMaps/CoordinateMap.tpp"

namespace grmhd::AnalyticData {
PolarMagnetizedFmDisk::PolarMagnetizedFmDisk(
    MagnetizedFmDisk fm_disk, domain::CoordinateMaps::SphericalTorus torus_map)
    : fm_disk_(std::move(fm_disk)),
      torus_map_(
          domain::make_coordinate_map<Frame::BlockLogical, Frame::Inertial>(
              torus_map)),
      bare_torus_map_(std::move(torus_map)) {}

void PolarMagnetizedFmDisk::pup(PUP::er& p) {
  p | fm_disk_;
  p | torus_map_;
  p | bare_torus_map_;
}

bool operator==(const PolarMagnetizedFmDisk& lhs,
                const PolarMagnetizedFmDisk& rhs) {
  return lhs.fm_disk_ == rhs.fm_disk_ and lhs.torus_map_ == rhs.torus_map_;
}

bool operator!=(const PolarMagnetizedFmDisk& lhs,
                const PolarMagnetizedFmDisk& rhs) {
  return not(lhs == rhs);
}
}  // namespace grmhd::AnalyticData
