// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "Utilities/PrettyType.hpp"

#define MAX_NUMBER_OF_NEUTRINO_ENERGY_BINS 12

/// Namespace for neutrino physics
namespace neutrinos {

template <size_t EnergyBin>
struct ElectronNeutrinos {
  static constexpr size_t energy_bin = EnergyBin;
};
template <size_t EnergyBin>
struct ElectronAntiNeutrinos {
  static constexpr size_t energy_bin = EnergyBin;
};
template <size_t EnergyBin>
struct HeavyLeptonNeutrinos {
  static constexpr size_t energy_bin = EnergyBin;
};

template <typename Species>
std::string get_name(const Species& /*species*/) {
  return pretty_type::short_name<Species>() +
         std::to_string(Species::energy_bin);
}

}  // namespace neutrinos
