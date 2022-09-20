// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Imex/Mode.hpp"
#include "Evolution/Tags.hpp"
#include "Options/Options.hpp"

namespace imex {
namespace OptionTags {
struct ImexMode {
  static constexpr Options::String help{"IMEX implementation to use"};
  using type = ::imex::Mode;
  using group = evolution::OptionTags::Group;
};
}  // namespace OptionTags

namespace Tags {
struct Mode : db::SimpleTag {
  using type = ::imex::Mode;
  using option_tags = tmpl::list<::imex::OptionTags::ImexMode>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& mode) { return mode; }
};
}  // namespace Tags
}  // namespace imex
