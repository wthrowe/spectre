// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Evolution/Imex/Tags/ImplicitHistory.hpp"
#include "Evolution/Imex/Tags/Mode.hpp"
#include "Time/History.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace imex {
template <typename System, typename = typename System::implicit_sectors>
struct Initialize;

template <typename System, typename... Sectors>
struct Initialize<System, tmpl::list<Sectors...>> {
  static_assert(tt::assert_conforms_to_v<System, protocols::ImexSystem>);

  using const_global_cache_tags = tmpl::list<imex::Tags::Mode>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<Tags::ImplicitHistory<Sectors>...>;
  using compute_tags = tmpl::list<>;

  using return_tags = simple_tags;
  using argument_tags = tmpl::list<::Tags::HistoryEvolvedVariables<>>;

  static void apply(
      const gsl::not_null<typename Tags::ImplicitHistory<Sectors>::type*>...
          histories,
      const TimeSteppers::History<typename System::variables_tag::type>&
          explicit_history) {
    const auto order = explicit_history.integration_order();
    expand_pack((histories->integration_order(order), 0)...);
  }
};
}  // namespace imex
