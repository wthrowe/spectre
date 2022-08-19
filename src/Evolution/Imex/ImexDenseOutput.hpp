// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace imex::Tags {
template <typename ImplicitSector>
struct ImplicitHistory;
}  // namespace imex::Tags
/// \endcond

namespace imex {
/// Mutator to apply the implicit portion of dense output, intended
/// for use in `RunEventsAndDenseTriggers`.
template <typename ImexSystem>
struct ImexDenseOutput {
 private:
  static_assert(tt::assert_conforms_to_v<ImexSystem, protocols::ImexSystem>);

 public:
  using return_tags = tmpl::list<typename ImexSystem::variables_tag>;
  using argument_tags = tmpl::push_front<
      tmpl::transform<typename ImexSystem::implicit_sectors,
                      tmpl::bind<Tags::ImplicitHistory, tmpl::_1>>,
      ::Tags::TimeStepper<>, ::Tags::Time>;

  template <typename... Histories>
  static void apply(
      const gsl::not_null<typename ImexSystem::variables_tag::type*> variables,
      const TimeStepper& time_stepper, const double time,
      const Histories&... histories) {
    const auto update_sector = [&](auto sector_v, auto history) {
      using sector = tmpl::type_from<decltype(sector_v)>;
      auto sector_variables =
          variables->template reference_subset<typename sector::tensors>();
      time_stepper.dense_update_u_implicit(make_not_null(&sector_variables),
                                           history, time);
      return 0;
    };
    tmpl::as_pack<typename ImexSystem::implicit_sectors>([&](auto... sectors_v) {
      expand_pack(update_sector(sectors_v, histories)...);
    });
  }
};
}  // namespace imex
