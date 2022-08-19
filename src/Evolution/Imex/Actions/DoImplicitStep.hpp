// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

//FIXME
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Evolution/Imex/SolveImplicitSector.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
// IWYU pragma: no_forward_declare TimeDelta
// IWYU pragma: no_forward_declare db::DataBox
/// \endcond

namespace imex::Actions {
/// \ingroup ActionsGroup
/// \brief Perform implicit variable updates for one substep
///
/// Uses://FIXME
/// - DataBox:
///   - variables_tag (either the provided `VariablesTag` or the
///   `system::variables_tag` if none is provided)
///   - Tags::HistoryEvolvedVariables<variables_tag>
///   - Tags::TimeStep
///   - Tags::TimeStepper<>
///   - Tags::IsUsingTimeSteppingErrorControlBase
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - variables_tag
///   - Tags::HistoryEvolvedVariables<variables_tag>
struct DoImplicitStep {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    using system = typename Metavariables::system;
    static_assert(tt::conforms_to_v<system, protocols::ImexSystem>);
    //FIXME skip if no-op.  Here or in s_i_s?
    tmpl::for_each<typename system::implicit_sectors>([&](auto sector_v) {
      using sector = tmpl::type_from<decltype(sector_v)>;
      solve_implicit_sector<sector>(make_not_null(&box));
    });
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

//FIXME
template <typename System>
struct InitializeImex {
  static_assert(tt::conforms_to_v<System, protocols::ImexSystem>);
  using initialization_tags = tmpl::list<>;
  using simple_tags = tmpl::transform<typename System::implicit_sectors, tmpl::bind<Tags::ImplicitHistory, tmpl::_1>>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    const auto order = get<::Tags::HistoryEvolvedVariables<typename System::variables_tag>>(box).integration_order();
    tmpl::for_each<simple_tags>([&](auto tag) {
      using Tag = tmpl::type_from<decltype(tag)>;
      Initialization::mutate_assign<tmpl::list<Tag>>(
          make_not_null(&box), typename Tag::type(order));
    });
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace imex::Actions
