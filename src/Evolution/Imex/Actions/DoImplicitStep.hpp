// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>

#include "Evolution/Imex/Protocols/ImexSystem.hpp"
#include "Evolution/Imex/SolveImplicitSector.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace db {
template <typename TagsList>
class DataBox;
}  // namespace db
namespace tuples {
template <class... Tags>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace imex::Actions {
/// \ingroup ActionsGroup
/// \brief Perform implicit variable updates for one substep
///
/// Uses:
/// - DataBox:
///   - Tags::TimeStep
///   - Tags::TimeStepper<>
///   - imex::Tags::Mode
///   - as required by system implicit sectors
///
/// DataBox changes:
/// - variables_tag
/// - imex::Tags::ImplicitHistory<sector> for each sector
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
    static_assert(tt::assert_conforms_to_v<system, protocols::ImexSystem>);
    tmpl::for_each<typename system::implicit_sectors>([&](auto sector_v) {
      using sector = tmpl::type_from<decltype(sector_v)>;
      solve_implicit_sector<sector>(make_not_null(&box));
    });
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace imex::Actions
