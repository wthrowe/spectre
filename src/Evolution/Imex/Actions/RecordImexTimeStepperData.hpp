// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

//FIXME
#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/NoSuchType.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
class TimeStepId;
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
// IWYU pragma: no_forward_declare db::DataBox
/// \endcond

namespace imex::Actions {
//FIXME doc test
/// \ingroup ActionsGroup
/// \ingroup TimeGroup
/// \brief Records the variables and their time derivatives in the
/// time stepper history.
///
/// With `dt_variables_tag = db::add_tag_prefix<Tags::dt, variables_tag>`:
///
/// Uses:
/// - GlobalCache: nothing
/// - DataBox:
///   - variables_tag (either the provided `VariablesTag` or the
///   `system::variables_tag` if none is provided)
///   - dt_variables_tag
///   - Tags::HistoryEvolvedVariables<variables_tag>
///   - Tags::TimeStepId
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - Tags::HistoryEvolvedVariables<variables_tag>
struct RecordImexTimeStepperData {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {  // NOLINT const
    using system = typename Metavariables::system;
    static_assert(tt::conforms_to_v<system, imex::protocols::ImexSystem>);
    static constexpr size_t volume_dim = system::volume_dim;

    tmpl::for_each<typename system::implicit_sectors>([&](auto sector_v) {
      using sector = tmpl::type_from<decltype(sector_v)>;
      using source = typename sector::source;
      using history_tag = imex::Tags::ImplicitHistory<sector>;
      db::mutate_apply<tmpl::list<history_tag>,
                       tmpl::push_front<typename source::argument_tags, domain::Tags::Mesh<volume_dim>, ::Tags::TimeStepId>>(
        [&](const gsl::not_null<typename history_tag::type*> history,
            const Mesh<volume_dim>& mesh, const TimeStepId& time_step_id,
            const auto&... source_arguments) {
          Variables<db::wrap_tags_in<::Tags::Source, typename sector::tensors>>
              implicit_source(mesh.number_of_grid_points(), 0.0);
          tmpl::as_pack<typename source::return_tags>(
              [&](auto... return_tags) {
                source::apply(
                    &get<tmpl::type_from<decltype(return_tags)>>(implicit_source)...,
                    source_arguments...);
              });
          history->insert(time_step_id, implicit_source.template reference_with_different_prefixes<Variables<db::wrap_tags_in<::Tags::dt, typename sector::tensors>>>());
        },
        make_not_null(&box));
    });

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace imex::Actions
