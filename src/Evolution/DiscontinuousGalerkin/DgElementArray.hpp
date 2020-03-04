// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <functional>
#include <queue>
#include <vector>

#include "AlgorithmArray.hpp"
#include "Domain/Block.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementId.hpp"
#include "Domain/ElementIndex.hpp"
#include "Domain/InitialElementIds.hpp"
#include "Domain/OptionTags.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Utilities/Numeric.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/*!
 * \brief The parallel component responsible for managing the DG elements that
 * compose the computational domain
 *
 * This parallel component will perform the actions specified by the
 * `PhaseDepActionList`.
 *
 */
template <class Metavariables, class PhaseDepActionList>
struct DgElementArray {
  static constexpr size_t volume_dim = Metavariables::volume_dim;

  using chare_type = Parallel::Algorithms::Array;
  using metavariables = Metavariables;
  using phase_dependent_action_list = PhaseDepActionList;
  using array_index = ElementIndex<volume_dim>;

  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<volume_dim>>;

  struct CostExponentOption {
    using type = double;
    static std::string name() noexcept { return "CostExponent"; }
    static constexpr OptionString help{
        "Exponent for estimating cost based on linear extents"};
  };

  struct CostExponent : db::SimpleTag {
    using type = double;
    using option_tags = tmpl::list<CostExponentOption>;

    static constexpr bool pass_metavariables = false;
    static double create_from_options(const double cost_exponent) noexcept {
      return cost_exponent;
    }
  };

  using array_allocation_tags =
      tmpl::list<domain::Tags::InitialRefinementLevels<volume_dim>,
                 domain::Tags::InitialExtents<volume_dim>,
                 CostExponent>;

  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>,
      array_allocation_tags>;

  static void allocate_array(
      Parallel::CProxy_ConstGlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<initialization_tags>&
          initialization_items) noexcept;

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_ConstGlobalCache<Metavariables>& global_cache) noexcept {
    auto& local_cache = *(global_cache.ckLocalBranch());
    Parallel::get_parallel_component<DgElementArray>(local_cache)
        .start_phase(next_phase);
  }
};

template <class Metavariables, class PhaseDepActionList>
void DgElementArray<Metavariables, PhaseDepActionList>::allocate_array(
    Parallel::CProxy_ConstGlobalCache<Metavariables>& global_cache,
    const tuples::tagged_tuple_from_typelist<initialization_tags>&
        initialization_items) noexcept {
  auto& local_cache = *(global_cache.ckLocalBranch());
  auto& dg_element_array =
      Parallel::get_parallel_component<DgElementArray>(local_cache);
  const auto& domain =
      Parallel::get<domain::Tags::Domain<volume_dim>>(local_cache);
  const auto& initial_refinement_levels =
      get<domain::Tags::InitialRefinementLevels<volume_dim>>(
          initialization_items);
  const auto& initial_extents =
      get<domain::Tags::InitialExtents<volume_dim>>(initialization_items);

  const int number_of_procs = Parallel::number_of_procs();
  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>> proc_loads;
  for (int i = 0; i < number_of_procs; ++i) {
    proc_loads.emplace(0.0, i);
  }

  const double cost_exponent = get<CostExponent>(initialization_items);
  const auto element_cost =
      [&cost_exponent](const double linear_size) noexcept {
    return pow(linear_size, cost_exponent);
  };

  std::priority_queue<std::pair<double, size_t>> block_costs;
  for (const auto& block : domain.blocks()) {
    const auto grid_points = initial_extents[block.id()];
    const double linear_size =
        pow(alg::accumulate(grid_points, 1.0, std::multiplies<double>{}),
            1.0 / volume_dim);
    block_costs.emplace(element_cost(linear_size), block.id());
  }

  while (not block_costs.empty()) {
    const auto cost = block_costs.top().first;
    const auto block_id = block_costs.top().second;
    block_costs.pop();
    const auto initial_ref_levs = initial_refinement_levels[block_id];
    const std::vector<ElementId<volume_dim>> element_ids =
        initial_element_ids(block_id, initial_ref_levs);
    for (size_t i = 0; i < element_ids.size(); ++i) {
      auto proc = proc_loads.top();
      proc_loads.pop();
      dg_element_array(ElementIndex<volume_dim>(element_ids[i]))
          .insert(global_cache, initialization_items, proc.second);
      proc.first += cost;
      proc_loads.push(proc);
    }
  }
  dg_element_array.doneInserting();
}
