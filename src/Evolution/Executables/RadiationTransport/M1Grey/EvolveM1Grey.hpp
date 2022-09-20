// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <vector>

#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "Evolution/ComputeTags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/DiscontinuousGalerkin/Limiters/LimiterActions.hpp"
#include "Evolution/DiscontinuousGalerkin/Limiters/Minmod.hpp"
#include "Evolution/DiscontinuousGalerkin/Limiters/Tags.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Factory.hpp"
#include "Evolution/Imex/Actions/DoImplicitStep.hpp"
#include "Evolution/Imex/ImexDenseOutput.hpp"
#include "Evolution/Imex/Actions/RecordImexTimeStepperData.hpp"
#include "Evolution/Initialization/ConservativeSystem.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Initialization/GrTagsForHydro.hpp"
#include "Evolution/Initialization/Limiter.hpp"
#include "Evolution/Initialization/SetVariables.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Initialize.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/M1Closure.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/M1HydroCoupling.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/System.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "Evolution/Systems/RadiationTransport/Tags.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "Options/Options.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/CheckpointAndExitAfterWallclock.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/PhaseControl/VisitAndReturn.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/Events/Factory.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/HomogeneousSphere.hpp"
#include "PointwiseFunctions/AnalyticData/RadiationTransport/M1Grey/SphericalGaussian.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/RadiationTransport/M1Grey/ConstantM1.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Actions/RecordTimeStepperData.hpp"
#include "Time/Actions/SelfStartActions.hpp"  // IWYU pragma: keep
#include "Time/Actions/UpdateU.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/StepControllers/Factory.hpp"
#include "Time/StepControllers/StepController.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Inertial;
}  // namespace Frame
namespace PUP {
class er;
}  // namespace PUP
namespace Parallel {
template <typename Metavariables>
class CProxy_GlobalCache;
}  // namespace Parallel
/// \endcond

#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"
struct FIXME_Flatten {
  using return_tags = tmpl::list<
      RadiationTransport::M1Grey::Tags::TildeE<Frame::Inertial, neutrinos::ElectronNeutrinos<1>>,
      RadiationTransport::M1Grey::Tags::TildeS<Frame::Inertial, neutrinos::ElectronNeutrinos<1>>>;
  using argument_tags = tmpl::list<
      gr::Tags::InverseSpatialMetric<3>,
      domain::Tags::DetInvJacobian<Frame::ElementLogical, Frame::Inertial>,
      domain::Tags::Mesh<3>, Tags::Time, domain::Tags::Element<3>>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> tilde_e,
                    const gsl::not_null<tnsr::i<DataVector, 3>*> tilde_s,
                    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
                    const Scalar<DataVector>& det_jacobian,
                    const Mesh<3>& mesh,

const double time, const Element<3>& element) {
    //FIXME allocations
    const auto energy_excess = tenex::evaluate(
        square((*tilde_e)()) - (*tilde_s)(ti::i) * inverse_spatial_metric(ti::I, ti::J) * (*tilde_s)(ti::j));
    if (/*min(get(*tilde_e)) >= 0.0 and */min(get(energy_excess)) >= 0.0) {
      return;
    }
#if 1  // Scale all
    const auto element_volume = definite_integral(get(det_jacobian), mesh);
    const auto average_e =
        definite_integral(get(det_jacobian) * get(*tilde_e), mesh) / element_volume;
    double average_energy_excess = square(average_e);
    tnsr::i<double, 3> average_s;
    for (size_t i = 0; i < 3; ++i) {
      average_s.get(i) =
          definite_integral(get(det_jacobian) * tilde_s->get(i), mesh) / element_volume;
      average_energy_excess -= square(average_s.get(i)); // FIXME metric
    }
    if (average_energy_excess < 0.0) {
      ERROR("Average values are nonphysical.\n" << *tilde_e << "\n" << *tilde_s << "\n" << average_e << "\n" << average_s << "\nt = " << time);
    }
    double flatten_factor = 1.0;
    // if (min(get(*tilde_e)) < 0.0) {
    //   flatten_factor = average_e / (average_e - min(get(*tilde_e)));
    // }
    //FIXME try rewriting with 0 = no change
    for (size_t i = 0; i < get(*tilde_e).size(); ++i) {
      if (get(*tilde_e)[i] < 0.0 or get(energy_excess)[i] < 0.0) {
        //FIXME guarantee single positive?
        const auto point_flatten_factor = real_roots(
            square(get(*tilde_e)[i] - average_e) -
            square(get<0>(*tilde_s)[i] - get<0>(average_s)) -
            square(get<1>(*tilde_s)[i] - get<1>(average_s)) -
            square(get<2>(*tilde_s)[i] - get<2>(average_s)),
            2.0 * ((get(*tilde_e)[i] - average_e) * average_e -
                   (get<0>(*tilde_s)[i] - get<0>(average_s)) * get<0>(average_s) -
                   (get<1>(*tilde_s)[i] - get<1>(average_s)) * get<1>(average_s) -
                   (get<2>(*tilde_s)[i] - get<2>(average_s)) * get<2>(average_s)),
            average_energy_excess);
        ASSERT(point_flatten_factor.has_value(),
               get(*tilde_e)[i] << " "
               << get<0>(*tilde_s)[i] << " "
               << get<1>(*tilde_s)[i] << " "
               << get<2>(*tilde_s)[i] << " "
               << average_e << " "
               << get<0>(average_s) << " "
               << get<1>(average_s) << " "
               << get<2>(average_s));
        if (point_flatten_factor.has_value()) {
          if ((*point_flatten_factor)[0] > 0.0 and
              (*point_flatten_factor)[0] < flatten_factor) {
            flatten_factor = (*point_flatten_factor)[0];
          } else if ((*point_flatten_factor)[1] < flatten_factor) {
            flatten_factor = (*point_flatten_factor)[1];
          }
        }
      }
    }
    ASSERT(flatten_factor >= 0.0 and flatten_factor <= 1.0,
           "Not a valid adjustment: " << flatten_factor
           << "\n" << *tilde_e << "\n" << *tilde_s);
    for (size_t i = 0; i < get(*tilde_e).size(); ++i) {
      get(*tilde_e)[i] =
          average_e + (get(*tilde_e)[i] - average_e) * flatten_factor;
      for (size_t j = 0; j < 3; ++j) {
        tilde_s->get(j)[i] =
            average_s.get(j) + (tilde_s->get(j)[i] -
                                average_s.get(j)) * flatten_factor;
      }
    }
#endif  // scale all
#if 0  // scale s
    const auto element_volume = definite_integral(get(det_jacobian), mesh);
    tnsr::i<double, 3> average_s;
    double average_s_norm = 0.0;
    for (size_t i = 0; i < 3; ++i) {
      average_s.get(i) =
          definite_integral(get(det_jacobian) * tilde_s->get(i), mesh) /
          element_volume;
      average_s_norm += square(average_s.get(i)); //FIXME metric
    }
    // 0 >= (s0^2 - e^2) + x 2 s0.(s - s0) + x^2 (s - s0)^2
    double flatten_factor = 1.0;
    //FIXME try rewriting with 0 = no change
    for (size_t i = 0; i < get(*tilde_e).size(); ++i) {
      if (get(energy_excess)[i] < 0.0) {
        const auto point_flatten_factor = real_roots(
            square(get<0>(*tilde_s)[i] - get<0>(average_s)) +
            square(get<1>(*tilde_s)[i] - get<1>(average_s)) +
            square(get<2>(*tilde_s)[i] - get<2>(average_s)),
            2.0 *
            ((get<0>(*tilde_s)[i] - get<0>(average_s)) * get<0>(average_s) +
             (get<1>(*tilde_s)[i] - get<1>(average_s)) * get<1>(average_s) +
             (get<2>(*tilde_s)[i] - get<2>(average_s)) * get<2>(average_s)),
            average_s_norm - square(get(*tilde_e)[i]));
        if (not point_flatten_factor.has_value()) {
          ERROR(get(*tilde_e)[i] << " "
                << get<0>(*tilde_s)[i] << " "
                << get<1>(*tilde_s)[i] << " "
                << get<2>(*tilde_s)[i] << " "
                << get<0>(average_s) << " "
                << get<1>(average_s) << " "
                << get<2>(average_s));
        }

        if ((*point_flatten_factor)[1] < 0.0 or
            (*point_flatten_factor)[1] > 1.0) {
          ERROR("Can't flatten: "
                << (*point_flatten_factor)[1] << " " << i << "\n"
                << get(*tilde_e)[i] << " "
                << get<0>(*tilde_s)[i] << " "
                << get<1>(*tilde_s)[i] << " "
                << get<2>(*tilde_s)[i] << " "
                << get<0>(average_s) << " "
                << get<1>(average_s) << " "
                << get<2>(average_s));
          if ((*point_flatten_factor)[1] < flatten_factor) {
            flatten_factor = (*point_flatten_factor)[1];
          }
        }
      }
    }

    for (size_t j = 0; j < 3; ++j) {
      tilde_s->get(j) =
          average_s.get(j) +
          (tilde_s->get(j) - average_s.get(j)) * flatten_factor;
    }

    if (min(square(get(*tilde_e)) - square(get<0>(*tilde_s))
            - square(get<1>(*tilde_s))
            - square(get<2>(*tilde_s))) < 1.0) {
      ERROR("Flattening failed: "
            << flatten_factor << " "
            << get(*tilde_e) << " "
            << get<0>(*tilde_s) << " "
            << get<1>(*tilde_s) << " "
            << get<2>(*tilde_s) << " "
            << get<0>(average_s) << " "
            << get<1>(average_s) << " "
            << get<2>(average_s));
    }
#endif  // scale s
    // if (flatten_factor < 0.9 and std::abs(time - 4.0) < 1e-5) {
    //   Parallel::printf("%f\t%f\t%s\n", time, flatten_factor, element.id());
    // }
  }
};

struct FIXME_Atmosphere {
  using return_tags = tmpl::list<
      RadiationTransport::M1Grey::Tags::TildeE<
          Frame::Inertial, neutrinos::ElectronNeutrinos<1>>,
      RadiationTransport::M1Grey::Tags::TildeS<
          Frame::Inertial, neutrinos::ElectronNeutrinos<1>>>;
  using argument_tags = tmpl::list<>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> tilde_e,
                    const gsl::not_null<tnsr::i<DataVector, 3>*> tilde_s) {
    // FIXME should there be factors of det(J)?
    const double energy_cutoff = 1.0e-8;
    const double momentum_cutoff = 1.0e-8;
    const double energy_atmosphere = 5.0e-9;
    for (size_t s = 0; s < get(*tilde_e).size(); ++s) {
      if (std::abs(get(*tilde_e)[s]) < energy_cutoff and
          //FIXME
          square(get<0>(*tilde_s)[s]) + square(get<1>(*tilde_s)[s]) +
          square(get<2>(*tilde_s)[s]) < square(momentum_cutoff)) {
        get(*tilde_e)[s] = energy_atmosphere;
        for (size_t i = 0; i < 3; ++i) {
          tilde_s->get(i)[s] = 0.0;
        }
      }
    }
  }
};

struct EvolutionMetavars {
  static constexpr size_t volume_dim = 3;
  static constexpr dg::Formulation dg_formulation =
      dg::Formulation::StrongInertial;

  // To switch which initial data is evolved you only need to change the
  // line `using initial_data = ...;` and include the header file for the
  // solution.
  using initial_data =
      RadiationTransport::M1Grey::AnalyticData::SphericalGaussian<2>;
  static_assert(
      is_analytic_data_v<initial_data> xor is_analytic_solution_v<initial_data>,
      "initial_data must be either an analytic_data or an analytic_solution");

  // Set list of neutrino species to be used by M1 code
  using neutrino_species = tmpl::list<neutrinos::ElectronNeutrinos<1>>;

  using system = RadiationTransport::M1Grey::System<neutrino_species>;
  using temporal_id = Tags::TimeStepId;
  static constexpr bool local_time_stepping = false;
  using initial_data_tag =
      tmpl::conditional_t<is_analytic_solution_v<initial_data>,
                          Tags::AnalyticSolution<initial_data>,
                          Tags::AnalyticData<initial_data>>;
  using analytic_variables_tags = typename system::variables_tag::tags_list;
  using limiter = Tags::Limiter<
      Limiters::Minmod<3, typename system::variables_tag::tags_list>>;

  using analytic_compute =
      evolution::Tags::AnalyticSolutionsCompute<volume_dim,
                                                analytic_variables_tags>;
  using error_compute = Tags::ErrorsCompute<analytic_variables_tags>;
  using error_tags = db::wrap_tags_in<Tags::Error, analytic_variables_tags>;
  using observe_fields = tmpl::push_back<
      tmpl::append<typename system::variables_tag::tags_list,
                   typename system::primitive_variables_tag::tags_list,
                   error_tags>,
      domain::Tags::Coordinates<volume_dim, Frame::Grid>,
      domain::Tags::Coordinates<volume_dim, Frame::Inertial>>;
  using non_tensor_compute_tags =
      tmpl::list<::Events::Tags::ObserverMeshCompute<volume_dim>,
                 analytic_compute, error_compute>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DenseTrigger, DenseTriggers::standard_dense_triggers>,
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<Event, tmpl::flatten<tmpl::list<
                              Events::Completion,
                              dg::Events::field_observations<
                                  volume_dim, Tags::Time, observe_fields,
                                  non_tensor_compute_tags>,
                              Events::time_events<system>>>>,
        tmpl::pair<LtsTimeStepper, TimeSteppers::lts_time_steppers>,
        tmpl::pair<PhaseChange,
                   tmpl::list<PhaseControl::VisitAndReturn<
                                  Parallel::Phase::LoadBalancing>,
                              PhaseControl::CheckpointAndExitAfterWallclock>>,
        tmpl::pair<
            RadiationTransport::M1Grey::BoundaryConditions::BoundaryCondition<
                metavariables::neutrino_species>,
            RadiationTransport::M1Grey::BoundaryConditions::
                standard_boundary_conditions<metavariables::neutrino_species>>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                   StepChoosers::standard_step_choosers<system, false>>,
        tmpl::pair<StepChooser<StepChooserUse::Slab>,
                   StepChoosers::standard_slab_choosers<
                       system, local_time_stepping, false>>,
        tmpl::pair<StepController, StepControllers::standard_step_controllers>,
        tmpl::pair<TimeSequence<double>,
                   TimeSequences::all_time_sequences<double>>,
        tmpl::pair<TimeSequence<std::uint64_t>,
                   TimeSequences::all_time_sequences<std::uint64_t>>,
        tmpl::pair<TimeStepper, TimeSteppers::time_steppers>,
        tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                         Triggers::time_triggers>>>;
  };

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
          tmpl::at<typename factory_creation::factory_classes, Event>>>>;

  static_assert(not local_time_stepping);
  using step_actions = tmpl::flatten<tmpl::list<
      evolution::dg::Actions::ComputeTimeDerivative<volume_dim, system,
                                                    AllStepChoosers>,
      tmpl::conditional_t<
          local_time_stepping,
          tmpl::list<evolution::Actions::RunEventsAndDenseTriggers<
                         tmpl::list<evolution::dg::ApplyBoundaryCorrections<
                             EvolutionMetavars, true>,
                                    imex::ImexDenseOutput<system>>>,
                     evolution::dg::Actions::ApplyLtsBoundaryCorrections<
                         EvolutionMetavars>>,
          tmpl::list<
              evolution::dg::Actions::ApplyBoundaryCorrectionsToTimeDerivative<
                  EvolutionMetavars>,
              Actions::RecordTimeStepperData<>,
              imex::Actions::RecordImexTimeStepperData,
              evolution::Actions::RunEventsAndDenseTriggers<tmpl::list<imex::ImexDenseOutput<system>>>,
              Actions::UpdateU<>>>,
      Actions::MutateApply<FIXME_Atmosphere>,
      Actions::MutateApply<FIXME_Flatten>,
      imex::Actions::DoImplicitStep,//FIXME at end?
      Limiters::Actions::SendData<EvolutionMetavars>,
      Limiters::Actions::Limit<EvolutionMetavars>,
      Actions::MutateApply<FIXME_Atmosphere>,
      Actions::MutateApply<FIXME_Flatten>,
      Actions::MutateApply<typename RadiationTransport::M1Grey::
                               ComputeM1Closure<neutrino_species>>/*,
      Actions::MutateApply<typename RadiationTransport::M1Grey::
                               ComputeM1HydroCoupling<neutrino_species>>*/>>;

  using dg_registration_list =
      tmpl::list<observers::Actions::RegisterEventsWithObservers>;

  using initialization_actions = tmpl::list<
      Initialization::Actions::TimeAndTimeStep<EvolutionMetavars>,
      evolution::dg::Initialization::Domain<volume_dim>,
      Initialization::Actions::GrTagsForHydro<system>,
      Initialization::Actions::ConservativeSystem<system>,
      evolution::Initialization::Actions::SetVariables<
          domain::Tags::Coordinates<volume_dim, Frame::ElementLogical>>,
      Initialization::Actions::TimeStepperHistory<EvolutionMetavars>,
      imex::Actions::InitializeImex<system>,
      RadiationTransport::M1Grey::Actions::InitializeM1Tags<system>,
      Actions::MutateApply<typename RadiationTransport::M1Grey::
                               ComputeM1Closure<neutrino_species>>,
      /*Actions::MutateApply<typename RadiationTransport::M1Grey::
                               ComputeM1HydroCoupling<neutrino_species>>,*/
      Initialization::Actions::AddComputeTags<
          StepChoosers::step_chooser_compute_tags<EvolutionMetavars>>,
      ::evolution::dg::Initialization::Mortars<volume_dim, system>,
      Initialization::Actions::Minmod<3>,
      evolution::Actions::InitializeRunEventsAndDenseTriggers,
      Initialization::Actions::RemoveOptionsAndTerminatePhase>;

  using dg_element_array = DgElementArray<
      EvolutionMetavars,
      tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,

          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions, system>>,

          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,

          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<Actions::RunEventsAndTriggers, Actions::ChangeSlabSize,
                         step_actions, Actions::AdvanceTime,
                         PhaseControl::Actions::ExecutePhaseChange>>>>;

  template <typename ParallelComponent>
  struct registration_list {
    using type =
        std::conditional_t<std::is_same_v<ParallelComponent, dg_element_array>,
                           dg_registration_list, tmpl::list<>>;
  };

  using component_list =
      tmpl::list<observers::Observer<EvolutionMetavars>,
                 observers::ObserverWriter<EvolutionMetavars>,
                 dg_element_array>;

  using const_global_cache_tags = tmpl::list<initial_data_tag>;

  static constexpr Options::String help{
      "Evolve the M1Grey system (without coupling to hydro).\n\n"};

  static constexpr std::array<Parallel::Phase, 5> default_phase_order{
      {Parallel::Phase::Initialization,
       Parallel::Phase::InitializeTimeStepperHistory, Parallel::Phase::Register,
       Parallel::Phase::Evolve, Parallel::Phase::Exit}};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling,
    &setup_memory_allocation_failure_reporting,
    &disable_openblas_multithreading,
    &domain::creators::register_derived_with_charm,
    &domain::creators::time_dependence::register_derived_with_charm,
    &domain::FunctionsOfTime::register_derived_with_charm,
    &Parallel::register_derived_classes_with_charm<
        RadiationTransport::M1Grey::BoundaryCorrections::BoundaryCorrection<
            metavariables::neutrino_species>>,
    &Parallel::register_factory_classes_with_charm<metavariables>};

static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
