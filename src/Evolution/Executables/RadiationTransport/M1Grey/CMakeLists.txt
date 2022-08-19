# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  GeneralRelativitySolutions
  Hydro
  IO
  Imex
  Informer
  Limiters
  LinearOperators
  MathFunctions
  M1Grey
  M1GreyAnalyticData
  M1GreySolutions
  Options
  Parallel
  PhaseControl
  Time
  Utilities
  )

add_spectre_parallel_executable(
  EvolveM1GreyConstantM1
  EvolveM1Grey
  Evolution/Executables/RadiationTransport/M1Grey
  EvolutionMetavars
  "${LIBS_TO_LINK}"
  )

add_spectre_parallel_executable(
  EvolveM1GreyHomogeneousSphere
  EvolveM1Grey
  Evolution/Executables/RadiationTransport/M1Grey
  EvolutionMetavars
  "${LIBS_TO_LINK}"
  )
