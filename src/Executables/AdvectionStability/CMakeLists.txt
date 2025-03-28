# Distributed under the MIT License.
# See LICENSE.txt for details.

set(EXECUTABLE AdvectionStability)

add_spectre_executable(
  ${EXECUTABLE}
  EXCLUDE_FROM_ALL
  AdvectionStability.cpp
  )

# FIXME
target_link_libraries(
  ${EXECUTABLE}
  PRIVATE
  #Boost::boost
  #Boost::program_options
  DataStructures
  DiscontinuousGalerkin
  ErrorHandling
  LinearAlgebra
  Printf
  Spectral
  Time
  Utilities
  )
