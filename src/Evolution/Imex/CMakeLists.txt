# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Imex)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Mode.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Mode.hpp
  NamespaceDocs.hpp
  SolveImplicitSector.hpp
  )

target_link_libraries(
  ${LIBRARY}
  INTERFACE
  DataStructures
  LinearSolver
  RootFinding
  Time
  Utilities
  PUBLIC
  Options
  )

add_subdirectory(Protocols)
add_subdirectory(Tags)
