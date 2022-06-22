# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_Imex")

set(LIBRARY_SOURCES
  Test_SolveImplicitSector.cpp
  )

add_subdirectory(Tags)

add_test_library(
  ${LIBRARY}
  "Evolution/Imex/"
  "${LIBRARY_SOURCES}"
  ""
  )

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DataStructures
  DataStructuresHelpers
  Framework
  Imex
  ImexHelpers
  Time
  Utilities
  )

add_dependencies(
  ${LIBRARY}
  module_GlobalCache
  )
