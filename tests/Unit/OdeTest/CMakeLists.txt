# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_Ode")

set(LIBRARY_SOURCES
  Test_Ode.cpp
  )

add_test_library(
  ${LIBRARY}
  "OdeTest"
  "${LIBRARY_SOURCES}"
  ""
  )

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  Time
  TimeStepperHelpers
  )

add_dependencies(
  ${LIBRARY}
  module_GlobalCache
  )
