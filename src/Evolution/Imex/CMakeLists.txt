# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Imex)

add_spectre_library(${LIBRARY} INTERFACE)

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  NamespaceDocs.hpp
  )

target_link_libraries(
  ${LIBRARY}
  INTERFACE
  DataStructures
  Time
  Utilities
  )

add_subdirectory(Protocols)
add_subdirectory(Tags)
