# Distributed under the MIT License.
# See LICENSE.txt for details.

add_subdirectory(DataStructures)
add_subdirectory(Domain)
add_subdirectory(ErrorHandling)
add_subdirectory(Executables)
add_subdirectory(Numerical)

# Propagate core labraries up to global scope.
set(SPECTRE_LIBRARIES
    "${SPECTRE_CORE_LIBRARIES};${SPECTRE_LIBRARIES}" PARENT_SCOPE)
