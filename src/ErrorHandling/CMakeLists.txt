# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ErrorHandling)

set(LIBRARY_SOURCES
    AbortWithErrorMessage.cpp
    FloatingPointExceptions.cpp)

add_library(${LIBRARY} ${LIBRARY_SOURCES})

# ${LIBRARY} is considered a core library on which all executables depend
set(SPECTRE_CORE_LIBRARIES "${SPECTRE_CORE_LIBRARIES};${LIBRARY}" PARENT_SCOPE)
