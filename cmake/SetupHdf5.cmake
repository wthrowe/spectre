# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(HDF5 REQUIRED C)

message(STATUS "HDF5 libs: " ${HDF5_LIBRARIES})
message(STATUS "HDF5 incl: " ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
list(APPEND SPECTRE_LIBRARIES ${HDF5_LIBRARIES})
