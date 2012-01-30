project (drainagebasin)

cmake_minimum_required (VERSION 2.6.2)

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

find_package (GSL REQUIRED)
find_package (NetCDF REQUIRED)
find_package (MPI REQUIRED)

set (DEFAULT_PEDANTIC_CFLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Wsign-compare -Wno-redundant-decls -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -funit-at-a-time")
set (DEFAULT_PEDANTIC_CXXFLAGS "${DEFAULT_PEDANTIC_CFLAGS} -Woverloaded-virtual")

set (CMAKE_BUILD_TYPE "Debug")
set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")

include_directories (${GSL_INCLUDES} ${NETCDF_INCLUDES} ${MPI_INCLUDE_PATH})
set (EXTERNAL_LIBS ${GSL_LIBRARIES})# ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})

add_executable (db drainagebasin.cc system.cc)
target_link_libraries(db ${EXTERNAL_LIBS})