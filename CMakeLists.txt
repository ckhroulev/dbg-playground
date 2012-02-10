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

set (DB_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
# Force the default cache entry to obey Pism_BUILD_TYPE
set (CMAKE_BUILD_TYPE ${DB_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
mark_as_advanced (CMAKE_BUILD_TYPE) # hide the default entry

set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math")

include_directories (src src/dem src/io ${GSL_INCLUDES} ${NETCDF_INCLUDES} ${MPI_INCLUDE_PATH})
set (EXTERNAL_LIBS ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_LIBRARIES})

# Tell GSL to use inline versions of some functions:
add_definitions(-DHAVE_INLINE -DGSL_RANGE_CHECK_OFF)

add_executable (db
  src/drainagebasin.cc
  src/init.cc
  src/streamline.cc
  src/dem/DEM.cc
  src/basins.cc
  src/io/IO.cc
  src/io/PISMNC3File.cc
  src/io/PISMNCFile.cc)
target_link_libraries(db ${EXTERNAL_LIBS})