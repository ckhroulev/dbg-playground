project(foobar)

cmake_minimum_required (VERSION 2.6.2)

# Put executables in the build directory:
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

find_package(OpenMP REQUIRED)

add_executable(foobar_par foobar.c)
set_target_properties (foobar_par PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})

add_executable(foobar foobar.c)
