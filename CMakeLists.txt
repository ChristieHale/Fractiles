# CMake project file for the fractile program (same compiler options as HAZ)

cmake_minimum_required (VERSION 3.10)
project (Fractiles Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif ()

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # Fixed form with long lines, and static storage of local variables
  set (CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -fno-automatic")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
endif ()

add_executable (Fractiles
  "Fractiles Files/fract_main.f"
  "Fractiles Files/rd_flt.f"
  "Fractiles Files/rd_input.f"
  "Fractiles Files/random.f")
