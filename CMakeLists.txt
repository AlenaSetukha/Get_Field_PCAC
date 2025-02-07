cmake_minimum_required(VERSION 3.20)

project(PCAC_LIB CXX)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "release")
endif()

#compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -Wall -Wextra -pedantic -fsanitize=leak, -fsanitize=undefinded, -fsanitize=memory")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -flto -march=native")

# packages
find_package(OpenMP QUIET)

include_directories(include)
include_directories(include/ED_Field)
include_directories(include/ED_Operators)
include_directories(include/grids)
include_directories(include/integrators)
include_directories(include/io)
include_directories(include/param)

set(SOURCE
    src/grids/Cell.cpp
    src/grids/Cell_Triangle.cpp
    src/grids/common_type_triangle.cpp
    src/grids/common_type.cpp
    src/ED_Field/E0.cpp
    src/ED_Field/get_einc.cpp
    src/ED_Field/get_hinc.cpp
    src/ED_Field/get_j.cpp
    src/io/write_to_file.cpp
    src/param/ED_Par.cpp
    src/param/Kernel_Par.cpp
    src/param/Integral_Par.cpp
    src/param/Num_Par.cpp
    src/element_geom.cpp
    src/integrals_analytic.cpp
    src/kernel_lib.cpp
    )

add_library(pcac_lib STATIC ${SOURCE})




add_executable(test_PCAC ./tests/field/main_field.cpp)

if(OPENMP_FOUND)
    target_link_libraries(test_PCAC PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(test_PCAC PUBLIC pcac_lib)