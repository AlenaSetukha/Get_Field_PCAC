cmake_minimum_required(VERSION 3.20)

project(Get_field_Ideal CXX)

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

set(SOURCE
    src/cell.cpp
    src/common_type.cpp
    src/e0.cpp
    src/ED_Par.cpp
    src/element_geom.cpp
    src/f_par.cpp
    src/get_einc.cpp
    src/get_j.cpp
    src/integral_par.cpp
    src/integrals_analytic.cpp
    src/K.cpp
    src/kernel_lib.cpp
    src/Num_Par.cpp
    src/write_to_file.cpp
    )

add_library(field_lib STATIC ${SOURCE})


add_executable(get_field main_field_ideal.cpp)

if(OPENMP_FOUND)
    target_link_libraries(get_field PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(get_field PUBLIC field_lib)