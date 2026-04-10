cmake_minimum_required(VERSION 3.16)
project(creamer_flow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(creamer_flow creamer_flow.cpp)
set_target_properties(creamer_flow PROPERTIES OUTPUT_NAME creamer_flow_plan)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
  target_link_libraries(creamer_flow PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(creamer_flow PRIVATE CREAMER_HAVE_OPENMP=1)
endif()

if (MSVC)
  target_compile_options(creamer_flow PRIVATE /O2)
else()
  target_compile_options(creamer_flow PRIVATE -O3 -march=native)
endif()
