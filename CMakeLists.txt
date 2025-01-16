cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_COMPILER mpic++)

project(MPIReduction)

find_package(MPI REQUIRED)

add_executable(reduce main.cpp)
find_package(OpenMP REQUIRED)

target_compile_options(reduce PRIVATE -fopenmp)
target_link_libraries(reduce OpenMP::OpenMP_CXX)
target_link_libraries(reduce MPI::MPI_CXX)