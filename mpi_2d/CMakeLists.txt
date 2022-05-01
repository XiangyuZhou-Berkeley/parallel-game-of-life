cmake_minimum_required(VERSION 3.22)
project(gameOfLife)

set (CMAKE_CXX_STANDARD 14)
find_package(MPI REQUIRED)
add_executable(mpi2d main.cpp mpi2d.cpp common.h)
# add_library(board serial.cpp)
target_link_libraries(mpi2d PRIVATE MPI::MPI_CXX)