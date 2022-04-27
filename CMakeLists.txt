project(gameOfLife)

set (CMAKE_CXX_STANDARD 14)

add_executable(serial main.cpp)
add_library(board serial.cpp)
target_link_libraries(serial board)