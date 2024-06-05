cmake_minimum_required(VERSION 3.28)
project(PF)

set(CMAKE_CXX_STANDARD 17)

add_executable(PF main.cpp)

include_directories(/opt/homebrew/opt/armadillo/include)

link_directories(/opt/homebrew/opt/armadillo/lib)

target_link_libraries(PF /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)
