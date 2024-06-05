cmake_minimum_required(VERSION 3.28)
project(UKF)

set(CMAKE_CXX_STANDARD 17)

add_executable(UKF main.cpp)

include_directories(/opt/homebrew/opt/armadillo/include)

link_directories(/opt/homebrew/opt/armadillo/lib)

target_link_libraries(UKF /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)