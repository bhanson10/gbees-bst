cmake_minimum_required(VERSION 3.28)
project(EKF)

set(CMAKE_CXX_STANDARD 17)

add_executable(EKF main.cpp)

include_directories(/opt/homebrew/opt/armadillo/include)

link_directories(/opt/homebrew/opt/armadillo/lib)

target_link_libraries(EKF /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)