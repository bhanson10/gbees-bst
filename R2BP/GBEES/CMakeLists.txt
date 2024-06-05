cmake_minimum_required(VERSION 3.28)
project(GBEES)

set(CMAKE_CXX_STANDARD 17)

add_executable(GBEES main.cpp BST.h)

include_directories(/opt/homebrew/opt/armadillo/include)

link_directories(/opt/homebrew/opt/armadillo/lib)

target_link_libraries(GBEES /opt/homebrew/opt/armadillo/lib/libarmadillo.dylib)
