#CMakeLists.txt to build project using cmake
#run $cmake CMakeLists.txt

cmake_minimum_required(VERSION 2.6)

project(pgm)
include_directories("projet/src/")

file(
GLOB_RECURSE
source_files
projet/src/*.[cht]pp
)

add_executable(pgm ${source_files})

SET(CMAKE_BUILD_TYPE Debug)
ADD_DEFINITIONS( -O2 -Wall -Wextra -std=c++11 -Wno-comment -Wno-unused-parameter -Wno-unused-function -DCORRECTION -fopenmp)

TARGET_LINK_LIBRARIES(pgm -lm -ldl -fopenmp)


