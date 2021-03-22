cmake_minimum_required(VERSION 3.8...3.15)

project (NeoPZExamples LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Finding the pz package
find_package(NeoPZ REQUIRED)
add_subdirectory(SpecialMaps)
add_subdirectory(F17DirectionalRefinement)

