cmake_minimum_required(VERSION 3.10)
project(Limdd) 


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CORE_INCLUDE_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DDPackage/mqt-core")
message(STATUS "add submodule")
find_package(xtensor REQUIRED)
if(NOT xtensor_FOUND)
    message(FATAL_ERROR "xtensor not found. Please make sure xtensor is installed.")
else()
    message(STATUS "xtensor include directories: ${xtensor_INCLUDE_DIRS}")
endif()

find_package(xtl REQUIRED)
if(NOT xtl_FOUND)
    message(FATAL_ERROR "xtl not found. Please make sure xtl is installed.")
else()
    message(STATUS "xtl include directories: ${xtl_INCLUDE_DIRS}")
endif()

include_directories(${PROJECT_SOURCE_DIR}/DDPackage)

add_subdirectory(DDPackage)
add_subdirectory(test)
