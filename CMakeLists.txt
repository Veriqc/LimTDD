cmake_minimum_required(VERSION 3.10)
project(Limtdd) 


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(MY_LOCAL_INSTALL_DIR "$ENV{HOME}/local")
include_directories(
    ${MY_LOCAL_INSTALL_DIR}/include
)


set(CORE_INCLUDE_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DDPackage/core")
set(CORE_NAME ${PROJECT_NAME}_core)

message(STATUS "Using local xtensor and xtl from: ${MY_LOCAL_INSTALL_DIR}")
# find_package(xtensor REQUIRED)
# if(NOT xtensor_FOUND)
#     message(FATAL_ERROR "xtensor not found.")
# else()
#     message(STATUS "xtensor include directories: ${xtensor_INCLUDE_DIRS}")
# endif()

# find_package(xtl REQUIRED)
# if(NOT xtl_FOUND)
#     message(FATAL_ERROR "xtl not found.")
# else()
#     message(STATUS "xtl include directories: ${xtl_INCLUDE_DIRS}")
# endif()

include_directories(${PROJECT_SOURCE_DIR}/DDPackage)

add_subdirectory(DDPackage)
add_subdirectory(test)