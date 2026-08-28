cmake_minimum_required(VERSION 3.10)
project(Limtdd) 


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(LIMTDD_ROOT_OF_UNIT 8 CACHE STRING "Order N of the root-of-unity map set")
set(LIMTDD_T_GATE_ORDER 8 CACHE STRING "Order of diag(1, exp(2*pi*i/order)) used for QASM t/tdg gates")
if(NOT LIMTDD_ROOT_OF_UNIT MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "LIMTDD_ROOT_OF_UNIT must be a positive integer")
endif()
if(NOT LIMTDD_T_GATE_ORDER MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "LIMTDD_T_GATE_ORDER must be a positive integer")
endif()
add_compile_definitions(
    LIMTDD_ROOT_OF_UNIT=${LIMTDD_ROOT_OF_UNIT}
    LIMTDD_T_GATE_ORDER=${LIMTDD_T_GATE_ORDER}
)
message(STATUS "LimTDD precision N: ${LIMTDD_ROOT_OF_UNIT}")
message(STATUS "LimTDD QASM t-gate order: ${LIMTDD_T_GATE_ORDER}")

set(MY_LOCAL_INSTALL_DIR "$ENV{HOME}/thu/LimTDDexpr")
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