# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dd)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/operations)

# File groups
file(GLOB SOURCES
    "dd/*.cpp"
    "operations/*.cpp"
)

# Create the library
add_library(Limdd_lib INTERFACE ${SOURCES})
target_include_directories(Limdd_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/graph ${XTENSOR_INCLUDE_DIRS} ${XTL_INCLUDE_DIRS})