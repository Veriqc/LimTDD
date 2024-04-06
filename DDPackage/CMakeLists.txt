add_subdirectory(src)
# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dd)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/operations)

# File groups
file(GLOB SOURCES
    "dd/*.cpp"
    "operations/*.cpp"
)

# Create the library
add_library(Limdd_lib ${SOURCES})
target_include_directories(Limdd_lib PUBLIC ${xtensor_INCLUDE_DIRS} ${xtl_INCLUDE_DIRS})
