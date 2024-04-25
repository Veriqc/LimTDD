add_subdirectory(src)
# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dd)

# File groups
file(GLOB SOURCES
    "dd/*.cpp"
)

# Create the library
add_library(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${xtensor_INCLUDE_DIRS} ${xtl_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${CORE_NAME})