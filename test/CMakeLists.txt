# List your test source files here
add_definitions(-DPROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
set(TEST_SOURCES
    # test_whole.cpp
    # test_tn.cpp
    # test_tdd.cpp
    # test_bv.cpp
    # test_qft.cpp
    test_cl.cpp
    test_io.cpp
    # test_clt.cpp
    # test_id.cpp
    # Add new test files here
)
message(STATUS "CORE_NAME is: ${CORE_NAME}")
message(STATUS "CORE_INCLUDE_BUILD_DIR is: ${CORE_INCLUDE_BUILD_DIR}")
if(NOT CORE_INCLUDE_DIRS)
  message(STATUS "No include directories found for target: ${CORE_NAME}")
endif()
message(STATUS "Include directories for ${CORE_NAME}: ${CORE_INCLUDE_DIRS}")

# Loop through each source file to create an executable and link libraries
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAME} ${CORE_NAME})
endforeach()
