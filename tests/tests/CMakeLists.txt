cmake_minimum_required(VERSION 3.21)

set(gtesthelper
    ${CMAKE_CURRENT_SOURCE_DIR}/c_api_test_helpers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/testhelpers.hpp
)

macro(cavc_add_test name)
    add_executable(${name} ${name}.cpp ${gtesthelper})

    target_include_directories(${name}
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include)

    target_link_libraries(${name}
    PRIVATE
        PolylineFactory
        ${CAVC_CPP_HEADER_ONLY_LIB}
        ${CAVC_C_API_LIB}
        GTest::gmock_main)

    gtest_discover_tests(${name}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DISCOVERY_MODE PRE_TEST)
endmacro()

cavc_add_test(TEST_sample)
cavc_add_test(TEST_cavc_pline)
cavc_add_test(TEST_cavc_pline_function)
cavc_add_test(TEST_cavc_parallel_offset)
cavc_add_test(TEST_cavc_combine_plines)
cavc_add_test(TEST_staticspatialindex)
cavc_add_test(TEST_cavc_api_regression)
cavc_add_test(TEST_cavc_offset_islands)
cavc_add_test(TEST_cavc_internal_slice_view)
cavc_add_test(TEST_cavc_parallel_offset_fuzz_regression)
cavc_add_test(TEST_cavc_compound_combine)
