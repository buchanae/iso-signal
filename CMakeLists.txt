cmake_minimum_required( VERSION 2.6 )
project( iso_signal )

include_directories(
    include
    /Users/abuchanan/chang/duncan/include
    /Users/abuchanan/bamtools/include
)
set( 
    SRCS
    src/JunctionIndex.cpp
    src/StackReader.cpp
    src/Alignment.cpp
    src/Coverage.cpp
    src/helpers.cpp
)
set(
    LIBS
    uriparser
    /Users/abuchanan/chang/duncan/build/libduncan.a
    /Users/abuchanan/bamtools/lib/libbamtools.dylib
)
add_executable( iso_signal ${SRCS} src/main.cpp)
target_link_libraries( iso_signal ${LIBS} )

option( build_tests "Build tests" ON )
if( build_tests )

    enable_testing()
    set(TEST_SRCS
      tests/JunctionIndex.cpp
      tests/StackReader.cpp
      tests/Alignment.cpp
      tests/Coverage.cpp
      tests/BamReader.cpp
      tests/JunctionFilter.cpp
      tests/helpers.cpp
    )
    add_custom_command(
        OUTPUT dummies
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/dummies"
            "${CMAKE_CURRENT_BINARY_DIR}/dummies"
            DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/dummies"
    )
    add_custom_target(TestDummies ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/dummies")
    include_directories(
        libs/gmock/include
        libs/gmock/gtest/include
        libs/gmock
        libs/gmock/gtest
    )
    add_library(
        gmock
        libs/gmock/src/gmock-all.cc
        libs/gmock/gtest/src/gtest-all.cc
        libs/gmock/gtest/src/gtest_main.cc
    )
    add_executable(
        tests
        ${SRCS}
        ${TEST_SRCS}
    )
    target_link_libraries(
        tests
        ${LIBS}
        gmock
    )
    add_test( all-tests tests )
endif()
