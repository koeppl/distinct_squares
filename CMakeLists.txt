cmake_minimum_required(VERSION 3.0.2)
project (distinct_squares)

set(CXX_STANDARD c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=${CXX_STANDARD}")


include(ExternalProject)

# Global dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

# Usage: find_or_download_package(PACKAGE_NAME, FOUND_FLAG_PREFIX, TARGET_PREFIX)
include(find_or_download_package)
find_or_download_package(Glog GLOG glog)
find_or_download_package(SDSL SDSL sdsl)
find_or_download_package(GFlags GFLAGS gflags)
find_or_download_package(GTest GTEST gtest)



add_executable(check_file check_file.cpp)
target_link_libraries(check_file
    glog
    sdsl
	)
add_executable(distinct_squares squares.cpp)
target_link_libraries(distinct_squares
	gflags
    glog
    sdsl
	)
add_executable(test_squares test_squares.cpp)
target_link_libraries(test_squares
	gtest
    glog
    sdsl
	${TEST_TARGET_DEPS}
	)

enable_testing()
add_test(NAME check COMMAND test_squares)


