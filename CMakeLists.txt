cmake_minimum_required(VERSION 3.0.2)
project (distinct_squares)

set(CXX_STANDARD c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=${CXX_STANDARD}")


include(ExternalProject)

# Global dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# Usage: find_or_download_package(PACKAGE_NAME, FOUND_FLAG_PREFIX, TARGET_PREFIX)
include(find_or_download_package)
find_or_download_package(Glog GLOG glog)
find_or_download_package(SDSL SDSL sdsl)
find_or_download_package(GFlags GFLAGS gflags)


add_executable(distinct_squares squares.cpp)

target_link_libraries(distinct_squares
    glog
    sdsl
	)

