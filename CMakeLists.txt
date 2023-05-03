cmake_minimum_required(VERSION 3.10)

project(pir_scheme LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
# Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_COMPILER "clang++")

#option(CMAKE_BUILD_TYPE "Building Type" Debug)

# Set build type to `Release` if non was specified:
# (cf. https://gitlab.kitware.com/cmake/community/wikis/FAQ#how-can-i-change-the-default-build-mode-and-see-it-reflected-in-the-gui)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)
#set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O2 -Wall -march=native")

option(EXTRA_TESTS "EXTRA_TESTS: enable extra tests" ON)

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "EXTRA_TESTS:      ${EXTRA_TESTS}")


# Using intel hexl
find_package(HEXL 1.2.5
    HINTS ${HEXL_HINT_DIR}
    REQUIRED)
if (NOT TARGET HEXL::hexl)
    message(FATAL_ERROR "TARGET HEXL::hexl not found")
endif()

message(STATUS "HEXL_VERSION:     ${HEXL_VERSION}")
message(STATUS "HEXL_HINT_DIR:    ${HEXL_HINT_DIR}")
message(STATUS "HEXL_INCLUDE_DIRS: ${HEXL_INCLUDE_DIRS}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

add_executable(test-first-dimension ${PROJECT_SOURCE_DIR}/tests/test-first-dimension.cpp)
target_include_directories(test-first-dimension PUBLIC ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(test-first-dimension PRIVATE pir HEXL::hexl)

if (EXTRA_TESTS)
    add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif()
