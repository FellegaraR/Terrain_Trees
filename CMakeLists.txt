cmake_minimum_required(VERSION 2.8)

project(Terrain_Trees)
set(terrain_trees_VERSION_MAJOR 1)
set(terrain_trees_VERSION_MINOR 0)

SET(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/Modules/)
MESSAGE(STATUS "local cmake root: " ${CMAKE_MODULE_PATH})

set(LANGUAGE C++)
#set(CMAKE_VERBOSE_MAKEFILE on)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

set(CMAKE_BUILD_TYPE Release)

find_package(Boost REQUIRED)
find_package(Doxygen REQUIRED)

#if(CMAKE_COMPILER_IS_GNUCXX) #do not work on mac
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11") ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s -lrt") ## Strip binary
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s") ## Strip binary
endif()

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()

################################################################## prova non funzionante
#add_custom_target(clean-cmake-files
#   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
#)
## clean-all.cmake
#set(cmake_generated ${CMAKE_BINARY_DIR}/CMakeCache.txt
#                    ${CMAKE_BINARY_DIR}/cmake_install.cmake
#                    ${CMAKE_BINARY_DIR}/Makefile
#                    ${CMAKE_BINARY_DIR}/CMakeFiles
#)
#foreach(file ${cmake_generated})
#  if (EXISTS ${file})
#     file(REMOVE_RECURSE ${file})
#  endif()
#endforeach(file)
#add_custom_target(clean-all
#   COMMAND ${CMAKE_BUILD_TOOL} clean
#   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
#)
############################################################### ===============================

include_directories(sources)
add_subdirectory(sources)

### MAIN OF THE LIBRARY
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
add_executable(terrain_trees targets/main.cpp)
target_link_libraries(terrain_trees _terrain_trees_apps)


### TEST THE LIBRARY
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/test )
add_executable(test_index_generation EXCLUDE_FROM_ALL tests/test_spatial_index_generation.cpp)
target_link_libraries(test_index_generation _terrain_trees_core)

add_executable(test_curvature EXCLUDE_FROM_ALL tests/test_curvature.cpp)
target_link_libraries(test_curvature _terrain_trees_apps)

add_executable(test_soup EXCLUDE_FROM_ALL tests/test_soup.cpp)
target_link_libraries(test_soup _terrain_trees_apps)

add_executable(test_batched_queries EXCLUDE_FROM_ALL tests/test_batched_queries.cpp)
target_link_libraries(test_batched_queries _terrain_trees_apps)

add_executable(test_terrain_feature_extraction EXCLUDE_FROM_ALL tests/test_terrain_feature_extraction.cpp)
target_link_libraries(test_terrain_feature_extraction _terrain_trees_apps)

add_executable(test_spatial_queries EXCLUDE_FROM_ALL tests/test_spatial_queries.cpp)
target_link_libraries(test_spatial_queries _terrain_trees_apps)

add_custom_target(tests DEPENDS test_index_generation test_curvature test_soup test_batched_queries test_terrain_feature_extraction test_spatial_queries)


